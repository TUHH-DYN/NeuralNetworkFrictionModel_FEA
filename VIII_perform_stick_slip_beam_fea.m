function [output] = performStickSlipBeamFEA(hmax, fmodelType, tfinal, dt, eps, beta)

%% Define friction model

if strcmp(fmodelType, 'neuralNet')

    % Load trained friction model
    regModel = load('optimizedRegNet.mat');
    fmodel.trainedModel.RegressionNeuralNetwork = regModel.regNet;

end

% Set further computational model parameters
vbelt = 0.02;               % Belt velocity (default: 0.02)
Fn = 100.0;                 % Normal force (default: 100.0)

% Define friction model
fmodel.type = fmodelType;
fmodel.a    = 5.0;          % Control parameter for the negative gradient of the friction curve (default: 5.0)
fmodel.v0   = 0.3;          % Reference velocity if exponential decay is used   (default: 0.3) 
fmodel.musd = 2.5;          % Ratio static to dynamic fric coeff, mu_st/mu_d    (default: 2.5)
fmodel.mud  = 0.2;          % Dynamic coeff of friction, mu_d                   (default: 0.2)
fmodel.muv  = 0.001;        % Linear strengthening parameter                    (default: 0.001)
fmodel.eps  = eps;          % Threshold value for sliding velocity              (default: 1e-3)

% Compute maximum static friction force
Fcrit = Fn * fmodel.mud * fmodel.musd;

%% Create geometry
L = 2.0;    % Beam length
h = 0.01;   % Beam height
I = h^3/12; % Moment of inertia

x = [0.; L; L; 0.];
y = [-h/2; -h/2; h/2; h/2];

g = decsg([3, 4, x' y']');

%% Modal analysis

% Create modal analysis model for 2-D plane-stress problem
modelModal = createpde('structural', 'modal-planestress');

geometryFromEdges(modelModal, g);

mesh = generateMesh(modelModal, "Hmax", hmax)

% Define material properties
E   = 2e11;
nu  = 0.3;
rho = 1000.;

structuralProperties(modelModal, "YoungsModulus", E, "PoissonsRatio", nu, "MassDensity", rho);

% Specify left edge of the beam as fixed boundary
structuralBC(modelModal, "Edge", 4, "Constraint", "fixed");

% Solve problem for given frequency range
modalRes = solve(modelModal, 'FrequencyRange', [-0.1, 1e4]');

%% Static analysis

% Create static 2-D plane-stress model
modelStatic = createpde('structural', 'static-planestress');

% Use same geometry and mesh
geometryFromEdges(modelStatic, g);
modelStatic.Mesh = mesh;

% Specify material properties
structuralProperties(modelStatic, "YoungsModulus", E, "PoissonsRatio", nu, "MassDensity", rho); 

% Specify same constraint on left end of the beam
structuralBC(modelStatic, 'Edge', 4, 'Constraint', 'fixed');

% Apply static vertical load on the right side of the beam
structuralBoundaryLoad(modelStatic, 'Vertex', 2, 'Force', [0; 0.98*Fcrit]);

% Solve static model
staticRes = solve(modelStatic);

%% Transient analysis

% Create transient 2-D plane-stress model
modelTransient = createpde("structural", "transient-planestress");

% Use same geometry and mesh
geometryFromEdges(modelTransient, g);
modelTransient.Mesh = mesh;

% Specify material properties
structuralProperties(modelTransient, "YoungsModulus", E, "PoissonsRatio", nu, "MassDensity", rho)

% Specify same constraint on left end of the beam
structuralBC(modelTransient, 'Edge', 4, 'Constraint', 'fixed');

% Define initial conditions (from static solution)
% structuralIC(model, "Displacement", zeros(2,1), "Velocity", zeros(2,1));
structuralIC(modelTransient, staticRes);

% Introduce stiffness proportional damping
structuralDamping(modelTransient, "Alpha", 0.0, "Beta", 0.0);

% Set up time stepping
tStart  = 0.0;
tEnd    = tfinal;
tStep   = dt; % Time step size (default: 1e-4)

nSteps = size((tStart : tStep : tEnd), 2);
recordSteps = 240;

history.t  = zeros(1, nSteps);
history.uy = zeros(1, nSteps);
history.Fy = zeros(1, nSteps);

% Initialize friction state to sticktion
slip = false;

model.SolverOptions.ReportStatistics  = 'off';
model.SolverOptions.AbsoluteTolerance = 1.e-6; % Default: 1e-6
model.SolverOptions.RelativeTolerance = 1.e-3; % Default: 1e-3

for step = 1 : nSteps-1

    if mod(step, recordSteps) == 0; fprintf('t = %f\n', history.t(step)); end

    tListTmp = [step-1 step] .* tStep;

    % Solve structural model for current time fraction
    if slip == false
        rhs = solve(modelTransient, tListTmp);
    else
        rhs = solve(modelTransient, tListTmp, "ModalResults", modalRes);
    end

    % Append partial solution to global solution vectors
    history.t (step+1) = rhs.SolutionTimes(end);
    history.uy(step+1) = interpolateDisplacement(rhs, [L; -h/2.]).uy(end);
    
    % Compute reaction forces due to current deflection at beam tip
    FreactBernoulli = -(history.uy(step+1)*3*E*I)/(L^3);
    Freact = FreactBernoulli;

    % Compute relative sliding velocity
    vy = (history.uy(step+1)-history.uy(step)) / tStep;
    vrel = vy - vbelt;

    % Update initial conditions
    structuralIC(modelTransient, rhs);

    if slip == false

        if abs(Freact) >= Fcrit

            % Stick -> Slip
            fprintf("Transition stick -> slip at t = " + num2str(history.t(step+1)) + " s\n");
            slip = true;
            structuralBC(modelTransient, "Vertex", 2, "YDisplacement", []);
            structuralDamping(modelTransient, "Alpha", 0.0, "Beta", beta);
          
        else

            % Stick
            slip = false;
            dy = vbelt*tStep;
            structuralBC(modelTransient, "Vertex", 2, "YDisplacement", history.uy(step+1)+dy);

        end

    elseif slip == true

        if (abs(vrel) < fmodel.eps) && (abs(Freact) < Fcrit)

            % Slip -> Stick
            fprintf("Transition slip -> stick at t = " + num2str(history.t(step+1)) + " s\n");
            slip = false;
            structuralDamping(modelTransient, "Alpha", 0.0, "Beta", 0.0);
            structuralBoundaryLoad(modelTransient, "Vertex", 2, "Force", []);
            structuralBC(modelTransient, "Vertex", 2, "YDisplacement", history.uy(step+1)+vbelt*tStep);

        else

            % Slip
            slip = true;

            % Compute sliding friction force at the beam tip
            Ff = getFrictionForce(Fn, vrel, fmodel);
            history.Fy(step+1) = Ff;

            structuralBoundaryLoad(modelTransient, "Vertex", 2, "Force", [0. Ff]);

        end

    end

end

%% Write results
save("output_hmax" + num2str(hmax) + ...
    "_" + fmodelType + ...
    "_tfinal" + num2str(tfinal) + ...
    "_dt" + num2str(dt) + ...
    "_eps" + num2str(eps) + ...
    "_b" + num2str(beta) + ".mat");

output = 1;

end
