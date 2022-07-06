% Time stepping convergence study for a 2D transient cantilever beam Finite
% Element model
%
% Note: requires function getFrictionForce.m

%% Initialization
clear
close all

% Load custom colormaps
load("colors\viridis.mat")
load("colors\Set1.mat")
set(0, "DefaultAxesColorOrder", Set1)

% Set some parameter values
hmax    = 0.02; % Target maximum element edge size  [m]
tfinal  = 1.30; % Final time                        [s]
beta    = 1e-2; % Damping

%% Define friction model

% Set further computational model parameters
vbelt       = 0.02;         % Belt velocity [m/s]
Fn          = 100.0;        % Normal force  [N]

% Define friction model
fmodel.type = 'exponential';
fmodel.a    = 5.0;          % Control parameter for the negative gradient of the friction curve
fmodel.v0   = 0.3;          % Reference velocity if exponential decay is used
fmodel.mus  = 0.5;          % Static coefficient of friction
fmodel.muk  = 0.2;          % Kinetic coefficient of friction
fmodel.eps  = 1e-4;         % Threshold value for sliding velocity [m/s]

% Compute maximum static friction force
Fcrit = Fn * fmodel.mus;

%% Define geometry

L = 2.0;    % Beam length       [m]
h = 0.01;   % Beam height       [m]
I = h^3/12; % Moment of inertia

x = [0.; L; L; 0.];
y = [-h/2; -h/2; h/2; h/2];

g = decsg([3, 4, x' y']');

%% Static analysis
%  Static results may be used as an initial condition for the Finite
%  Element model according to a state just before the first stick-slip
%  transition occurs. This procedure saves time by skipping the 'boring'
%  first stick phase. It also avoids the breakaway point being a multiple
%  of the time step size.

fprintf('\nPerforming static analysis of the cantilever beam Finite Element model...');

% Create static analysis model for 2-D plane-stress problem
modelStatic = createpde('structural', 'static-planestress');

% Create geometry
geometryFromEdges(modelStatic, g);

% Create mesh
mesh = generateMesh(modelStatic, "Hmax", hmax);

% Define material properties
E   = 2e11;     % Young's modulus   [Pa]
nu  = 0.3;      % Poisson's ratio   [-]
rho = 1000.;    % Mass density      [kg/m^3]

% Specify material properties
structuralProperties(modelStatic, "YoungsModulus", E, "PoissonsRatio", nu, "MassDensity", rho); 

% Specify same constraint at left end of the beam
structuralBC(modelStatic, 'Edge', 4, 'Constraint', 'fixed');

% Apply static vertical load at the right side of the beam
structuralBoundaryLoad(modelStatic, 'Vertex', 2, 'Force', [0; 0.01*Fcrit]);

% Solve static model
staticRes = solve(modelStatic);

% Compute initial deflection
deltaInit = interpolateDisplacement(staticRes, [L; 0.]).uy;

fprintf(' done.\nThe initial deflection of the beam''s free end is %6.4e m.\n', deltaInit);

%% Transient analysis
%  This is the actual transient analysis of the stick-slip behaviour of a
%  cantilever beam subjected to frictional contact with a moving conveyor
%  belt at the free end.

% Compute analytical solution
deltaBreak = (Fcrit * L^3) / (3 * E * I);
tBreak = deltaBreak/vbelt - deltaInit/vbelt;
fprintf('The analytical breakaway point is %6.4e s, when the deflection is %6.4e m.\n\n', tBreak, deltaBreak);

modelTransient.SolverOptions.ReportStatistics  = 'off';
modelTransient.SolverOptions.AbsoluteTolerance = 1.e-6; % Default: 1e-6
modelTransient.SolverOptions.RelativeTolerance = 1.e-3; % Default: 1e-3

figure
sgtitle("Cantilever beam 2D transient Finite Element analysis")
set(gcf, "WindowState", "maximized")
drawnow

fprintf('\nStarting time stepping convergence study...\n');

transitionPoints = [];
dtVec = [1e-3, 5e-4, 4e-4, 3e-4, 2e-4, 1e-4];

for dt = dtVec

    fprintf('\nPerforming transient analysis using time step size dt = %6.4e s...\n\n', dt);

    % Create transient analysis model for 2-D plane-stress problem
    modelTransient = createpde("structural", "transient-planestress");
    
    % Use same geometry and mesh
    geometryFromEdges(modelTransient, g);
    modelTransient.Mesh = mesh;
    
    % Specify material properties
    structuralProperties(modelTransient, "YoungsModulus", E, "PoissonsRatio", nu, "MassDensity", rho);
    
    % Specify same constraint at left end of the beam
    structuralBC(modelTransient, 'Edge', 4, 'Constraint', 'fixed');

    % Define initial conditions
%     structuralIC(modelTransient, "Displacement", zeros(2,1), "Velocity", zeros(2,1)); % Zero ICs
    structuralIC(modelTransient, staticRes); % Initial deflection from static results
    
    % Set up time stepping
    tStart  = 0.0;
    tEnd    = tfinal;
    tStep   = dt;
    
    nSteps = size((tStart : tStep : tEnd), 2);
    
    % Initialize container to store results
    history.t  = zeros(1, nSteps);
    history.uy = zeros(1, nSteps);
    history.Fy = zeros(1, nSteps);

    history.uy(1) = interpolateDisplacement(staticRes, [L; 0.]).uy;
    
    % Initialize friction state to sticktion
    slip = false;
    
    for step = 1 : nSteps-1
    
        if mod(step, 20) == 0; fprintf('t = %6.4f\n', history.t(step)); end
    
        tListTmp = [step-1 step] .* tStep;
    
        % Solve structural model for current time fraction
        rhs = solve(modelTransient, tListTmp);
    
        % Append partial solution to global solution vectors
        history.t (step+1) = rhs.SolutionTimes(end);
        history.uy(step+1) = interpolateDisplacement(rhs, [L; 0.]).uy(end);
        
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
                fprintf("t = %6.4f s: transition stick -> slip\n", history.t(step+1));
                transitionPoints = [transitionPoints, Freact];
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

%            if (vrel >= 0) && (abs(Freact) < Fcrit)
            if (abs(vrel) < fmodel.eps) && (abs(Freact) < Fcrit)
    
                % Slip -> Stick
                fprintf("t = %6.4f s: transition slip -> stick\n", history.t(step+1));
                transitionPoints = [transitionPoints, vrel];
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
    
                structuralBoundaryLoad(modelTransient, "Vertex", 2, "Force", [0. -Ff]);
    
            end
    
        end
    
    end

    fprintf('\nThe analysis for time step size dt = %6.4e s has been completed.\n', dt);
    
    % Plot the time-series of the deflection of the beam's free end
    ax = subplot(2,2,1);
    plot(history.t, history.uy)
    hold on
    xlabel("t [s]")
    ylabel("u_y(x=L) [m]")
    title({'y-displacement at the beam tip'})
    grid on
    drawnow
    
    % Plot the phase portrait at the beam's free end
    ax = subplot(2,2,2);
    plot(history.uy(2:end), diff(history.uy)./diff(history.t)-vbelt)
    hold on
    xlabel("u_y(x=L) [m]")
    ylabel("v_{rel}(x=L) [m/s]")
    title({'Phase portrait at the beam tip'})
    grid on
    drawnow

end

%% Add legends to plots
ax = subplot(2,2,1);
legend("dt = " + string(dtVec) + " s", "Location", "southeast");

ax = subplot(2,2,2);
legend("dt = " + string(dtVec) + " s");
xlim([0.004 0.0096])
ylim([-0.12 0.04])

% Plot convergence behavior for stick-slip transition points
ax = subplot(2,2,3);
FreactTrans = transitionPoints(transitionPoints < -45);
FreactErr = abs(abs(FreactTrans) - Fcrit) ./ abs(Fcrit) * 100;
plot(ceil(tfinal ./ dtVec), FreactErr);
xlabel("n_{steps} [-]")
ylabel("e_{rel} [%]")
title({'Stick-slip transition points error'})
grid on
legend("1st trans.", "2nd trans.", "3rd trans.", "4th trans.")

% Plot convergence behavior for slip-stick transition points
ax = subplot(2,2,4);
vrelTrans = transitionPoints(transitionPoints > -1);
vrelErr = abs(vrelTrans) * 100;
plot(ceil(tfinal ./ dtVec), vrelErr);
xlabel("n_{steps} [-]")
ylabel("e_{rel}=|v_s^t| [%]")
title({'Slip-stick transition points error'})
grid on
legend("1st trans.", "2nd trans.", "3rd trans.", "4th trans.")