% Mesh convergence study for a 2D static cantilever beam Finite Element model

%% Initialization
clear
close all

% Create static analysis model for 2-D plane-stress problem
modelStatic = createpde("structural", "static-planestress")

%% Create Geometry

L = 2.0;    % Beam length
h = 0.01;   % Beam height
I = h^3/12; % Moment of inertia

x = [0.; L; L; 0.];
y = [-h/2; -h/2; h/2; h/2];

geometryFromEdges(modelStatic, decsg([3, 4, x' y']'))

% Plot geometry
figure
sgtitle("Cantilever beam 2D static Finite Element analysis")
set(gcf, "WindowState", "maximized")
ax = subplot(2,2,1)
pdegplot(modelStatic, "VertexLabels", "on", "EdgeLabels", "on");
title('Geometry with labels')
xlabel("x")
ylabel("y")
xlim([0. 2.2])
ylim([-0.06 0.06])
grid on
ax.DataAspectRatio = [1 0.1 0.1];

%% Preprocess

% Define material properties
E   = 2e11;     % Young's modulus
nu  = 0.3;      % Poisson's ratio
rho = 1000.;    % Mass density

structuralProperties(modelStatic, "YoungsModulus", E, "PoissonsRatio", nu, "MassDensity", rho)

% Specify left edge of the beam as fixed boundary
structuralBC(modelStatic, "Edge", 4, "Constraint", "fixed");

% Apply static vertical load corresponding to maximum static friction force
% at the right side of the beam
Fcrit = 50;
structuralBoundaryLoad(modelStatic, 'Vertex', 2, 'Force', [0; Fcrit]);

% Compute analytical solution
deltaRef = (Fcrit * L^3) / (3 * E * I);
% Note: Euler-Bernoulli beam theory
% https://en.wikipedia.org/wiki/Euler%E2%80%93Bernoulli_beam_theory (accessed 2022-06-30)
% https://en.wikipedia.org/wiki/List_of_second_moments_of_area (accessed 2022-06-30)
fprintf('Expected max. deflection: delta = %6.4f m.\n', deltaRef);

%% Solve

% Perform mesh convergence study using increasingly fine meshes until error
% limit is reached
errorLimit = 7e-3;

hmax = h .* 200; % Initial value for target maximum element edge size
i = 1;

ax = subplot(2,2,2)
while true

    % Create mesh
    generateMesh(modelStatic, "Hmax", hmax, "GeometricOrder", "quadratic");

    % Plot mesh
    pdeplot(modelStatic)
    xlabel("x")
    ylabel("y")
    title("Finite Element mesh")
    xlim([0. 2.2])
    ylim([-0.06 0.06])
    grid on
    ax.DataAspectRatio = [1 0.1 0.1];
    drawnow

    % Solve the structural analysis problem
    staticRes = solve(modelStatic);

    % Interpolate the y-displacement at the geometric center of the beam
    delta = interpolateDisplacement(staticRes, [L; 0]).uy;

    numNodes(i) = size(modelStatic.Mesh.Nodes, 2);
    err(i) = abs(deltaRef-delta) / abs(deltaRef) * 100;
    fprintf('h = %6.4f, nodes = %8d, e = %6.4f %%.\n', hmax, numNodes(i), err(i));

    if err(i) < errorLimit; break; end

    i = i + 1;
    hmax = hmax .* 0.7; % Reduce element size by 30%

end

%% Postprocess

% Plot mesh convergence behavior
subplot(2,2,3)
loglog(numNodes, err, "Marker", "o", "MarkerSize", 5.0)
hold on
loglog(numNodes, (1./numNodes.^1).*2e2, "LineStyle", "--")
loglog(numNodes, (1./numNodes.^2).*1e3, "LineStyle", "--")
xlabel("n_{nodes} [-]")
ylabel("e_{rel} [%]")
title("Mesh convergence")
xlim([10 inf])
ylim([0.2*errorLimit 10.0])
legend('Rel. error', '1/n^1', '1/n^2')
grid off

% Plot deformed cantilever beam
ax = subplot(2,2,4)
pdeplot(modelStatic, "XYData", staticRes.Displacement.uy, ...
    "Deformation", staticRes.Displacement, 'DeformationScaleFactor', 1e0, ...
    "ColorMap", "Jet");
xlabel("x")
ylabel("y")
title('y-displacements')
xlim([0. 2.2])
ylim([-0.06 0.06])
grid on
ax.DataAspectRatio = [1 0.1 0.1];
