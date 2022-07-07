% Analytical friction model parametrization and visualization
%
% Note: requires function getFrictionForce.m

clear
close all

% Load custom colormaps
load("colors\viridis.mat")
load("colors\Set1.mat")
set(0, "DefaultAxesColorOrder", Set1)

Fn = 1.0;                       % Normal force [N]
vs = linspace(-1.0, 1.0, 1e4)'; % Sample sliding velocities [m/s]

% Set common friction parameters
fmodel.mus = 0.5; % Static coefficient of friction
fmodel.muk = 0.2; % Kinetic coefficient of friction

% Compute maximum static friction force [N]
Fcrit = Fn * fmodel.mus;

% Set up figure
figure
sgtitle("Analytical friction models")
set(gcf, "WindowState", "maximized")

%% Exponential friction model

fmodel.type = 'exponential';
fmodel.a = 5.0; % Control parameter for the negative gradient of the friction curve

% Sample friction coefficient and friction force
[Ff, mu] = getFrictionForce(Fn, vs, fmodel);

% Plot Stribeck curve
ax = subplot(2,2,1);
plot(vs, mu)
hold on
plot(vs, +fmodel.muk*ones(size(vs)), '--', 'Color', Set1(end,:))
plot(vs, -fmodel.muk*ones(size(vs)), '--', 'Color', Set1(end,:))
title('Exponential friction model')
xlabel('v_{s} [m/s]')
ylabel('\mu [-]')
xlim([0.0, max(vs)])
ylim([0.0, 1.2*fmodel.mus])
grid on

yyaxis right
plot(vs, +fmodel.muk*ones(size(vs)), '--')
plot(vs, -fmodel.muk*ones(size(vs)), '--')
plot(vs, +fmodel.mus*ones(size(vs)), '--')
plot(vs, -fmodel.mus*ones(size(vs)), '--')
ylim([0.0, 1.2*fmodel.mus])
yticks([-fmodel.mus, -fmodel.muk, fmodel.muk, fmodel.mus])
yticklabels({'-\mu_s', '-\mu_k', '\mu_k', '\mu_s'})

% Plot friction force
ax = subplot(2,2,3);
plot(vs, Ff)
hold on
plot(vs, +Fn*fmodel.muk*ones(size(vs)), '--', 'Color', Set1(end,:))
plot(vs, -Fn*fmodel.muk*ones(size(vs)), '--', 'Color', Set1(end,:))
xlabel('v_{s} [m/s]')
ylabel('F_{f} [N]')
xlim([min(vs), max(vs)])
ylim([-1.2*Fcrit, 1.2*Fcrit])
grid on

yyaxis right
plot(vs, +Fn*fmodel.muk*ones(size(vs)), '--')
plot(vs, -Fn*fmodel.muk*ones(size(vs)), '--')
plot(vs, +Fn*fmodel.mus*ones(size(vs)), '--')
plot(vs, -Fn*fmodel.mus*ones(size(vs)), '--')
ylim([-1.2*Fcrit, 1.2*Fcrit])
yticks([-Fn*fmodel.mus, -Fn*fmodel.muk, Fn*fmodel.muk, Fn*fmodel.mus])
yticklabels({'-F_n \mu_s', '-F_n \mu_k', 'F_n \mu_k', 'F_n \mu_s'})


%% Polynomial friction model

fmodel.type = 'polynomial';
fmodel.v0 = 0.3; % Reference velocity if exponential decay is used

% Sample friction coefficient and friction force
[Ff, mu] = getFrictionForce(Fn, vs, fmodel);

% Plot Stribeck curve
ax = subplot(2,2,2);
plot(vs, mu)
hold on
plot([+fmodel.v0 +fmodel.v0], [-Fn +Fn], '--', 'Color', Set1(end,:))
plot([-fmodel.v0 -fmodel.v0], [-Fn +Fn], '--', 'Color', Set1(end,:))
xticks([-0.9, -0.6, -fmodel.v0, 0, fmodel.v0, 0.6, 0.9])
xticklabels({-0.9, -0.6, '-v_0 = -0.3', 0, 'v_0 = 0.3', 0.6, 0.9})
title('Polynomial friction model')
xlabel('v_{s} [m/s]')
ylabel('\mu [-]')
xlim([0.0, max(vs)])
ylim([0.0, 1.2*fmodel.mus])
grid on

yyaxis right
plot(vs, +fmodel.muk*ones(size(vs)), '--')
plot(vs, -fmodel.muk*ones(size(vs)), '--')
plot(vs, +fmodel.mus*ones(size(vs)), '--')
plot(vs, -fmodel.mus*ones(size(vs)), '--')
ylim([0.0, 1.2*fmodel.mus])
yticks([-fmodel.mus, -fmodel.muk, fmodel.muk, fmodel.mus])
yticklabels({'-\mu_s', '-\mu_k', '\mu_k', '\mu_s'})

% Plot friction force
ax = subplot(2,2,4);
plot(vs, Ff)
hold on
plot([+fmodel.v0 +fmodel.v0], [-Fn +Fn], '--', 'Color', Set1(end,:))
plot([-fmodel.v0 -fmodel.v0], [-Fn +Fn], '--', 'Color', Set1(end,:))
xticks([-0.9, -0.6, -fmodel.v0, 0, fmodel.v0, 0.6, 0.9])
xticklabels({-0.9, -0.6, '-v_0 = -0.3', 0, 'v_0 = 0.3', 0.6, 0.9})
xlabel('v_{s} [m/s]')
ylabel('F_{f} [N]')
xlim([min(vs), max(vs)])
ylim([-1.2*Fcrit, 1.2*Fcrit])
grid on

yyaxis right
plot(vs, +Fn*fmodel.muk*ones(size(vs)), '--')
plot(vs, -Fn*fmodel.muk*ones(size(vs)), '--')
plot(vs, +Fn*fmodel.mus*ones(size(vs)), '--')
plot(vs, -Fn*fmodel.mus*ones(size(vs)), '--')
ylim([-1.2*Fcrit, 1.2*Fcrit])
yticks([-Fn*fmodel.mus, -Fn*fmodel.muk, Fn*fmodel.muk, Fn*fmodel.mus])
yticklabels({'-F_n \mu_s', '-F_n \mu_k', 'F_n \mu_k', 'F_n \mu_s'})