% Analytical friction model parametrization and visualization

clear
close all

% Load custom colormaps
load("colors\viridis.mat")
load("colors\Set1.mat")
set(0, "DefaultAxesColorOrder", Set1)

% Set up figure
figure
sgtitle("Analytical friction models")
set(gcf, "WindowState", "maximized")

% Set common friction parameters
fmodel.musd = 2.5; % Ratio static to dynamic fric coeff, mu_st/mu_d
fmodel.mud  = 0.2; % Dynamic coeff of friction, mu_d

%% Exponential friction model

fmodel.type = 'exponential';
fmodel.a = 5.0; % Control parameter for the negative gradient of the friction curve

% Sample friction coefficient and friction force
vs = linspace(-1.0, 1.0, 1e4)'; % Sample sliding velocities
Fn = 100.0; % Normal force [N]
[Ff, mu] = getFrictionForce(Fn, vs, fmodel);

% Plot Stribeck curve
ax = subplot(2,2,1);
plot(vs, mu)
hold on
plot(vs, +fmodel.mud*ones(size(vs)), '--', 'Color', Set1(end,:))
plot(vs, -fmodel.mud*ones(size(vs)), '--', 'Color', Set1(end,:))
title('Exponential friction model')
xlabel('v_{s} [m/s]')
ylabel('\mu [-]')
xlim([-1.0, 1.0])
ylim([-0.60, 0.60])
grid on

yyaxis right
plot(vs, +fmodel.mud*ones(size(vs)), '--')
plot(vs, -fmodel.mud*ones(size(vs)), '--')
plot(vs, +fmodel.mud*fmodel.musd*ones(size(vs)), '--')
plot(vs, -fmodel.mud*fmodel.musd*ones(size(vs)), '--')
ylim([-0.60, 0.60])
yticks([-fmodel.mud*fmodel.musd, -fmodel.mud, fmodel.mud, fmodel.mud*fmodel.musd])
yticklabels({'-\mu_s', '-\mu_k', '\mu_k', '\mu_s'})

% Plot friction force
ax = subplot(2,2,3);
plot(vs, Ff)
hold on
plot(vs, +Fn*fmodel.mud*ones(size(vs)), '--', 'Color', Set1(end,:))
plot(vs, -Fn*fmodel.mud*ones(size(vs)), '--', 'Color', Set1(end,:))
xlabel('v_{s} [m/s]')
ylabel('F_{f} [N]')
xlim([-1.0, 1.0])
ylim([-60, 60])
grid on

yyaxis right
plot(vs, +Fn*fmodel.mud*ones(size(vs)), '--')
plot(vs, -Fn*fmodel.mud*ones(size(vs)), '--')
plot(vs, +Fn*fmodel.mud*fmodel.musd*ones(size(vs)), '--')
plot(vs, -Fn*fmodel.mud*fmodel.musd*ones(size(vs)), '--')
ylim([-60, 60])
yticks([-Fn*fmodel.mud*fmodel.musd, -Fn*fmodel.mud, Fn*fmodel.mud, Fn*fmodel.mud*fmodel.musd])
yticklabels({'-F_n \mu_s', '-F_n \mu_k', 'F_n \mu_k', 'F_n \mu_s'})


%% Polynomial friction model

fmodel.type = 'polynomial';
fmodel.v0 = 0.3; % Reference velocity if exponential decay is used

% Sample friction coefficient and friction force
vs = linspace(-1.0, 1.0, 1e4)';
Fn = 100.0; % Normal force [N]
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
xlim([-1.0, 1.0])
ylim([-0.60, 0.60])
grid on

yyaxis right
plot(vs, +fmodel.mud*ones(size(vs)), '--')
plot(vs, -fmodel.mud*ones(size(vs)), '--')
plot(vs, +fmodel.mud*fmodel.musd*ones(size(vs)), '--')
plot(vs, -fmodel.mud*fmodel.musd*ones(size(vs)), '--')
ylim([-0.60, 0.60])
yticks([-fmodel.mud*fmodel.musd, -fmodel.mud, fmodel.mud, fmodel.mud*fmodel.musd])
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
xlim([-1.0, 1.0])
ylim([-60, 60])
grid on

yyaxis right
plot(vs, +Fn*fmodel.mud*ones(size(vs)), '--')
plot(vs, -Fn*fmodel.mud*ones(size(vs)), '--')
plot(vs, +Fn*fmodel.mud*fmodel.musd*ones(size(vs)), '--')
plot(vs, -Fn*fmodel.mud*fmodel.musd*ones(size(vs)), '--')
ylim([-60, 60])
yticks([-Fn*fmodel.mud*fmodel.musd, -Fn*fmodel.mud, Fn*fmodel.mud, Fn*fmodel.mud*fmodel.musd])
yticklabels({'-F_n \mu_s', '-F_n \mu_k', 'F_n \mu_k', 'F_n \mu_s'})