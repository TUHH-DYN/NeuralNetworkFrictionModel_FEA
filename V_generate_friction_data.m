clear
close all

% Define analytical friction model
fmodel.type = 'exponential';    % Choose between 'exponential' and 'polynomial'
fmodel.a    = 5.0;              % Control parameter for the negative gradient of the friction curve (default: 5.0)
fmodel.v0   = 0.3;              % Reference velocity if exponential decay is used   (default: 0.3) 
fmodel.musd = 2.5;              % Ratio static to dynamic fric coeff, mu_st/mu_d    (default: 2.5)
fmodel.mud  = 0.2;              % Dynamic coeff of friction, mu_d                   (default: 0.2)
fmodel.muv  = 0.001;            % Linear strengthening parameter                    (default: 0.001)
fmodel.eps  = 1e-3;             % Threshold value for sliding velocity              (default: 1e-3)

vrel = (-1.0 : fmodel.eps : 1.0);   % Relative sliding velocity
Fn   = 100;                         % Normal force (default: 100)
[Fn_, vrel_] = meshgrid(Fn, vrel');

% Unroll matrices to column vectors
vrel_ = reshape(vrel_, [], 1);
Fn_   = reshape(Fn_  , [], 1);

% Add noise to the input data
sigma = 0.02; % Standard deviation set to 2%
vrel_ = vrel_ + (sigma .* vrel_) .* randn(size(vrel_));
Fn_ = Fn_ + (sigma .* Fn_) .* randn(size(Fn_));

% Compute friction force
Ff = getFrictionForce(Fn_, vrel_, fmodel);
% Note: Noise should not be added to the output data in order to represent
% the analytical friction model properly

% Assemble data table
data = table(Fn_, vrel_, Ff, VariableNames=["Fn", "vrel", "Ff"]);

data.Properties.Description = 'Friction data samples';
data.Properties.VariableUnits = {'N', 'm/s', 'N'};
data.Properties.VariableDescriptions{'Fn'  } = 'Normal force';
data.Properties.VariableDescriptions{'vrel'} = 'Relative sliding velocity';
data.Properties.VariableDescriptions{'Ff'  } = 'Kinetic friction force';

% Plot the data
figure
scatter3(data, "Fn", "vrel", "Ff", "Marker", ".")
xlabel("F_n [N]");
ylabel("v_{rel} [m/s]");
zlabel("F_f [N]");
title("Friction data samples")
subtitle("from " + fmodel.type + " friction model")
% view([-90 0])
hold on

% Add surface plot as reference
Fn = linspace(min(Fn_), max(Fn_), 10);
vrel = (min(vrel_) : fmodel.eps : max(vrel_));
[Fn_, vrel_] = meshgrid(Fn, vrel');
Ff = getFrictionForce(Fn_, vrel_, fmodel);
surf(Fn_, vrel_, Ff, "EdgeColor", "none", "FaceAlpha", 0.5);

% Preview and summarize table
head(data)
summary(data)

% Write data to file
writetable(data, fmodel.type + "FrictionModelSamples.csv")
