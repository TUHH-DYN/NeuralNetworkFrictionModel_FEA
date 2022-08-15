% Generation of training and test data for neural network regression
%
% This script generates training and test data sets based on analytical
% friction models. The data is used to fit a regression neural network
% model and assess its performance.
%
% Requires: function getFrictionForce.m
%
% Author: Kerstin Vater, MSc
% Machine Learning Dynamics Group (M-14)
% Hamburg University of Technology
% Am Schwarzenberg-Campus 1
% 21073 Hamburg, Germany
% E-mail: kerstin.vater@tuhh.de  
% URL: https://www.tuhh.de/dyn

%------------- BEGIN CODE --------------

% clear
% close all

% Load custom colormaps
load("colors\viridis.mat")
load("colors\Set1.mat")
set(0, "DefaultAxesColorOrder", Set1)

% Set up figure
figure
sgtitle("Friction data sampling and partitioning")
set(gcf, "WindowState", "maximized")

% Set target normal force [N]
Fn = 100.0;

% Set (approximate) number of data samples
nsamp = 5000;

% Set holdout value for data partitioning training/test
holdout = 0.3;

% Generate sample points for the relative sliding velocity according to
% Gaussian distribution with zero mean in order to capture the friction
% force discontinuity at the origin
rng("default") % Restore the random number generator seed for reproducibility
vs = normrnd(0, 0.1, [1, nsamp])';
vs = vs(abs(vs)>1e-4); % Exclude velocity values less than eps

% Generate sample points for the normal force according to uniform
% distribution around the target value of the noraml force
rng("default")
Fn = unifrnd(0.98*Fn, 1.02*Fn, size(vs));

% Set common friction parameters
fmodel.mus = 0.5; % Static coefficient of friction
fmodel.muk = 0.2; % Kinetic coefficient of friction

%% Exponential friction model

fmodel.type = 'exponential';
fmodel.a = 5.0; % Control parameter for the negative gradient of the friction curve

% Compute friction force
Ff = getFrictionForce(Fn, vs, fmodel);

% Add some (normal distributed, zero mean) noise
sigma = 0.02; % Standard deviation
rng("default")
Ff = Ff + (sigma .* Ff) .* randn(size(Ff));

% Assemble data tables
dataFull = table(Fn, vs, Ff, 'VariableNames', ["Fn", "vs", "Ff"]);

dataFull.Properties.Description = 'Friction data samples';
dataFull.Properties.VariableUnits = {'N', 'm/s', 'N'};
dataFull.Properties.VariableDescriptions{'Fn'} = 'Normal force';
dataFull.Properties.VariableDescriptions{'vs'} = 'Relative sliding velocity';
dataFull.Properties.VariableDescriptions{'Ff'} = 'Kinetic friction force';

% Split the data into training and test set
rng("default") % For reproducibility of the data partition
c = cvpartition(height(dataFull), "HoldOut", holdout);
idxTrain  = training(c); % Training set indices
idxTest   = test    (c); % Test set indices
dataTrain = dataFull(idxTrain, :);
dataTest  = dataFull(idxTest , :);

fprintf("Friction data generation completed.\n");

% Plot the data samples
figure
scatter3(dataTrain.Fn, dataTrain.vs, dataTrain.Ff, "Marker", ".")
hold on
scatter3(dataTest.Fn, dataTest.vs, dataTest.Ff, "Marker", ".")
xlabel("F_n [N]");
ylabel("v_{s} [m/s]");
zlabel("F_f [N]");
title("Exponential Friction Model")
view([110 15])

% Add analytical surface plot as reference
Fn = linspace(min(Fn), max(Fn), 10);
[Fn_, vrel_] = meshgrid(Fn, linspace(min(vs), max(vs), 1000)');
Ff = getFrictionForce(Fn_, vrel_, fmodel);
surf(Fn_, vrel_, Ff, "EdgeColor", "none", "FaceAlpha", 0.5);
legend("Training data (" + num2str((1-holdout)*100) + " %)", ...
           "Test data (" + num2str(   holdout *100) + " %)", ...
           "Analytical friction model", "Location", "northeast")

% Plot histogram
figure
ax = subplot(1,2,1);

% -- ecpdf
cdfplot(dataFull.vs); hold on; 
cdfplot(dataTrain.vs); 
cdfplot(dataTest.vs);
ylabel('Empirical cumulative probability [-]');
title('Inputs')

% -- envelope of histogram
% [counts, edges] = histcounts(dataFull.vs, 50, 'Normalization', "probability"); locs = movmean(edges, 2, 'Endpoints', 'discard');
% plot(locs, counts, 'LineWidth', 2); hold on;
% 
% [counts, edges] = histcounts(dataTrain.vs, 50, 'Normalization', "probability"); locs = movmean(edges, 2, 'Endpoints', 'discard');
% plot(locs, counts, 'LineWidth', 2);
% 
% [counts, edges] = histcounts(dataTest.vs, 50, 'Normalization', "probability"); locs = movmean(edges, 2, 'Endpoints', 'discard');
% plot(locs, counts, 'LineWidth', 2);

% -- histogram
% histogram(dataFull.vs, 'BinWidth',0.1, 'Normalization', "probability", 'faceAlpha', 0.2, 'EdgeCOlor', 'k'); %'DisplayStyle','stairs');
% hold on
% histogram(dataTrain.vs, 'BinWidth', 0.1, 'Normalization',"probability", 'faceAlpha', 0.2); %, 'DisplayStyle','stairs');
% histogram(dataTest.vs, 'BinWidth', 0.1, 'Normalization', "probability", 'faceAlpha', 0.2); %, 'DisplayStyle','stairs');

xlabel("v_{s} [m/s]");
legend("Entire data set (100 %)", ...
    "Training data (" + num2str((1-holdout)*100) + " %)", ...
        "Test data (" + num2str(   holdout *100) + " %)", ...
        "Location", "northwest")

ax = subplot(1,2,2);

% -- ecpdf
cdfplot(dataFull.Ff); hold on; 
cdfplot(dataTrain.Ff); 
cdfplot(dataTest.Ff);
ylabel('Empirical cumulative probability [-]');
title('Outputs')
xlabel("F_{f} [N]");
legend("Entire data set (100 %)", ...
    "Training data (" + num2str((1-holdout)*100) + " %)", ...
        "Test data (" + num2str(   holdout *100) + " %)", ...
        "Location", "northwest")

% Preview and summarize table
head   (dataFull)
summary(dataFull)

% Write data to file
writetable(dataFull,  fmodel.type + "_friction_model_samples_full.csv")
writetable(dataTrain, fmodel.type + "_friction_model_samples_training.csv")
writetable(dataTest,  fmodel.type + "_friction_model_samples_test.csv")

%------------- END OF CODE --------------