% Create a regression neural network friction model with low error by using
% the OptimizeHyperparameters argument. This argument causes fitrnet to
% minimize cross-validation loss over some problem hyperparameters by using
% Bayesian optimization.
clear
close all

%% Prepare data

% Read friction data into a table
data = readtable("exponentialFrictionModelSamples.csv");

% Partition the data into training and test sets
rng("default") % For reproducibility of the data partition
c = cvpartition(height(data), "HoldOut", 0.30);

trainingIdx = training(c); % Training set indices
dataTrain = data(trainingIdx,:);
testIdx = test(c); % Test set indices
dataTest = data(testIdx,:);

%% Perform hyperparameter optimization
rng("default") % For reproducibility
regNet = fitrnet(...
    dataTrain, "Ff", ...
    "OptimizeHyperparameters", "auto");
%     "HyperparameterOptimizationOptions", ...
%     struct("AcquisitionFunctionName", "expected-improvement-plus"));

% Find the mean squared error of the resulting model on the test data set
testMSE = loss(regNet, dataTest, "Ff")

%% Determine network size
numNeurons = sum(regNet.LayerSizes)
numLearnables = 2 * regNet.LayerSizes(1); % Count weights
for l = 2 : numel(regNet.LayerSizes)
    numLearnables = numLearnables + ...
        regNet.LayerSizes(l) * regNet.LayerSizes(l-1);
end
numLearnables = numLearnables + 1 * regNet.LayerSizes(end);
numLearnables = numLearnables + numNeurons + 1 % Add bias

%% Perform cross-validation
partitionedModel = crossval(regNet, 'KFold', 5);

% Compute regression loss for observations not used for training
kfLoss = kfoldLoss(partitionedModel, "Mode", "individual", "LossFun", "mse");

% Compute max value, mean value and standard deviation of loss values
kflMax  = max(kfLoss)
kflMean = mean(kfLoss)
kflStd  = std(kfLoss)

% Compute validation RMSE
validationRMSE = sqrt(kflMean);

%% Train a regression neural network model using the full data set
% This code specifies all the model options and trains the model.
fullRegNet = fitrnet(...
    data, "Ff", ...
    'LayerSizes',       regNet.ModelParameters.LayerSizes , ...
    'Activations',      regNet.ModelParameters.Activations, ...
    'Lambda',           regNet.ModelParameters.Lambda, ...
    'IterationLimit',   regNet.ModelParameters.IterationLimit, ...
    'Standardize',      regNet.ModelParameters.StandardizeData, ...
    'StoreHistory',     true);

%% Plot predicted friction force
Fn = 100.0;  % Normal force                         (default: 100)
mud  = 0.2;  % Dynamic coeff of friction            (default: 0.2)
eps  = 1e-3; % Threshold value for sliding velocity (default: 1e-3)

% Sample friction coefficient and friction force
vrel = linspace(-1.0, 1.0, 1e5)';
Xsample = table(Fn.*ones(size(vrel)), vrel, 'VariableNames', {'Fn','vrel'});
Ff = predict(fullRegNet, Xsample);

% Plot the friction force as function of relative sliding velocity
figure
plot(vrel, Ff)
hold on
plot(vrel, +Fn*mud*ones(size(vrel)), '--', 'Color', [.5 .5 .5])
plot(vrel, -Fn*mud*ones(size(vrel)), '--', 'Color', [.5 .5 .5])
plot([+eps +eps], [-Fn +Fn], '--', 'Color', [.5 .5 .5])
plot([-eps -eps], [-Fn +Fn], '--', 'Color', [.5 .5 .5])
title('Friction force curve', 'F_n = 100N, \mu_k = 0.2')
xlabel('v_{rel} [m/s]')
ylabel('F_f [N]')
xlim([-10*eps, 10*eps])
ylim([-Fn, Fn])
grid on
hold on

%% Plot reference curve
fmodel.type = 'exponential';    % Choose between 'exponential' and 'polynomial'
fmodel.a    = 5.0;              % Control parameter for the negative gradient of the friction curve (default: 5.0)
fmodel.v0   = 0.3;              % Reference velocity if exponential decay is used   (default: 0.3) 
fmodel.musd = 2.5;              % Ratio static to dynamic fric coeff, mu_st/mu_d    (default: 2.5)
fmodel.mud  = 0.2;              % Dynamic coeff of friction, mu_d                   (default: 0.2)
fmodel.muv  = 0.001;            % Linear strengthening parameter                    (default: 0.001)
fmodel.eps  = 1e-3;             % Threshold value for sliding velocity              (default: 1e-3)

FfRef = getFrictionForce(Fn, vrel, fmodel);
plot(vrel, FfRef)
legend('Neural network prediction', 'Analytical friction model')