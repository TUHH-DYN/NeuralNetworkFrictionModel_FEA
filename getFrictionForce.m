function [Ff, mu] = getFrictionForce(Fn, vs, fmodel)
% GETFRICTIONFORCE returns the kinetic friction force as a function of the
% normal force and the relative sliding velocity, depending on the
% governing analytical friction model.
%
%   Input parameters:
%   Fn      normal force                    [N]
%   vs      relative sliding velocity       [m/s]
%   fmodel  friction model with parameters
%           including
%   .type   'exponential', 'polynomial', 'nn', ...
%   .mus    static coefficient of friction  [-]
%   .muk    kinetic coefficient of friction [-]
%   
%   Output parameters:
%   Ff      kinetic friction force          [N]
%   mu      Coefficient of friction         [-]

    % Extract common friction parameters
    muS = fmodel.mus;   % Static coefficient of friction
    muK = fmodel.muk;   % Kinetic coefficient of friction
    
    switch fmodel.type

        case 'coulomb'

            % Simple Coulomb law
            Ff = sign(vs) .* muK .* ones(size(vs)) .* Fn;

        case 'exponential'

            % Exponential decaying
            alpha = fmodel.a;

            mu = muK + (muS - muK) .* exp(-alpha.*abs(vs));
            Ff = sign(vs) .* mu .* Fn;

        case 'polynomial'

            % Polynomial, kind of Stribeck
            v0 = fmodel.v0; % Reference velocity

            mu = muS - 3/2 * (muS - muK) * ((vs./v0) - 1/3 * (vs./v0).^3) .*sign(vs);
            Ff = sign(vs) .* mu .* Fn;

        case 'nn'

            % Predict the friction force using a trained linear Regression model
            Ff = predict(fmodel.trainedModel.RegressionNeuralNetwork, ...
                table(Fn.*ones(size(vs)), vs, 'VariableNames', {'Fn', 'vs'}));

        otherwise

            Ff = 0.0;
    end

    % Compute coefficient of friction
    mu = abs(Ff ./ Fn);

end
