function [Ff, mu] = getFrictionForce(Fn, vrel, fmodel)
% GETFRICTIONFORCE returns the kinetic friction force as a function of the
% normal force and the relative sliding velocity, depending on the
% governing analytical friction model.
%
%   Input parameters:
%   Fn      normal force                    [N]
%   vrel    relative sliding velocity       [m/s]
%   fmodel  friction model with parameters
%           including
%   .type   'exponential', 'polynomial', 'nn', ...
%   .mud    dynamic coeff of friction       [-]
%   .musd   ratio static to dynamic fric coeff
%
%   Output parameters:
%   Ff      kinetic friction force          [N]
%   mu      Coefficient of friction         [-]

    % Extract common parameters
    musd = fmodel.musd;  % Ratio static to dynamic fric coeff, mu_st/mu_d
    mud  = fmodel.mud;   % Dynamic coeff of friction
    
    switch fmodel.type

        case 'coulomb'

            % Simple Coulomb law
            Ff = -sign(vrel) .* mud .* ones(size(vrel)) .* Fn;

        case 'exponential'

            % Exponential decaying
            alpha = fmodel.a;

            mu = mud + mud .* (musd-1) .* exp(-alpha.*abs(vrel));
            Ff = -sign(vrel) .* mu .* Fn;

        case 'polynomial'

            % Polynomial, kind of Stribeck
            v0 = fmodel.v0; % Reference velocity

            mu = mud*musd - 3/2*mud*(musd - 1)*((vrel./v0) - 1/3*(vrel./v0).^3).*sign(vrel);
            Ff = -sign(vrel) .* mu .* Fn;

        case 'nn'

            % Predict the friction force using a trained linear Regression model
            Ff = predict(fmodel.trainedModel.RegressionNeuralNetwork, ...
                table(Fn.*ones(size(vrel)), vrel, 'VariableNames', {'Fn','vrel'}));

        otherwise

            Ff = 0.0;
    end

    % Compute coefficient of friction
    mu = Ff ./ Fn;

end
