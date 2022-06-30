function [F, mu] = getFrictionForce(Fn, vrel, fmodel)

    % Extract common parameters
    musd    = fmodel.musd;  % Ratio static to dynamic fric coeff, mu_st/mu_d
    mud     = fmodel.mud;   % Dynamic coeff of friction
    
    switch fmodel.type

        case 'coulomb'

            % Simple Coulomb law
            F = -sign(vrel) .* mud .* ones(size(vrel)) .* Fn;

        case 'leine1998eq11'

            % Friction force according to [Leine et al., 1998] eq. (11)
            delta = 3.0; % [s\m]

            F = -(Fn .* mud .* musd .* sign(vrel)) ./ (1. + delta .* abs(vrel));

        case 'smoothed'

            % Smoothed friction force according to [Leine et al., 1998] eq. (12)
            eps = 1e6; % Smoothing control parameter
            delta = 3.0; % [s\m]

            F = -(Fn .* mud .* musd .* (2./pi) .* atan(eps .* vrel)) ./ (1. + delta .* abs(vrel));

        case 'exponential'

            % Exponential decaying (plus linear strengthening) as used by [Won and Chung, 2018]
%             muv = fmodel.muv; % Linear strengthening parameter
            alpha = fmodel.a;

            mu = mud + mud .* (musd-1) .* exp(-alpha.*abs(vrel)); % + muv .* abs(vrel)./v0;
            F = -sign(vrel) .* mu .* Fn;

        case 'polynomial'

            % Polynomial, kind of Stribeck
            v0 = fmodel.v0; % Reference velocity

            mu = mud*musd - 3/2*mud*(musd - 1)*((vrel./v0) - 1/3*(vrel./v0).^3).*sign(vrel);
            F = -sign(vrel) .* mu .* Fn;

        case 'neuralNet'

            % Predict the friction force using a trained linear Regression model
            F = predict(fmodel.trainedModel.RegressionNeuralNetwork, ...
                table(Fn.*ones(size(vrel)), vrel, 'VariableNames', {'Fn','vrel'}));

        otherwise

            F = 0.0;
    end

    % Compute coefficient of friction
    mu = F ./ Fn;

end
