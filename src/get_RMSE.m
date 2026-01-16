% function BVP = get_RMSE(BVP, iter, model_type)
% % Compute RMSE between statFEM predicted and observed displacements
% % model_type: 'LE' or 'HE'
% 
% %% === Check input ===
% if nargin < 3
%     error('get_RMSE requires 3 inputs: BVP, iter, model_type');
% end
% 
% if ~strcmp(model_type, 'LE') && ~strcmp(model_type, 'HE')
%     error('model_type must be either "LE" or "HE"');
% end
% 
% %% === Assign fields based on model type ===
% switch model_type
%     case 'LE'
%         mu_zx = BVP.statFEM.LE.mu_zx;
%         mu_zy = BVP.statFEM.LE.mu_zy;
%     case 'HE'
%         mu_zx = BVP.euclied.HE.mu_zx;
%         mu_zy = BVP.euclied.HE.mu_zy;
% end
% 
% Y_exp_x = BVP.obs.Y_exp_x(:, iter);
% Y_exp_y = BVP.obs.Y_exp_y(:, iter);
% nSen = BVP.obs.nSen;
% 
% %% === Compute RMSE ===
% rmse_x = norm(mu_zx - Y_exp_x)^2;
% rmse_y = norm(mu_zy - Y_exp_y)^2;
% 
% RMSE_model = sqrt((rmse_x + rmse_y) / nSen);
% 
% %% === Display ===
% fprintf('%s RMSE at iteration %d: %.6f\n', model_type, iter, RMSE_model);
% 
% % %% === Store in BVP ===
% BVP.discovered.RMSE(end+1) = RMSE_model;
% 
% end
% 
% 



function BVP = get_RMSE(BVP, iter, model_type)
%GET_RMSE Computes the RMSE between statFEM predictions and observed displacements.
%
% Inputs:
%   BVP        - Problem structure containing statFEM results and observations
%   iter       - Current iteration index for accessing experimental data
%   model_type - 'LE' (Linear Elastic) or 'HE' (Hyperelastic) model
%
% Output:
%   BVP        - Updated BVP structure with appended RMSE

    % --- Check input validity ---
    if nargin < 3
        error('get_RMSE requires 3 inputs: BVP, iter, model_type');
    end

    valid_models = {'LE', 'HE'};
    if ~ismember(model_type, valid_models)
        error('model_type must be either "LE" or "HE"');
    end

%% === Assign fields based on model type ===
switch model_type
    case 'LE'
        mu_zx = BVP.statFEM.LE.mu_zx;
        mu_zy = BVP.statFEM.LE.mu_zy;
    case 'HE'
        mu_zx = BVP.euclied.HE.mu_zx;
        mu_zy = BVP.euclied.HE.mu_zy;
end

    % --- Observed experimental data ---
    Y_exp_x = BVP.obs.Y_exp_x(:, iter);
    Y_exp_y = BVP.obs.Y_exp_y(:, iter);

    % --- Number of sensors ---
    nSen = BVP.obs.nSen;

    % --- Compute RMSE (root-mean-square error) ---
    rmse_x = norm(mu_zx - Y_exp_x)^2;
    rmse_y = norm(mu_zy - Y_exp_y)^2;
    RMSE_model = sqrt((rmse_x + rmse_y) / (2 * nSen));

    % --- Display result ---
    fprintf('[%s] RMSE at iteration %d: %.6e\n', model_type, iter, RMSE_model);

    % --- Append RMSE to BVP ---
    if ~isfield(BVP.discovered, 'RMSE') || isempty(BVP.discovered.RMSE)
        BVP.discovered.RMSE = RMSE_model;
    else
        BVP.discovered.RMSE(end+1) = RMSE_model;
    end
end
