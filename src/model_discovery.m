function BVP = model_discovery(BVP, MR_order, material_model, criterion)
%% Assign from BVP
nElm = BVP.preProc.msh.nElm;
nElNode = BVP.preProc.msh.nElNode;
dim = BVP.preProc.msh.DIM;
elementNodes = BVP.preProc.msh.elementNodes;
nodeCoordinates = BVP.preProc.msh.nodeCoordinates;
GDOFs = BVP.preProc.msh.GDOFs;
Q_pts = BVP.preProc.Q_pts;
Q_wts = BVP.preProc.Q_wts;
dN_dxi = BVP.preProc.dN_dxi;
dN_deta = BVP.preProc.dN_deta;

switch material_model
    case 'LE'
        mu_ux_y = BVP.statFEM.LE.mu_ux_y;
        mu_uy_y = BVP.statFEM.LE.mu_uy_y;
    case 'HE'
        mu_ux_y = BVP.euclied.HE.mu_ux_y;
        mu_uy_y = BVP.euclied.HE.mu_uy_y;
    case 'ref'
        mu_ux_y = BVP.proc.NH.ux;
        mu_uy_y = BVP.proc.NH.uy;
    otherwise
        error('Unsupported material model. Use "LE" or "HE".');
end

%
activeDOFs = BVP.preProc.BC.activeDOFs;
force = BVP.preProc.BC.force;
prop0 = BVP.preProc.material.prop0;
W_total0 = BVP.proc.NH.W_total0;

u_true = BVP.proc.NH.u;
hyperelastic_matparameters = BVP.discovered.hyperelastic_matparameters;
epsilon_u_iter = BVP.discovered.epsilon_u_iter;
epsilon_kappa_iter = BVP.discovered.epsilon_kappa_iter;
epsilon_I_iter = BVP.discovered.epsilon_I_iter;
w_latex_full_iter = BVP.discovered.w_latex_full_iter;
W_true = BVP.discovered.W_true;
J1_samples = BVP.proc.NH.J1_samples;
J2_samples = BVP.proc.NH.J2_samples;
J3_samples = BVP.proc.NH.J3_samples;
relative_rmse_iter = BVP.discovered.relative_rmse_iter;
%% Calculation
voigtMap = [0, 1; 2, 3];

gradNa = cell(nElNode, 1);

for e = 1:nElm % loop over elements
    indx = elementNodes(e, :);

    for i = 1:size(Q_wts, 1)
        xi = Q_pts(i, 1); eta = Q_pts(i, 2);
        naturalDerivatives = [dN_dxi(xi, eta); dN_deta(xi, eta)]';
        %
        [JacobianMatrix, ~, XYDerivatives] = Jacobian2D(nodeCoordinates(indx, :), naturalDerivatives);

        B = zeros(4, dim * nElNode);
        B(1, 1:2:end) = XYDerivatives(:, 1)';
        B(2, 1:2:end) = XYDerivatives(:, 2)';
        B(3, 2:2:end) = XYDerivatives(:, 1)';
        B(4, 2:2:end) = XYDerivatives(:, 2)';
        B_element{e, 1} = B;
        %
        gradNa{1, 1}(e, :) = XYDerivatives(1, :);
        gradNa{2, 1}(e, :) = XYDerivatives(2, :);
        gradNa{3, 1}(e, :) = XYDerivatives(3, :);
        qpWeights(e, 1) = Q_wts(i) * det(JacobianMatrix);
    end

end

u = cell(nElNode, 1);
u_nodes = [mu_ux_y mu_uy_y];
connectivity{1, 1} = elementNodes(:, 1);
connectivity{2, 1} = elementNodes(:, 2);
connectivity{3, 1} = elementNodes(:, 3);

for i = 1:nElNode
    u{i, 1} = u_nodes(connectivity{i, 1}, :);
end

F = zeros(nElm, 4);

for a = 1:nElNode

    for i = 1:dim

        for j = 1:dim
            F(:, voigtMap(i, j) + 1) = F(:, voigtMap(i, j) + 1) + u{a}(:, i) .* gradNa{a}(:, j);
        end

    end

end

F(:, 1) = F(:, 1) + 1.0;
F(:, 4) = F(:, 4) + 1.0;

% Compute strain and invariants
J = computeJacobian(F);
C = computeCauchyGreenStrain(F);
[I1, I2, I3] = computeStrainInvariants(C);
dI1dF = computeStrainInvariantDerivatives(F, 1);
dI2dF = computeStrainInvariantDerivatives(F, 2);
dI3dF = computeStrainInvariantDerivatives(F, 3);

if MR_order == 3
    featuresFunc_I1 = @(I) computeFeatures_MR3(I, I2, I3);
    featuresFunc_I2 = @(I) computeFeatures_MR3(I1, I, I3);
    featuresFunc_I3 = @(I) computeFeatures_MR3(I1, I2, I);
elseif MR_order == 2
    featuresFunc_I1 = @(I) computeFeatures_MR2(I, I2, I3);
    featuresFunc_I2 = @(I) computeFeatures_MR2(I1, I, I3);
    featuresFunc_I3 = @(I) computeFeatures_MR2(I1, I2, I);
elseif MR_order == 1
    featuresFunc_I1 = @(I) computeFeatures_MR1(I, I2, I3);
    featuresFunc_I2 = @(I) computeFeatures_MR1(I1, I, I3);
    featuresFunc_I3 = @(I) computeFeatures_MR1(I1, I2, I);
else
    error('Invalid MR_order. Choose 1, 2, or 3.');
end

% Compute derivatives with finite difference
[d_features_dI1, numFeatures] = differentiateFeaturesWithInvariantsFD(featuresFunc_I1, I1);
[d_features_dI2, ~] = differentiateFeaturesWithInvariantsFD(featuresFunc_I2, I2);
[d_features_dI3, ~] = differentiateFeaturesWithInvariantsFD(featuresFunc_I3, I3);

LHS = zeros(GDOFs, numFeatures);

for ele = 1:nElm
    d_features_dF_element = assembleFeatureDerivative(d_features_dI1, d_features_dI3, dI1dF, dI3dF, ele);
    % B_element0 = assembleB(gradNa, ele, dim, nElNode);
    % LHS_element = (B_element0' * d_features_dF_element') * qpWeights(ele);
    LHS_element = (B_element{ele, 1}' * d_features_dF_element') * qpWeights(ele);
    LHS = assembleGlobalMatrix(LHS, LHS_element, connectivity, ele, nElNode);
end

LHS_active = LHS(activeDOFs, :);
RHS_active = force(activeDOFs);
A = LHS_active;
b = RHS_active;
%% Matlab lasso regression
%% Non-negative Lasso regression (with sparsity criterion)
% Inputs:
%   - A : (n x p) matrix
%   - b : (n x 1) vector

[n, p] = size(A);

% Quadratic programming setup
optsQP = optimoptions('quadprog','Display','off');
H_base = 2 * (A' * A);
f0 = -2 * (A' * b);
lb = zeros(p,1);   % enforce B ≥ 0

% Lambda values (logarithmic scale)
lambda_list = logspace(-2, 4.5, 500);
num_k = zeros(1,1000);
% Initialize result storage
errors = zeros(size(lambda_list));
AIC = zeros(size(lambda_list));
BIC = zeros(size(lambda_list));
B_all = zeros(p, length(lambda_list));
nonzeros_all = zeros(size(lambda_list));

r = 3; % volumetric penalty multiplier
Aineq = [3*ones(1, p-1), -1];  % sum(B(1:end-1)) - B(end)
% Aineq = [3*ones(1, 2), -1];  % sum(B(1:end-1)) - B(end)
bineq = 0;

% Solve non-negative lasso for each lambda
for i = 1:length(lambda_list)
    lambda = lambda_list(i);

    % Solve via quadratic programming (QP)
    f = f0 + lambda * ones(p,1);
    % B = quadprog(H_base, f, [], [], [], [], lb, [], [], optsQP);
    B = quadprog(H_base, f, [], [], Aineq, bineq, lb, [], [], optsQP);

    % Save coefficients
    B_all(:, i) = B;

    % Evaluate performance metrics
    residuals = b - A * B;
    mse = mean(residuals.^2);
    k = nnz(B > 1e-13); % active coefficients (> threshold)
    % k = nnz(B > 1e-6); % active coefficients (> threshold)
    % if i ==1 || i ==2
    %     comlexity_k(i) = k;
    % end
    %
    % if i>2
    %     if comlexity_k(i-1) < k && comlexity_k(i-2) < k
    %         comlexity_k(i) = k-1;
    %     else
    %         comlexity_k(i) = k;
    %     end
    % end
    num_k(i) = k;



    errors(i) = sqrt(mse);
    nonzeros_all(i) = k;
    AIC(i) = 2*k + n*log(mse + eps);
    BIC(i) = k*log(n) + n*log(mse + eps);
end
% close all
% semilogx(lambda_list,errors)
% hold on
% semilogx(lambda_list,rmoutliers(comlexity_k,"mean"))

% set(gca, 'XScale','log');
% plot(lambda_list, num_k, '--')

% stop

%% Choose model explicitly based on sparsity
% criterion = 'sparse'; % Options: 'RMSE', 'AIC', 'BIC', 'sparse'

switch criterion
    case 'RMSE'
        [min_rmse, ~] = min(errors);
        [max_rmse, ~] = max(errors);
        r_l = 1e-5;
        m_th = min_rmse + r_l*(max_rmse - min_rmse);

        best_idx = find(errors >= m_th, 1);

    case 'AIC'
        [~, best_idx] = min(AIC);
    case 'BIC'
        [~, best_idx] = min(BIC);
    case 'sparse'
        min_nz = min(nonzeros_all);
        sparse_candidates = find(nonzeros_all == min_nz);
        [~, best_in_group] = min(errors(sparse_candidates));
        best_idx = sparse_candidates(best_in_group);
end
%% Final selected solution
lambda_best = lambda_list(best_idx);
B = B_all(:, best_idx);
max(errors)
% stop
% Threshold tiny coefficients for clean output
B(abs(B) < 1e-15) = 0;

%% Display Results
fprintf('Selected by %s:\n', criterion);
fprintf('  lambda   = %.4e\n', lambda_best);
fprintf('  RMSE     = %.4f\n', errors(best_idx));
fprintf('  Nonzeros = %d\n', nonzeros_all(best_idx));

%% ---------- Pareto-style plot with BOTH selected lambdas ----------
% Inputs assumed from your script:
%   lambda_list, errors, B_all, nonzeros_all
%   best_idx (from your switch/criterion code)
%   lambda_best = lambda_list(best_idx)

close all

% Red curve (complexity) — L1 norm. If you prefer kappa: use nonzeros_all.
complexity_L1 = sum(abs(B_all), 1);

% Choose an MSE threshold (like the figure). Adjust as you wish.
MSE_th = prctile(errors, 20);
MSE_th = 70;
% Pareto-style pick: first lambda (from left to right) whose MSE <= threshold
idx_feas = find(errors <= MSE_th, 1, 'first');
if ~isempty(idx_feas)
    idx_star = idx_feas;
else
    [~, idx_star] = min(errors);     % fallback
end
lambda_star = lambda_list(idx_star);

% For convenience
mse_star   = errors(idx_star);
l1_star    = complexity_L1(idx_star);
mse_best   = errors(best_idx);
l1_best    = complexity_L1(best_idx);

figure('Color','w'); hold on; box on;
set(gca, 'XScale','log');           % if you want decreasing λ to the left: set(gca,'XDir','reverse')

% --- Left axis: MSE (blue)
% yyaxis left
% set(gca,'YTick',[])    % remove tick values
pMSE = plot(lambda_list, errors, '-', 'LineWidth', 1, ...
    'MarkerIndices', round(linspace(1,numel(lambda_list),12)));

% pMSE = plot(lambda_list, num_k, '-o', 'LineWidth', 1.8, ...
%     'MarkerIndices', round(linspace(1,numel(lambda_list),12)));

ylabel('$\text{RMSE}^{(\lambda)}$', 'Interpreter','latex')

pYline = plot(lambda_list,MSE_th*ones(1,500), 'k--');
xPos = lambda_list(100);             % x-position for the label (end of line)
yPos = MSE_th;                       % y-position for the label (on the line)
text(xPos, yPos, '$\tau$', ...
    'Interpreter', 'latex', ...
    'VerticalAlignment', 'bottom', ... % position relative to the line
    'HorizontalAlignment', 'right');   
% yline(MSE_th, 'k--', 'MSE_{th}', 'Interpreter','latex');
hold on

% Markers for both chosen lambdas on MSE
plot(lambda_best, mse_best, 'ro', 'MarkerSize', 2, 'LineWidth', 1.6);



% --- Right axis: complexity (red)
% yyaxis right
% % % set(gca,'YTick',[])    % remove tick values
% pL1 = plot(lambda_list, complexity_L1, '-', 'LineWidth', 0.5, ...
%     'MarkerIndices', round(linspace(1,numel(lambda_list),12)));
% ylabel('||\kappa||_1')

% Markers for both chosen lambdas on complexity
% plot(lambda_best, l1_best, 'ro', 'MarkerSize', 8, 'LineWidth', 1.6);
% xlim padded
% ylim padded

% --- Vertical lines and labels
% xline(lambda_star, 'k--', '\lambda^*_{Pareto}', ...
%     'Interpreter','tex','LabelVerticalAlignment','bottom');

hold on
% xline(lambda_best, 'r--', '\lambda^*', ...
%     'Interpreter','tex','LabelVerticalAlignment','top');


% Draw vertical line at lambda_best
pXline = plot([lambda_best lambda_best], ylim, 'r--'); 
hold on

% Add label near the top of the line
xPos = lambda_best;
yPos = ylim;        % returns [ymin ymax]
yPos = yPos(2);     % pick top of the y-axis

text(xPos, yPos, '$\lambda^*$', ...
    'Interpreter', 'latex', ...
    'VerticalAlignment', 'top', ... % just below top edge
    'HorizontalAlignment', 'left', 'Color','red');  % centered on the line



ylim padded
xlim padded

% --- Cosmetics
% xlabel('\lambda')
% title('Pareto analysis with criterion-based and Pareto-selected \lambda^*', 'Interpreter','tex')
% legend(pMSE, {'$\text{RMSE}^{(\lambda)}$'}, 'Location','northwest', 'Interpreter','latex')
% legend([pMSE pL1], {'MSE','||\kappa||_1'}, 'Location','best')


% convPDFV4(BVP.plotSettings, sprintf('pareto_fullpath_tau'), 'none'); close;


%%
disp('LASSO coefficients:');
disp(B);
B(end) = 2 * B(end);
prop = zeros(10, 1);
prop(1:size(B, 1) - 1, :) = B(1:end - 1);
prop(end) = B(end);

% stop
A_mn = prop(:) / 10 ^ 6;
A_mn(end) = A_mn(end) / 2;
%
%
% 1. Define symbolic variables
syms J1 J2 J3 real

% 2. Define symbolic features explicitly (same order as in your code)
Q{1} = (J1 - 3);
Q{2} = (J2 - 3);
Q{3} = (J1 - 3) ^ 2;
Q{4} = (J1 - 3) * (J2 - 3);
Q{5} = (J2 - 3) ^ 2;
Q{6} = (J1 - 3) ^ 3;
Q{7} = (J1 - 3) ^ 2 * (J2 - 3);
Q{8} = (J1 - 3) * (J2 - 3) ^ 2;
Q{9} = (J2 - 3) ^ 3;
Q{10} = (J3 - 1) ^ 2;
%
latex_terms = strings(1, 10); % initialize array of strings

for i = 1:10

    if abs(A_mn(i)) > 1e-12
        coeff_str = sprintf('%.3f', A_mn(i));

        if i == 1 || i == 2
            term_str = coeff_str + "\,(" + latex(Q{i}) + ")";
        else
            term_str = coeff_str + "\," + latex(Q{i});
        end

        latex_terms(i) = term_str;
    end

end

% Concatenate terms with ' + '
latex_expression = strjoin(latex_terms(latex_terms ~= ""), ' + ');


% Define W_discovered function handle
W_discovered = @(J1, J2, J3) ...
    prop(1) * (J1 - 3) + ...
    prop(2) * (J2 - 3) + ...
    prop(3) * (J1 - 3).^2 + ...
    prop(4) * (J1 - 3).*(J2 - 3) + ...
    prop(5) * (J2 - 3).^2 + ...
    prop(6) * (J1 - 3).^3 + ...
    prop(7) * (J1 - 3).^2 .* (J2 - 3) + ...
    prop(8) * (J1 - 3) .* (J2 - 3).^2 + ...
    prop(9) * (J2 - 3).^3 + ...
    0.5 * prop(10) * (J3 - 1).^2;

% Number of samples
nJ1 = length(J1_samples);
nJ2 = length(J2_samples);
nJ3 = length(J3_samples);

% Initialize total error
total_error = 0;

% Loop through all combinations
for i = 1:nJ1
    for j = 1:nJ2
        for k = 1:nJ3
            J1 = J1_samples(i);
            J2 = J2_samples(j);
            J3 = J3_samples(k);

            % Evaluate W_true and W_pred
            W_true_val = W_true(J1, J2, J3);  % Define this function
            W_pred_val = W_discovered(J1, J2, J3);  % Define this function

            % Avoid division by zero
            if abs(W_true_val) > 1e-10
                err = (W_true_val - W_pred_val)^2 / W_true_val^2;
                total_error = total_error + err;
            end
        end
    end
end

% Compute average relative error
relative_rmse = total_error / (nJ1 * nJ2 * nJ3);


fprintf('Relative RMSE of statFEM-Euclied = %.5f\n', relative_rmse);



% Wrap in equation
w_latex_full = latex_expression;

epsilon_kappa = norm(prop - prop0) / norm(prop0);
epsilon_kappa_iter = [epsilon_kappa_iter, epsilon_kappa];
hyperelastic_matparameters = [hyperelastic_matparameters, prop];

%%
[u_euclied, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, W_total, ~] = solver_HE(BVP.preProc, BVP.preProc.BC.force, prop);
% epsilon_u = norm(P_full * u_euclied - True_exp) / norm(True_exp);
epsilon_u = norm(u_euclied - u_true) / norm(u_true);
epsilon_u_iter = [epsilon_u_iter, epsilon_u];
epsilon_I = (norm(W_total - W_total0) / norm(W_total0));
epsilon_I_iter = [epsilon_I_iter, epsilon_I];
w_latex_full_iter = [w_latex_full_iter, w_latex_full];
relative_rmse_iter = [relative_rmse_iter, relative_rmse];
%% Assign back to BVP

BVP.discovered.material.prop = prop;
BVP.discovered.hyperelastic_matparameters = hyperelastic_matparameters;
BVP.discovered.epsilon_kappa_iter = epsilon_kappa_iter;
BVP.discovered.epsilon_I_iter = epsilon_I_iter;
BVP.discovered.epsilon_u_iter = epsilon_u_iter;
BVP.discovered.w_latex_full = w_latex_full;
BVP.discovered.w_latex_full_iter = w_latex_full_iter;
BVP.discovered.W_total = W_total;
BVP.discovered.W_discovered = W_discovered;
BVP.discovered.relative_rmse_iter = relative_rmse_iter;
end