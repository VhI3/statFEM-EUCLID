function BVP = model_discovery_Euclied(BVP, MR_order, criterion)
%% Assign from BVP
nElm = BVP.obs.nuSenElm;
nElNode = BVP.preProc.msh.nElNode;
dim = BVP.preProc.msh.DIM;
elementNodes = BVP.obs.sensor_connectivity;
nodeCoordinates = BVP.obs.sensorCoordinates;
GDOFs = BVP.obs.GSensorDOFs;
Q_pts = BVP.preProc.Q_pts;
Q_wts = BVP.preProc.Q_wts;
dN_dxi = BVP.preProc.dN_dxi;
dN_deta = BVP.preProc.dN_deta;

mu_ux_y = BVP.obs.Y_exp_x(:, 1);
mu_uy_y = BVP.obs.Y_exp_y(:, 1);

% mu_ux_y = BVP.obs.Y_exp_x;
% mu_uy_y = BVP.obs.Y_exp_y;

% mu_ux_y = BVP.proc.NH.ux;
% mu_uy_y = BVP.proc.NH.uy;
%
activeDOFs = BVP.obs.activeSensorDOFs;
force = BVP.obs.sensor_force;
prop0 = BVP.preProc.material.prop0;
P_full = BVP.obs.P_full;
True_exp = BVP.obs.True_exp;
W_total0 = BVP.proc.NH.W_total0;
W_true = BVP.discovered.W_true;

J1_samples = BVP.proc.NH.J1_samples;
J2_samples = BVP.proc.NH.J2_samples;
J3_samples = BVP.proc.NH.J3_samples;
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

%%
% A : (n x p) matrix
% b : (n x 1) vector
[n, p] = size(A);

% Precompute for speed
H_base = 2 * (A' * A);
f_base = -2 * (A' * b);

% Lambda values (logarithmic scale)
lambda_list = logspace(-5, 5, 100);

% Initialize
errors = zeros(size(lambda_list)); % RMSE or residual norms
AIC = zeros(size(lambda_list));
BIC = zeros(size(lambda_list));
B_all = zeros(p, length(lambda_list));
nonzeros_all = zeros(size(lambda_list));

% quadprog options
opts = optimoptions('quadprog', 'Display', 'off');

for i = 1:length(lambda_list)
    lambda = lambda_list(i);
    
    H = H_base;
    f = f_base + lambda * ones(p, 1); % Add L1 penalty approximation
    
    % Solve with non-negativity constraint
    lb = zeros(p, 1);
    [B, ~, exitflag] = quadprog(H, f, [], [], [], [], lb, [], [], opts);
    
    % Handle solver failure
    if exitflag ~= 1
        warning('quadprog failed at lambda = %.4f', lambda);
        B = nan(p, 1);
    end
    
    % Save solution
    B_all(:, i) = B;
    pred = A * B;
    residuals = b - pred;
    mse = mean(residuals .^ 2);
    k = sum(B > 1e-6); % Number of active features
    
    errors(i) = sqrt(mse);
    nonzeros_all(i) = k;
    
    AIC(i) = 2 * k + n * log(mse + eps); % AIC
    BIC(i) = k * log(n) + n * log(mse + eps); % BIC
end

% === Choose selection criterion ===
% criterion = 'BIC';  % Options: 'RMSE' | 'AIC' | 'BIC' | 'sparse'

switch criterion
    case 'RMSE'
        [~, best_idx] = min(errors);
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

% Final outputs
lambda_best = lambda_list(best_idx);
B = B_all(:, best_idx);

fprintf('Selected by %s:\n  lambda = %.4e\n  RMSE = %.4f\n  Nonzeros = %d\n', ...
    criterion, lambda_best, errors(best_idx), nonzeros_all(best_idx));

B_Raw = B;
% Optional: threshold tiny coefficients
B(abs(B) < 1e-6) = 0;
%%

% B(abs(B) < 1e1) = 0;
disp('LASSO coefficients:');
disp(B);
B(end) = 2 * B(end);
prop = zeros(10, 1);
prop(1:size(B, 1) - 1, :) = B(1:end - 1);
prop(end) = B(end);

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

% Wrap in equation
w_latex_full = latex_expression;

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
    prop(10) * (J3 - 1).^2;


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


fprintf('Relative RMSE of Euclied = %.5f\n', relative_rmse);

epsilon_kappa = norm(prop - prop0) / norm(prop0);
% hyperelastic_matparameters = [hyperelastic_matparameters, prop];
%%
[u_euclied, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, W_total, ~] = solver_HE(BVP.preProc, BVP.preProc.BC.force, prop);
epsilon_u = norm(P_full * u_euclied - True_exp) / norm(True_exp);
epsilon_I = (norm(W_total - W_total0) / norm(W_total0));
%%
BVP.discovered.Euclied.prop = prop;
BVP.discovered.Euclied.u_euclied = u_euclied;
BVP.discovered.Euclied.epsilon_kappa = epsilon_kappa;
BVP.discovered.Euclied.epsilon_u = epsilon_u;
BVP.discovered.Euclied.w_latex_full = w_latex_full;
BVP.discovered.Euclied.B_Raw = B_Raw;
BVP.discovered.Euclied.epsilon_I = epsilon_I;
BVP.discovered.Euclied.W_discovered = W_discovered;
BVP.discovered.Euclied.relative_rmse = relative_rmse;
end


