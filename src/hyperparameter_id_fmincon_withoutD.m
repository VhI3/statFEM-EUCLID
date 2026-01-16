function BVP = hyperparameter_id_fmincon_withoutD(BVP, nrep, material_model, sigmad_d_ub)

    %% Assign from BVP
    sensorCoordinates = {BVP.obs.sensorCoordinates, BVP.obs.sensorCoordinates};
    Y_exp = {BVP.obs.Y_exp_x(:, nrep), BVP.obs.Y_exp_y(:, nrep)};
    Pmatrix = {BVP.obs.Px, BVP.obs.Py};

    switch material_model
        case 'LE'
            mean_u = {BVP.lsNIPCE.LE.mean_ux, BVP.lsNIPCE.LE.mean_uy};
            cov_u = {nearestSPD(BVP.lsNIPCE.LE.cov_ux), nearestSPD(BVP.lsNIPCE.LE.cov_uy)};
        case 'HE'
            mean_u = {BVP.euclied.HE.mean_ux, BVP.euclied.HE.mean_uy};
            cov_u = {nearestSPD(BVP.euclied.HE.cov_ux), nearestSPD(BVP.euclied.HE.cov_uy)};
        otherwise
            error('Unsupported material model. Use "LE" or "HE".');
    end

    DOFs = BVP.preProc.msh.DOFs;
    nsen = BVP.obs.nSen;
    % C_e = BVP.obs.C_e;
    GSensorDOFs = BVP.obs.GSensorDOFs;
    fixedPrescribedSensorDOFs = BVP.obs.fixedPrescribedSensorDOFs;
    identified_hyperparameters = BVP.discovered.identified_hyperparameters;

    %% Initialize
    rho_est = zeros(DOFs, 1);
    sigd_est = zeros(DOFs, 1);
    ld_est = zeros(DOFs, 1);
    C_d0 = cell(DOFs, 1);

    %% Optimization parameters
    % nTrials = 1000; % You can adjust this
    % rng(0); % For reproducibility

    % Define log-space bounds
    % lb = log([0.5, BVP.preProc.statFEM.sig_e/10, 1]);
    % ub = log([1.0, 10^-(sigmad_d_ub), 1.0e02]);

    % lb = log([0.5, BVP.preProc.statFEM.sig_e / 10, 1]);
    % ub = log([2.0, 1, 1.0e02]);

    % LHS sampling in log-space
    % paramRanges = [lb; ub]'; % for easier indexing
    % lhsPoints = lhsdesign(nTrials, 3);
    % init_guesses = zeros(nTrials, 3);
    % 
    % for j = 1:3
    %     init_guesses(:, j) = paramRanges(j, 1) + lhsPoints(:, j) * (paramRanges(j, 2) - paramRanges(j, 1));
    % end

    % p = sobolset(3);
    % p = scramble(p, 'MatousekAffineOwen');
    % init_guesses = net(p, nTrials);  % Normalized [0, 1]
    % % Scale to log-bounds
    % for j = 1:3
    %     init_guesses(:, j) = lb(j) + init_guesses(:, j) .* (ub(j) - lb(j));
    % end

    % Optimization options for fmincon
    % options = optimoptions('fmincon', ...
    %     'Algorithm', 'interior-point', ...
    %     'SpecifyObjectiveGradient', true, ...
    %     'Display', 'none', ...
    %     'MaxIterations', 1000, ...
    %     'OptimalityTolerance', 1e-6, ...
    %     'StepTolerance', 1e-8, ...
    %     'FunctionTolerance', 1e-6);

    %% Loop over DOFs
    % nrep = 1;

    for ii = 1:DOFs

        % MultiStart setup
        % ms = MultiStart('UseParallel', true, 'Display', 'off');
        % startPoints = CustomStartPointSet(init_guesses);
        % 
        % % Define the negative log-likelihood function
        % negLogLikelihoodFunc = @(w) negativeLogLikelihood(w, ...
        %     cov_u{ii}, Pmatrix{ii}, C_e{ii}, Y_exp{ii}, ...
        %     mean_u{ii}, sensorCoordinates{ii}, nrep, nsen);
        % 
        % % Create optimization problem for this DOF
        % problem = createOptimProblem('fmincon', ...
        %     'objective', negLogLikelihoodFunc, ...
        %     'x0', init_guesses(1, :), ...
        %     'lb', lb, ...
        %     'ub', ub, ...
        %     'options', options);
        % 
        % % Run MultiStart optimization
        % [best_w, best_loss] = run(ms, problem, startPoints);

        % Store estimated parameters (convert from log-space)
        rho_est(ii) = 1;
        sigd_est(ii) = 0;
        ld_est(ii) = 100;

        % Covariance matrix
        % C_d0{ii} = sqexp(sensorCoordinates{ii}, sensorCoordinates{ii}, ...
        %     sigd_est(ii), ld_est(ii)); % Still in log-space here

        C_d0{ii} = zeros(nsen); % Still in log-space here

        % Output log
        fprintf('DOF %d Estimated: rho = %.4f, sig_d = %.4f, l_d = %.4f\n', ...
            ii, rho_est(ii), sigd_est(ii), ld_est(ii));

        % Save to BVP
        identified_hyperparameters{ii} = [identified_hyperparameters{ii}; ...
                                              rho_est(ii), sigd_est(ii), ld_est(ii)];
    end

    %% Assemble full covariance matrix
    C_d_full = zeros(GSensorDOFs);
    C_d_full(1:2:end, 1:2:end) = C_d0{1};
    C_d_full(2:2:end, 2:2:end) = C_d0{2};
    C_d_full(fixedPrescribedSensorDOFs, fixedPrescribedSensorDOFs) = 0;

    C_d = {C_d_full(1:2:end, 1:2:end), C_d_full(2:2:end, 2:2:end)};

    %% Store final results
    switch material_model
        case 'LE'
            BVP.hyp.LE.rho = rho_est;
            BVP.hyp.LE.sigd = sigd_est;
            BVP.hyp.LE.ld = ld_est;
            BVP.hyp.LE.C_d = C_d;
            BVP.hyp.LE.C_d_full = C_d_full;
        case 'HE'
            BVP.hyp.HE.rho = rho_est;
            BVP.hyp.HE.sigd = sigd_est;
            BVP.hyp.HE.ld = ld_est;
            BVP.hyp.HE.C_d = C_d;
            BVP.hyp.HE.C_d_full = C_d_full;
    end

    BVP.discovered.identified_hyperparameters = identified_hyperparameters;

end

% function [f, g, H] = negativeLogLikelihood(w, C_u_active_pc, P_active, C_e, y_obs, mu_u_active_pc, senCoor, nrep, nsen)
function [f, g] = negativeLogLikelihood(w, C_u_active_pc, P_active, C_e, y_obs, mu_u_active_pc, senCoor, nrep, nsen)
    % function [f] = negativeLogLikelihood(w, C_u_active_pc, P_active, C_e, y_obs, mu_u_active_pc, senCoor, nrep, nsen)
    % NEGATIVELOGLIKELIHOOD - Computes the negative log-likelihood, its gradient, and Hessian.
    %
    % Inputs:
    %   w               - Log-space hyperparameter vector: [log_rho, log_sigma, log_ell]
    %   C_u_active_pc   - Prior covariance from PCE
    %   P_active        - Projection matrix
    %   C_e             - Measurement noise covariance
    %   y_obs           - Observations (nsen x nrep)
    %   mu_u_active_pc  - Mean displacement from PCE
    %   senCoor         - Sensor coordinates
    %   nrep            - Number of experiment repetitions
    %   nsen            - Number of sensors
    %
    % Outputs:
    %   f - Scalar value of the negative log-likelihood (cost function)
    %   g - (3x1) Gradient vector of the negative log-likelihood
    %   H - (3x3) Hessian matrix of the negative log-likelihood

    %% Step 1: Negative log-posterior
    f = logposterior1D(w, C_u_active_pc, P_active, C_e, y_obs, mu_u_active_pc, senCoor, nrep, nsen);

    %% Step 2: Gradient
    g = logpost_deriv1D(w, C_u_active_pc, P_active, C_e, y_obs, mu_u_active_pc, senCoor, nrep, nsen);

    %% Step 3: Hessian
    % H = hessian_logposterior1D(w, C_u_active_pc, P_active, C_e, y_obs, mu_u_active_pc, senCoor, nrep, nsen);
end

function log_posterior = logposterior1D(w, C_u_active_pc, P_active, C_e, y_obs, mu_u_active_pc, senCoor, nrep, nsen)

    log_rho = w(1); % Now in log-space
    log_sigma = w(2);
    log_ell = w(3);

    rho = exp(log_rho); % Transform to real space

    A = P_active * C_u_active_pc * P_active';
    mu_proj = P_active * mu_u_active_pc;

    K = sqexp(senCoor, senCoor, log_sigma, log_ell);
    Sigma = nearestSPD(rho ^ 2 * A + K + C_e);

    L = chol(Sigma, 'lower');
    invSigma = chol_solve(L, eye(nsen));

    residual_term = 0;

    for i = 1:nrep
        r = y_obs(:, i) - rho * mu_proj;
        residual_term = residual_term + r' * invSigma * r;
    end

    log_posterior = 0.5 * ( ...
        nrep * nsen * log(2 * pi) + ...
        2 * nrep * sum(log(diag(L))) + ...
        residual_term);
end

function deriv = logpost_deriv1D(w, C_u_active_pc, P_active, C_e, y_obs, mu_u_active_pc, senCoor, nrep, nsen)

    log_rho = w(1); % log-space
    log_sigma = w(2);
    log_ell = w(3);

    rho = exp(log_rho); % transform to real space

    A = P_active * C_u_active_pc * P_active';
    mu_proj = P_active * mu_u_active_pc;

    K = sqexp(senCoor, senCoor, log_sigma, log_ell);
    Sigma = nearestSPD(rho ^ 2 * A + K + C_e);

    L = chol(Sigma, 'lower');
    invSigma = chol_solve(L, eye(nsen));

    residuals = y_obs - rho * mu_proj;

    % First derivatives of K
    K_deriv = sqexp_deriv(senCoor, senCoor, log_sigma, log_ell);
    dK_dlogsigma = K_deriv(:, :, 1);
    dK_dlogell = K_deriv(:, :, 2);

    % ---- Derivative w.r.t. log(ρ) ----
    dSigma_drho = 2 * rho * A;
    tr_term = trace(invSigma * dSigma_drho);
    quad_term = 0;
    mean_term = 0;

    for i = 1:nrep
        r = residuals(:, i);
        v = invSigma * r;
        quad_term = quad_term + v' * dSigma_drho * v;
        mean_term = mean_term + mu_proj' * v;
    end

    dlogrho = 0.5 * (nrep * tr_term - quad_term) - mean_term;
    deriv(1) = dlogrho * rho; % Chain rule: ∂ϑ/∂log(ρ) = ∂ϑ/∂ρ × ∂ρ/∂log(ρ)

    % ---- Derivative w.r.t. log(σ) ----
    tr_sigma = trace(invSigma * dK_dlogsigma);
    quad_sigma = sum(arrayfun(@(i) ...
        (invSigma * residuals(:, i))' * dK_dlogsigma * (invSigma * residuals(:, i)), 1:nrep));
    deriv(2) = 0.5 * (nrep * tr_sigma - quad_sigma);

    % ---- Derivative w.r.t. log(ℓ) ----
    tr_ell = trace(invSigma * dK_dlogell);
    quad_ell = sum(arrayfun(@(i) ...
        (invSigma * residuals(:, i))' * dK_dlogell * (invSigma * residuals(:, i)), 1:nrep));
    deriv(3) = 0.5 * (nrep * tr_ell - quad_ell);
end

% options = optimoptions('fminunc', ...
%     'Algorithm', 'trust-region', ...
%     'SpecifyObjectiveGradient', true, ...
%     'FiniteDifferenceType', 'central', ...
%     'CheckGradients', true, ...
%     'Display', 'iter-detailed', ...
%     'MaxIterations', 1000, ...
%     'OptimalityTolerance', 1e-6, ...
%     'StepTolerance', 1e-10, ...
%     'FunctionTolerance', 1e-10);

% function Ahat = nearestSPD(A)
% % nearestSPD - Fast nearest Symmetric Positive Definite matrix to A
%
%     if nargin ~= 1
%         error('Exactly one input argument is required.');
%     end
%
%     % Check square
%     [r, c] = size(A);
%     if r ~= c
%         error('Input matrix must be square.');
%     elseif r == 1 && A <= 0
%         Ahat = eps;
%         return;
%     end
%
%     % Step 1: Symmetrize
%     B = (A + A') / 2;
%
%     % Step 2: Use eig instead of svd for speed (B is symmetric)
%     [V, D] = eig(B);
%     D = max(D, eps(max(diag(D))));  % Ensure non-negative eigenvalues
%     H = V * D * V';
%
%     % Step 3: Average to get nearest SPD
%     Ahat = (B + H) / 2;
%
%     % Step 4: Ensure symmetry
%     Ahat = (Ahat + Ahat') / 2;
%
%     % Step 5: Ensure positive-definiteness via modified Cholesky
%     [R, p] = chol(Ahat);
%     k = 0;
%     I = eye(r, 'like', A);  % GPU-compatible identity
%     while p ~= 0
%         k = k + 1;
%         % Add jitter using scaled identity matrix
%         jitter = eps(norm(Ahat)) + k * 1e-12;
%         Ahat = Ahat + jitter * I;
%         [R, p] = chol(Ahat);
%     end
% end
