function BVP = statFEM_update_HE(BVP)
    %% Assign from BVP
    % Degrees of freedom
    GDOFs = BVP.preProc.msh.GDOFs;

    % Observation-related data
    elementNodes = BVP.preProc.msh.elementNodes;
    Px = BVP.obs.Px;
    Py = BVP.obs.Py;


    mu_ux = BVP.euclied.HE.mean_ux;
    mu_uy = BVP.euclied.HE.mean_uy;
    C_ux = BVP.euclied.HE.cov_ux;
    C_uy = BVP.euclied.HE.cov_uy;

    Y_exp_x = BVP.obs.Y_exp_x;
    Y_exp_y = BVP.obs.Y_exp_y;

    % Hyperparameters and covariance matrices
    rho_x = BVP.hyp.HE.rho(1); % rho for X dimension
    rho_y = BVP.hyp.HE.rho(2); % rho for Y dimension
    C_ex = BVP.obs.C_e{1};
    C_ey = BVP.obs.C_e{2};
    C_dx = BVP.hyp.HE.C_d{1}; % Covariance for X dimension
    C_dy = BVP.hyp.HE.C_d{2}; % Covariance for Y dimension

    nrep = BVP.obs.nrep;

    %% Initialize Posterior Mean and Covariance
    mu_u_y = zeros(GDOFs, 1);
    C_u_y = zeros(GDOFs, GDOFs);

    L_C_ux = chol(nearestSPD(C_ux), 'lower');
    inv_C_ux = chol_solve(L_C_ux, eye(size(C_ux)));

    L_C_uy = chol(nearestSPD(C_uy), 'lower');
    inv_C_uy = chol_solve(L_C_uy, eye(size(C_uy)));

    %% Compute Inverse of (C_dx + C_ex)
    L_C_dx_C_ex = chol(nearestSPD(C_dx + C_ex), 'lower');
    inv_C_dx_C_ex = chol_solve(L_C_dx_C_ex, eye(size(C_dx + C_ex))); %#ok<ELARLOG> % Efficient inversion: inv(C_d + C_e)

    %% Compute Inverse of (C_dy + C_ey) using Cholesky Decomposition
    L_C_dy_C_ey = chol(nearestSPD(C_dy + C_ey), 'lower');
    inv_C_dy_C_ey = chol_solve(L_C_dy_C_ey, eye(size(C_dy + C_ey))); %#ok<ELARLOG> % Efficient inversion: inv(C_d + C_e)

    %% Compute the Posterior Covariance (C_ux_y)
    % Original: C_ux_y = inv(rho_x^2 * nrep * Px' * inv(C_dx + C_ex) * Px + inv(C_ux));
    M_x = rho_x ^ 2 * nrep * (Px' * inv_C_dx_C_ex * Px) + inv_C_ux;
    L_M_x = chol(M_x, 'lower');
    inv_M_x = chol_solve(L_M_x, eye(size(M_x)));
    C_ux_y0 = inv_M_x;

    %% Compute the Posterior Covariance (C_uy_y)
    % Original: C_uy_y = inv(rho_y^2 * nrep * Py' * inv(C_dy + C_ey) * Py + inv(C_uy));
    M_y = rho_y ^ 2 * nrep * (Py' * inv_C_dy_C_ey * Py) + inv_C_uy;
    L_M_y = chol(M_y, 'lower');
    inv_M_y = chol_solve(L_M_y, eye(size(M_y)));
    C_uy_y0 = inv_M_y;

    %%  Compute the Posterior Mean (mu_ux_y)
    % Original: mu_ux_y = C_ux_y * (rho_x * Px' * inv(C_dx + C_ex) * sum(y_obs, 2) + inv(C_ux) * mu_ux);
    right_hand_side_x = rho_x * (Px' * inv_C_dx_C_ex * sum(Y_exp_x(:, 1:nrep), 2)) + (inv_C_ux * mu_ux);
    mu_ux_y0 = C_ux_y0 * right_hand_side_x;

    %%  Compute the Posterior Mean (mu_uy_y)
    right_hand_side_y = rho_y * (Py' * inv_C_dy_C_ey * sum(Y_exp_y(:, 1:nrep), 2)) + (inv_C_uy * mu_uy);
    mu_uy_y0 = C_uy_y0 * right_hand_side_y;

    mu_ux_y = rho_x * mu_ux_y0;
    mu_uy_y = rho_y * mu_uy_y0;

    C_ux_y = rho_x ^ 2 * C_ux_y0;
    C_uy_y = rho_y ^ 2 * C_uy_y0;

    mu_u_y(1:2:GDOFs) = mu_ux_y;
    mu_u_y(2:2:GDOFs) = mu_ux_y;

    C_u_y(1:2:end, 1:2:end) = C_ux_y;
    C_u_y(2:2:end, 2:2:end) = C_uy_y;

    std_ux_y = sqrt(diag(C_ux_y));
    std_uy_y = sqrt(diag(C_uy_y));
    %
    mu_ux_y_Surf = makeSurf(elementNodes, mu_ux_y);
    mu_uy_y_Surf = makeSurf(elementNodes, mu_uy_y);

    std_ux_y_Surf = makeSurf(elementNodes, std_ux_y);
    std_uy_y_Surf = makeSurf(elementNodes, std_uy_y);

    % True response in x and y directions
    mu_zx = Px * mu_ux_y;
    mu_zy = Py * mu_uy_y;
    cov_zx = Px * C_ux_y * Px' + C_dx;
    cov_zy = Py * C_uy_y * Py' + C_dy;
    std_zx = sqrt(diag(cov_zx));
    std_zy = sqrt(diag(cov_zy));

    %% Assign Results Back to BVP
    BVP.euclied.HE.mu_u_y = mu_u_y;
    BVP.euclied.HE.mu_ux_y = mu_ux_y;
    BVP.euclied.HE.mu_uy_y = mu_uy_y;
    BVP.euclied.HE.mu_ux_y_Surf = mu_ux_y_Surf;
    BVP.euclied.HE.mu_uy_y_Surf = mu_uy_y_Surf;
    BVP.euclied.HE.cov_u_y = C_u_y;
    BVP.euclied.HE.cov_ux_y = C_ux_y;
    BVP.euclied.HE.cov_uy_y = C_uy_y;
    BVP.euclied.HE.std_ux_y = std_ux_y;
    BVP.euclied.HE.std_uy_y = std_uy_y;
    BVP.euclied.HE.std_ux_y_Surf = std_ux_y_Surf;
    BVP.euclied.HE.std_uy_y_Surf = std_uy_y_Surf;

    BVP.euclied.HE.mu_zx = mu_zx;
    BVP.euclied.HE.mu_zy = mu_zy;
    BVP.euclied.HE.cov_zx = cov_zx;
    BVP.euclied.HE.cov_zy = cov_zy;
    BVP.euclied.HE.std_zx = std_zx;
    BVP.euclied.HE.std_zy = std_zy;
end
