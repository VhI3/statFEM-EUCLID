function BVP = process_LE_PC(BVP)
    %% Assign from BVP
    % Mesh and boundary conditions
    GDOFs = BVP.preProc.msh.GDOFs; % Global degrees of freedom
    activeDOFsX = BVP.preProc.BC.activeDOFsX; % Active DOFs in X direction
    activeDOFsY = BVP.preProc.BC.activeDOFsY; % Active DOFs in Y direction
    activeDOFs = BVP.preProc.BC.activeDOFs; % Active DOFs
    elementNodes = BVP.preProc.msh.elementNodes; % Element connectivity
    nuNodes = BVP.preProc.msh.nuNodes;
    top_right_node = BVP.preProc.msh.top_right_node;

    xi = BVP.preProc.UQ.xi;
    nSampled_force = BVP.preProc.UQ.nSampled_force;
    sampled_force = BVP.preProc.UQ.sampled_force;
    P_u = BVP.preProc.UQ.P_u;
    Psi_s = BVP.preProc.UQ.Psi_s;
    PsiSqNorm = BVP.preProc.UQ.PsiSqNorm;
    xi_MC = BVP.preProc.UQ.xi_MC;
    nMC = BVP.preProc.UQ.nMC;
    %% Initialize Storage for Displacement and Strain
    bigPsi = zeros(nSampled_force, P_u); % Basis matrix
    u_realization = zeros(nSampled_force, GDOFs); % Displacement realizations
    cov_u_pc = zeros(GDOFs); % Covariance matrix for displacement
    %% Process Using Least-Squares NIPCE
    for i = 1:nSampled_force
        % Compute basis functions for the sample
        for j = 1:P_u
            xi_sym = sym(sprintf('%.17g', xi(i))); % Convert to a symbolic rational value
            bigPsi(i, j) = (1 / sqrt(factorial(j))) * double(subs(Psi_s{j, 1}, xi_sym));
        end

        % Solve FEM problem for the current realization
        [u_realization(i, :), ~, ~, ~, ~, ~, ~, ~, ~, ~] = solver_LE(BVP.preProc, sampled_force(:, i));
    end

    %% Compute Statistical Quantities
    % Mean displacement coefficients
    u_NIPC = (bigPsi' * bigPsi) \ (bigPsi' * u_realization);
    u_NIPC = u_NIPC';
    mean_u_pc = u_NIPC(:, 1);
    mean_ux_pc = mean_u_pc(1:2:end);
    mean_uy_pc = mean_u_pc(2:2:end);
    mean_ux_active_pc = mean_u_pc(activeDOFsX);
    mean_uy_active_pc = mean_u_pc(activeDOFsY);
    mean_u_active_pc = mean_u_pc(activeDOFs);
    % Covariance of displacement
    for j = 2:P_u
        cov_u_pc = cov_u_pc + (1 / sqrt(factorial(j))) ^ 2 * PsiSqNorm(j) * u_NIPC(:, j) * u_NIPC(:, j)';
    end

    cov_u_pc = 0.5 * (cov_u_pc + cov_u_pc'); % Ensure symmetry
    ci_cov_u_pc = sqrt(diag(cov_u_pc)) * 1.96;
    cov_ux_pc = cov_u_pc(1:2:end, 1:2:end);
    cov_uy_pc = cov_u_pc(2:2:end, 2:2:end);
    cov_ux_active_pc = cov_u_pc(activeDOFsX, activeDOFsX);
    cov_uy_active_pc = cov_u_pc(activeDOFsY, activeDOFsY);
    cov_u_active_pc = cov_u_pc(activeDOFs, activeDOFs);
    % Surface representation
    mean_ux_pc_Surf = makeSurf(elementNodes, mean_ux_pc);
    std_ux_pc = sqrt(diag(cov_u_pc(1:2:end, 1:2:end)));
    std_ux_pc_Surf = makeSurf(elementNodes, std_ux_pc);
    mean_uy_pc_Surf = makeSurf(elementNodes, mean_uy_pc);
    std_uy_pc = sqrt(diag(cov_u_pc(2:2:end, 2:2:end)));
    std_uy_pc_Surf = makeSurf(elementNodes, std_uy_pc);

    %% Representation of the Response Displacement
    % Initialize displacement storage
    u_samples = zeros(GDOFs, nMC);

    % Define symbolic variable for PCE evaluation
    xi_v = sym('xi_1'); % Symbolic variable
    xi_n = xi_MC; % Numeric samples

    % Evaluate PCE Expansion for Displacement
    for j = 0:(P_u - 1)
        fprintf('Evaluating PC expansion %g/%g\n', j, P_u - 1);

        % Extract PCE polynomial basis function
        psi_j = double(subs(Psi_s{j + 1}, xi_v, xi_n));

        % Extract corresponding PCE coefficient
        d_j = u_NIPC(:, j + 1);

        % Compute displacement realization
        for i = 1:nuNodes
            u_samples(i, :) = u_samples(i, :) + d_j(i) * psi_j';
        end

    end

    ux_samples = u_samples(1:2:end, :);
    uy_samples = u_samples(2:2:end, :);

    %% **Extract & Analyze Displacement at Top Right Corner**
    ux_samples_top_right = ux_samples(top_right_node, :);
    uy_samples_top_right = uy_samples(top_right_node, :);

    % % Compute Probability Density Function (PDF)
    [pdf_ux_top_right, x_pdf_ux_top_right] = ksdensity(ux_samples_top_right); % pdf in x direction
    [pdf_uy_top_right, x_pdf_uy_top_right] = ksdensity(uy_samples_top_right); % pdf in y direction

    % % Compute Cumulative Distribution Function (CDF)
    [cdf_ux_top_right, x_cdf_ux_top_right] = ecdf(ux_samples_top_right); % cdf in x direction
    [cdf_uy_top_right, x_cdf_uy_top_right] = ecdf(uy_samples_top_right); % cdf in y direction

    %% Assign Results Back to BVP
    BVP.lsNIPCE.LE.mean_u = mean_u_pc;
    BVP.lsNIPCE.LE.mean_ux = mean_ux_pc;
    BVP.lsNIPCE.LE.mean_uy = mean_uy_pc;
    BVP.lsNIPCE.LE.mean_ux_Surf = mean_ux_pc_Surf;
    BVP.lsNIPCE.LE.mean_uy_Surf = mean_uy_pc_Surf;
    BVP.lsNIPCE.LE.mean_u_active = mean_u_active_pc;
    BVP.lsNIPCE.LE.mean_ux_active = mean_ux_active_pc;
    BVP.lsNIPCE.LE.mean_uy_active = mean_uy_active_pc;
    BVP.lsNIPCE.LE.std_ux_Surf = std_ux_pc_Surf;
    BVP.lsNIPCE.LE.std_uy_Surf = std_uy_pc_Surf;
    BVP.lsNIPCE.LE.cov_u = cov_u_pc;
    BVP.lsNIPCE.LE.ci_cov_u = ci_cov_u_pc;
    BVP.lsNIPCE.LE.cov_ux_active = cov_ux_active_pc;
    BVP.lsNIPCE.LE.cov_uy_active = cov_uy_active_pc;
    BVP.lsNIPCE.LE.cov_u_active = cov_u_active_pc;
    BVP.lsNIPCE.LE.cov_ux = cov_ux_pc;
    BVP.lsNIPCE.LE.cov_uy = cov_uy_pc;
    BVP.lsNIPCE.LE.pdf_ux_top_right = pdf_ux_top_right;
    BVP.lsNIPCE.LE.x_pdf_ux_top_right = x_pdf_ux_top_right;
    BVP.lsNIPCE.LE.cdf_ux_top_right = cdf_ux_top_right;
    BVP.lsNIPCE.LE.x_cdf_ux_top_right = x_cdf_ux_top_right;
    BVP.lsNIPCE.LE.ux_samples_top_right = ux_samples_top_right;
    BVP.lsNIPCE.LE.uy_samples_top_right = uy_samples_top_right;
    BVP.lsNIPCE.LE.pdf_uy_top_right = pdf_uy_top_right;
    BVP.lsNIPCE.LE.x_pdf_uy_top_right = x_pdf_uy_top_right;
    BVP.lsNIPCE.LE.cdf_uy_top_right = cdf_uy_top_right;
    BVP.lsNIPCE.LE.x_cdf_uy_top_right = x_cdf_uy_top_right;
end
