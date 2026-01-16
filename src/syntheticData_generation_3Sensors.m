function BVP = syntheticData_generation_3Sensors(BVP)
    %% Assign from BVP

    GSensorDOFs = BVP.obs.GSensorDOFs;
    nrep = BVP.obs.nrep;
    nsen = BVP.obs.nSen;
    e_sampled = BVP.obs.e_sampled;
    ex_sampled = BVP.obs.ex_sampled;
    ey_sampled = BVP.obs.ey_sampled;

    u_samples = BVP.proc.NH.u;
    ux_samples = BVP.proc.NH.ux;
    uy_samples = BVP.proc.NH.uy;

    P_full = BVP.obs.P_full;
    Px = BVP.obs.Px;
    Py = BVP.obs.Py;
    % sensor_connectivity = BVP.obs.sensor_connectivity;
    activeSensorDOFsX = BVP.obs.activeSensorDOFsX;
    activeSensorDOFsY = BVP.obs.activeSensorDOFsY;
    % BVP_sensor = BVP.obs.BVP_sensor;
    % force = BVP.obs.sensor_force;
    % prop0 = BVP.preProc.material.prop0;
    %% Calculation
    % To be sure, if the observation data are generated in the correct way, we
    % generate the observation data in 3 ways. Once based on the full form of the
    % observation matrix and once based on the decomposed version in x and y directions
    % and once based on the active form of the observation matrix
    % The first one is the full form of the observation matrix
    % Full form of the observation matrix

    Y_exp_full = zeros(GSensorDOFs, nrep);
    True_exp = zeros(GSensorDOFs, nrep);

    % [Y_simulated, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = solver_HE(BVP_sensor, force, prop0);

    for irep = 1:nrep
        % Generation of obsrtvatin
        True_exp(:, irep) = P_full * u_samples(:, :);
        Y_exp_full(:, irep) = P_full * u_samples(:, :) + e_sampled(:, irep);

        % True_exp(:, irep) = Y_simulated(:, :);
        % Y_exp_full(:, irep) = Y_simulated(:, :) + e_sampled(:, irep);

    end

    % The second one is the decomposed version in x and y directions
    % Generating of observation data in x and y directions
    Y_exp_x = zeros(nsen, nrep);
    Y_exp_y = zeros(nsen, nrep);

    for i = 1:nrep
        Y_exp_x(:, i) = Px * ux_samples(:, :) + ex_sampled(:, i);
        Y_exp_y(:, i) = Py * uy_samples(:, :) + ey_sampled(:, i);
    end

    Y_exp_full0 = zeros(GSensorDOFs, nrep);
    Y_exp_full0(1:2:end, :) = Y_exp_x;
    Y_exp_full0(2:2:end, :) = Y_exp_y;

    norm_full = norm(Y_exp_full - Y_exp_full0);

    if norm_full < 1e-10
        disp('The full form of the observation matrix is correct!');
    else
        disp('Error detected: The full form of the observation matrix may be incorrect.');
    end

    % The third one is the active form of the observation matrix
    Y_exp_x_active = Y_exp_full(activeSensorDOFsX, :);
    Y_exp_y_active = Y_exp_full(activeSensorDOFsY, :);

    % Compute the sum and mean of the generated observations across all realizations
    sum_yObs = {};
    sum_yObs{1} = sum(Y_exp_x, 2);
    sum_yObs{2} = sum(Y_exp_y, 2);

    % mean of observations
    mean_Y_exp = {};
    mean_Y_exp{1} = sum_yObs{1} / nrep;
    mean_Y_exp{2} = sum_yObs{2} / nrep;

    %% Assign back to BVP
    BVP.obs.Y_exp_x = Y_exp_x;
    BVP.obs.Y_exp_y = Y_exp_y;
    BVP.obs.Y_exp = Y_exp_full;
    BVP.obs.True_exp = True_exp;
    BVP.obs.Y_exp_x_active = Y_exp_x_active;
    BVP.obs.Y_exp_y_active = Y_exp_y_active;
    BVP.obs.sum_yObs = sum_yObs; % Sum of observations
    BVP.obs.mean_Y_exp = mean_Y_exp; % Mean of observations
    % BVP.obs.mean_Y_exp_x_Surf = mean_Y_exp_x_Surf;
    % BVP.obs.mean_Y_exp_y_Surf = mean_Y_exp_y_Surf;
end
