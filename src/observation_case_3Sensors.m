function BVP = observation_case_3Sensors(BVP, nSen, nrep)
    %% Assign from BVP
    DOFs = BVP.preProc.msh.DOFs;
    sig_e = BVP.preProc.statFEM.sig_e;
    L = BVP.preProc.msh.L;
    H = BVP.preProc.msh.H;
    tractionX = BVP.preProc.BC.tractionX;
    tractionY = BVP.preProc.BC.tractionY;
    %% Calculation
    switch nSen
        case 7
            quarterPlate_L100_esf100
        case 13
            quarterPlate_L100_esf90
        case 15000
            quarterPlate_L100_esf50
        case 38
            quarterPlate_L100_esf25
        case 144
            quarterPlate_L100_esf10
        case 292
            quarterPlate_L100_esf7
        case 8
            quarterPlate_L100_adapt_esf500
        case 15
            quarterPlate_L100_adapt_esf200
        case 28
            quarterPlate_L100_adapt_esf120
        case 33
            quarterPlate_L100_adapt_esf110
        case 39
            quarterPlate_L100_adapt_esf100
        case 54
            quarterPlate_L100_adapt_esf90
        case 57
            quarterPlate_L100_adapt_esf70
        case 3
            quarterPlate_L100_adapt_threeNodes
        otherwise
            error('Wrong number of sensors')
    end

    % Sensor coordinates
    sensorCoordinates = msh.POS(:, 1:end - 1);
    x_sensor = sensorCoordinates(:, 1);
    y_sensor = sensorCoordinates(:, 2);
    % Global degree of freedoms on sensors
    GSensorDOFs = DOFs * nSen;
    u_sensors = zeros(GSensorDOFs, 1);
    sensorCoordinates_vector = zeros(GSensorDOFs, 1);
    sensorCoordinates_vector(1:2:end) = x_sensor;
    sensorCoordinates_vector(2:2:end) = y_sensor;
    % Active sensors and prescribed sensors
    leftSensors = find(abs(x_sensor) < eps * 10);
    bottomSensors = find(abs(y_sensor) < eps * 10);
    rightSensors = find(abs(x_sensor - L) < eps * 10);
    topSensors = find(abs(y_sensor - H) < eps * 10);

    fixedPrescribedSensorDOFsX = DOFs * leftSensors - 1;
    fixedPrescribedSensorDOFsY = DOFs * bottomSensors;
    fixedPrescribedSensorDOFs = sort([fixedPrescribedSensorDOFsX; fixedPrescribedSensorDOFsY]);

    activeSensorDOFsX = setdiff(1:2:GSensorDOFs, fixedPrescribedSensorDOFsX);
    activeSensorDOFsY = setdiff(2:2:GSensorDOFs, fixedPrescribedSensorDOFsY);
    activeSensorDOFs = setdiff(1:GSensorDOFs, fixedPrescribedSensorDOFs);
    prescribedDOFs = sort([fixedPrescribedSensorDOFs]);

    sensorCoordinates_active_vector = sensorCoordinates_vector(activeSensorDOFs);
    x_sensor_active = sensorCoordinates_vector(activeSensorDOFsX);
    y_sensor_active = sensorCoordinates_vector(activeSensorDOFsY);

    % Number of active sensors in x and y directions
    nSen_active{1} = length(x_sensor_active);
    nSen_active{2} = length(y_sensor_active);

    C_e0{1} = sig_e * eye(nSen);
    C_e0{2} = sig_e * eye(nSen);

    C_e_full = zeros(GSensorDOFs, GSensorDOFs);
    C_e_full(1:2:end, 1:2:end) = C_e0{1};
    C_e_full(2:2:end, 2:2:end) = C_e0{2};

    C_e_full(fixedPrescribedSensorDOFs, fixedPrescribedSensorDOFs) = zeros(length(fixedPrescribedSensorDOFs), length(fixedPrescribedSensorDOFs));

    C_e{1} = C_e_full(1:2:end, 1:2:end);
    C_e{2} = C_e_full(2:2:end, 2:2:end);

    C_e_active{1} = C_e_full(activeSensorDOFsX, activeSensorDOFsX);
    C_e_active{2} = C_e_full(activeSensorDOFsY, activeSensorDOFsY);

    rng(0);
    e_x = mvnrnd(zeros(nSen, 1), C_e{1}, nrep)';
    rng(0);
    e_y = mvnrnd(zeros(nSen, 1), C_e{2}, nrep)';

    e_sampled = zeros(GSensorDOFs, nrep);
    e_sampled(1:2:end, :) = e_x;
    e_sampled(2:2:end, :) = e_y;


    %%
    sensor_force = zeros(GSensorDOFs, 1);

    % Sort rightSesnsors by their Y coordinate (2nd column of nodeCoordinates)
    [~, sortIdx] = sort(sensorCoordinates(rightSensors, 2), 'ascend');
    rightSensorCorr = rightSensors(sortIdx);
    rightEdge_y = sensorCoordinates(rightSensorCorr, 2); % Y-coordinates
    rightEdge_y_pairs = [rightEdge_y(1:end - 1) rightEdge_y(2:end)];

    rightEdge_dofX = DOFs * rightSensorCorr - 1; % X DOFs only
    rightEdge_dofX_pairs = [rightEdge_dofX(1:end - 1) rightEdge_dofX(2:end)];

    for ey = 1:size(rightEdge_y_pairs, 1)
        ed = rightEdge_dofX_pairs(ey, :);
        dy = abs(rightEdge_y_pairs(ey, 2) - rightEdge_y_pairs(ey, 1));
        sensor_force(ed, 1) = sensor_force(ed, 1) + dy / 2 * tractionX * ones(2, 1);
    end

    % Top nodes
    [~, sortIdy] = sort(sensorCoordinates(topSensors, 1), 'ascend');
    topSensorCorr = topSensors(sortIdy);
    topEdge_x = sensorCoordinates(topSensorCorr, 1); % X-coordinates
    topEdge_x_pairs = [topEdge_x(1:end - 1) topEdge_x(2:end)];
    %
    topEdge_dofY = DOFs * topSensorCorr; % Y DOFs only
    topEdge_dofY_pairs = [topEdge_dofY(1:end - 1) topEdge_dofY(2:end)];

    % Deterministic force
    for ex = 1:size(topEdge_x_pairs, 1)
        ed = topEdge_dofY_pairs(ex, :);
        dx = abs(topEdge_x_pairs(ex, 2) - topEdge_x_pairs(ex, 1));
        sensor_force(ed, 1) = sensor_force(ed, 1) + dx / 2 * tractionY * ones(2, 1);
    end

    tbc_sensor = [];
    tbc_sensor = [tbc_sensor rightSensors];
    % An indicator if it is forced control
    ntbc_sensor = size(tbc_sensor, 1);

    ndbc_sensor = size(prescribedDOFs, 1);
    dbc_sensor = u_sensors(prescribedDOFs);

    % % To calculate the experimental directly with given sensor coordinates
    BVP_sensor.msh.GDOFs = GSensorDOFs;
    BVP_sensor.BC.ndbc = ndbc_sensor;
    BVP_sensor.BC.ntbc = ntbc_sensor;
    BVP_sensor.BC.dbc = dbc_sensor;
    BVP_sensor.BC.prescribedDOFs = prescribedDOFs;
    BVP_sensor.tol = BVP.preProc.tol;
    BVP_sensor.maxit = BVP.preProc.maxit;
    BVP_sensor.maxreit = BVP.preProc.maxreit;
    BVP_sensor.timeInterval = BVP.preProc.timeInterval;

    BVP_sensor.Q_pts = BVP.preProc.Q_pts;
    BVP_sensor.Q_wts = BVP.preProc.Q_wts;
    BVP_sensor.Nshape = BVP.preProc.Nshape;
    BVP_sensor.dN_dxi = BVP.preProc.dN_dxi;
    BVP_sensor.dN_deta = BVP.preProc.dN_deta;
    BVP_sensor.msh.nodeCoordinates = sensorCoordinates;
    BVP_sensor.msh.GDOFs = GSensorDOFs;
    BVP_sensor.msh.T = BVP.preProc.msh.T;
    BVP_sensor.msh.nElNode = 3;
    BVP_sensor.msh.DOFs = 2;

    %% Assign back to BVP
    BVP.obs.sensorCoordinates = sensorCoordinates;
    BVP.obs.sensor_force = sensor_force;
    BVP.obs.nSen = nSen;
    BVP.obs.nrep = nrep;
    BVP.obs.tbc_sensor = tbc_sensor;
    BVP.obs.ntbc_sensor = ntbc_sensor;
    BVP.obs.ndbc_sensor = ndbc_sensor;
    BVP.obs.dbc_sensor = dbc_sensor;
    BVP.obs.nSen_active = nSen_active;
    BVP.obs.x_sensor = x_sensor;
    BVP.obs.y_sensor = y_sensor;
    BVP.obs.fixedPrescribedSensorDOFsX = fixedPrescribedSensorDOFsX;
    BVP.obs.fixedPrescribedSensorDOFsY = fixedPrescribedSensorDOFsY;
    BVP.obs.fixedPrescribedSensorDOFs = fixedPrescribedSensorDOFs;
    BVP.obs.prescribedDOFs = prescribedDOFs;
    BVP.obs.activeSensorDOFsX = activeSensorDOFsX;
    BVP.obs.activeSensorDOFsY = activeSensorDOFsY;
    BVP.obs.activeSensorDOFs = activeSensorDOFs;
    BVP.obs.GSensorDOFs = GSensorDOFs;
    BVP.obs.sensorCoordinates_vector = sensorCoordinates_vector;
    BVP.obs.sensorCoordinates_active_vector = sensorCoordinates_active_vector;
    BVP.obs.x_sensor_active = x_sensor_active;
    BVP.obs.y_sensor_active = y_sensor_active;

    BVP.obs.C_e = C_e;
    BVP.obs.C_e0 = C_e0;
    BVP.obs.C_e_active = C_e_active;
    BVP.obs.e_sampled = e_sampled;
    BVP.obs.ex_sampled = e_x;
    BVP.obs.ey_sampled = e_y;
    BVP.obs.C_e_full = C_e_full;

    BVP.obs.BVP_sensor = BVP_sensor;
end
