%% statFEM-EUCLID: Main Driver Script
% Author: Vahab Narouie, TU-Braunschweig, 2024
% Description: Executes preprocessing, model discovery, hyperparameter estimation,
%              and iterative statFEM refinement for a quarter-plate with hole.
% License: GNU GPL v3.0

clearvars; close all; clc;

%% Add Required Paths
addpath('../lib');
addpath('../Meshes/');

BVP = [];

%% Step 1: Preprocessing
BVP = preprocess(BVP);

%% Step 2: Linear Elastic FEM (Deterministic)
BVP = process_LE(BVP);
fprintf('Linear Elastic - Top-right node: ux = %.6f, uy = %.6f\n', ...
    BVP.proc.LE.ux(BVP.preProc.msh.top_right_node), ...
    BVP.proc.LE.uy(BVP.preProc.msh.top_right_node));

%% Step 3: Hyperelastic FEM (Deterministic)
BVP = process_HE(BVP);
fprintf('Hyperelastic - Top-right node: ux = %.6f, uy = %.6f\n', ...
    BVP.proc.NH.ux(BVP.preProc.msh.top_right_node), ...
    BVP.proc.NH.uy(BVP.preProc.msh.top_right_node));

%% Step 4: Polynomial Chaos Expansion (Linear Elastic)
BVP = process_LE_PC(BVP);

%% Step 5: Define Sensor and Noise Settings
noise_vector   = [1e-4];   % Noise levels
sensor_vector  = [38];     % Sensor counts
criterion      = 'RMSE';   % Model selection: 'RMSE', 'AIC', 'BIC', 'sparse'
nrep           = 1;        % Number of synthetic repetitions
MR_order       = 3;        % Maximum order for model discovery

%% Step 6: Loop Over Configurations
for inoise = 1:length(noise_vector)
    BVP.preProc.statFEM.sig_e = noise_vector(inoise);
    sigma_e     = BVP.preProc.statFEM.sig_e;
    noise_floor = sqrt(2) * sigma_e;

    for isensor = 1:length(sensor_vector)
        nSen = sensor_vector(isensor);

        % Observation Setup
        if nSen == 3
            BVP = observation_case_3Sensors(BVP, nSen, nrep);
        else
            BVP = observation_case(BVP, nSen, nrep);
        end

        BVP.discovered.RMSE = [];

        %% Step 7: Generate Projection Matrix
        BVP = projection_tri(BVP);

        %% Step 8: Synthetic Data Generation and Model Discovery
        if nSen == 3
            BVP = syntheticData_generation_3Sensors(BVP);
        else
            BVP = syntheticData_generation(BVP);
            BVP = model_discovery_Euclied(BVP, MR_order, criterion);
        end

        %% Step 9: Hyperparameter Estimation (Linear Elastic)
        BVP = hyperparameter_id_fmincon_withoutD(BVP, nrep, 'LE', 0);

        %% Step 10: StatFEM Update (Linear Elastic)
        BVP = statFEM_update_LE(BVP);
        fprintf('StatFEM LE - Top-right node: ux = %.6f, uy = %.6f\n', ...
            BVP.statFEM.LE.mu_ux_y(BVP.preProc.msh.top_right_node), ...
            BVP.statFEM.LE.mu_uy_y(BVP.preProc.msh.top_right_node));

        BVP = get_RMSE(BVP, 1, 'LE');
        BVP = model_discovery(BVP, MR_order, 'LE', criterion);

        %% Step 11: Hyperelastic Model Refinement Loop
        max_iter    = 10;
        rmse_tol    = 1e-1;
        sigmad_d_ub = 2;

        for iter = 1:max_iter
            fprintf('\nIteration %d: Hyperelastic refinement\n', iter);

            BVP = process_HE_PC(BVP);
            fprintf('Hyperelastic PC - Top-right node: ux = %.6f, uy = %.6f\n', ...
                BVP.euclied.HE.mean_ux(BVP.preProc.msh.top_right_node), ...
                BVP.euclied.HE.mean_uy(BVP.preProc.msh.top_right_node));

            BVP = hyperparameter_id_fmincon_withoutD(BVP, nrep, 'HE', sigmad_d_ub);
            BVP = statFEM_update_HE(BVP);

            fprintf('StatFEM HE - Top-right node: ux = %.6f, uy = %.6f\n', ...
                BVP.euclied.HE.mean_ux(BVP.preProc.msh.top_right_node), ...
                BVP.euclied.HE.mean_uy(BVP.preProc.msh.top_right_node));

            BVP = get_RMSE(BVP, 1, 'HE');
            current_rmse = BVP.discovered.RMSE(end);

            % Check convergence with noise floor
            if current_rmse <= noise_floor
                fprintf('Stopping: RMSE %.3e is below noise floor %.3e\n', current_rmse, noise_floor);
                break;
            end

            % Check convergence with tolerance
            if iter > 1
                rmse_prev = BVP.discovered.RMSE(end - 1);
                delta_rmse = abs(current_rmse - rmse_prev);

                if delta_rmse < 1e-8
                    fprintf('Stopping: RMSE change below tolerance (%.2e)\n', delta_rmse);
                    break;
                elseif current_rmse > rmse_prev + rmse_tol
                    fprintf('Warning: RMSE increased. Performing recovery.\n');
                    improvement = true;

                    while improvement
                        sigmad_d_ub = sigmad_d_ub + 0.25;
                        BVP = hyperparameter_id_fmincon_withoutD(BVP, nrep, 'HE', sigmad_d_ub);
                        BVP = statFEM_update_HE(BVP);
                        recovery_rmse = get_RMSE0(BVP, 1, 'HE');

                        delta = rmse_prev - recovery_rmse;
                        improvement = (delta < 0);

                        if improvement
                            BVP.discovered.RMSE(end) = recovery_rmse;
                        end

                        if abs(delta) < rmse_tol
                            fprintf('Recovery: Improvement below tolerance. Exiting.\n');
                            break;
                        end
                    end

                    current_rmse = recovery_rmse;
                end
            end

            % Perform further model discovery if needed
            if current_rmse > noise_floor
                BVP = model_discovery(BVP, MR_order, 'HE', criterion);
            else
                fprintf('Stopping: RMSE below noise floor (%.6f)\n', current_rmse);
                break;
            end

            BVP.iteration = iter;
            disp('Material parameters after refinement:');
            disp(BVP.discovered.hyperelastic_matparameters);
        end

        %% Step 12: Finalize and Save Results
        BVP = chooseBestModel(BVP);
        BVP = displacement_stress_discovered(BVP);
        BVP = displacement_stress_error(BVP);

        % save(sprintf('BVP_nsen%d_%.0e.mat', nSen, sigma_e), 'BVP');

        % Clean up results before next run
        BVP.discovered.identified_hyperparameters{1} = [];
        BVP.discovered.identified_hyperparameters{2} = [];
        BVP.discovered.hyperelastic_matparameters = BVP.preProc.material.prop0;
        BVP.discovered.RMSE = [];
        BVP.discovered.epsilon_u_iter = [];
        BVP.discovered.epsilon_kappa_iter = [];
        BVP.discovered.epsilon_I_iter = [];
        BVP.discovered.w_latex_full_iter = [];
        BVP.discovered.relative_rmse_iter = [];
    end
end
