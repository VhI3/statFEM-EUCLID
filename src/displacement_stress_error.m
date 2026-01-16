% function BVP = displacement_stress_error(BVP)
%     %% Assign from BVP
%     ux_Surf_discovered = BVP.discovered.best_model.ux_Surf;
%     uy_Surf_discovered = BVP.discovered.best_model.uy_Surf;
% 
%     ux_Surf_True = BVP.proc.NH.ux_Surf;
%     uy_Surf_True = BVP.proc.NH.uy_Surf;
% 
%     cauchyStress_xx_Surf_True = BVP.discovered.best_model.cauchyStress_xx_Surf;
%     cauchyStress_yy_Surf_True = BVP.discovered.best_model.cauchyStress_yy_Surf;
%     cauchyStress_xy_Surf_True = BVP.discovered.best_model.cauchyStress_xy_Surf;
% 
%     cauchyStress_xx_Surf = BVP.proc.NH.cauchyStress_xx_Surf;
%     cauchyStress_yy_Surf = BVP.proc.NH.cauchyStress_yy_Surf;
%     cauchyStress_xy_Surf = BVP.proc.NH.cauchyStress_xy_Surf;
% 
%     %% Calculation
%     % Calculate the norm error in displacement
%     ux_error = sqrt((ux_Surf_discovered - ux_Surf_True) .^ 2); %ux_Surf_discovered - ux_Surf_True;
%     uy_error = sqrt((uy_Surf_discovered - uy_Surf_True) .^ 2); %uy_Surf_discovered - uy_Surf_True;
% 
%     %% Calculate the norm error in stress
%     cauchyStress_xx_error = sqrt((cauchyStress_xx_Surf - cauchyStress_xx_Surf_True) .^ 2); %cauchyStress_xx_Surf - cauchyStress_xx_Surf_True;
%     cauchyStress_yy_error = sqrt((cauchyStress_yy_Surf - cauchyStress_yy_Surf_True) .^ 2); %cauchyStress_yy_Surf - cauchyStress_yy_Surf_True;
%     cauchyStress_xy_error = sqrt((cauchyStress_xy_Surf - cauchyStress_xy_Surf_True) .^ 2); %cauchyStress_xy_Surf - cauchyStress_xy_Surf_True;
% 
%     %%
%     %% Assign back to BVP
%     BVP.discovered.best_model.ux_error = ux_error;
%     BVP.discovered.best_model.uy_error = uy_error;
%     BVP.discovered.best_model.cauchyStress_xx_error = cauchyStress_xx_error;
%     BVP.discovered.best_model.cauchyStress_yy_error = cauchyStress_yy_error;
%     BVP.discovered.best_model.cauchyStress_xy_error = cauchyStress_xy_error;
% end



function BVP = displacement_stress_error(BVP)
    %% Assign from BVP
    % Displacement
    ux_discovered = BVP.discovered.best_model.ux_Surf;
    uy_discovered = BVP.discovered.best_model.uy_Surf;
    ux_true = BVP.proc.NH.ux_Surf;
    uy_true = BVP.proc.NH.uy_Surf;

    % Stress (components)
    sxx_discovered = BVP.discovered.best_model.cauchyStress_xx_Surf;
    syy_discovered = BVP.discovered.best_model.cauchyStress_yy_Surf;
    sxy_discovered = BVP.discovered.best_model.cauchyStress_xy_Surf;

    sxx_true = BVP.proc.NH.cauchyStress_xx_Surf;
    syy_true = BVP.proc.NH.cauchyStress_yy_Surf;
    sxy_true = BVP.proc.NH.cauchyStress_xy_Surf;

    %% Displacement magnitude and error
    disp_discovered = sqrt(ux_discovered.^2 + uy_discovered.^2);
    disp_true = sqrt(ux_true.^2 + uy_true.^2);
    disp_error = abs(disp_discovered - disp_true);

    %% Von Mises stress computation
    vonMises_discovered = sqrt(sxx_discovered.^2 - sxx_discovered .* syy_discovered + syy_discovered.^2 + 3 * sxy_discovered.^2);
    vonMises_true = sqrt(sxx_true.^2 - sxx_true .* syy_true + syy_true.^2 + 3 * sxy_true.^2);
    vonMises_error = abs(vonMises_discovered - vonMises_true);

    %% Assign results back to BVP structure
    BVP.discovered.best_model.disp_error = disp_error;
    BVP.discovered.best_model.vonMises_error = vonMises_error;

    % Optionally store full magnitudes
    BVP.discovered.best_model.disp_discovered = disp_discovered;
    BVP.discovered.best_model.disp_true = disp_true;
    BVP.discovered.best_model.vonMises_discovered = vonMises_discovered;
    BVP.discovered.best_model.vonMises_true = vonMises_true;
end

