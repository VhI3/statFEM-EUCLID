function BVP = displacement_stress_discovered(BVP)

    %% Assign from BVP
    elementNodes = BVP.preProc.msh.elementNodes;
    force = BVP.preProc.BC.force;
    %% Calculation
    % This is the calculation of displacement field
    [u, force, Krr, Kru, Kur, Kuu, Ur, Uu, Rr, Ru, stretch_history, energy_history, W_total, J_min_max] = solver_HE(BVP.preProc, force, BVP.discovered.best_model.prop);
    [cauchyStress_xx, cauchyStress_yy, cauchyStress_xy] = stress_HE(BVP.preProc, u, BVP.discovered.best_model.prop);
    %%
    ux = u(1:2:end, :);
    uy = u(2:2:end, :);
    ux_Surf = makeSurf(elementNodes, ux);
    uy_Surf = makeSurf(elementNodes, uy);
    cauchyStress_xx_Surf = makeSurf(elementNodes, cauchyStress_xx);
    cauchyStress_yy_Surf = makeSurf(elementNodes, cauchyStress_yy);
    cauchyStress_xy_Surf = makeSurf(elementNodes, cauchyStress_xy);
    %% Assign back to BVP
    BVP.discovered.best_model.u = u;
    BVP.discovered.best_model.force = force;
    BVP.discovered.best_model.ux = ux;
    BVP.discovered.best_model.uy = uy;
    BVP.discovered.best_model.ux_Surf = ux_Surf;
    BVP.discovered.best_model.uy_Surf = uy_Surf;
    BVP.discovered.best_model.Krr = Krr;
    BVP.discovered.best_model.Kru = Kru;
    BVP.discovered.best_model.Kur = Kur;
    BVP.discovered.best_model.Kuu = Kuu;
    BVP.discovered.best_model.Ur = Ur;
    BVP.discovered.best_model.Uu = Uu;
    BVP.discovered.best_model.Ru = Ru;
    BVP.discovered.best_model.Rr = Rr;
    BVP.discovered.best_model.W_total = W_total;
    BVP.discovered.best_model.J_min_max = J_min_max;
    BVP.discovered.best_model.stretch_history = stretch_history;
    BVP.discovered.best_model.energy_history = energy_history;
    BVP.discovered.best_model.cauchyStress_xx = cauchyStress_xx;
    BVP.discovered.best_model.cauchyStress_yy = cauchyStress_yy;
    BVP.discovered.best_model.cauchyStress_xy = cauchyStress_xy;
    BVP.discovered.best_model.cauchyStress_xx_Surf = cauchyStress_xx_Surf;
    BVP.discovered.best_model.cauchyStress_yy_Surf = cauchyStress_yy_Surf;
    BVP.discovered.best_model.cauchyStress_xy_Surf = cauchyStress_xy_Surf;
end
