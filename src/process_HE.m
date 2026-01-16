function BVP = process_HE(BVP)
disp('5. FEM analysis with Hyperelastic material...')
%% Assign from BVP
elementNodes = BVP.preProc.msh.elementNodes;
force = BVP.preProc.BC.force;
%% Calculation
% This is the calculation of displacement field
[u, force, Krr, Kru, Kur, Kuu, Ur, Uu, Rr, Ru, stretch_history, energy_history, W_total0, J_min_max] = solver_HE(BVP.preProc, force, BVP.preProc.material.prop0);
[cauchyStress_xx, cauchyStress_yy, cauchyStress_xy] = stress_HE(BVP.preProc, u, BVP.preProc.material.prop0);
%%
ux = u(1:2:end, :);
uy = u(2:2:end, :);
ux_Surf = makeSurf(elementNodes, ux);
uy_Surf = makeSurf(elementNodes, uy);
cauchyStress_xx_Surf = makeSurf(elementNodes, cauchyStress_xx);
cauchyStress_yy_Surf = makeSurf(elementNodes, cauchyStress_yy);
cauchyStress_xy_Surf = makeSurf(elementNodes, cauchyStress_xy);
% Sample the J_min_max values based on the linspace
J1_samples = linspace(J_min_max(1,1), J_min_max(1,2), 100);
J2_samples = linspace(J_min_max(2,1), J_min_max(2,2), 100);
J3_samples = linspace(J_min_max(3,1), J_min_max(3,2), 100);

%% Assign back to BVP
BVP.proc.NH.u = u;
BVP.proc.NH.force = force;
BVP.proc.NH.ux = ux;
BVP.proc.NH.uy = uy;
BVP.proc.NH.ux_Surf = ux_Surf;
BVP.proc.NH.uy_Surf = uy_Surf;
BVP.proc.NH.Krr = Krr;
BVP.proc.NH.Kru = Kru;
BVP.proc.NH.Kur = Kur;
BVP.proc.NH.Kuu = Kuu;
BVP.proc.NH.Ur = Ur;
BVP.proc.NH.Uu = Uu;
BVP.proc.NH.Ru = Ru;
BVP.proc.NH.Rr = Rr;
BVP.proc.NH.W_total0 = W_total0;
BVP.proc.NH.stretch_history = stretch_history;
BVP.proc.NH.energy_history = energy_history;
BVP.proc.NH.cauchyStress_xx = cauchyStress_xx;
BVP.proc.NH.cauchyStress_yy = cauchyStress_yy;
BVP.proc.NH.cauchyStress_xy = cauchyStress_xy;
BVP.proc.NH.cauchyStress_xx_Surf = cauchyStress_xx_Surf;
BVP.proc.NH.cauchyStress_yy_Surf = cauchyStress_yy_Surf;
BVP.proc.NH.cauchyStress_xy_Surf = cauchyStress_xy_Surf;
BVP.proc.NH.J_min_max = J_min_max;
BVP.proc.NH.J1_samples = J1_samples;
BVP.proc.NH.J2_samples = J2_samples;
BVP.proc.NH.J3_samples = J3_samples;
end
