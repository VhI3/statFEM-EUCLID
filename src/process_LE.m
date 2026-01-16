function BVP = process_LE(BVP)
    disp('2. FEM Linear Elastic Processing')
    %% Assign from BVP
    elementNodes = BVP.preProc.msh.elementNodes;
    force = BVP.preProc.BC.force;
    %% Calculation
    % This is the calculation of displacement field
    [u, force, Krr, Kru, Kur, Kuu, Ur, Uu, Rr, Ru] = solver_LE(BVP.preProc, force);
    %%
    ux = u(1:2:end, :);
    uy = u(2:2:end, :);
    ux_Surf = makeSurf(elementNodes, ux);
    uy_Surf = makeSurf(elementNodes, uy);
    %% Assign back to BVP
    BVP.proc.LE.u = u;
    BVP.proc.LE.force = force;
    BVP.proc.LE.ux = ux;
    BVP.proc.LE.uy = uy;
    BVP.proc.LE.ux_Surf = ux_Surf;
    BVP.proc.LE.uy_Surf = uy_Surf;
    BVP.proc.LE.Krr = Krr;
    BVP.proc.LE.Kru = Kru;
    BVP.proc.LE.Kur = Kur;
    BVP.proc.LE.Kuu = Kuu;
    BVP.proc.LE.Ur = Ur;
    BVP.proc.LE.Uu = Uu;
    BVP.proc.LE.Ru = Ru;
    BVP.proc.LE.Rr = Rr;
end
