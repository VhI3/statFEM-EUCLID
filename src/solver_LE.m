function [u, force, Krr, Kru, Kur, Kuu, Ur, Uu, Rr, Ru] = solver_LE(preProc_Variables, force)
    %% Assign from BVP
    activeDOFs = preProc_Variables.BC.activeDOFs;
    prescribedDOFs = preProc_Variables.BC.prescribedDOFs;
    GDOFs = preProc_Variables.msh.GDOFs;
    nElm = preProc_Variables.msh.nElm;
    nodeCoordinates = preProc_Variables.msh.nodeCoordinates;
    et = preProc_Variables.msh.elementNodes;
    ed = preProc_Variables.msh.elementDOFs;
    tickness = preProc_Variables.msh.T;
    C_Constitutive = preProc_Variables.material.C_Constitutive;
    Q_pts = preProc_Variables.Q_pts;
    Q_wts = preProc_Variables.Q_wts;
    dN_dxi = preProc_Variables.dN_dxi;
    dN_deta = preProc_Variables.dN_deta;
    E_linear = preProc_Variables.material.E_linear;
    %% Calcualtion
    % calculation of the system stiffness matrix
    dtan = E_linear * C_Constitutive;
    u = zeros(GDOFs, 1);
    stiffness = zeros(GDOFs, GDOFs); % reserve stiffness matrix

    for e = 1:nElm % loop over elements
        indice = et(e, :);
        elementDof = ed(e, :);
        % E_e = E_vector(indice);
        ke = zeros(size(elementDof, 2));

        for i = 1:size(Q_wts, 1)
            xi = Q_pts(i, 1); eta = Q_pts(i, 2);
            naturalDerivatives = [dN_dxi(xi, eta); dN_deta(xi, eta)]';
            %
            [JacobianMatrix, ~, XYDerivatives] = Jacobian2D(nodeCoordinates(indice, :), naturalDerivatives);
            % B matrix
            B = zeros(3, size(elementDof, 2));
            B(1, 1:2:end) = XYDerivatives(:, 1)';
            B(2, 2:2:end) = XYDerivatives(:, 2)';
            B(3, 1:2:end) = XYDerivatives(:, 2)';
            B(3, 2:2:end) = XYDerivatives(:, 1)';
            %
            ke = ke + Q_wts(i) * (B' * dtan * B) * det(JacobianMatrix) * tickness;
        end

        stiffness(elementDof, elementDof) = stiffness(elementDof, elementDof) + ke;
    end

    %
    Krr = sparse(stiffness(activeDOFs, activeDOFs));
    Kru = sparse(stiffness(activeDOFs, prescribedDOFs));
    Kur = sparse(stiffness(prescribedDOFs, activeDOFs));
    Kuu = sparse(stiffness(prescribedDOFs, prescribedDOFs));
    %
    % Solution
    % static analysis
    Rr = sparse(force(activeDOFs));
    Uu = sparse(u(prescribedDOFs));
    Ur = Krr \ (Rr - Kru * Uu);
    Ru = Kuu * Uu + Kur * Ur;
    %
    u(activeDOFs) = Ur;
    force(prescribedDOFs) = Ru;
end
