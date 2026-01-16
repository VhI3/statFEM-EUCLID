function [cauchyStress_xx, cauchyStress_yy, cauchyStress_xy] = stress_HE(preProc, u, prop)
    % This function calculates the stress fields based on the displacement field u
    % preProc: preprocessed data containing mesh and material properties
    % u: displacement field
    DOFs = preProc.msh.DOFs;
    nElm = preProc.msh.nElm;
    nElNode = preProc.msh.nElNode;
    nodeCoordinates = preProc.msh.nodeCoordinates;
    et = preProc.msh.elementNodes;
    ed = preProc.msh.elementDOFs;
    Q_pts = preProc.Q_pts;
    Q_wts = preProc.Q_wts;
    dN_dxi = preProc.dN_dxi;
    dN_deta = preProc.dN_deta;


    cauchyStress_all = zeros(nElm, size(Q_wts, 1), 3);

    for e = 1:nElm
        indice = et(e, :);
        elementDof = ed(e, :);
        elDisp = u(elementDof);
        elDisp = reshape(elDisp, DOFs, nElNode);
        c = 1; % Counter for Gauss points

        for i = 1:size(Q_wts, 1)
            xi = Q_pts(i, 1); eta = Q_pts(i, 2);
            naturalDerivatives = [dN_dxi(xi, eta); dN_deta(xi, eta)]';
            [~, ~, XYDerivatives] = Jacobian2D(nodeCoordinates(indice, :), naturalDerivatives);

            F = eye(3);
            F(1:2, 1:2) = elDisp * XYDerivatives + eye(2);
            F(3, 3) = 1;

            C11 = F(1, 1) .^ 2 + F(2, 1) .^ 2;
            C12 = F(1, 1) .* F(1, 2) + F(2, 1) .* F(2, 2);
            C21 = F(1, 2) .* F(1, 1) + F(2, 2) .* F(2, 1);
            C22 = F(1, 2) .^ 2 + F(2, 2) .^ 2;

            C = [C11, C12; C21, C22];
            I3 = C(1, 1) .* C(2, 2) - C(1, 2) .* C(2, 1);
            J = sqrt(I3);

            [S, ~, ~] = stressTangent3D_MR3(C, prop);

            sigma = (1 / J) * F(1:2, 1:2) * S * F(1:2, 1:2)';

            cauchyStress_all(e, c, :) = voigt(sigma);

            c = c + 1;
        end

    end

    nNodes = size(nodeCoordinates, 1);
    stressSum = zeros(nNodes, 3);
    counter = zeros(nNodes, 1);

    for e = 1:nElm
        nodes = et(e, :);

        for c = 1:size(Q_wts, 1)
            % Assign the same stress to all element nodes (lumped projection)
            stress_gp = squeeze(cauchyStress_all(e, c, :));

            for ni = 1:length(nodes)
                nID = nodes(ni);
                stressSum(nID, :) = stressSum(nID, :) + stress_gp';
                counter(nID) = counter(nID) + 1;
            end

        end

    end

    cauchyStress_xx = stressSum(:, 1) ./ counter;
    cauchyStress_yy = stressSum(:, 2) ./ counter;
    cauchyStress_xy = stressSum(:, 3) ./ counter;
end
