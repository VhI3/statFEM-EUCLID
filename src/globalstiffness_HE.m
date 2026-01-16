function [stiffness, R, lambda_out, W_out, W_total, J_min_max] = globalstiffness_HE(BVP_tmp, u, prop)
%% Assign from BVP
GDOFs = BVP_tmp.GDOFs;
DOFs = BVP_tmp.DOFs;
nElm = BVP_tmp.nElm;
nElNode = BVP_tmp.nElNode;
nodeCoordinates = BVP_tmp.nodeCoordinates;
et = BVP_tmp.elementNodes;
ed = BVP_tmp.elementDOFs;
Q_pts = BVP_tmp.Q_pts;
Q_wts = BVP_tmp.Q_wts;
dN_dxi = BVP_tmp.dN_dxi;
dN_deta = BVP_tmp.dN_deta;
thickness = BVP_tmp.thickness;
top_right_element = BVP_tmp.top_right_element;
bottom_left_element = BVP_tmp.bottom_left_element;
top_left_element = BVP_tmp.top_left_element;
bottom_right_element = BVP_tmp.bottom_right_element;

%% Initialization
stiffness = zeros(GDOFs, GDOFs);
R = zeros(GDOFs, 1);
W_out = [0; 0; 0; 0]; % Initialize W_out
lambda_out = [1; 1; 1; 1]; % Initialize lambda_out
W_total = 0;
% Store the minimum and maximum values of J1, J2, and J3
J1_min = Inf;
J1_max = -Inf;

J2_min = Inf;
J2_max = -Inf;

J3_min = Inf;
J3_max = -Inf;


for e = 1:nElm
    indice = et(e, :);
    elementDof = ed(e, :);
    elDisp = u(elementDof);
    elDisp = reshape(elDisp, DOFs, nElNode);
    
    ke = zeros(size(elementDof, 2));
    re = zeros(size(elementDof, 2), 1);
    
    for i = 1:size(Q_wts, 1)
        xi = Q_pts(i, 1); eta = Q_pts(i, 2);
        naturalDerivatives = [dN_dxi(xi, eta); dN_deta(xi, eta)]';
        [JacobianMatrix, ~, XYDerivatives] = Jacobian2D(nodeCoordinates(indice, :), naturalDerivatives);
        
        F = eye(3);
        F(1:2, 1:2) = elDisp * XYDerivatives + eye(2);
        
        C11 = F(1, 1) .^ 2 + F(2, 1) .^ 2;
        C12 = F(1, 1) .* F(1, 2) + F(2, 1) .* F(2, 2);
        C21 = F(1, 2) .* F(1, 1) + F(2, 2) .* F(2, 1);
        C22 = F(1, 2) .^ 2 + F(2, 2) .^ 2;
        
        % Concatenate results along columns (similar to PyTorch torch.cat)
        C = [C11, C12; C21, C22];
        
        [S3D, D3D, W_gp, current_J1, current_J2, current_J3] = stressTangent3D_MR3(C, prop);
        W_total = W_total + W_gp;
        
        J1_min = min(J1_min, current_J1);
        J1_max = max(J1_max, current_J1);
        
        J2_min = min(J2_min, current_J2);
        J2_max = max(J2_max, current_J2);
        
        J3_min = min(J3_min, current_J3);
        J3_max = max(J3_max, current_J3);
        
        
        S_voigt = [S3D(1, 1); S3D(2, 2); S3D(1, 2)];
        D_voigt = D3D([1 2 3], [1 2 3]);
        
        dNdx = XYDerivatives(:, 1);
        dNdy = XYDerivatives(:, 2);
        
        BN = zeros(3, size(elementDof, 2));
        BG = zeros(4, size(elementDof, 2));
        
        BN(1, 1:2:end) = F(1, 1) * dNdx;
        BN(1, 2:2:end) = F(2, 1) * dNdx;
        BN(2, 1:2:end) = F(1, 2) * dNdy;
        BN(2, 2:2:end) = F(2, 2) * dNdy;
        BN(3, 1:2:end) = F(1, 1) * dNdy + F(1, 2) * dNdx;
        BN(3, 2:2:end) = F(2, 1) * dNdy + F(2, 2) * dNdx;
        
        BG(1, 1:2:end) = dNdx;
        BG(2, 1:2:end) = dNdy;
        BG(3, 2:2:end) = dNdx;
        BG(4, 2:2:end) = dNdy;
        
        sigma = [S_voigt(1), S_voigt(3); S_voigt(3), S_voigt(2)];
        stan = zeros(4);
        stan(1:2, 1:2) = sigma;
        stan(3:4, 3:4) = sigma;
        
        ke = ke + Q_wts(i) * (BN' * D_voigt * BN + BG' * stan * BG) * det(JacobianMatrix) * thickness;
        re = re + Q_wts(i) * BN' * S_voigt * det(JacobianMatrix) * thickness;
    end
    
    stiffness(elementDof, elementDof) = stiffness(elementDof, elementDof) + ke;
    R(elementDof) = R(elementDof) + re;
end
J_min_max = [J1_min, J1_max; J2_min, J2_max; J3_min, J3_max];
end
