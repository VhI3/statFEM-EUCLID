function BVP = preprocess(BVP)
% disp('1. Pre pocessing the plate with hole ...')
% %% Assign from BVP
% quarterPlate_L100_esf100
% quarterPlate_L100_esf90
% quarterPlate_L100_esf50
% quarterPlate_L100_esf25
% quarterPlate_L100_esf10
quarterPlate_L100_esf7
% quarterPlate_L100_esf5
% quarterPlate_L100_esf3
% quarterPlate_L100_esf2

% quarterPlate_L100_adapt_esf20

%% Assig from BVP
% Node coordinates of FEM mesh
nodeCoordinates = msh.POS(:, 1:end - 1);
% Element connectivity of FEM mesh
elementNodes = msh.TRIANGLES(:, 1:end - 1);
% Number of nodes per element
nElNode = 3;
% Number of elements
nuElm = size(elementNodes, 1);
% Define orders for Jacobian calculation and for 2D shape functions
order_jacobi = [1 3 2];
% Dimension of the problem
DIM = 2;
% Number of Degree of Freedom per Nodes
DOFs = 2;
% Number of DOFs per element
eDOFs = DOFs * nElNode;
% Number of nodes
nuNodes = size(nodeCoordinates, 1);
% Number of Global DOFs
GDOFs = nuNodes * DOFs;
% node coordinate in vector format
nodeCoordinates_vector = zeros(GDOFs, 1);
nodeCoordinates_vector(1:2:end) = nodeCoordinates(:, 1);
nodeCoordinates_vector(2:2:end) = nodeCoordinates(:, 2);
% Generate every element DOFs
elementDOFs = zeros(nuElm, eDOFs);
elementDOFs(:, 1:2:eDOFs) = 2 * elementNodes - 1;
elementDOFs(:, 2:2:eDOFs) = 2 * elementNodes;
%
% Length of the plate
L = max(nodeCoordinates(:, 1));
% Height of the plate
H = max(nodeCoordinates(:, 2));
T = 1;
leftNodes = find(abs(nodeCoordinates(:, 1)) < 1e-13);
rightNodes = find(abs(nodeCoordinates(:, 1) - L) < 1e-13);
topNodes = find(abs(nodeCoordinates(:, 2) - H) < 1e-13);
bottomNodes = find(abs(nodeCoordinates(:, 2)) < 1e-13);

bcs_nodes = zeros(nuNodes, DOFs);
bcs_nodes(leftNodes, 1) = 1;
bcs_nodes(rightNodes, 1) = 2;
bcs_nodes(topNodes, 2) = 4;
bcs_nodes(bottomNodes, 2) = 3;

% Define orders for Jacobian calculation and for 2D shape functions
NGT = 1;
[Q_tri, W_tri] = Gauss_int_tri(DIM, NGT);

% Initialization of shape functions and their derivatives
Nshape = @(xi, eta) [1 - xi - eta, xi, eta];
dN_dxi = @(xi, eta) [-1, 1, 0];
dN_deta = @(xi, eta) [-1, 0, 1];

% Combine the natural derivatives into a single matrix
naturalDerivatives = [dN_dxi(Q_tri(1), Q_tri(2)); dN_deta(Q_tri(1), Q_tri(2))]';

% Initialize a counter for elements with negative Jacobian determinants
count = 0;
% Loop through each element in the mesh
for e = 1:size(elementNodes, 1)
    % Extract node indices for the current element
    indice = elementNodes(e, :);
    
    % Get the coordinates of the nodes for the current element
    coords = nodeCoordinates(elementNodes(e, :), :);
    
    % Calculate the Jacobian matrix using the node coordinates and natural derivatives
    [JacobianMatrix, ~, ~] = Jacobian2D(coords, naturalDerivatives);
    
    % Calculate the determinant of the Jacobian matrix
    detJ = det(JacobianMatrix);
    
    % If the determinant is negative, adjust the node ordering
    if detJ < 0
        count = count + 1; % Increment the counter
        elementNodes(e, :) = indice(order_jacobi); % Reorder the nodes
    end
    
end

top_right_node = find(nodeCoordinates(:, 1) == max(nodeCoordinates(:, 1)) & nodeCoordinates(:, 2) == max(nodeCoordinates(:, 2)), 1);
bottom_left_node = find(nodeCoordinates(:, 1) == min(nodeCoordinates(:, 1)) & nodeCoordinates(:, 2) - 10 < 1e-13, 1);
top_left_node = find(nodeCoordinates(:, 1) == min(nodeCoordinates(:, 1)) & nodeCoordinates(:, 2) == max(nodeCoordinates(:, 2)), 1);
bottom_right_node = find(nodeCoordinates(:, 1) == max(nodeCoordinates(:, 1)) & nodeCoordinates(:, 2) - L < 1e-13, 1);

% Print the number of elements with negative Jacobian determinants
disp(['Number of elements with negative Jacobian determinants: ', num2str(count)]);
% Create a subplot which plots the mesh, and sensor locations
xSurf = makeSurf(elementNodes, nodeCoordinates(:, 1));
ySurf = makeSurf(elementNodes, nodeCoordinates(:, 2));

% finding element index of the top-right element which contains the top-right node
% This is done by checking which element contains the top-right node
top_right_element = find(any(elementNodes == top_right_node, 2), 1);
bottom_left_element = find(any(elementNodes == bottom_left_node, 2), 1);
top_left_element = find(any(elementNodes == top_left_node, 2), 1);
bottom_right_element = find(any(elementNodes == bottom_right_node, 2), 1);

% Create a subplot which plots the mesh for top-right element.
xSurf_top_right_element = makeSurf(elementNodes(top_right_element, :), nodeCoordinates(:, 1));
ySurf_top_right_element = makeSurf(elementNodes(top_right_element, :), nodeCoordinates(:, 2));

% Create a subplot which plots the mesh for bottom-left element.
xSurf_bottom_left_element = makeSurf(elementNodes(bottom_left_element, :), nodeCoordinates(:, 1));
ySurf_bottom_left_element = makeSurf(elementNodes(bottom_left_element, :), nodeCoordinates(:, 2));

% Create a subplot which plots the mesh for top-left element.
xSurf_top_left_element = makeSurf(elementNodes(top_left_element, :), nodeCoordinates(:, 1));
ySurf_top_left_element = makeSurf(elementNodes(top_left_element, :), nodeCoordinates(:, 2));

% Create a subplot which plots the mesh for bottom-right element.
xSurf_bottom_right_element = makeSurf(elementNodes(bottom_right_element, :), nodeCoordinates(:, 1));
ySurf_bottom_right_element = makeSurf(elementNodes(bottom_right_element, :), nodeCoordinates(:, 2));

% -------------------- Boundary conditions -----------------------------
% Dirichlet
% Initialize the displacement vector
u = zeros(GDOFs, 1);
% The prescribed displacements are set to zero on the left in the x-direction and on the bottom in the y-direction
fixedPrescribedDOFsX = DOFs * leftNodes - 1;
% fixedPrescribedDOFsY = sort(DOFs * [topNodes; bottomNodes]);
fixedPrescribedDOFsY = DOFs * bottomNodes;

fixedPrescribedDOFs = sort([fixedPrescribedDOFsX; fixedPrescribedDOFsY]);
u(fixedPrescribedDOFs) = 0;
% The prescribed are not set to zero on the right in the x-direction and on the top in the y-direction
nonFixedPrescribedDOFsX = [];
nonFixedPrescribedDOFsY = [];
% nonFixedPrescribedDOFsX = DOFs*rightNodes-1;
% nonFixedPrescribedDOFsY = DOFs*topNodes;
nonFixedPrescribedDOFs = sort([nonFixedPrescribedDOFsX; nonFixedPrescribedDOFsY]);

% -------------------- Load -----------------------------
% This is the load combination for the plate with hole
loadCase = 'UT'; % Uniaxial tension (UT), Biaxial tension (BT), Uniaxial compression (UC), Biaxial compression (BC)
% loadCase = 'BT'; % Biaxial tension (BT)
% loadCase = 'UC'; % Uniaxial compression (UC)
% loadCase = 'BC'; % Biaxial compression (BC)
switch loadCase
    case 'UT' % Uniaxial tension (UT)
        lstp_x = 1; lstp_y = 0;
        tractionX = lstp_x * 0.5e6;
        tractionY = lstp_y * 0.5e6;
    case 'BT' % Biaxial tension (BT)
        lstp_x = 1; lstp_y = 1;
        tractionX = lstp_x * 0.5e6;
        tractionY = lstp_y * 0.5e6;
    case 'UC' % Uniaxial compression (UC)
        lstp_x = -1; lstp_y = 0;
        tractionX = lstp_x * 0.5e6;
        tractionY = lstp_y * 0.5e6;
    case 'BC' % Biaxial compression (BC)
        lstp_x = -1; lstp_y = -1;
        tractionX = lstp_x * 0.5e6;
        tractionY = lstp_y * 0.5e6;
    otherwise
        error('Load case not defined')
end

mean_tractionX = tractionX;
mean_tractionY = tractionY;
std_tractionX = 0.05 * mean_tractionX;
std_tractionY = 0.05 * mean_tractionY;

% The uncertainty in prior is orginating from the uncertain force, which is Gaussian distributed with defined mean and standard deviation
% Sample Gaussain of force with mean and standard deviation
nSampled_force = 5;
nMC = 1000;
% Generate the standard normal distribution
rng(3);
xi = randn(nSampled_force, 1);
xi_MC = randn(nMC, 1);
% Generate the force based in xi
sampled_tractionX = mean_tractionX + std_tractionX * xi;
sampled_tractionY = mean_tractionY + std_tractionY * xi;

P_PC = 1;
[~, Psi_s, ~, PsiSqNorm, P_u] = Hermite_PC(1, P_PC); % Hermite basis

% Sort rightNodes by their Y coordinate (2nd column of nodeCoordinates)
[~, sortIdx] = sort(nodeCoordinates(rightNodes, 2), 'ascend');
rightNodesCorr = rightNodes(sortIdx);
rightEdge_y = nodeCoordinates(rightNodesCorr, 2); % Y-coordinates
rightEdge_y_pairs = [rightEdge_y(1:end - 1) rightEdge_y(2:end)];

rightEdge_dofX = DOFs * rightNodesCorr - 1; % X DOFs only
rightEdge_dofX_pairs = [rightEdge_dofX(1:end - 1) rightEdge_dofX(2:end)];

sampled_force = zeros(GDOFs, nSampled_force);

for i = 1:nSampled_force
    
    for ey = 1:size(rightEdge_y_pairs, 1)
        ed = rightEdge_dofX_pairs(ey, :);
        dy = abs(rightEdge_y_pairs(ey, 2) - rightEdge_y_pairs(ey, 1));
        sampled_force(ed, i) = sampled_force(ed, i) + dy / 2 * sampled_tractionX(i) * ones(2, 1);
    end
    
end

force = zeros(GDOFs, 1);

for ey = 1:size(rightEdge_y_pairs, 1)
    ed = rightEdge_dofX_pairs(ey, :);
    dy = abs(rightEdge_y_pairs(ey, 2) - rightEdge_y_pairs(ey, 1));
    force(ed, 1) = force(ed, 1) + dy / 2 * tractionX * ones(2, 1);
end

% Top nodes
[~, sortIdy] = sort(nodeCoordinates(topNodes, 1), 'ascend');
topNodesCorr = topNodes(sortIdy);
topEdge_x = nodeCoordinates(topNodesCorr, 1); % X-coordinates
topEdge_x_pairs = [topEdge_x(1:end - 1) topEdge_x(2:end)];
%
topEdge_dofY = DOFs * topNodesCorr; % Y DOFs only
topEdge_dofY_pairs = [topEdge_dofY(1:end - 1) topEdge_dofY(2:end)];

for i = 1:nSampled_force
    
    for ex = 1:size(topEdge_x_pairs, 1)
        ed = topEdge_dofY_pairs(ex, :);
        dx = abs(topEdge_x_pairs(ex, 2) - topEdge_x_pairs(ex, 1));
        sampled_force(ed, i) = sampled_force(ed, i) + dx / 2 * sampled_tractionY(i) * ones(2, 1);
    end
    
end

% Deterministic force
for ex = 1:size(topEdge_x_pairs, 1)
    ed = topEdge_dofY_pairs(ex, :);
    dx = abs(topEdge_x_pairs(ex, 2) - topEdge_x_pairs(ex, 1));
    force(ed, 1) = force(ed, 1) + dx / 2 * tractionY * ones(2, 1);
end

tbc = [];

switch loadCase
    case 'UT' % Uniaxial tension (UT)
        tbc = [tbc rightNodes]; % the number of neumann boundary condition e.g. forces.
    case 'BT' % Biaxial tension (BT)
        tbc = [tbc rightNodes topNodes]; % the number of neumann boundary condition e.g. forces.
    case 'UC' % Uniaxial compression (UC)
        tbc = [tbc rightNodes]; % the number of neumann boundary condition e.g. forces.
    case 'BC' % Biaxial compression (BC)
        tbc = [tbc rightNodes topNodes]; % the number of neumann boundary condition e.g. forces.
    otherwise
        error('Load case not defined')
end

% An indicator if it is forced control
ntbc = size(tbc, 1);
%
prescribedDOFsX = sort([fixedPrescribedDOFsX; nonFixedPrescribedDOFsX]);
prescribedDOFsY = sort([fixedPrescribedDOFsY; nonFixedPrescribedDOFsY]);
prescribedDOFs = sort([fixedPrescribedDOFs; nonFixedPrescribedDOFs]);
ndbc = size(prescribedDOFs, 1);
dbc = u(prescribedDOFs);
activeDOFsX = setdiff(1:2:GDOFs, prescribedDOFsX);
activeDOFsY = setdiff(2:2:GDOFs, prescribedDOFsY);
activeDOFs = setdiff(1:GDOFs, prescribedDOFs);

% -------------------- 3. Material properties -----------------------------
% for the purpose of hyperelastic materila, I use the following values
% A10 = 0.5e6;
% A01 = 0;

A10 = 0.3e6;
A01 = 0.2e6;

A20 = 0;
A11 = 0;
A02 = 0;
A30 = 0;
A21 = 0;
A12 = 0;
A03 = 0;
K_bulk = 6*(A10+ A01+ A20+ A11+ A02+ A30+ A21+ A12+ A03);

% Define W_true function handle
W_true = @(J1, J2, J3) ...
    A10 * (J1 - 3) + ...
    A01 * (J2 - 3) + ...
    A20 * (J1 - 3) .^ 2 + ...
    A11 * (J1 - 3) .* (J2 - 3) + ...
    A02 * (J2 - 3) .^ 2 + ...
    A30 * (J1 - 3) .^ 3 + ...
    A21 * (J1 - 3) .^ 2 .* (J2 - 3) + ...
    A12 * (J1 - 3) .* (J2 - 3) .^ 2 + ...
    A03 * (J2 - 3) .^ 3 + ...
    0.5 * K_bulk * (J3 - 1) .^ 2;

prop0 = [
    A10; A01; A20; A11; A02; A30; A21; A12; A03; K_bulk
    ];

mu = 2 * (A10 + A01);
E = (9 * K_bulk * mu) / (3 * K_bulk + mu);
nu = (3 * K_bulk - 2 * mu) / (2 * (3 * K_bulk + mu));

E_linear = E / 2;


%%
% 1. Define symbolic variables
syms J1 J2 J3 real
A_mn = prop0(:) / 10 ^ 6 ;
A_mn(end) = A_mn(end) / 2;

% 2. Define symbolic features explicitly (same order as in your code)
Q{1} = (J1 - 3);
Q{2} = (J2 - 3);
Q{3} = (J1 - 3) ^ 2;
Q{4} = (J1 - 3) * (J2 - 3);
Q{5} = (J2 - 3) ^ 2;
Q{6} = (J1 - 3) ^ 3;
Q{7} = (J1 - 3) ^ 2 * (J2 - 3);
Q{8} = (J1 - 3) * (J2 - 3) ^ 2;
Q{9} = (J2 - 3) ^ 3;
Q{10} = (J3 - 1) ^ 2;
%
latex_terms = strings(1, 10); % initialize array of strings

for i = 1:10
    
    if abs(A_mn(i)) > 1e-12
        coeff_str = sprintf('%.3f', A_mn(i));
        
        if i == 1 || i == 2
            term_str = coeff_str + "\,(" + latex(Q{i}) + ")";
        else
            term_str = coeff_str + "\," + latex(Q{i});
        end
        
        latex_terms(i) = term_str;
    end
    
end

% Concatenate terms with ' + '
latex_expression = strjoin(latex_terms(latex_terms ~= ""), ' + ');

% Wrap in equation
w_latex_true = latex_expression;

%%
Constitutive_matrix = 'plane stress'; % 'plane stress' or 'plane strain'

switch Constitutive_matrix
    case 'plane strain'
        C_Constitutive = ... % Constitutive matrix for plain strain
            1 / (1 - nu ^ 2) * [1 nu 0;
            nu 1 0;
            0 0 (1 - nu) / 2];
    case 'plane stress'
        C_Constitutive = ... % Constitutive matrix for plain stress
            1 / (1 + nu) / (1 - 2 * nu) * [1 - nu nu 0;
            nu 1 - nu 0;
            0 0 1/2 - nu];
    otherwise
        error('Constitutive matrix not defined')
end

tol = 1e-8; % tolerance of convergence
maxit = 50; % maximum iterative steps
maxreit = 6; % maximum load/displacement step reduction
timeInterval = 0.1; % initial time interval

%% Save the identified Hyperparameters and Hyperelastic model parameters for every iteration
identified_hyperparameters{1} = [];
identified_hyperparameters{2} = [];
hyperelastic_matparameters = prop0;
RMSE = [];
epsilon_u_iter = [];
epsilon_kappa_iter = [];
epsilon_I_iter = [];
w_latex_full_iter = [];
relative_rmse_iter = [];
%% Assign back to BVP
BVP.preProc.msh.nodeCoordinates = nodeCoordinates;
BVP.preProc.msh.elementNodes = elementNodes;
BVP.preProc.msh.nElm = nuElm;
BVP.preProc.msh.xSurf = xSurf;
BVP.preProc.msh.ySurf = ySurf;
BVP.preProc.msh.nElNode = nElNode;
BVP.preProc.msh.nuNodes = nuNodes;
BVP.preProc.msh.DIM = DIM;
BVP.preProc.msh.DOFs = DOFs;
BVP.preProc.msh.eDOFs = eDOFs;
BVP.preProc.msh.GDOFs = GDOFs;
BVP.preProc.msh.L = L;
BVP.preProc.msh.H = H;
BVP.preProc.msh.T = T;
BVP.preProc.msh.top_right_node = top_right_node;
BVP.preProc.msh.elementDOFs = elementDOFs;
BVP.preProc.msh.leftNodes = leftNodes;
BVP.preProc.msh.rightNodes = rightNodes;
BVP.preProc.msh.topNodes = topNodes;
BVP.preProc.msh.bottomNodes = bottomNodes;
BVP.preProc.msh.top_right_element = top_right_element;
BVP.preProc.msh.bottom_left_element = bottom_left_element;
BVP.preProc.msh.top_left_element = top_left_element;
BVP.preProc.msh.bottom_right_element = bottom_right_element;
BVP.preProc.msh.xSurf_top_right_element = xSurf_top_right_element;
BVP.preProc.msh.ySurf_top_right_element = ySurf_top_right_element;
BVP.preProc.msh.xSurf_bottom_left_element = xSurf_bottom_left_element;
BVP.preProc.msh.ySurf_bottom_left_element = ySurf_bottom_left_element;
BVP.preProc.msh.xSurf_top_left_element = xSurf_top_left_element;
BVP.preProc.msh.ySurf_top_left_element = ySurf_top_left_element;
BVP.preProc.msh.xSurf_bottom_right_element = xSurf_bottom_right_element;
BVP.preProc.msh.ySurf_bottom_right_element = ySurf_bottom_right_element;
BVP.preProc.msh.nodeCoordinates_vector = nodeCoordinates_vector;
% Gauss quadrature points
BVP.preProc.Q_pts = Q_tri;
BVP.preProc.Q_wts = W_tri;
BVP.preProc.Nshape = Nshape;
BVP.preProc.dN_dxi = dN_dxi;
BVP.preProc.dN_deta = dN_deta;
% Boundary condition
BVP.preProc.BC.fixedPrescribedDOFs = fixedPrescribedDOFs;
BVP.preProc.BC.fixedPrescribedDOFsX = fixedPrescribedDOFsX;
BVP.preProc.BC.fixedPrescribedDOFsY = fixedPrescribedDOFsY;
BVP.preProc.BC.prescribedDOFs = prescribedDOFs;
BVP.preProc.BC.prescribedDOFsX = prescribedDOFsX;
BVP.preProc.BC.prescribedDOFsY = prescribedDOFsY;
BVP.preProc.BC.activeDOFs = activeDOFs;
BVP.preProc.BC.activeDOFsX = activeDOFsX;
BVP.preProc.BC.activeDOFsY = activeDOFsY;
BVP.preProc.BC.dbc = dbc;
BVP.preProc.BC.ntbc = ntbc;
BVP.preProc.BC.ndbc = ndbc;
BVP.preProc.BC.bcs_nodes = bcs_nodes;
BVP.preProc.BC.force = force;
BVP.preProc.BC.tractionX = tractionX;
BVP.preProc.BC.tractionY = tractionY;
% -------------------- 3. Material properties -----------------------------
BVP.preProc.material.E = E;
BVP.preProc.material.E_linear = E_linear;
BVP.preProc.material.C_Constitutive = C_Constitutive;
BVP.preProc.material.prop0 = prop0;
% -------------------- 9. non-intrusive least square PCE ------------------
BVP.preProc.UQ.sampled_force = sampled_force;
BVP.preProc.UQ.mean_tractionX = mean_tractionX;
BVP.preProc.UQ.mean_tractionY = mean_tractionY;
BVP.preProc.UQ.std_tractionX = std_tractionX;
BVP.preProc.UQ.std_tractionY = std_tractionY;
BVP.preProc.UQ.xi = xi;
BVP.preProc.UQ.nMC = nMC;
BVP.preProc.UQ.xi_MC = xi_MC;
BVP.preProc.UQ.nSampled_force = nSampled_force;
BVP.preProc.UQ.P_u = P_u;
BVP.preProc.UQ.Psi_s = Psi_s;
BVP.preProc.UQ.PsiSqNorm = PsiSqNorm;
% -------------------- 10. statFEM -------------------------------
% BVP.preProc.statFEM.sig_e = sig_e;
% -------------------- 11. Newton-Raphson ---------------------
BVP.preProc.tol = tol;
BVP.preProc.maxit = maxit;
BVP.preProc.maxreit = maxreit;
BVP.preProc.timeInterval = timeInterval;
%% --- discovered model parameters and hyperparameters
BVP.discovered.identified_hyperparameters = identified_hyperparameters;
BVP.discovered.hyperelastic_matparameters = hyperelastic_matparameters;
BVP.discovered.RMSE = RMSE;
BVP.discovered.epsilon_u_iter = epsilon_u_iter;
BVP.discovered.epsilon_kappa_iter = epsilon_kappa_iter;
BVP.discovered.epsilon_I_iter = epsilon_I_iter;
BVP.discovered.w_latex_full_iter = w_latex_full_iter;
BVP.discovered.W_true = W_true;
BVP.discovered.relative_rmse_iter = relative_rmse_iter;
BVP.discovered.w_latex_true = w_latex_true;
end
