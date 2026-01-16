%% Gauss_int_tri.m
% -------------------------------------------------------------------------
% Function: Gauss_int_tri
%
% Purpose:
%   Provides quadrature points and weights for triangular and tetrahedral
%   elements for numerical integration.
%
% Inputs:
%   sdim      - Space dimension (2 for triangles, 3 for tetrahedra)
%   quadorder - Order of quadrature (1, 2, 3 for tetrahedra; 1, 2, 5, 7 for triangles)
%
% Outputs:
%   Q - Quadrature points (each row is a point)
%   W - Quadrature weights
%
% Author: Vahab Narouie, TU-Braunschweig, 2025
% License: GNU GPL v3.0
% -------------------------------------------------------------------------

function [Q, W] = Gauss_int_tri(sdim, quadorder)

% ---------------------- TETRAHEDRAL ELEMENTS (sdim = 3) -------------------
if (sdim == 3)
    % Ensure valid quadrature order
    if ~ismember(quadorder, [1, 2, 3])
        warning('Incorrect quadrature order for tetrahedral quadrature, defaulting to 1.');
        quadorder = 1;
    end
    
    % Quadrature points and weights for tetrahedra
    switch quadorder
        case 1
            Q = [0.25, 0.25, 0.25]; 
            W = 1/6;
            
        case 2
            Q = [0.58541020, 0.13819660, 0.13819660;
                 0.13819660, 0.58541020, 0.13819660;
                 0.13819660, 0.13819660, 0.58541020;
                 0.13819660, 0.13819660, 0.13819660];
            W = [1; 1; 1; 1] / 24;
            
        case 3
            Q = [0.25, 0.25, 0.25;
                 1/2, 1/6, 1/6;
                 1/6, 1/2, 1/6;
                 1/6, 1/6, 1/2;
                 1/6, 1/6, 1/6];
            W = [-4/5, 9/20, 9/20, 9/20, 9/20] / 6;
    end

% ---------------------- TRIANGULAR ELEMENTS (sdim = 2) -------------------
else
    % Ensure valid quadrature order
    if quadorder > 7
        warning('Quadrature order too high for triangular quadrature, defaulting to 1.');
        quadorder = 1;
    end
    
    % Quadrature points and weights for triangles
    switch quadorder
        case 1
            Q = [1/3, 1/3];
            W = 1/2;
            
        case 2
            Q = [1/6, 1/6;
                 2/3, 1/6;
                 1/6, 2/3];
            W = [1/3; 1/3; 1/3] / 2;
            
        case 5
            Q = [0.1012865073235, 0.1012865073235;
                 0.7974269853531, 0.1012865073235;
                 0.1012865073235, 0.7974269853531;
                 0.4701420641051, 0.0597158717898;
                 0.4701420641051, 0.4701420641051;
                 0.0597158717898, 0.4701420641051;
                 1/3, 1/3];
            W = [0.1259391805448;
                 0.1259391805448;
                 0.1259391805448;
                 0.1323941527885;
                 0.1323941527885;
                 0.1323941527885;
                 0.2250000000000] / 2;
            
        case 7
            Q = [0.0651301029022, 0.0651301029022;
                 0.8697397941956, 0.0651301029022;
                 0.0651301029022, 0.8697397941956;
                 0.3128654960049, 0.0486903154253;
                 0.6384441885698, 0.3128654960049;
                 0.0486903154253, 0.6384441885698;
                 0.6384441885698, 0.0486903154253;
                 0.3128654960049, 0.6384441885698;
                 0.0486903154253, 0.3128654960049;
                 0.2603459660790, 0.2603459660790;
                 0.4793080678419, 0.2603459660790;
                 0.2603459660790, 0.4793080678419;
                 1/3, 1/3];
            W = [0.0533472356088;
                 0.0533472356088;
                 0.0533472356088;
                 0.0771137608903;
                 0.0771137608903;
                 0.0771137608903;
                 0.0771137608903;
                 0.0771137608903;
                 0.0771137608903;
                 0.1756152576332;
                 0.1756152576332;
                 0.1756152576332;
                 -0.1495700444677] / 2;
    end
end

end
