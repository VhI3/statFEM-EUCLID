function dIdF = computeStrainInvariantDerivatives(F, i)
% Extract deformation gradient components
F11 = F(:, 1);
F12 = F(:, 2);
F21 = F(:, 3);
F22 = F(:, 4);

% Initialize dIdF as a zero matrix
dIdF = zeros(size(F));

if i == 1
    % dI1/dF:
    dIdF = 2.0 * F;
    
elseif i == 2
    % dI2/dF:
    dIdF11 = 2.0 * F11 - 2.0 .* F12 .* F21 .* F22 + 2.0 .* F11 .* (F22 .^ 2);
    dIdF12 = 2.0 .* F12 + 2.0 .* F12 .* (F21 .^ 2) - 2.0 .* F11 .* F21 .* F22;
    dIdF21 = 2.0 .* F21 + 2.0 .* (F12 .^ 2) .* F21 - 2.0 .* F11 .* F12 .* F22;
    dIdF22 = 2.0 .* F22 - 2.0 .* F11 .* F12 .* F21 + 2.0 .* (F11 .^ 2) .* F22;
    
    % Concatenating results to match the shape of F
    dIdF = [dIdF11, dIdF12, dIdF21, dIdF22];
    
elseif i == 3
    % dI3/dF:
    J = F11 .* F22 - F12 .* F21; % Compute determinant (Jacobian)
    dIdF11 = 2.0 .* F22 .* J;
    dIdF12 = -2.0 .* F21 .* J;
    dIdF21 = -2.0 .* F12 .* J;
    dIdF22 = 2.0 .* F11 .* J;
    
    % Concatenating results to match the shape of F
    dIdF = [dIdF11, dIdF12, dIdF21, dIdF22];
    
else
    error('Incorrect invariant index');
end

end