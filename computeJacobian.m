function J = computeJacobian(F)
% Extract deformation gradient components
F11 = F(:, 1);
F12 = F(:, 2);
F21 = F(:, 3);
F22 = F(:, 4);

% Compute the determinant of F (Jacobian)
J = F11 .* F22 - F12 .* F21;
end
