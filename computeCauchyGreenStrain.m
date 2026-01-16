function C = computeCauchyGreenStrain(F)
% Extract deformation gradient components
F11 = F(:, 1);
F12 = F(:, 2);
F21 = F(:, 3);
F22 = F(:, 4);

% Compute Cauchy-Green strain tensor components
C11 = F11 .^ 2 + F21 .^ 2;
C12 = F11 .* F12 + F21 .* F22;
C21 = F12 .* F11 + F22 .* F21;
C22 = F12 .^ 2 + F22 .^ 2;

% Concatenate results along columns (similar to PyTorch torch.cat)
C = [C11, C12, C21, C22];
end