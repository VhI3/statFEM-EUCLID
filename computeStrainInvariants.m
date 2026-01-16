function [I1, I2, I3] = computeStrainInvariants(C)
% Extract components of the Right Cauchy-Green strain tensor
C11 = C(:, 1);
C12 = C(:, 2);
C21 = C(:, 3);
C22 = C(:, 4);

% Compute strain invariants
I1 = C11 + C22 + 1.0;
% I2 = C11 + C22 - C12 .* C21 + C11 .* C22;
I2 = C11 .* C22 - C12 .* C21 + C11 + C22; % or simply: I3 + C11 + C22
I3 = C11 .* C22 - C12 .* C21;
end