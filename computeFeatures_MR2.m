

function Q = computeFeatures_MR2(I1, I2, I3)
% Compute feature values
J1 = I1 .* (I3 .^ (-1/3));
J2 = (I1 + I3 - 1) .* (I3 .^ (-2/3));
J3 = sqrt(I3);

Q(:, 1) = J1 - 3.0;
Q(:, 2) = J2 - 3.0;
Q(:, 3) = (J1 - 3.0) .^ 2;
Q(:, 4) = (J1 - 3.0) .* (J2 - 3.0);
Q(:, 5) = (J2 - 3.0) .^ 2;
Q(:, 6) = (J3 - 1.0) .^ 2;

end