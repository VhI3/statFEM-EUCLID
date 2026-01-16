function [S, D_voigt, W, J1, J2, J3] = stressTangent3D_MR3(C, prop)
% C: 3x3 right Cauchy-Green tensor
% K: bulk modulus

A10 = prop(1);
A01 = prop(2);
A20 = prop(3);
A11 = prop(4);
A02 = prop(5);
A30 = prop(6);
A21 = prop(7);
A12 = prop(8);
A03 = prop(9);
K   = prop(10);

I = eye(2);

I1 = C(1, 1) + C(2, 2) + 1.0;
I2 = C(1, 1) .* C(2, 2) - C(1, 2) .* C(2, 1) + C(1, 1) + C(2, 2); % or simply: I3 + C11 + C22
I3 = C(1, 1) .* C(2, 2) - C(1, 2) .* C(2, 1);
J = sqrt(I3);
invC = inv(C);

% First derivatives
I1_E = 2 * I;
I2_E = 2 * (I1 * I - C);
I3_E = 2 * I3 * invC;

% Reduced invariants
J1 = I1 * I3 ^ (-1/3);
J2 = I2 * I3 ^ (-2/3);
J3 = J;

% First derivatives of reduced invariants
J1_E = I1_E * I3 ^ (-1/3) - (1/3) * I1 * I3 ^ (-4/3) * I3_E;
J2_E = I2_E * I3 ^ (-2/3) - (2/3) * I2 * I3 ^ (-5/3) * I3_E;
J3_E = 0.5 * I3 ^ (-0.5) * I3_E;

% Stress: Second Piola-Kirchhoff (3x3)
S = A10 * J1_E + A01 * J2_E + 2 * A20 * (J1 - 3) * J1_E + A11 * (J2 - 3) * J1_E ...
    + A11 * (J1 - 3) * J2_E + 2 * A02 * (J2 - 3) * J2_E + 3 * A30 * (J1 - 3) ^ 2 * J1_E ...
    + A21 * (J1 - 3) ^ 2 * J2_E + 2 * A21 * (J1 - 3) * (J2 - 3) * J1_E ...
    + 2 * A12 * (J1 - 3) * (J2 - 3) * J2_E + A12 * (J2 - 3) ^ 2 * J1_E ...
    + 3 * A03 * (J2 - 3) ^ 2 * J2_E + K * (J3 - 1) * J3_E;

% --- Tangent Modulus ---

% Second derivative terms
invCdyad = dyad2D(invC, invC); % invC ⊗ invC
I2_EE = 4 * dyad2D(I, I) - symFourthOrder2D();
I3_EE = 4 * I3 * (invCdyad - dyadContract2D(invC, invC));

outer_J1J3 = outer2D(J1_E, I3_E);
outer_J2J3 = outer2D(J2_E, I3_E);
outer_J3J3 = outer2D(J3_E, J3_E);

J1_EE =- (1/3) * I1 * I3 ^ (-4/3) * tensor4ToVoigt2D(I3_EE) ...
    - (1/3) * I3 ^ (-4/3) * tensor4ToVoigt2D(outer_J1J3 + permute(outer_J1J3, [3 4 1 2])) ...
    + (4/9) * I1 * I3 ^ (-7/3) * tensor4ToVoigt2D(outer_J3J3);

J2_EE = I3 ^ (-2/3) * tensor4ToVoigt2D(I2_EE) ...
    - (2/3) * I2 * I3 ^ (-5/3) * tensor4ToVoigt2D(I3_EE) ...
    - (2/3) * I3 ^ (-5/3) * tensor4ToVoigt2D(outer_J2J3 + permute(outer_J2J3, [3 4 1 2])) ...
    + (10/9) * I2 * I3 ^ (-8/3) * tensor4ToVoigt2D(outer_J3J3);

J3_EE =- (1/4) * I3 ^ (-3/2) * tensor4ToVoigt2D(outer_J3J3) ...
    + 0.5 * I3 ^ (-1/2) * tensor4ToVoigt2D(I3_EE);

% Total tangent in Voigt notation (6x6)


% Tangent tensor
D_voigt = ...
    A10 * J1_EE + A01 * J2_EE ...
    + 2 * A20 * tensor4ToVoigt2D(outer2D(J1_E, J1_E)) ...
    + 2 * A20 * (J1 - 3) * J1_EE ...
    + A11 * tensor4ToVoigt2D(outer2D(J2_E, J1_E) + outer2D(J1_E, J2_E)) ...
    + A11 * (J2 - 3) * J1_EE + A11 * (J1 - 3) * J2_EE ...
    + 2 * A02 * tensor4ToVoigt2D(outer2D(J2_E, J2_E)) ...
    + 2 * A02 * (J2 - 3) * J2_EE ...
    + 6 * A30 * (J1 - 3) * tensor4ToVoigt2D(outer2D(J1_E, J1_E)) + 3 * A30 * (J1 - 3) ^ 2 * J1_EE ...
    + A21 * (J1 - 3) ^ 2 * J2_EE + 2 * A21 * (J1 - 3) * (J2 - 3) * J1_EE ...
    + 2 * A21 * (J1 - 3) * tensor4ToVoigt2D(outer2D(J1_E, J2_E)) + 2 * A21 * (J2 - 3) * tensor4ToVoigt2D(outer2D(J1_E, J1_E)) + 2 * A21 * (J1 - 3) * tensor4ToVoigt2D(outer2D(J2_E, J1_E)) ...
    + 2 * A12 * (J2 - 3) * tensor4ToVoigt2D(outer2D(J1_E, J2_E)) + 2 * A12 * (J1 - 3) * tensor4ToVoigt2D(outer2D(J2_E, J2_E)) ...
    + A12 * (J2 - 3) ^ 2 * J1_EE + 2 * A12 * (J1 - 3) * (J2 - 3) * J2_EE ...
    + 2 * A12 * (J2 - 3) * tensor4ToVoigt2D(outer2D(J2_E, J1_E)) ...
    + 6 * A03 * (J2 - 3) * tensor4ToVoigt2D(outer2D(J2_E, J2_E)) + 3 * A03 * (J2 - 3) ^ 2 * J2_EE ...
    + K * (J3 - 1) * J3_EE + K * tensor4ToVoigt2D(outer_J3J3);

% Strain energy
W = A10 * (J1 - 3) + A01 * (J2 - 3) + A20 * (J1 - 3) ^ 2 + A11 * (J1 - 3) * (J2 - 3) + A02 * (J2 - 3) ^ 2 ...
    + A30 * (J1 - 3) ^ 3 + A21 * (J1 - 3) ^ 2 * (J2 - 3) + A12 * (J1 - 3) * (J2 - 3) ^ 2 + A03 * (J2 - 3) ^ 3 ...
    + 0.5 * K * (J3 - 1) ^ 2;

end

function T = dyad2D(A, B)
T = zeros(2, 2, 2, 2);

for i = 1:2, for j = 1:2, for k = 1:2, for l = 1:2
                T(i, j, k, l) = A(i, j) * B(k, l);
end, end, end, end
end

function T = dyadContract2D(A, B)
T = zeros(2, 2, 2, 2);

for i = 1:2, for j = 1:2, for k = 1:2, for l = 1:2
                T(i, j, k, l) = A(i, k) * B(j, l);
end, end, end, end
end

function T = outer2D(A, B)
T = zeros(2, 2, 2, 2);

for i = 1:2, for j = 1:2, for k = 1:2, for l = 1:2
                T(i, j, k, l) = A(i, j) * B(k, l);
end, end, end, end
end

function V = tensor4ToVoigt2D(D)
% Converts 4th-order tensor to 3x3 Voigt notation (2D symmetric)
V = zeros(3, 3);

voigtPairs = [1 1; 2 2; 1 2];

for I = 1:3
    i = voigtPairs(I, 1);
    j = voigtPairs(I, 2);

    for J = 1:3
        k = voigtPairs(J, 1);
        l = voigtPairs(J, 2);
        V(I, J) = D(i, j, k, l);
    end

end

end

function I4 = symFourthOrder2D()
I4 = zeros(2, 2, 2, 2);

for i = 1:2

    for j = 1:2

        for k = 1:2

            for l = 1:2
                I4(i, j, k, l) = 0.5 * (delta(i, k) * delta(j, l) + delta(i, l) * delta(j, k));
            end

        end

    end

end

end

function val = delta(i, j)
val = double(i == j);
end
