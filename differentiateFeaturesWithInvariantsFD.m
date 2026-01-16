
function [d_feature_dI, numFeatures] = differentiateFeaturesWithInvariantsFD(featuresFunc, I, h)
% featuresFunc: function handle @(I) computeFeatures(I, ...)
% I: input vector of invariants (n x 1)
% h: finite difference step size (e.g. 1e-6)

if nargin < 3
    h = 1e-6; % default step size
end

% Number of samples and features
n = size(I, 1);
features = featuresFunc(I); % compute base features
numFeatures = size(features, 2);

% Initialize output
d_feature_dI = zeros(n, numFeatures);

% Loop over samples and apply perturbation
for i = 1:n
    I_perturbed = I;
    I_perturbed(i) = I_perturbed(i) + h;
    features_perturbed = featuresFunc(I_perturbed);
    d_feature_dI(i, :) = (features_perturbed(i, :) - features(i, :)) / h;
end

end