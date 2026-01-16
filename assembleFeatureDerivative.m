
function d_features_dF_element = assembleFeatureDerivative(d_features_dI1, d_features_dI3, dI1dF, dI3dF, ele)
% Extract element-specific rows
d_features_dI1_element = d_features_dI1(ele, :);
d_features_dI3_element = d_features_dI3(ele, :);
dI1dF_element = dI1dF(ele, :);
dI3dF_element = dI3dF(ele, :);

% Compute the outer product (equivalent to np.outer in Python)
d_features_dF_element = (d_features_dI1_element' * dI1dF_element) + ...
    (d_features_dI3_element' * dI3dF_element);
end
