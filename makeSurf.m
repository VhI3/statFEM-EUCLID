function [vectSurf] = makeSurf(elementNodes, vect)
% makeSurf creates a surface vector from element nodes and a vector.
%
% Inputs:
%   elementNodes: A matrix where each row defines nodes for an element.
%   vect: A vector that will be associated with the surface elements.
%
% Outputs:
%   vectSurf: A reshaped vector representing surface values corresponding to elementNodes.

% Check if vect is a column vector, transpose if not.
if size(vect, 2) ~= 1
    vect = vect';
end

% Get the number of elements.
e = 1:size(elementNodes, 1);
numNodes = size(elementNodes, 2);

% Extract node indices from elementNodes for each element.
indice2 = elementNodes(e, :);

% Map the values in vect to the corresponding surface nodes.
vectSurf = vect(indice2', 1);

% Reshape vectSurf to have numNodes rows corresponding to nodes of each element.
% The number of columns will be equal to the number of elements.
vectSurf = reshape(vectSurf, [numNodes, size(elementNodes, 1)]);

end
