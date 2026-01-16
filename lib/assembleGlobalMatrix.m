function matrix_global = assembleGlobalMatrix(matrix_global, matrix_element, connectivity, ele, numNodesPerElement)
% Assembles the element-level matrix into the global matrix

for a = 1:numNodesPerElement
    node = connectivity{a}(ele); % node index from connectivity
    global_row1 = 2 * node - 1; % MATLAB index for DOF x
    global_row2 = 2 * node; % MATLAB index for DOF y
    
    matrix_global(global_row1, :) = matrix_global(global_row1, :) + matrix_element(2 * a - 1, :);
    matrix_global(global_row2, :) = matrix_global(global_row2, :) + matrix_element(2 * a, :);
end

end