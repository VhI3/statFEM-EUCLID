function BVP = projection_tri(BVP)
    %   PROJECTION_TRI constructs the observation matrix for FEM analysis.
    %   This function creates a mapping matrix P that projects nodes of FE mesh
    %   onto sensor locations. The mapping is done by checking if the sensor
    %   location is inside an element in the mesh. If the sensor is inside an element,
    %   the function computes the local coordinates (ξ, η) and shape functions N(ξ, η).
    %   The function then constructs the observation matrices Px and Py for the x and y
    %   directions, respectively. The function also checks if the sensor is outside the domain.
    %   If the sensor is outside the domain, it marks it as such and does not compute
    %   the local coordinates or shape functions. The are some numerical tests to ensure
    %   the correctness of the projection matrices. 
    %
    % Inputs:
    %   BVP - Boundary value problem structure containing:
    %       Sensor coordinates, node coordinates, DOFs, and active DOFs.
    %
    % Outputs:
    %   BVP - Updated BVP structure with computed projection matrices.
    %
    % Project:
    % Author: Vahab Narouie, TU-Braunschweig, 2025
    % License: GNU GPL v3.0 (see LICENSE file for details)

    %% Assign from BVP
    % Nummber of global DOFs of FEM nodes
    GDOFs = BVP.preProc.msh.GDOFs;
    % FEM element connectivity
    elementNodes = BVP.preProc.msh.elementNodes; % Element connectivity
    % FEM node coordinates in a matrix form
    nodeCoordinates = BVP.preProc.msh.nodeCoordinates; % Node coordinates
    % FEM node coordinates in a vector form
    nodeCoordinates_vector = BVP.preProc.msh.nodeCoordinates_vector;
    % number of nodes
    numberNodes = BVP.preProc.msh.nuNodes;
    % Number of elements
    nElm = BVP.preProc.msh.nElm;
    % Active DOFs of FEM nodes in x-direction
    activeDOFsX = BVP.preProc.BC.activeDOFsX;
    % Active DOFs of FEM nodes in y-direction
    activeDOFsY = BVP.preProc.BC.activeDOFsY;
    % Active DOFs of FEM nodes in x and y-direction together
    activeDOFs = BVP.preProc.BC.activeDOFs;

    % sensor coordinates in x-direction
    x_sensor = BVP.obs.x_sensor;
    % sensor coordinates in y-direction
    y_sensor = BVP.obs.y_sensor;
    % Active Dofs of sensors in x-direction
    x_sensor_active = BVP.obs.x_sensor_active;
    % Active Dofs of sensors in y-direction
    y_sensor_active = BVP.obs.y_sensor_active;
    % Sensor coordinates in a vector
    sensorCoordinates_vector = BVP.obs.sensorCoordinates_vector;
    % Number of sensors
    nsen = BVP.obs.nSen;
    % Active Dofs of sensors in x-direction
    activeSensorDOFsX = BVP.obs.activeSensorDOFsX;
    % Active Dofs of sensors in y-direction
    activeSensorDOFsY = BVP.obs.activeSensorDOFsY;
    % Active Dofs of sensors in x and y-direction together
    activeSensorDOFs = BVP.obs.activeSensorDOFs;
    % Nummber of global DOFs of sensors
    GSensorDOFs = BVP.obs.GSensorDOFs;
    % Sensor coordinates in a vector
    sensorCoordinates_active_vector = BVP.obs.sensorCoordinates_active_vector;

    %% Initialize Observation Matrix and Sensor Nodes
    tolerance = 1e-6; % Define a small numerical tolerance
    Px = zeros(nsen, numberNodes);
    Py = zeros(nsen, numberNodes);
    %% Compute Projection Matrix for Each Sensor
    % Step 1: Loop over all sensors
    for i = 1:nsen
        sensor_x = x_sensor(i);
        sensor_y = y_sensor(i);
        elementFound = 0; % Reset for each sensor

        % Step 2: Try to find the sensor in an element in mesh
        for j = 1:nElm
            % Get element nodes (4-node quadrilateral assumed)
            nodes = elementNodes(j, :);
            x_elem = nodeCoordinates(nodes, 1); % X-coordinates of element nodes
            y_elem = nodeCoordinates(nodes, 2); % Y-coordinates of element nodes

            % Step 3: Check if Sensor is Inside the Element
            if isInsideTriangle(sensor_x, sensor_y, x_elem, y_elem)
                % Step 4: Compute Local Coordinates (ξ, η)
                [xi, eta] = directMappingTri(sensor_x, sensor_y, x_elem, y_elem);

                % If (ξ, η) are slightly outside [-1,1] due to numerical errors, correct them
                if xi < -1 || xi > 1 || eta < -1 || eta > 1
                    % Clamping to boundary
                    xi = max(-1, min(1, xi));
                    eta = max(-1, min(1, eta));

                    % Highlight in red (indicating minor out-of-bounds correction)
                    scatter(sensor_x, sensor_y, ...
                        'SizeData', 50, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', ...
                        'MarkerFaceAlpha', .5, 'MarkerEdgeAlpha', 1);
                end

                % Step 5: Compute Shape Functions at (ξ, η)
                N = [1 - xi - eta, xi, eta];

                Px(i, nodes) = N;
                Py(i, nodes) = N;
                elementFound = 1; % Mark as found
                break; % Exit loop early if element is found
            end

        end

        % Step 6: If sensor is not found inside any element, mark it as "out of domain"
        if elementFound == 0
            fprintf('Sensor at (%.3f, %.3f) is OUTSIDE the domain.\n', sensor_x, sensor_y);

            % Scatter plot in black (indicating out-of-domain sensor)
            scatter(sensor_x, sensor_y, ...
                'SizeData', 10, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', ...
                'MarkerFaceAlpha', .5, 'MarkerEdgeAlpha', 1);
        end

    end

    % Contacate the Px and Py matrices
    P_full = zeros(GSensorDOFs, GDOFs);
    P_full(1:2:GSensorDOFs, 1:2:GDOFs) = Px;
    P_full(2:2:GSensorDOFs, 2:2:GDOFs) = Py;

    % Extract the Px_active and Py_active matrices
    Px_active = P_full(activeSensorDOFsX, activeDOFsX);
    Py_active = P_full(activeSensorDOFsY, activeDOFsY);
    Px_active0 = P_full(1:2:end, activeDOFsX);
    Py_active0 = P_full(2:2:end, activeDOFsY);
    %% Test the projection matrices
    % Test the P_full matrix
    sensorCoordinates_vector_proj = P_full * nodeCoordinates_vector;
    differen_full = sensorCoordinates_vector_proj - sensorCoordinates_vector;

    if norm(differen_full) < tolerance
        disp('Projection matrix P_full is correctly computed!');
    else
        disp('Error detected: Projection matrix P_full may be incorrect.');
    end

    % Test the Px and Py matrices
    x_sensor_proj = Px * nodeCoordinates(:, 1);
    y_sensor_proj = Py * nodeCoordinates(:, 2);
    differen_x = x_sensor_proj - x_sensor;
    differen_y = y_sensor_proj - y_sensor;

    if norm(differen_x) < tolerance && norm(differen_y) < tolerance
        disp('Projection matrices Px and Py are correctly computed!');
    else
        disp('Error detected: Projection matrices Px and Py may be incorrect.');
    end

    % Test the P_active matrix
    P_active = P_full(activeSensorDOFs, activeDOFs);
    sensorCoordinates_active_vector_proj = P_active * nodeCoordinates_vector(activeDOFs);
    differen_active = sensorCoordinates_active_vector_proj - sensorCoordinates_active_vector;

    if norm(differen_active) < tolerance
        disp('Projection matrix P_active is correctly computed!');
    else
        disp('Error detected: Projection matrix P_active may be incorrect.');
    end

    % Test the Px_active and Py_active matrices
    x_sensor_active_proj = Px_active * nodeCoordinates_vector(activeDOFsX);
    y_sensor_active_proj = Py_active * nodeCoordinates_vector(activeDOFsY);
    differen_x_active = x_sensor_active_proj - x_sensor_active;
    differen_y_active = y_sensor_active_proj - y_sensor_active;
    error_x_active = norm(differen_x_active);
    error_y_active = norm(differen_y_active);

    if error_x_active < tolerance && error_y_active < tolerance
        disp('Projection matrices Px_active and Py_active are correctly computed!');
    else
        disp('Error detected: Projection matrices Px_active and Py_active may be incorrect.');
    end

    %% Step 7: Assign Back to BVP
    BVP.obs.Px = Px;
    BVP.obs.Py = Py;
    % BVP.obs.Px_active = Px_active;
    % BVP.obs.Py_active = Py_active;
    BVP.obs.Px_active = Px_active0;
    BVP.obs.Py_active = Py_active0;
    BVP.obs.P_full = P_full;
    disp('5. Projection Matrix Computation is Completed.');
end

function inside = isInsideTriangle(x, y, x_elem, y_elem)
    % ISINSIDETRIANGLE Uses inpolygon to check if (x, y) is inside a triangle.
    % Inputs:
    %   x, y       - Coordinates of the point
    %   x_elem     - X-coordinates of triangle nodes (1x3)
    %   y_elem     - Y-coordinates of triangle nodes (1x3)
    % Output:
    %   inside     - Logical true/false

    % Ensure it's a closed polygon (optional but safer for numerical checks)
    x_poly = [x_elem(:); x_elem(1)];
    y_poly = [y_elem(:); y_elem(1)];

    % Use MATLAB's built-in inpolygon
    inside = inpolygon(x, y, x_poly, y_poly);
end

function [xi, eta] = directMappingTri(x, y, x_elem, y_elem)
    % Map (x,y) into reference triangle (ξ,η) using affine transformation
    x1 = x_elem(1); y1 = y_elem(1);
    x2 = x_elem(2); y2 = y_elem(2);
    x3 = x_elem(3); y3 = y_elem(3);

    A = [x2 - x1, x3 - x1;
         y2 - y1, y3 - y1];
    b = [x - x1; y - y1];

    xi_eta = A \ b;
    xi = xi_eta(1);
    eta = xi_eta(2);
end
