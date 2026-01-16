function [u, force, Krr, Kru, Kur, Kuu, Ur, Uu, Rr, Ru, stretch_history, energy_history, W_total, J_min_max] = solver_HE(preProc_Variables, force, prop)
%% Assign from BVP
% global energy_history stretch_history
GDOFs = preProc_Variables.msh.GDOFs;
ndbc = preProc_Variables.BC.ndbc;
ntbc = preProc_Variables.BC.ntbc;
dbc = preProc_Variables.BC.dbc;
prescribedDOFs = preProc_Variables.BC.prescribedDOFs;
tol = preProc_Variables.tol;
maxit = preProc_Variables.maxit;
maxreit = preProc_Variables.maxreit;
timeInterval = preProc_Variables.timeInterval;

BVP_tmp.Q_pts = preProc_Variables.Q_pts;
BVP_tmp.Q_wts = preProc_Variables.Q_wts;
BVP_tmp.Nshape = preProc_Variables.Nshape;
BVP_tmp.dN_dxi = preProc_Variables.dN_dxi;
BVP_tmp.dN_deta = preProc_Variables.dN_deta;
BVP_tmp.nElm = preProc_Variables.msh.nElm;
BVP_tmp.nodeCoordinates = preProc_Variables.msh.nodeCoordinates;
BVP_tmp.elementNodes = preProc_Variables.msh.elementNodes;
BVP_tmp.elementDOFs = preProc_Variables.msh.elementDOFs;
BVP_tmp.GDOFs = preProc_Variables.msh.GDOFs;
BVP_tmp.thickness = preProc_Variables.msh.T;
BVP_tmp.nElNode = preProc_Variables.msh.nElNode;
BVP_tmp.DOFs = preProc_Variables.msh.DOFs;
BVP_tmp.top_right_element = preProc_Variables.msh.top_right_element;
BVP_tmp.bottom_left_element = preProc_Variables.msh.bottom_left_element;
BVP_tmp.top_left_element = preProc_Variables.msh.top_left_element;
BVP_tmp.bottom_right_element = preProc_Variables.msh.bottom_right_element;
tracking_enabled = true;
%% Calculation
energy_history = [0; 0; 0; 0]; % Initialize W_out
stretch_history = [1; 1; 1; 1]; % Initialize lambda_out

u = zeros(GDOFs, 1); % Nodal displacements
u_steps = u;
cu = zeros(GDOFs, 1); % Converged nodal displacements
step = 0; % Load/displacement step index
curtime = 0; % Current time
reit = 0; % Reduction index for time stepping
cnit = []; % Record number of Newton iterations
counter = 0;

while curtime ~= 1 % Load/time stepping loop
    counter = counter + 1;
    curtime = curtime + timeInterval;
    
    if curtime > 1
        timeInterval = 1 - curtime + timeInterval;
        curtime = 1;
    end
    
    err = 1e6; % Initial error
    perr = err;
    nit = 0; % Newton iteration count
    iu = zeros(GDOFs, 1); % Iterative displacement
    
    while (err > tol) && (nit <= maxit)
        nit = nit + 1;
        
        [k, r, lambda_gp, W_gp, W_total, J_min_max] = globalstiffness_HE(BVP_tmp, u, prop);
        
        % External force setup
        f = zeros(GDOFs, 1);
        
        if ntbc ~= 0
            f = force;
        end
        
        % Apply displacement boundary conditions
        if ndbc ~= 0
            k(prescribedDOFs, :) = 0;
            k(prescribedDOFs, prescribedDOFs) = eye(ndbc);
            f(prescribedDOFs) = 0;
            
            if nit == 1
                f(prescribedDOFs) = dbc;
            end
            
        end
        
        % Residual vector
        b = curtime * f - r;
        
        if ndbc ~= 0
            b(prescribedDOFs) = curtime * dbc - u(prescribedDOFs);
        end
        
        % Solve linear system
        
        
        %%
        du = k \ b;
        u = u + du;
        iu = iu + du;
        % norm(du0 - du)
        % Compute relative error
        alldof = 1:GDOFs;
        freedof = setdiff(alldof, prescribedDOFs);
        
        if nit > 1
            num = b(freedof)' * b(freedof);
            denom = 1 + f(freedof)' * f(freedof);
            err = num / denom;
        end
        
        % Divergence detection
        if err / perr > 1e6 && nit > 2
            nit = maxit + 1;
        else
            perr = err;
        end
        
    end
    
    % Handle convergence or failure
    if nit <= maxit
        reit = 0;
        step = step + 1;
        cu = u;
        cnit = [cnit, nit];
        
        if length(cnit) >= 2 && all(cnit(end - 1:end) <= 5)
            timeInterval = timeInterval * 1.5;
        end
        
        if tracking_enabled && ~isempty(lambda_gp)
            stretch_history(:, step + 1) = lambda_gp;
            energy_history(: , step + 1) = W_gp;
        end
        
        Krr = 0;
        Kru = 0;
        Kur = 0;
        Kuu = 0;
        
        Rr = 0;
        Uu = 0;
        Ur = 0;
        Ru = 0;
        
        u_steps = [u_steps u];
    else
        
        if reit <= maxreit
            curtime = curtime - timeInterval;
            timeInterval = timeInterval / 4;
            reit = reit + 1;
            u = cu;
        else
            return; % Terminate analysis due to non-convergence
        end
        
    end
    
end

% Round the small values to zero
u_steps(abs(u_steps) < 1e-12) = 0;
u = u_steps(:, end);
end
