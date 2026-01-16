function BVP = chooseBestModel(BVP)

    % Choose the best model based on some criteria (e.g., Relative RMSE)
    [relative_rmse, index_min_relative_rmse] = min(BVP.discovered.relative_rmse_iter);
    eps_u = BVP.discovered.epsilon_u_iter(index_min_relative_rmse);
    w = BVP.discovered.w_latex_full_iter(index_min_relative_rmse);
    prop = BVP.discovered.hyperelastic_matparameters(:,index_min_relative_rmse+1);


       % [a,b] = min(BVP.discovered.relative_rmse_iter);
        % BVP.discovered.hyperelastic_matparameters(:,b+1)

    % Example: Select the model with the lowest RMSE

    % Return the updated BVP structure
    BVP.discovered.best_model.eps_u = eps_u;
    BVP.discovered.best_model.w = w;
    BVP.discovered.best_model.relative_rmse = relative_rmse;
    BVP.discovered.best_model.prop = prop;
end
