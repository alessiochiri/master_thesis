% CONVERGENCE CURVE - LOO error vs Degree D
% Constraints N=2D and N=3D, k-means selection
% 2 distributions, 10 selected outputs

setpath
uqlab

% --- Output map of interest ---
output_indices = [30, 60, 90, 159, 189, 219, 259, 260, 261, 262];
output_names   = {'$u$ - node 30', '$u$ - node 60', '$u$ - node 90', ...
                  '$v$ - node 30', '$v$ - node 60', '$v$ - node 90', ...
                  '$\mathrm{Re}(\lambda)$', '$\mathrm{Im}(\lambda)$', '$E$', '$\mathcal{E}$'};

% --- Distribution definitions ---
% Variable/file names match the PCE-fitting script:
%   dati_LOG_4000_1000.mat   -> X_LHS_1, Y_LHS_1
%   dati_NORM_4000_1000.mat  -> X_LHS_6, Y_LHS_6
dist_info = {
    'dati_LOG_4000_1000.mat',  'X_LHS_1', 'Y_LHS_1', 'Lognormal', [4000 1000], 'Log $\mu$=4000 $\sigma$=1000';
    'dati_NORM_4000_1000.mat', 'X_LHS_6', 'Y_LHS_6', 'Gaussian',  [4000 1000], 'Norm $\mu$=4000 $\sigma$=1000';
};

% --- Degrees to explore ---
D_2N = 1:20;
D_3N = 1:13;

n_out  = length(output_indices);
n_dist = size(dist_info, 1);

% --- Storage for results ---
LOO_2D = nan(n_dist, length(D_2N), n_out);
LOO_3D = nan(n_dist, length(D_3N), n_out);

% =========================================================
% MAIN LOOP
% =========================================================
for k = 1:n_dist
    rng(100, 'twister')
    fprintf('\n=== Distribution %d/%d: %s ===\n', k, n_dist, dist_info{k,6});

    % Load data
    data = load(dist_info{k,1});
    X    = data.(dist_info{k,2});   % 40x1
    Y    = data.(dist_info{k,3});   % 40x262

    % Keep only the 10 outputs of interest
    Y_sel = Y(:, output_indices);   % 40x10

    % Create UQLab input
    InputOpts.Marginals.Type = dist_info{k,4};
    if strcmp(dist_info{k,4}, 'Lognormal')
        InputOpts.Marginals.Moments    = dist_info{k,5};
    else
        InputOpts.Marginals.Parameters = dist_info{k,5};
    end
    InputOpts.Marginals.Bounds = [1 8000];
    myInput = uq_createInput(InputOpts);
    uq_selectInput(myInput);

    % --- N = 2D ---
    fprintf('  Constraint N=2D...\n')
    for di = 1:length(D_2N)
        D = D_2N(di);
        N = 2 * D;
        fprintf('    D=%d, N=%d\n', D, N);

        [X_sub, idx_kmeans] = uq_subsample(X, N, 'k-means');
        Y_sub = Y_sel(idx_kmeans, :);   % Nx10

        MetaOpts.Type        = 'Metamodel';
        MetaOpts.MetaType    = 'PCE';
        MetaOpts.Degree      = D;
        MetaOpts.Method      = 'OLS';
        MetaOpts.ExpDesign.X = X_sub;
        MetaOpts.ExpDesign.Y = Y_sub;

        try
            myPCE = uq_createModel(MetaOpts);
            for s = 1:n_out
                LOO_2D(k, di, s) = myPCE.Error(s).LOO;  % index 1:10
            end
        catch ME
            fprintf('    ERROR D=%d N=%d: %s\n', D, N, ME.message);
        end

        clear MetaOpts myPCE X_sub Y_sub idx_kmeans
    end

    % --- N = 3D ---
    fprintf('  Constraint N=3D...\n')
    for di = 1:length(D_3N)
        D = D_3N(di);
        N = 3 * D;
        fprintf('    D=%d, N=%d\n', D, N);

        [X_sub, idx_kmeans] = uq_subsample(X, N, 'k-means');
        Y_sub = Y_sel(idx_kmeans, :);   % Nx10

        MetaOpts.Type        = 'Metamodel';
        MetaOpts.MetaType    = 'PCE';
        MetaOpts.Degree      = D;
        MetaOpts.Method      = 'OLS';
        MetaOpts.ExpDesign.X = X_sub;
        MetaOpts.ExpDesign.Y = Y_sub;

        try
            myPCE = uq_createModel(MetaOpts);
            for s = 1:n_out
                LOO_3D(k, di, s) = myPCE.Error(s).LOO;  % index 1:10
            end
        catch ME
            fprintf('    ERROR D=%d N=%d: %s\n', D, N, ME.message);
        end

        clear MetaOpts myPCE X_sub Y_sub idx_kmeans
    end

    clear data X Y Y_sel InputOpts myInput
    uqlab
end

% Save results
if exist('convergence_LOO.mat', 'file')
    delete('convergence_LOO.mat')
end
save('convergence_LOO.mat', 'LOO_2D', 'LOO_3D', ...
     'D_2N', 'D_3N', 'output_indices', 'output_names', 'dist_info');
disp('=== Results saved to convergence_LOO.mat ===')

% =========================================================
% PLOT
% =========================================================
for k = 1:n_dist

    fig = figure('Name', sprintf('Convergence - Dist %d', k), ...
                 'NumberTitle', 'off', ...
                 'Position', [100 100 1600 900]);

    for s = 1:n_out
        subplot(2, 5, s);
        hold on

        Q2_2D = 1 - squeeze(LOO_2D(k, :, s));
        Q2_3D = 1 - squeeze(LOO_3D(k, :, s));

        plot(D_2N, Q2_2D, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 5);
        plot(D_3N, Q2_3D, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 5);

        yline(0.95, 'k--', 'LineWidth', 1);
        yline(0.90, 'k:',  'LineWidth', 1);

        hold off
        xlabel('Degree $D$', 'Interpreter', 'latex')
        ylabel('$Q^2$', 'Interpreter', 'latex')
        title(output_names{s}, 'Interpreter', 'latex')
        legend('N=2D', 'N=3D', '$Q^2$=0.95', '$Q^2$=0.90', ...
               'Interpreter', 'latex', 'Location', 'southeast')
        ylim([-0.1 1.05])
        grid on
        box on
    end

    sgtitle(sprintf('$Q^2$ Convergence vs Degree - %s', dist_info{k,6}), ...
            'Interpreter', 'latex', 'FontSize', 13);

    fname = sprintf('Convergence_dist_%02d.png', k);
    if exist(fname, 'file'), delete(fname); end
    exportgraphics(fig, fname, 'Resolution', 300);
    close(fig)
    fprintf('Figure saved: %s\n', fname);
end

disp('=== All plots generated ===')