% PCE BOOTSTRAP UNCERTAINTY — Distributions
% Log(4000,1000), Norm(4000,1000)

clearvars
setpath
uqlab

output_indices = [10, 65, 120, 139, 195, 248, 259, 260, 261, 262];

% LaTeX versions, used for axis labels / figure titles
output_labels  = {'$u_{inner,start}$', '$u_{center}$', '$u_{inner,end}$', ...
                  '$v_{inner,start}$', '$v_{center}$', '$v_{inner,end}$', ...
                  '$\mathrm{Re}(\lambda)$', '$\mathrm{Im}(\lambda)$', ...
                  '$E$', '$\mathcal{E}$'};

% Plain-text versions, used for filenames (no LaTeX markup, filesystem-safe)
output_labels_plain = {'u_inner_start', 'u_center', 'u_inner_end', ...
                  'v_inner_start', 'v_center', 'v_inner_end', ...
                  'Re_lambda', 'Im_lambda', ...
                  'E', 'Enstrophy'};

% --- Distribution configurations ---
configs = {
    1,  'LOG',  'PCE_OLS_LOG_1_bootstrap.mat',  'PCE_OLS_1b',  'dati_LOG_4000_1000.mat', 'X_LHS_1', 'Y_LHS_1', 'Lognormal', [4000 1000], 'Moments';
    6,  'NORM', 'PCE_OLS_NORM_6_bootstrap.mat', 'PCE_OLS_6b',  'dati_NORM_4000_1000.mat','X_LHS_6', 'Y_LHS_6', 'Gaussian',  [4000 1000], 'Parameters';
};

Xval  = linspace(1, 8000, 1000)';
cmap  = lines(2);
n_out = length(output_indices);

% =========================================================
% LOOP OVER DISTRIBUTIONS
% =========================================================
for d = 1:size(configs, 1)

   idx         = configs{d,1};
   dist_label  = configs{d,2};
   pce_file    = configs{d,3};
   pce_varname = configs{d,4};
   % configs{d,5} = LHS file
   x_varname   = configs{d,6};
   y_varname   = configs{d,7};
   dist_type   = configs{d,8};
   params      = configs{d,9};
   param_field = configs{d,10};

    fprintf('\n=== Distribution %d — %s mu=%d sigma=%d ===\n', ...
            idx, dist_type, params(1), params(2))

   % --- Load PCE bootstrap ---
    data_pce = load(pce_file);
    PCE_var  = data_pce.(pce_varname);

    % --- Load training data ---
    data_lhs = load(configs{d,5});
    X_data   = data_lhs.(x_varname);
    Y_data   = data_lhs.(y_varname);

    % --- Recreate input ---
    InputVal.Marginals.Type          = dist_type;
    InputVal.Marginals.(param_field) = params;
    InputVal.Marginals.Bounds        = [1, 8000];
    myInputVal = uq_createInput(InputVal);
    uq_selectInput(myInputVal);

    % --- Evaluate PCE and bootstrap ---
    fprintf('  Evaluating PCE bootstrap on validation grid...\n')
    [YPCval_all, ~, YPCval_Bootstrap_all] = uq_evalModel(PCE_var, Xval);
    fprintf('  PCE output size:       %s\n', mat2str(size(YPCval_all)))
    fprintf('  Bootstrap output size: %s\n', mat2str(size(YPCval_Bootstrap_all)))

    % --- Title string ---
    if strcmp(dist_type, 'Lognormal')
        title_str = sprintf('Lognormal $\\mu$=%d $\\sigma$=%d', params(1), params(2));
    else
        title_str = sprintf('Normal $\\mu$=%d $\\sigma$=%d', params(1), params(2));
    end

    % =====================================================
    % FIGURES
    % =====================================================
    for s = 1:n_out
        out_idx    = output_indices(s);
        out_name   = output_labels{s};
        out_plain  = output_labels_plain{s};

        % Extract output s
        Ymean = YPCval_all(:, out_idx);
        Yboot = squeeze(YPCval_Bootstrap_all(:, :, out_idx));

        % Bootstrap confidence bounds
        Ylo = quantile(Yboot, 0.025, 2);
        Yhi = quantile(Yboot, 0.975, 2);

        % Training points
        Xed = X_data;
        Yed = Y_data(:, out_idx);

        % Plot
        fig = figure('Name', sprintf('Bootstrap dist %d — %s', idx, out_name), ...
                     'NumberTitle', 'off', 'Position', [100 100 900 500]);

        % Bootstrap replications in background
        pb = plot(Xval, Yboot, 'LineWidth', 0.5, ...
                  'Color', [cmap(2,:), 0.15]);
        hold on

        % 95% CI band
        fill([Xval; flipud(Xval)], [Ylo; flipud(Yhi)], cmap(1,:), ...
             'FaceAlpha', 0.2, 'EdgeColor', 'none')

        % PCE mean
        p2 = plot(Xval, Ymean, 'Color', cmap(1,:), 'LineWidth', 2, ...
                  'DisplayName', 'PCE mean');

        % CI bounds
        p3 = plot(Xval, Ylo, '--', 'Color', cmap(1,:), 'LineWidth', 1.2, ...
                  'DisplayName', '95\% CI');
             plot(Xval, Yhi, '--', 'Color', cmap(1,:), 'LineWidth', 1.2, ...
                  'HandleVisibility', 'off');

        % Training points on top
        p5 = plot(Xed, Yed, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k', ...
                  'DisplayName', 'Training points');

        % Reorder layers
        h = get(gca, 'Children');
        set(gca, 'Children', h(end:-1:1))

        hold off
        xlabel('$Re$', 'Interpreter', 'latex')
        ylabel(out_name, 'Interpreter', 'latex')
        title(sprintf('PCE Bootstrap - %s - %s', title_str, out_name), ...
              'Interpreter', 'latex')
        legend([p5, p2, p3, pb(1)], ...
               {'Training points', 'PCE mean', '95\% CI', 'Bootstrap replications'}, ...
               'Interpreter', 'latex', 'Location', 'best')
        grid on; box on

        % Save figure
        fname = sprintf('Bootstrap_PCE_%s_%d_%s.png', dist_label, idx, out_plain);
        if exist(fname, 'file'), delete(fname); end
        exportgraphics(fig, fname, 'Resolution', 300);
        close(fig)
        fprintf('  Saved: %s\n', fname)

    end

    fprintf('  All figures saved for distribution %d\n', idx)

    clear data_pce data_lhs PCE_var X_data Y_data InputVal myInputVal ...
      YPCval_all YPCval_Bootstrap_all
    uqlab

end

fprintf('\n=== Bootstrap uncertainty analysis completed ===\n')