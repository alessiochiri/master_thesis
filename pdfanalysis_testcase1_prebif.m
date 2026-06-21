% EMPIRICAL PDFs AND CONFIDENCE BANDS - PCE OLS
% 2 distributions: Lognormal(4000,1000), Gaussian(4000,1000)

setpath
uqlab

N_MC = 1e5;

% --- Distribution definitions ---
% Variable/file names match exactly what the PCE-fitting script saves:
%   PCE_OLS_LOG_4000_1000.mat   -> PCE_OLS_1
%   PCE_OLS_NORM_4000_1000.mat  -> PCE_OLS_6
dist_info = {
    'PCE_OLS_LOG_4000_1000.mat',  'PCE_OLS_1', 'Lognormal', [4000 1000], 'Log $\mu$=4000 $\sigma$=1000';
    'PCE_OLS_NORM_4000_1000.mat', 'PCE_OLS_6', 'Gaussian',  [4000 1000], 'Norm $\mu$=4000 $\sigma$=1000';
};

n_dist = size(dist_info, 1);

% --- Output indices ---
idx_u    = 1:129;       % u(y) profile
idx_v    = 130:258;     % v(x) profile
idx_ReL  = 259;
idx_ImL  = 260;
idx_KE   = 261;
idx_Enst = 262;

scalar_indices = [idx_ReL, idx_ImL, idx_KE, idx_Enst];
scalar_names   = {'$\mathrm{Re}(\lambda)$', '$\mathrm{Im}(\lambda)$', '$E$', '$\mathcal{E}$'};

% --- Physical coordinates ---
grid_data  = load('square_stokes_nobc.mat');
xy         = grid_data.xy;
tol        = 1e-6;
idx_vert   = find(abs(xy(:,1)) < tol);
y_nodes    = sort(xy(idx_vert, 2));          % 129x1, u(y)
idx_horiz  = find(abs(xy(:,2)) < tol);
x_nodes    = sort(xy(idx_horiz, 1));         % 129x1, v(x)

% --- Percentiles ---
prct_levels = [5, 25, 50, 75, 95];

% --- Storage ---
Y_MC       = cell(n_dist, 1);   % MC samples for each distribution
Y_pct_u    = zeros(n_dist, 129, length(prct_levels));
Y_pct_v    = zeros(n_dist, 129, length(prct_levels));
Y_scalar   = zeros(n_dist, 4, N_MC);

% =========================================================
% MAIN LOOP - Monte Carlo sampling
% =========================================================
for k = 1:n_dist

    fprintf('\n=== Distribution %d/%d: %s ===\n', k, n_dist, dist_info{k,5});

    % Load PCE
    data = load(dist_info{k,1});
    myPCE = data.(dist_info{k,2});

    % Create input
    InputOpts.Marginals.Type = dist_info{k,3};
    if strcmp(dist_info{k,3}, 'Lognormal')
        InputOpts.Marginals.Moments    = dist_info{k,4};
    else
        InputOpts.Marginals.Parameters = dist_info{k,4};
    end
    InputOpts.Marginals.Bounds = [1 8000];
    myInput = uq_createInput(InputOpts);

    % Sample input
    X_MC = uq_getSample(myInput, N_MC);

    % Evaluate PCE
    fprintf('  Evaluating PCE on %d samples...\n', N_MC);
    Y_full = uq_evalModel(myPCE, X_MC);   % N_MC x 262

    % Store scalar samples
    Y_scalar(k, :, :) = Y_full(:, scalar_indices)';

    % Profile percentiles
    Y_pct_u(k, :, :) = prctile(Y_full(:, idx_u), prct_levels, 1)';
    Y_pct_v(k, :, :) = prctile(Y_full(:, idx_v), prct_levels, 1)';

    % Store everything for KDE
    Y_MC{k} = Y_full(:, scalar_indices);

    clear data myPCE InputOpts myInput X_MC Y_full
    uqlab

end

fprintf('\n=== Sampling completed ===\n')

% =========================================================
% FIGURES - Scalar PDFs for each distribution
% =========================================================
for k = 1:n_dist

    fig = figure('Name', sprintf('PDF - Dist %d', k), ...
                 'NumberTitle', 'off', ...
                 'Position', [100 100 1400 400]);

    for s = 1:4
        subplot(1, 4, s)
        samples = squeeze(Y_scalar(k, s, :));
        [f, xi] = ksdensity(samples);
        plot(xi, f, 'LineWidth', 2)
        xlabel(scalar_names{s}, 'Interpreter', 'latex')
        ylabel('PDF')
        title(scalar_names{s}, 'Interpreter', 'latex')
        grid on
        box on
    end

    sgtitle(sprintf('Dist %d - Empirical PDF - %s', k, dist_info{k,5}), ...
            'Interpreter', 'latex', 'FontSize', 12);

    exportgraphics(fig, sprintf('PDF_scalars_dist_%02d.png', k), 'Resolution', 300);
    close(fig)
    fprintf('Figure saved: PDF_scalars_dist_%02d.png\n', k);

end

% =========================================================
% FIGURES - Confidence bands for u(y) and v(x)
% =========================================================
band_color   = [0.2 0.6 1.0];
median_color = [0.0 0.2 0.8];

for k = 1:n_dist

    fig = figure('Name', sprintf('Bands - Dist %d', k), ...
                 'NumberTitle', 'off', ...
                 'Position', [100 100 900 500]);

    % --- u(y) ---
    subplot(1, 2, 1)
    hold on

    u_p5  = squeeze(Y_pct_u(k,:,1))';   % 129x1
    u_p25 = squeeze(Y_pct_u(k,:,2))';
    u_p50 = squeeze(Y_pct_u(k,:,3))';
    u_p75 = squeeze(Y_pct_u(k,:,4))';
    u_p95 = squeeze(Y_pct_u(k,:,5))';

    % 5th-95th percentile band
    fill([u_p5;  flipud(u_p95)], ...
         [y_nodes; flipud(y_nodes)], ...
         band_color, 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
         'DisplayName', '5th-95th pct')

    % 25th-75th percentile band
    fill([u_p25; flipud(u_p75)], ...
         [y_nodes; flipud(y_nodes)], ...
         band_color, 'FaceAlpha', 0.4, 'EdgeColor', 'none', ...
         'DisplayName', '25th-75th pct')

    % Median
    plot(u_p50, y_nodes, '-', 'Color', median_color, 'LineWidth', 2, ...
         'DisplayName', 'Median')

    hold off
    xlabel('$u$', 'Interpreter', 'latex')
    ylabel('$y$', 'Interpreter', 'latex')
    title('u(y) confidence bands')
    legend('Location', 'northwest')
    grid on
    box on

    % --- v(x) ---
    subplot(1, 2, 2)
    hold on

    v_p5  = squeeze(Y_pct_v(k,:,1))';
    v_p25 = squeeze(Y_pct_v(k,:,2))';
    v_p50 = squeeze(Y_pct_v(k,:,3))';
    v_p75 = squeeze(Y_pct_v(k,:,4))';
    v_p95 = squeeze(Y_pct_v(k,:,5))';

    % 5th-95th percentile band
    fill([x_nodes;  flipud(x_nodes)], ...
         [v_p5; flipud(v_p95)], ...
         band_color, 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
         'DisplayName', '5th-95th pct')

    % 25th-75th percentile band
    fill([x_nodes; flipud(x_nodes)], ...
         [v_p25; flipud(v_p75)], ...
         band_color, 'FaceAlpha', 0.4, 'EdgeColor', 'none', ...
         'DisplayName', '25th-75th pct')

    % Median
    plot(x_nodes, v_p50, '-', 'Color', median_color, 'LineWidth', 2, ...
         'DisplayName', 'Median')

    hold off
    xlabel('$x$', 'Interpreter', 'latex')
    ylabel('$v$', 'Interpreter', 'latex')
    title('v(x) confidence bands')
    legend('Location', 'northeast')
    grid on
    box on

    sgtitle(sprintf('Dist %d - Velocity Confidence Bands - %s', k, dist_info{k,5}), ...
            'Interpreter', 'latex', 'FontSize', 12);

    fname = sprintf('Bands_dist_%02d.png', k);
    if exist(fname, 'file'), delete(fname); end
    exportgraphics(fig, fname, 'Resolution', 300);
    close(fig)
    fprintf('Figure saved: %s\n', fname);

end

% =========================================================
% FINAL FIGURE - Overlaid PDFs per scalar output (both distributions)
% =========================================================
dist_colors = lines(n_dist);
dist_styles = {'-', '--'};   % solid for Lognormal, dashed for Gaussian

fig = figure('Name', 'PDF comparison - all distributions', ...
             'NumberTitle', 'off', ...
             'Position', [100 100 1400 900]);

for s = 1:4
    subplot(2, 2, s)
    hold on

    for k = 1:n_dist
        samples = squeeze(Y_scalar(k, s, :));
        [f, xi] = ksdensity(samples);
        plot(xi, f, dist_styles{k}, 'Color', dist_colors(k,:), 'LineWidth', 1.5, ...
             'DisplayName', dist_info{k,5})
    end

    hold off
    xlabel(scalar_names{s}, 'Interpreter', 'latex')
    ylabel('PDF')
    title(scalar_names{s}, 'Interpreter', 'latex')
    legend('Interpreter', 'latex', 'Location', 'best', 'FontSize', 8)
    grid on
    box on
end

sgtitle('Empirical PDF Comparison - All Distributions', 'FontSize', 13)

exportgraphics(fig, 'PDF_comparison_all.png', 'Resolution', 300);
close(fig)
fprintf('Figure saved: PDF_comparison_all.png\n')

disp('=== All figures generated ===')

% =========================================================
% COMPARISON FIGURE - Medians and bands for u(y) and v(x)
% Lognormal vs Normal, side by side
% =========================================================

dist_labels = dist_info(:,5);   % {'Log $\mu$=4000 $\sigma$=1000', 'Norm $\mu$=4000 $\sigma$=1000'}

% =========================================================
% FIGURE u(y)
% =========================================================
fig = figure('Name', 'u(y) - All distributions', ...
             'NumberTitle', 'off', ...
             'Position', [100 100 900 600]);

hold on
for k = 1:n_dist
    u_p5  = squeeze(Y_pct_u(k,:,1))';
    u_p25 = squeeze(Y_pct_u(k,:,2))';
    u_p50 = squeeze(Y_pct_u(k,:,3))';
    u_p75 = squeeze(Y_pct_u(k,:,4))';
    u_p95 = squeeze(Y_pct_u(k,:,5))';
    c = dist_colors(k,:);

    % 5th-95th percentile band
    fill([u_p5; flipud(u_p95)], [y_nodes; flipud(y_nodes)], ...
         c, 'FaceAlpha', 0.08, 'EdgeColor', 'none', 'HandleVisibility', 'off')
    % 25th-75th percentile band
    fill([u_p25; flipud(u_p75)], [y_nodes; flipud(y_nodes)], ...
         c, 'FaceAlpha', 0.20, 'EdgeColor', 'none', 'HandleVisibility', 'off')
    % Median
    plot(u_p50, y_nodes, '-', 'Color', c, 'LineWidth', 2, ...
         'DisplayName', dist_labels{k})
end
hold off
xlabel('$u$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')
title('Lognormal vs Normal')
legend('Interpreter', 'latex', 'Location', 'northwest', 'FontSize', 9)
xlim auto
ylim([-1 1])
grid on
box on

sgtitle('u(y) - Median and Confidence Bands', 'FontSize', 13)

fname = 'Comparison_u_all.png';
if exist(fname, 'file'), delete(fname); end
exportgraphics(fig, fname, 'Resolution', 300);
close(fig)
fprintf('Figure saved: %s\n', fname)

% =========================================================
% FIGURE v(x)
% =========================================================
fig = figure('Name', 'v(x) - All distributions', ...
             'NumberTitle', 'off', ...
             'Position', [100 100 900 600]);

hold on
for k = 1:n_dist
    v_p5  = squeeze(Y_pct_v(k,:,1))';
    v_p25 = squeeze(Y_pct_v(k,:,2))';
    v_p50 = squeeze(Y_pct_v(k,:,3))';
    v_p75 = squeeze(Y_pct_v(k,:,4))';
    v_p95 = squeeze(Y_pct_v(k,:,5))';
    c = dist_colors(k,:);

    fill([x_nodes; flipud(x_nodes)], [v_p5; flipud(v_p95)], ...
         c, 'FaceAlpha', 0.08, 'EdgeColor', 'none', 'HandleVisibility', 'off')
    fill([x_nodes; flipud(x_nodes)], [v_p25; flipud(v_p75)], ...
         c, 'FaceAlpha', 0.20, 'EdgeColor', 'none', 'HandleVisibility', 'off')
    plot(x_nodes, v_p50, '-', 'Color', c, 'LineWidth', 2, ...
         'DisplayName', dist_labels{k})
end
hold off
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$v$', 'Interpreter', 'latex')
title('Lognormal vs Normal')
legend('Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 9)
xlim([-1 1])
ylim auto
grid on
box on

sgtitle('v(x) - Median and Confidence Bands', 'FontSize', 13)

fname = 'Comparison_v_all.png';
if exist(fname, 'file'), delete(fname); end
exportgraphics(fig, fname, 'Resolution', 300);
close(fig)
fprintf('Figure saved: %s\n', fname)