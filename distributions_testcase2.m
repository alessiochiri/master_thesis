% CONDITIONAL AND UNCONDITIONAL PDFs
% Bivariate case (nu, V) — 9 outputs

N_cond = 1e4;

% QoI set: 3 scalars + 6 velocity nodes
idx_ReL  = 259;
idx_KE   = 261;
idx_Enst = 262;
idx_nodes = [10, 65, 120, 139, 195, 248];
node_names = {'u - node 10', 'u - node 65', 'u - node 120', 'v - node 139', 'v - node 195', 'v - node 248'};

idx_all = [idx_ReL, idx_KE, idx_Enst, idx_nodes];
labels  = [{'Re(\lambda)', 'KE', 'Enstrophy'}, node_names];
n_qoi   = length(idx_all);

n_cols = 3;
n_rows = ceil(n_qoi / n_cols);


% 3a. CONDITIONAL DISTRIBUTION — fixed U_lid = 1,2,3,4,5

U_fixed = [1, 2, 3, 4, 5];
colors_U = lines(5);

X_full = uq_getSample(myInputA, N_cond, 'MC');
nu_samples = X_full(:, 1);

figure('Name', 'Conditional PDF - fixed U', 'Position', [100 100 1500 350*n_rows]);
for i = 1:n_qoi
    subplot(n_rows, n_cols, i); hold on
    for k = 1:5
        X_cond = [nu_samples, U_fixed(k) * ones(N_cond, 1)];
        Y_cond = uq_evalModel(PCE_A, X_cond);
        y_k = Y_cond(:, idx_all(i));
        [f, x] = ksdensity(y_k);
        plot(x, f, 'Color', colors_U(k,:), 'LineWidth', 1.6,'DisplayName', sprintf('U = %d', U_fixed(k)))
    end
    xlabel(labels{i}, 'FontSize', 9); ylabel('PDF', 'FontSize', 9)
    title(labels{i}, 'FontSize', 10)
    if i == 1
        legend('Location', 'best', 'FontSize', 7)
    end
    grid on
end
sgtitle('Conditional distributions: fixed U_{lid}', 'FontSize', 13)
exportgraphics(gcf, 'Conditional_PDF_fixedU.png', 'Resolution', 300);


% 3b. CONDITIONAL DISTRIBUTION — fixed nu at percentiles 10, 50, 90

pct = [0.10, 0.50, 0.90];
nu_pct = uq_all_invcdf(pct(:), myInputA.Marginals(1));
labels_pct = {'10th pct', '50th pct', '90th pct'};
colors_nu = [0.2 0.4 0.8; 0.1 0.7 0.3; 0.9 0.2 0.2];

X_full2 = uq_getSample(myInputA, N_cond, 'MC');
U_samples = X_full2(:, 2);

figure('Name', 'Conditional PDF - fixed nu', 'Position', [100 100 1500 350*n_rows]);
for i = 1:n_qoi
    subplot(n_rows, n_cols, i); hold on
    for k = 1:3
        X_cond = [nu_pct(k) * ones(N_cond, 1), U_samples];
        Y_cond = uq_evalModel(PCE_A, X_cond);
        y_k = Y_cond(:, idx_all(i));
        [f, x] = ksdensity(y_k);
        plot(x, f, 'Color', colors_nu(k,:), 'LineWidth', 1.6,'DisplayName', sprintf('\\nu = %.4f (%s)', nu_pct(k), labels_pct{k}))
    end
    xlabel(labels{i}, 'FontSize', 9); ylabel('PDF', 'FontSize', 9)
    title(labels{i}, 'FontSize', 10)
    if i == 1
        legend('Location', 'best', 'FontSize', 7)
    end
    grid on
end
sgtitle('Conditional distributions: fixed \nu (percentiles)', 'FontSize', 13)
exportgraphics(gcf, 'Conditional_PDF_fixedNu.png', 'Resolution', 300);


% 3c. MARGINAL) PDF — full uncertainty propagation
% Both nu and V vary jointly according to their true distributions

N_uncond = 5e4;   % larger sample for a smoother, more reliable PDF

X_uncond = uq_getSample(myInputA, N_uncond, 'MC');
Y_uncond = uq_evalModel(PCE_A, X_uncond);

figure('Name', 'Unconditional PDF', 'Position', [100 100 1500 350*n_rows]);
for i = 1:n_qoi
    subplot(n_rows, n_cols, i); hold on

    y_i = Y_uncond(:, idx_all(i));
    [f, x] = ksdensity(y_i);
    plot(x, f, 'k-', 'LineWidth', 2);

    % Mark mean and 90% credible interval for reference
    mu_i = mean(y_i);
    ci_lo = prctile(y_i, 5);
    ci_hi = prctile(y_i, 95);

    xline(mu_i, 'r--', 'LineWidth', 1.3, 'DisplayName', sprintf('mean = %.4g', mu_i));
    xline(ci_lo, 'b:', 'LineWidth', 1);
    xline(ci_hi, 'b:', 'LineWidth', 1);

    xlabel(labels{i}, 'FontSize', 9); ylabel('PDF', 'FontSize', 9)
    title(labels{i}, 'FontSize', 10)
    if i == 1
        legend('Location', 'best', 'FontSize', 7)
    end
    grid on
end
sgtitle('Marginal PDF — full propagation of (\nu, V) uncertainty', 'FontSize', 13)
exportgraphics(gcf, 'Unconditional_PDF.png', 'Resolution', 300);

% -------------------------------------------------------
% Summary statistics table (mean, std, 90% CI) 
% -------------------------------------------------------
fprintf('\n%-16s %-12s %-12s %-12s %-12s\n', 'Output', 'Mean', 'Std', 'P5', 'P95');
fprintf('----------------------------------------------------------------\n');
for i = 1:n_qoi
    y_i = Y_uncond(:, idx_all(i));
    fprintf('%-16s %-12.5g %-12.5g %-12.5g %-12.5g\n', ...
        labels{i}, mean(y_i), std(y_i), prctile(y_i,5), prctile(y_i,95));
end

save('PDF_analysis_results.mat', 'Y_uncond', 'X_uncond', 'idx_all', 'labels', '-v7.3');