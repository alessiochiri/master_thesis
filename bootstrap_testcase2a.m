% =========================================================
% BOOTSTRAP CI WIDTH — 2D heatmap over (nu, V) space
% Re(lambda), KE, Enstrophy
% =========================================================

clearvars
rng(100,'twister')
uqlab

fprintf('Loading data and bootstrap PCE...\n');
load('dati_due_variabili_1.mat', 'XA1', 'YA1');
load('PCE_due_variabili_bootstrap.mat', 'PCE_Bootstrap');

output_indices = [259, 261, 262];
output_labels  = {'$\mathrm{Re}(\lambda)$', '$E$', '$\mathcal{E}$'};
output_labels_plain = {'Re_lambda', 'KE', 'Enstrophy'};
n_out = length(output_indices);

nu_vec  = XA1(:,1);
vel_vec = XA1(:,2);

% -------------------------------------------------------
% Grid over the training domain, same resolution as A1
% -------------------------------------------------------
N_grid  = 150;
nu_grid  = linspace(min(nu_vec),  max(nu_vec),  N_grid);
vel_grid = linspace(min(vel_vec), max(vel_vec), N_grid);
[Nu_mesh, Vel_mesh] = meshgrid(nu_grid, vel_grid);
X_grid  = [Nu_mesh(:), Vel_mesh(:)];

fprintf('Evaluating bootstrap PCE on %d x %d grid...\n', N_grid, N_grid);
[YPC_grid, ~, YPC_Bootstrap_grid] = uq_evalModel(PCE_Bootstrap, X_grid);
fprintf('  PCE output size:       %s\n', mat2str(size(YPC_grid)));
fprintf('  Bootstrap output size: %s\n\n', mat2str(size(YPC_Bootstrap_grid)));

% =========================================================
% COMPUTE CI WIDTH ON THE GRID, FOR EACH OUTPUT
% =========================================================
CI_width_all = cell(n_out,1);

for s = 1:n_out
    out_idx = output_indices(s);

    Yboot_grid = squeeze(YPC_Bootstrap_grid(:, :, out_idx));  % [N_grid^2 x n_reps]

    Ylo_grid = quantile(Yboot_grid, 0.025, 2);
    Yhi_grid = quantile(Yboot_grid, 0.975, 2);

    CI_width = reshape(Yhi_grid - Ylo_grid, size(Nu_mesh));
    CI_width_all{s} = CI_width;

    fprintf('%s: CI width range [%.5f, %.5f], mean %.5f\n', ...
        output_labels{s}, min(CI_width(:)), max(CI_width(:)), mean(CI_width(:)));
end
fprintf('\n');

save('Bootstrap_CIwidth_grid.mat', 'CI_width_all', 'output_indices', ...
     'output_labels', 'Nu_mesh', 'Vel_mesh', '-v7.3');

% =========================================================
% FIGURE — one heatmap per output, side by side
% =========================================================
fig = figure('Position', [50 50 1500 450]);

for s = 1:n_out
    subplot(1, n_out, s)

    contourf(Nu_mesh, Vel_mesh, CI_width_all{s}, 30, 'LineColor', 'none');
    colormap(gca, 'parula');
    colorbar;

    hold on
    scatter(nu_vec, vel_vec, 15, 'k.', 'HandleVisibility', 'off');
    hold off

    xlabel('$\nu$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$V$', 'Interpreter', 'latex', 'FontSize', 12);
    title(sprintf('%s', output_labels{s}), 'Interpreter', 'latex', 'FontSize', 13);
    set(gca, 'FontSize', 10)
end

sgtitle('Bootstrap 95\% CI width over $(\nu, V)$ space', 'Interpreter', 'latex', 'FontSize', 14)

fname = 'Bootstrap_CIwidth_heatmap.png';
if exist(fname, 'file'), delete(fname); end
exportgraphics(fig, fname, 'Resolution', 300);
fprintf('Saved: %s\n', fname);