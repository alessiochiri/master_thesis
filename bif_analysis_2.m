% BIFURCATION ANALYSIS IN PARAMETER SPACE
%   A1 — Critical curve Re(lambda) = 0
%   A2 — Effective Re_c along the critical curve
%   A3 — Probability of instability (correct sampling)
%   A4 — Classification accuracy (PCE vs CFD)

clearvars
rng(100,'twister')
uqlab

load('dati_due_variabili_1.mat', 'XA1', 'YA1');
load('PCE_due_variabili.mat',    'PCE_A');

XA       = XA1;
YA       = YA1;
ReL_data = YA(:, 259);

nu_vec  = XA(:,1);
vel_vec = XA(:,2);

InputA.Marginals(1).Name    = 'nu';
InputA.Marginals(1).Type    = 'Lognormal';
InputA.Marginals(1).Moments = [0.002, 0.002];
InputA.Marginals(1).Bounds  = [0.0007, 0.01];

InputA.Marginals(2).Name       = 'V';
InputA.Marginals(2).Type       = 'Uniform';
InputA.Marginals(2).Parameters = [1, 5];

myInput = uq_createInput(InputA);
uq_selectInput(myInput);


N_grid  = 150;
nu_grid  = linspace(min(nu_vec),  max(nu_vec),  N_grid);
vel_grid = linspace(min(vel_vec), max(vel_vec), N_grid);
[Nu_mesh, Vel_mesh] = meshgrid(nu_grid, vel_grid);
X_grid  = [Nu_mesh(:), Vel_mesh(:)];

fprintf('Computing PCE on %d x %d grid...\n', N_grid, N_grid);
Y_grid   = uq_evalModel(PCE_A, X_grid);
ReL_grid = reshape(Y_grid(:,259), size(Nu_mesh));
fprintf('Done.\n\n');


% ANALYSIS 1 — Critical curve Re(lambda) = 0

fprintf('=== ANALYSIS 1: Critical curve Re(lambda)=0 ===\n');

fig1 = figure('Position',[100 100 700 550]);
contourf(Nu_mesh, Vel_mesh, ReL_grid, 30, 'LineColor','none');
colorbar; colormap('jet'); caxis([-0.3 0.05]);
hold on
[C1, ~] = contour(Nu_mesh, Vel_mesh, ReL_grid, [0 0], 'k-', 'LineWidth', 3);
scatter(nu_vec, vel_vec, 25, 'k.', 'HandleVisibility','off');
hold off
xlabel('\nu','FontSize',13); ylabel('V','FontSize',13);
title('Critical Curve Re(\lambda)=0','FontSize',13);
exportgraphics(fig1,'A1_critical_curve.png','Resolution',300);


% Extract critical curve coordinates
contour_data = C1;
idx = 1;
nu_crit_pts  = [];
vel_crit_pts = [];
while idx < size(contour_data,2)
    n_pts = contour_data(2, idx);
    nu_crit_pts  = [nu_crit_pts,  contour_data(1, idx+1:idx+n_pts)];
    vel_crit_pts = [vel_crit_pts, contour_data(2, idx+1:idx+n_pts)];
    idx = idx + n_pts + 1;
end
fprintf('  Critical curve: %d points extracted\n\n', length(nu_crit_pts));

% ANALYSIS 2 — Effective Re_c along the critical curve

fprintf('=== ANALYSIS 2: Effective Re_c along critical curve ===\n');

L    = 2.0;  % characteristic length, lid-driven cavity
Re_c = vel_crit_pts * L ./ nu_crit_pts;

fprintf('  Re_c range: [%.1f, %.1f]\n', min(Re_c), max(Re_c));
fprintf('  Re_c mean:  %.1f\n',  mean(Re_c));
fprintf('  Re_c std:   %.1f\n\n', std(Re_c));

fig2 = figure('Position',[100 100 700 400]);
[nu_sort, sort_idx] = sort(nu_crit_pts);
plot(nu_sort, Re_c(sort_idx), 'b-o', 'LineWidth', 2, 'MarkerSize', 4);
xlabel('\nu','FontSize',13); ylabel('Re_c = VL/\nu','FontSize',13);
title('Effective Re_c along bifurcation curve','FontSize',13);
yline(mean(Re_c), 'r--', sprintf('mean = %.0f', mean(Re_c)), 'LineWidth',1.5,'FontSize',11);
grid on; box on;
exportgraphics(fig2,'A2_Re_critical.png','Resolution',300);

% ANALYSIS 3 — Probability of instability
% FIXED: sample (nu, V) from their TRUE joint distribution

fprintf('=== ANALYSIS 3: Probability of instability ===\n');

N_MC = 10000;
X_MC = uq_getSample(myInput, N_MC, 'MC');   % correct joint sampling
nu_MC  = X_MC(:,1);
vel_MC = X_MC(:,2);

Y_MC   = uq_evalModel(PCE_A, X_MC);
ReL_MC = Y_MC(:,259);

P_unstable = mean(ReL_MC > 0);
fprintf('  P(unstable) = %.4f  (%.2f%%)\n\n', P_unstable, P_unstable*100);

fig4 = figure('Position',[100 100 700 550]);
scatter(nu_MC(ReL_MC<=0),  vel_MC(ReL_MC<=0),  10, 'b.', 'DisplayName','Stable');
hold on
scatter(nu_MC(ReL_MC>0),   vel_MC(ReL_MC>0),   10, 'r.', 'DisplayName','Unstable');
hold off
xlabel('\nu','FontSize',13); ylabel('V','FontSize',13);
title(sprintf('P(unstable) = %.2f%% (correct sampling: nu ~ Lognormal, V ~ Uniform)', P_unstable*100),'FontSize',12);
legend('FontSize',11,'Location','northeast');
grid on; box on;
exportgraphics(fig4,'A4_probability.png','Resolution',300);


% ANALYSIS 4 — Classification accuracy (PCE vs CFD)

fprintf('=== ANALYSIS 9: Classification accuracy (PCE vs CFD) ===\n');

ReL_PCE_data = uq_evalModel(PCE_A, XA);
ReL_PCE_ReL  = ReL_PCE_data(:,259);

true_unstable = ReL_data     > 0;
pred_unstable = ReL_PCE_ReL  > 0;

TP = sum( pred_unstable &  true_unstable);   % correctly flagged unstable
TN = sum(~pred_unstable & ~true_unstable);   % correctly flagged stable
FP = sum( pred_unstable & ~true_unstable);   % false alarm
FN = sum(~pred_unstable &  true_unstable);   % missed instability (most costly)

accuracy  = (TP + TN) / length(ReL_data);
precision = TP / (TP + FP + eps);
recall    = TP / (TP + FN + eps);
F1        = 2 * precision * recall / (precision + recall + eps);

fprintf('  Confusion matrix (positive class = Unstable):\n');
fprintf('              Pred Unstable  Pred Stable\n');
fprintf('  True Unstable    %3d           %3d\n', TP, FN);
fprintf('  True Stable      %3d           %3d\n', FP, TN);
fprintf('  Accuracy:  %.4f\n', accuracy);
fprintf('  Precision: %.4f\n', precision);
fprintf('  Recall:    %.4f\n', recall);
fprintf('  F1 score:  %.4f\n\n', F1);

fig9 = figure('Position',[100 100 500 450]);
conf_mat = [TP, FN; FP, TN];
imagesc(conf_mat);
colormap([1 1 1; 0.27 0.51 0.71]);
textStrings = {'TP','FN';'FP','TN'};
for ii = 1:2
    for jj = 1:2
        text(jj, ii, sprintf('%s\n%d', textStrings{ii,jj}, conf_mat(ii,jj)), ...
            'HorizontalAlignment','center','FontSize',14,'FontWeight','bold');
    end
end
set(gca,'XTick',[1 2],'XTickLabel',{'Pred Unstable','Pred Stable'},...
        'YTick',[1 2],'YTickLabel',{'True Unstable','True Stable'},'FontSize',11);
title('Confusion Matrix (PCE vs CFD, positive class = Unstable)','FontSize',13);
exportgraphics(fig9,'A9_confusion.png','Resolution',300);



% SUMMARY

fprintf('====================================================\n');
fprintf('           BIFURCATION ANALYSIS SUMMARY\n');
fprintf('====================================================\n');
fprintf('  A1 - Critical curve extracted (%d pts)\n',   length(nu_crit_pts));
fprintf('  A2 - Re_c in [%.0f, %.0f], mean = %.0f\n',   min(Re_c), max(Re_c), mean(Re_c));
fprintf('  A4 - P(unstable) = %.2f%%\n',                P_unstable*100);
fprintf('  A9 - Accuracy = %.2f%%, F1 = %.4f\n',        accuracy*100, F1);

save('Bifurcation_analysis_results.mat', 'nu_crit_pts', 'vel_crit_pts', 'Re_c', 'P_unstable','accuracy', 'precision', 'recall', 'F1', 'conf_mat', '-v7.3');