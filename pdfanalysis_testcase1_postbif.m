% PCE VISUALIZATION - PERIODIC BRANCH
load('dati_periodic_hb2_corrected.mat', 'X_LHS_HB2b', 'Y_LHS_HB2b');
load('PCE_OLS_HB2.mat',       'PCE_OLS_HB2');
uq_selectModel(PCE_OLS_HB2);

idx_QoI = [1, 2, 3, 4:12];
QoI_names = {'E', 'EN', '\omega', 'A_{max}(-0.8,-0.8)', 'A_{max}(0.8,-0.8)', 'A_{max}(-0.8, 0.8)', 'A_{max}(0.8, 0.8)', 'A_{max}(0,0)', 'A_{max}(0,0.5)', 'A_{max}(0,-0.5)',     'A_{max}(-0.5,0)',   'A_{max}(0.5,0)'};
n_QoI = length(idx_QoI);
N_MC  = 1e5;
n_LHS = size(X_LHS_HB2b, 1); 

% MC sampling on the PCE

fprintf('MC sampling: %d samples...\n', N_MC);
X_MC = uq_getSample(N_MC, 'MC');
Y_MC_all = uq_evalModel(PCE_OLS_HB2, X_MC);
Y_MC = Y_MC_all(:, idx_QoI);
Y_pct = prctile(Y_MC, [5 25 50 75 95], 1);

% OUTPUT PDF — KDE for each QoI
figure('Name','Output QoI PDF','Position',[100 100 1400 900]);
for q = 1:n_QoI
    subplot(3, 4, q);
    [f_q, xi_q] = ksdensity(Y_MC(:,q));
    fill([xi_q, fliplr(xi_q)], [f_q, zeros(1,length(f_q))], ...
        [0.4 0.7 0.5], 'FaceAlpha',0.4, 'EdgeColor',[0.1 0.4 0.2], 'LineWidth',1.5);
    hold on;
    xline(mean(Y_MC(:,q)),   '--b', '\mu', 'LabelVerticalAlignment','bottom', 'LineWidth',1.2);
    xline(median(Y_MC(:,q)), '--r', 'med.', 'LabelVerticalAlignment','bottom', 'LineWidth',1.2);
    xlabel(QoI_names{q}, 'Interpreter','tex');
    ylabel('PDF');
    title(QoI_names{q}, 'Interpreter','tex');
    grid on; box on;
end
sgtitle('Output PDFs', 'FontSize',12);

% SAMPLE TREND vs Re
[Re_sort, idx_s] = sort(X_LHS_HB2b(:,1));
Y_sort = Y_LHS_HB2b(idx_s, idx_QoI);
Re_grid    = linspace(8000, 8800, 500)';
Y_grid_all = uq_evalModel(PCE_OLS_HB2, Re_grid);
Y_grid     = Y_grid_all(:, idx_QoI);

figure('Name','QoI Trend vs Re','Position',[100 100 1400 900]);
for q = 1:n_QoI
    subplot(3, 4, q);
    hold on;
    % percentile bands
    yline(Y_pct(1,q), ':', 'Color',[0.6 0.6 0.9], 'LineWidth',1.0);   % P5
    yline(Y_pct(5,q), ':', 'Color',[0.6 0.6 0.9], 'LineWidth',1.0);   % P95
    yline(Y_pct(2,q), '--','Color',[0.3 0.3 0.8], 'LineWidth',1.0);   % P25
    yline(Y_pct(4,q), '--','Color',[0.3 0.3 0.8], 'LineWidth',1.0);   % P75
    % PCE curve and samples
    plot(Re_grid, Y_grid(:,q), 'b-', 'LineWidth',2);
    scatter(Re_sort, Y_sort(:,q), 30, 'k', 'filled', 'MarkerFaceAlpha',0.7);
    xlabel('Re');
    ylabel(QoI_names{q}, 'Interpreter','tex');
    title(QoI_names{q},  'Interpreter','tex');
    grid on; box on;
end
sgtitle('QoI vs Re — PCE curve + LHS samples','FontSize',12);

% Statistics
fprintf('\n%-25s  %12s  %12s  %12s  %12s\n', 'QoI','Mean','Std','P5','P95');
fprintf('%s\n', repmat('-',1,75));
for q = 1:n_QoI
    fprintf('%-25s  %12.4e  %12.4e  %12.4e  %12.4e\n', QoI_names{q}, mean(Y_MC(:,q)), std(Y_MC(:,q)), Y_pct(1,q), Y_pct(5,q)); 
end