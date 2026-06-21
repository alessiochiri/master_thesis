% PCE BOOTSTRAP — Post-bifurcation Regime
load('dati_periodic_hb2_corrected.mat', 'X_LHS_HB2b', 'Y_LHS_HB2b');

% Input
Input_per.Marginals.Type    = 'Lognormal';
Input_per.Marginals.Moments = [8500 200];
Input_per.Marginals.Bounds  = [8000 8800];
myInput_per = uq_createInput(Input_per);
uq_selectInput(myInput_per);

% Bootstrap PCE
MetaOpts_boot.Type                   = 'Metamodel';
MetaOpts_boot.MetaType               = 'PCE';
MetaOpts_boot.Degree                 = 1:15;
MetaOpts_boot.Method                 = 'OLS';
MetaOpts_boot.Bootstrap.Replications = 100;
MetaOpts_boot.ExpDesign.X            = X_LHS_HB2b;
MetaOpts_boot.ExpDesign.Y            = Y_LHS_HB2b;
PCE_boot = uq_createModel(MetaOpts_boot);
save('PCE_boot_HB2.mat', 'PCE_boot', '-v7.3');

% Visualizationn with confidence bands
Re_val = linspace(8000, 8800, 500)';
[Y_mean, Y_var, Y_boot] = uq_evalModel(PCE_boot, Re_val);

QoI_names = {'E', 'EN', '\omega', 'A_{max}(-0.8,-0.8)', 'A_{max}(0.8,-0.8)','A_{max}(-0.8,0.8)',  'A_{max}(0.8,0.8)','A_{max}(0,0)', 'A_{max}(0,0.5)', 'A_{max}(0,-0.5)',    'A_{max}(-0.5,0)',   'A_{max}(0.5,0)'};

[Re_sort, idx_s] = sort(X_LHS_HB2b(:,1));
Y_sort = Y_LHS_HB2b(idx_s, :);

figure('Position', [100 100 1400 900]);
for q = 1:12
    subplot(4,3,q); hold on

    % Bootstrap replications
    plot(Re_val, squeeze(Y_boot(:,:,q))', 'Color', [0.7 0.85 1.0 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off')

    % 95% bands — single legend entry
    q025 = quantile(squeeze(Y_boot(:,:,q)), 0.025, 2);
    q975 = quantile(squeeze(Y_boot(:,:,q)), 0.975, 2);
    fill([Re_val; flipud(Re_val)], [q025; flipud(q975)], [0.4 0.6 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', '95% CI')

    % PCE mean curve
    plot(Re_val, Y_mean(:,q), 'b-', 'LineWidth', 2, 'DisplayName', 'PCE mean')

    % LHS samples
    scatter(Re_sort, Y_sort(:,q), 30, 'k', 'filled', 'MarkerFaceAlpha', 0.7, 'DisplayName', 'LHS samples')

    xlabel('Re'); ylabel(QoI_names{q}, 'Interpreter', 'tex')
    title(QoI_names{q}, 'Interpreter', 'tex')
    legend('Location', 'best'); grid on; box on;
end
sgtitle('Bootstrap PCE', 'FontSize', 12);