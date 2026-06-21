% LOO CONVERGENCE CURVE — Error vs Degree D
% Post-bifurcation case H

load('dati_periodic_hb2_corrected.mat', 'X_LHS_HB2b', 'Y_LHS_HB2b');

QoI_names = {'E', 'EN', '\omega', 'A_{max}(-0.8,-0.8)', 'A_{max}(0.8,-0.8)', 'A_{max}(-0.8,0.8)',  'A_{max}(0.8,0.8)', 'A_{max}(0,0)', 'A_{max}(0,0.5)', 'A_{max}(0,-0.5)',    'A_{max}(-0.5,0)',   'A_{max}(0.5,0)'};
n_out = 12;

% Input
Input_per.Marginals.Type    = 'Lognormal';
Input_per.Marginals.Moments = [8500 200];
Input_per.Marginals.Bounds  = [8000 8800];
myInput_per = uq_createInput(Input_per);
uq_selectInput(myInput_per);

% Degrees to explore
% N=24 samples → N=2D → Dmax=12, N=3D → Dmax=8
D_2N = 1:12;
D_3N = 1:8;

LOO_2D = nan(length(D_2N), n_out);
LOO_3D = nan(length(D_3N), n_out);


% N = 2D
for di = 1:length(D_2N)
    D = D_2N(di);
    N = 2 * D;
    fprintf('  D=%d, N=%d\n', D, N);

    [X_sub, idx_k] = uq_subsample(X_LHS_HB2b, N, 'k-means');
    Y_sub = Y_LHS_HB2b(idx_k, :);

    MetaOpts.Type        = 'Metamodel';
    MetaOpts.MetaType    = 'PCE';
    MetaOpts.Degree      = D;
    MetaOpts.Method      = 'OLS';
    MetaOpts.ExpDesign.X = X_sub;
    MetaOpts.ExpDesign.Y = Y_sub;

    try
        myPCE = uq_createModel(MetaOpts);
        for s = 1:n_out
            LOO_2D(di, s) = myPCE.Error(s).LOO;
        end
    catch ME
        fprintf('  ERROR D=%d N=%d: %s\n', D, N, ME.message);
    end

    clear MetaOpts myPCE X_sub Y_sub idx_k
end

% N = 3D
for di = 1:length(D_3N)
    D = D_3N(di);
    N = 3 * D;
    fprintf('  D=%d, N=%d\n', D, N);

    [X_sub, idx_k] = uq_subsample(X_LHS_HB2b, N, 'k-means');
    Y_sub = Y_LHS_HB2b(idx_k, :);

    MetaOpts.Type        = 'Metamodel';
    MetaOpts.MetaType    = 'PCE';
    MetaOpts.Degree      = D;
    MetaOpts.Method      = 'OLS';
    MetaOpts.ExpDesign.X = X_sub;
    MetaOpts.ExpDesign.Y = Y_sub;

    try
        myPCE = uq_createModel(MetaOpts);
        for s = 1:n_out
            LOO_3D(di, s) = myPCE.Error(s).LOO;
        end
    catch ME
        fprintf('  ERROR D=%d N=%d: %s\n', D, N, ME.message);
    end

    clear MetaOpts myPCE X_sub Y_sub idx_k
end

save('convergence_LOO_HB2.mat', 'LOO_2D', 'LOO_3D', 'D_2N', 'D_3N', 'QoI_names');
fprintf('Saved to convergence_LOO_HB2.mat\n');

% PLOT
figure('Position', [100 100 1600 700]);
for s = 1:n_out
    subplot(3, 4, s); hold on
    Q2_2D = 1 - LOO_2D(:, s);
    Q2_3D = 1 - LOO_3D(:, s);
    plot(D_2N, Q2_2D, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'N=2D');
    plot(D_3N, Q2_3D, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'N=3D');
    yline(0.95, 'k--', 'LineWidth', 1, 'DisplayName', 'Q^2=0.95');
    yline(0.90, 'k:',  'LineWidth', 1, 'DisplayName', 'Q^2=0.90');
    xlabel('Degree D'); ylabel('Q^2')
    title(QoI_names{s}, 'Interpreter', 'tex')
    legend('Location', 'southeast'); grid on; box on;
    ylim([-0.1 1.05])
end
sgtitle('Q^2 Convergence vs Degree', 'FontSize', 13);