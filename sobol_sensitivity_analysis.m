% SOBOL' SENSITIVITY ANALYSIS
% Bivariate case (nu, V) — 3 scalar outputs + 6 velocity nodes

clearvars
rng(100,'twister')
uqlab

load('PCE_due_variabili.mat', 'PCE_A');
load('dati_due_variabili_1.mat', 'XA1', 'YA1');

XA = XA1;
YA = YA1;

% -------------------------------------------------------
% Output indices

idx_ReL  = 259;
idx_KE   = 261;
idx_Enst = 262;

idx_nodes = [10, 65, 120, 139, 195, 248];
node_names = {'u - node 10', 'u - node 65', 'u - node 120', 'v - node 139', 'v - node 195', 'v - node 248'};

output_idx = [idx_ReL, idx_KE, idx_Enst, idx_nodes];
names = [{'Re(\lambda)', 'KE', 'Enstrophy'}, node_names];
n_out = length(output_idx);

fprintf('Creating input distribution...\n');

InputA.Marginals(1).Name = 'nu';
InputA.Marginals(1).Type = 'Lognormal';
InputA.Marginals(1).Moments = [0.002, 0.002];
InputA.Marginals(1).Bounds = [0.0007, 0.01];

InputA.Marginals(2).Name = 'V';
InputA.Marginals(2).Type = 'Uniform';
InputA.Marginals(2).Parameters = [1, 5];

myInput = uq_createInput(InputA);
uq_selectInput(myInput);

% -------------------------------------------------------
% PCE creation 

YA_reduced = YA(:, output_idx);

MetaOpts = struct();
MetaOpts.Type = 'Metamodel';
MetaOpts.MetaType = 'PCE';
MetaOpts.Method = 'LARS';
MetaOpts.Degree = 1:40;
MetaOpts.TruncOptions.qNorm = 0.1:0.1:1;
MetaOpts.ExpDesign.X = XA;
MetaOpts.ExpDesign.Y = YA_reduced;

myPCE_reduced = uq_createModel(MetaOpts);
uq_selectModel(myPCE_reduced);

% SOBOL ANALYSIS
SobolOpts.Type = 'Sensitivity';
SobolOpts.Method = 'Sobol';

try
    mySobolAnalysis = uq_createAnalysis(SobolOpts);
    fprintf('Sobol analysis completed successfully.\n\n');
    success = true;
catch ME
    fprintf('Error: %s\n\n', ME.message);
    success = false;
end

if success
    mySobolResults = mySobolAnalysis.Results;

    FirstOrder = mySobolResults.FirstOrder;
    Total = mySobolResults.Total;
    VariableNames = mySobolResults.VariableNames;

    if size(FirstOrder, 1) < size(FirstOrder, 2)
        FirstOrder = FirstOrder';
        Total = Total';
    end

    % ---- Print results ----
    fprintf('SOBOL INDICES - FIRST ORDER (S_i)\n');
    fprintf('%-16s %-15s %-15s\n', 'Output', 'nu', 'V');
    for i = 1:n_out
        fprintf('%-16s %-15.6f %-15.6f\n', names{i}, FirstOrder(i,1), FirstOrder(i,2));
    end

    fprintf('\nSOBOL INDICES - TOTAL ORDER (S_i^T)\n');
    fprintf('%-16s %-15s %-15s\n', 'Output', 'nu', 'V');
    for i = 1:n_out
        fprintf('%-16s %-15.6f %-15.6f\n', names{i}, Total(i,1), Total(i,2));
    end

    fprintf('\nINTERACTIONS (Total - FirstOrder)\n');
    fprintf('%-16s %-15s %-15s\n', 'Output', 'nu', 'V');
    for i = 1:n_out
        inter_nu = Total(i,1) - FirstOrder(i,1);
        inter_v  = Total(i,2) - FirstOrder(i,2);
        fprintf('%-16s %-15.6f %-15.6f\n', names{i}, inter_nu, inter_v);
    end
    fprintf('\n');

    fprintf('OFFICIAL UQLab OUTPUT\n');
    uq_print(mySobolAnalysis);

    % FIGURE 1 — First-Order indices
    
    fig = figure('Position', [100 100 1200 650]);
    x_pos = [1, 2];
    width = 0.12;
    colors_bar = lines(n_out);

    hold on
    for i = 1:n_out
        bar(x_pos(1) + (i-1)*width*2, FirstOrder(i,1), width, 'FaceColor', colors_bar(i,:), 'DisplayName', names{i});
        bar(x_pos(2) + (i-1)*width*2, FirstOrder(i,2), width, 'FaceColor', colors_bar(i,:));
    end
    hold off

    xlabel('Input Variables', 'FontSize', 12)
    ylabel('First-Order Sobol Index (S_i)', 'FontSize', 12)
    title('First-Order Sobol Indices (S_i)', 'FontSize', 13)
    set(gca, 'XTick', [1, 2], 'XTickLabel', VariableNames, 'FontSize', 11)
    legend(names, 'FontSize', 9, 'Location', 'eastoutside')
    grid on; ylim([0, 1]); box on

    exportgraphics(fig, 'Sobol_FirstOrder.png', 'Resolution', 300);
    fprintf('Figure saved: Sobol_FirstOrder.png\n');
    close(fig);


    % FIGURE 2 — Total-Order indices

    fig = figure('Position', [100 100 1200 650]);
    hold on
    for i = 1:n_out
        bar(x_pos(1) + (i-1)*width*2, Total(i,1), width, 'FaceColor', colors_bar(i,:), 'DisplayName', names{i});
        bar(x_pos(2) + (i-1)*width*2, Total(i,2), width, 'FaceColor', colors_bar(i,:));
    end
    hold off

    xlabel('Input Variables', 'FontSize', 12)
    ylabel('Total-Order Sobol Index (S_i^T)', 'FontSize', 12)
    title('Total-Order Sobol Indices (S_i^T)', 'FontSize', 13)
    set(gca, 'XTick', [1, 2], 'XTickLabel', VariableNames, 'FontSize', 11)
    legend(names, 'FontSize', 9, 'Location', 'eastoutside')
    grid on; ylim([0, 1]); box on

    exportgraphics(fig, 'Sobol_TotalOrder.png', 'Resolution', 300);
    fprintf('Figure saved: Sobol_TotalOrder.png\n');
    close(fig);

    % FIGURE 3 — Si vs Si^T comparison (n_out subplots)

    n_cols = 3;
    n_rows = ceil(n_out / n_cols);
    fig = figure('Position', [100 100 1600 350*n_rows]);

    for i = 1:n_out
        subplot(n_rows, n_cols, i)

        x_compare = [1, 2, 3.5, 4.5];
        values = [FirstOrder(i,1), FirstOrder(i,2), Total(i,1), Total(i,2)];
        colors_compare = [0.2 0.6 1; 0.2 0.8 1; 1 0.5 0.2; 1 0.7 0.4];

        hold on
        for j = 1:4
            bar(x_compare(j), values(j), 0.6, 'FaceColor', colors_compare(j,:), ...
                'EdgeColor', 'black', 'LineWidth', 1.5);
        end
        hold off

        xline(2.75, 'k--', 'LineWidth', 1.5);
        set(gca, 'XTickLabel', {'S_i(\nu)', 'S_i(V)', 'S_i^T(\nu)', 'S_i^T(V)'}, ...
                 'XTick', x_compare, 'FontSize', 9)
        ylabel('Sensitivity Index', 'FontSize', 10)
        title(names{i}, 'FontSize', 10, 'FontWeight', 'bold')
        ylim([0, 1]); grid on; box on

        inter_nu = Total(i,1) - FirstOrder(i,1);
        inter_v  = Total(i,2) - FirstOrder(i,2);
        text(0.05, 0.95, sprintf('Int. \\nu: %.3f\nInt. V: %.3f', inter_nu, inter_v), ...
            'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'BackgroundColor', 'white', 'FontSize', 8, 'EdgeColor', 'black', 'Margin', 3);
    end

    sgtitle('Sobol Indices: First Order (S_i) vs Total Order (S_i^T)', 'FontSize', 13)
    exportgraphics(fig, 'Sobol_Comparison.png', 'Resolution', 300);
    fprintf('Figure saved: Sobol_Comparison.png\n');
    close(fig);



    % FIGURE 4 — Variance decomposition pie charts
     fig = figure('Position', [100 100 1600 350*n_rows]);
    colors_pie = [0.2 0.6 1; 0.8 0.2 0.2; 0.8 0.8 0.8];

    for i = 1:n_out
        subplot(n_rows, n_cols, i)

        S1_nu = FirstOrder(i,1);
        S1_v  = FirstOrder(i,2);
        S_residual = max(0, 1 - S1_nu - S1_v);

        values = [S1_nu, S1_v, S_residual];
        labels = {sprintf('\\nu\n(%.3f)', S1_nu), ...
                  sprintf('V\n(%.3f)', S1_v), ...
                  sprintf('Higher\n(%.3f)', S_residual)};

        eps_min = 1e-6;
        values_safe = max(values, eps_min);

        h = pie(values_safe, labels);


        patchHandles = findobj(h, 'Type', 'Patch');
        patchHandles = flipud(patchHandles(:));

        for j = 1:min(3, numel(patchHandles))
            set(patchHandles(j), 'FaceColor', colors_pie(j,:));
        end

        title(names{i}, 'FontSize', 10, 'FontWeight', 'bold')
    end

    sgtitle('Variance Decomposition (First-Order Sobol Indices)', 'FontSize', 13)
    exportgraphics(fig, 'Sobol_VarianceDecomposition.png', 'Resolution', 300);
    fprintf('Figure saved: Sobol_VarianceDecomposition.png\n') 

    % Save results, including reduced-PCE diagnostics
    save('Sobol_analysis_results.mat', 'mySobolAnalysis', 'FirstOrder', 'Total','VariableNames', 'names', 'output_idx','Q2_reduced', 'Degree_reduced', 'qNorm_reduced', '-v7.3');

else
    fprintf('Sobol analysis failed. Check the PCE reduced model.\n');
end