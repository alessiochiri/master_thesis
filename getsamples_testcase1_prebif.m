% =========================================================================
% Case 1: Lognormal Distribution, μ=4000, σ=1000
% =========================================================================
Input_1.Marginals.Type = 'Lognormal';
Input_1.Marginals.Moments = [4000 1000];
Input_1.Marginals.Bounds = [1 8000];
myInput_1 = uq_createInput(Input_1);
uq_selectInput(myInput_1);

fprintf('Generating 40 LHS samples...\n');
X_LHS_1 = uq_getSample(40, 'LHS');

fprintf('Evaluating solver...\n');
try
    Y_LHS_1 = ifiss_model(X_LHS_1);
    save('dati_LOG_4000_1000.mat', 'X_LHS_1', 'Y_LHS_1', '-v7.3');
    fprintf('✓ Solver evaluations saved\n');
catch ME
    fprintf('✗ Error: %s\n', ME.message);
    rethrow(ME);

% =========================================================================
% Case 2: Normal Distribution, μ=4000, σ=1000
% =========================================================================
Input_6.Marginals.Type = 'Gaussian';
Input_6.Marginals.Parameters = [4000 1000];
Input_6.Marginals.Bounds = [1 8000];
myInput_6 = uq_createInput(Input_6);
uq_selectInput(myInput_6);

fprintf('Generating 40 LHS samples...\n');
X_LHS_6 = uq_getSample(40, 'LHS');

fprintf('Evaluating solver...\n');
try
    Y_LHS_6 = ifiss_model(X_LHS_6);
    save('dati_NORM_4000_1000.mat', 'X_LHS_6', 'Y_LHS_6', '-v7.3');
    fprintf('✓ Solver evaluations saved\n');
catch ME
    fprintf('✗ Error: %s\n', ME.message);
    rethrow(ME);
end


