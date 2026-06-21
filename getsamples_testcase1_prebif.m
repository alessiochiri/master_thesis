% =========================================================================
% SCENARIO 1: Lognormale μ ALTO - μ=4000, σ=1000
% =========================================================================
fprintf('\n========== SCENARIO 1: Lognormale μ=4000, σ=1000 ==========\n');

Input_per.Marginals.Type = 'Lognormal';
Input_per.Marginals.Moments = [4000 1000];
Input_per.Marginals.Bounds = [1 8000];
myInput_per = uq_createInput(Input_per);
uq_selectInput(myInput_per);

fprintf('Generating 40 LHS samples...\n');
X_LHS_PER = uq_getSample(40, 'LHS');

fprintf('Evaluating solver...\n');
try
    Y_LHS_PER = ifiss_model(X_LHS_PER);
    save('dati_periodic.mat', 'X_LHS_PER', 'Y_LHS_PER', '-v7.3');
    fprintf('✓ Solver evaluations saved\n');
catch ME
    fprintf('✗ Error: %s\n', ME.message);
    rethrow(ME);
end

% =========================================================================
% SCENARIO 2: Normale BASE - μ=4000, σ=1000
% =========================================================================
fprintf('\n========== SCENARIO 6: Normale μ=4000, σ=1000 ==========\n');

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


