% Post-bifurcation case: Lognormal Distribution, μ=8500, σ=200

Input_per.Marginals.Type = 'Lognormal';
Input_per.Marginals.Moments = [8500 200];
Input_per.Marginals.Bounds = [8000 9000];
myInput_per = uq_createInput(Input_per);
uq_selectInput(myInput_per);

fprintf('Generating 30 LHS samples...\n');
X_LHS_HB2 = uq_getSample(30, 'LHS');

fprintf('Evaluating solver...\n');
try
    [Y_LHS_HB2, results_HB2] = ifiss_model_periodic(X_LHS_HB2);
    save('dati_periodic_hb2.mat', 'X_LHS_HB2', 'Y_LHS_HB2', 'results_HB2', '-v7.3');
    fprintf('✓ Solver evaluations saved\n');
catch ME
    fprintf('✗ Error: %s\n', ME.message);
    rethrow(ME);
end