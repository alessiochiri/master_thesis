% Test Case 2: bivariate
clearvars
rng(100,'twister')
uqlab

InputA.Marginals(1).Type = 'Lognormal'; % viscosity
InputA.Marginals(1).Moments = [0.002,0.002];
InputA.Marginals(1).Bounds = [0.0007,0.01];
InputA.Marginals(2).Type = 'Uniform'; % velocity
InputA.Marginals(2).Parameters = [1,5];


myInputA = uq_createInput(InputA);
uq_selectInput(myInputA);
uq_display(myInputA);

fprintf('Generating 100 LHS samples...\n');
XA = uq_getSample(100, 'LHS');

L = 2;
Re_samples = (XA(:,2) * L) ./ XA(:,1);
fprintf('\n=== RE COVERAGE ===\n');
fprintf('Re range: [%.2f, %.2f]\n', min(Re_samples), max(Re_samples));

try
    YA1 = ifiss_model_parametric(XA1);
    save('dati_due_variabili_1.mat', 'XA1', 'YA1', '-v7.3');
    fprintf('✓ Solver evaluations saved\n');
catch ME
    fprintf('✗ Error: %s\n', ME.message);
    rethrow(ME);
end


