% Test Case 2 — PCE with Bootstrap
% Bivariate case (nu, V)
clearvars
rng(100,'twister')
uqlab

InputA.Marginals(1).Type    = 'Lognormal'; % viscosity
InputA.Marginals(1).Moments = [0.002, 0.002];
InputA.Marginals(1).Bounds  = [0.0007, 0.01];

InputA.Marginals(2).Type       = 'Uniform'; % velocity
InputA.Marginals(2).Parameters = [1, 5];

load('dati_due_variabili_1.mat', 'XA1', 'YA1');
myInputA = uq_createInput(InputA);
uq_selectInput(myInputA);

MetaOpts_B.Type       = 'Metamodel';
MetaOpts_B.MetaType   = 'PCE';
MetaOpts_B.TruncOptions.qNorm = 0.1:0.1:1;
MetaOpts_B.qNormEarlyStop      = false;
MetaOpts_B.Degree     = 1:40;
MetaOpts_B.Method     = 'LARS';
MetaOpts_B.LARS.LarsEarlyStop = false;
MetaOpts_B.ExpDesign.X = XA1;
MetaOpts_B.ExpDesign.Y = YA1;

MetaOpts_B.Bootstrap.Replications = 100;

PCE_Bootstrap = uq_createModel(MetaOpts_B);
save('PCE_due_variabili_bootstrap.mat', 'PCE_Bootstrap', '-v7.3');

fprintf('Bootstrap PCE saved: PCE_due_variabili_bootstrap.mat\n');