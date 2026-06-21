setpath
rng(100,'twister')
uqlab

% Bootstrap PCE 1
% Input: Lognormal Distribution

load dati_LOG_4000_1000.mat

Input_1.Marginals.Type = 'Lognormal';
Input_1.Marginals.Moments = [4000 1000];
Input_1.Marginals.Bounds = [1 8000];
myInput_1 = uq_createInput(Input_1);
uq_selectInput(myInput_1);

MetaOpts_1b.Type = 'Metamodel';
MetaOpts_1b.MetaType = 'PCE';
MetaOpts_1b.Degree = 1:30;
MetaOpts_1b.Method = 'OLS';
MetaOpts_1b.ExpDesign.X = X_LHS_1;
MetaOpts_1b.ExpDesign.Y = Y_LHS_1;
MetaOpts_1b.Bootstrap.Replications = 100;
PCE_OLS_1b = uq_createModel(MetaOpts_1b);

save('PCE_OLS_LOG_1_bootstrap.mat', 'PCE_OLS_1b', '-v7.3');
clearvars

% Bootstrap PCE 2
% Input: Normal Distribution

load dati_NORM_4000_1000.mat

Input_6.Marginals.Type = 'Gaussian';
Input_6.Marginals.Parameters = [4000 1000];
Input_6.Marginals.Bounds = [1 8000];
myInput_6 = uq_createInput(Input_6);
uq_selectInput(myInput_6);

MetaOpts_6b.Type = 'Metamodel';
MetaOpts_6b.MetaType = 'PCE';
MetaOpts_6b.Degree = 1:30;
MetaOpts_6b.Method = 'OLS';
MetaOpts_6b.ExpDesign.X = X_LHS_6;
MetaOpts_6b.ExpDesign.Y = Y_LHS_6;
MetaOpts_6b.Bootstrap.Replications = 100;
PCE_OLS_6b = uq_createModel(MetaOpts_6b);

save('PCE_OLS_NORM_6_bootstrap.mat', 'PCE_OLS_6b', '-v7.3');
clearvars

