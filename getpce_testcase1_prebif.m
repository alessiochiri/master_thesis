rng(100,'twister')
uqlab

% PCE 1
% Input: Lognormal Distribution

load dati_LOG_4000_1000.mat

Input_1.Marginals.Type = 'Lognormal';
Input_1.Marginals.Moments = [4000 1000];
Input_1.Marginals.Bounds = [1 8000];
myInput_1 = uq_createInput(Input_1);
uq_selectInput(myInput_1);

MetaOpts_1.Type = 'Metamodel';
MetaOpts_1.MetaType = 'PCE';
MetaOpts_1.Degree = 1:30;
MetaOpts_1.Method = 'OLS';
MetaOpts_1.ExpDesign.X = X_LHS_1;
MetaOpts_1.ExpDesign.Y = Y_LHS_1;
PCE_OLS_1 = uq_createModel(MetaOpts_1);

save('PCE_OLS_LOG_4000_1000.mat', 'PCE_OLS_1', '-v7.3');

clearvars

% PCE 2
% Input: Normal Distribution

load dati_NORM_4000_1000.mat

Input_6.Marginals.Type = 'Gaussian';
Input_6.Marginals.Parameters = [4000 1000];
Input_6.Marginals.Bounds = [1 8000];
myInput_6 = uq_createInput(Input_6);
uq_selectInput(myInput_6);

MetaOpts_6.Type = 'Metamodel';
MetaOpts_6.MetaType = 'PCE';
MetaOpts_6.Degree = 1:30;
MetaOpts_6.Method = 'OLS';
MetaOpts_6.ExpDesign.X = X_LHS_6;
MetaOpts_6.ExpDesign.Y = Y_LHS_6;
PCE_OLS_6 = uq_createModel(MetaOpts_6);

save('PCE_OLS_NORM_4000_1000.mat', 'PCE_OLS_6', '-v7.3');

clearvars

