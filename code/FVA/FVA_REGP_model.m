%% This script runs FVA for REGP-constrained MCF10A model 
% written by: Cyriel Huijer
% Last updated: 1 April 2025


clc
clear

% Set up working directory:
root = ''; % This should be the base directory of the github repo
cd(root);
%% Init Cobra and set solver:

initCobraToolbox;
changeCobraSolver('gurobi');
setRavenSolver('gurobi');

%%

% Load REGP constrained model:
load('models/GP_model_20perc_bounds.mat')

% GP model = REGP_constrained mmodel
GP_model = GP_model_20perc_bounds;

% add back ATP maintenance rxn:
rxnsToAdd.rxns = {'MAR03964'};  % Reaction ID
%rxnsToAdd.equations = {'MAM01371c + MAM02040c => MAM01285c + MAM02039c + MAM02751c'};  % Reaction equation
rxnsToAdd.equations = {'ATP[c] + H2O[c] => ADP[c] + H+[c] + Pi[c]'};  % Reaction equation
rxnsToAdd.lb = 0;  
rxnsToAdd.ub = 1000;  
rxnsToAdd.rxnNames = {'ATP_maintenance'};  
rxnsToAdd.subSystems = {'Artificial reactions'}; 
GP_model = addRxns(GP_model,rxnsToAdd,3);

% Set model to experimental GR (20% flexibility allowed)
GP_model = setParam(GP_model,'lb','MAR13082',0.0486*0.8);
GP_model = setParam(GP_model,'ub','MAR13082',0.0486*1.2);

% Run FVA
rxn_index = 1:length(GP_model.rxns);
rxn = GP_model.rxns;
GP_min = [];
GP_max = [];

for i = 1:length(GP_model.rxns)
    disp(i);
    GP_model.c(:) = 0;
    GP_model.c(i) = 1;
    GP_sol_min = optimizeCbModel(GP_model,'min');
    GP_sol_max = optimizeCbModel(GP_model,'max'); 
    GP_min = [GP_min,GP_sol_min.f];
    GP_max = [GP_max,GP_sol_max.f];
end

% Save FVA table
fva_table_GP_simulated_exch_fluxes_20_percent_bounds = table(rxn_index',rxn, GP_min',GP_max',...
    'VariableNames',{'rxn_index','rxn','min','max'});
% Filter out exchange rxns
[exch_rxns, ~] = getExchangeRxns(CORE_model);
is_exchange = ismember(fva_table_CORE_simulated_exch_fluxes_700_percent_bounds.rxn, exch_rxns);
fva_table_CORE_simulated_exch_fluxes_700_percent_bounds =fva_table_CORE_simulated_exch_fluxes_700_percent_bounds(~is_exchange, :);
fva_table_GP_simulated_exch_fluxes_20_percent_bounds = fva_table_GP_simulated_exch_fluxes_20_percent_bounds(~is_exchange, :);
% Also filter out transport rxns
fva_table_CORE_simulated_exch_fluxes_700_percent_bounds = fva_remove_subsystems(fva_table_CORE_simulated_exch_fluxes_700_percent_bounds,CORE_model,'Transport reactions')
fva_table_GP_simulated_exch_fluxes_20_percent_bounds = fva_remove_subsystems(fva_table_GP_simulated_exch_fluxes_20_percent_bounds,GP_model,'Transport reactions')

% Write into table for further analysis (used for heatmap Fig 4)
writetable(fva_table_GP_simulated_exch_fluxes_20_percent_bounds,'data/fva_heatmap/fva_table_GP_simulated_exch_fluxes_20_percent_bounds_noTsp_rxns.csv');
% Also save updated REGP model (now includes ATPM model)
save('models/GP_model_20perc_bounds_ATPM.mat');



