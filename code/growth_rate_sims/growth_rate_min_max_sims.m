%% Growth rate simulation script
% Written by: Cyriel Huijer
% Last updated: 1 April 2025

% Used for determining the simulated growth rate for MCF10A REGP/CORE
% & NCI-60 models. Used for Figure 3.

% Set up working directory:
root = ''; % This should be the base directory of the github repo
cd(root);
%% Init Cobra and set solver:

initCobraToolbox;
changeCobraSolver('gurobi');
setRavenSolver('gurobi');

clear
clc

% New:
sim_exch_fluxes = readtable("data/growth_rate_simulations/simulated_exchange_fluxes.csv");
meas_exch_fluxes = readtable("data/exchange_fluxes/output/exch_fluxes_MCF10A_NCI60.csv");

% Calculate % eror rom 
sim_GP = sim_exch_fluxes.GP;
meas_GP = meas_exch_fluxes.REGP;

percent_error = abs(sim_GP - meas_GP) ./ abs(meas_GP) * 100;
percent_error_table = table(sim_exch_fluxes.Row, percent_error, ...
                            'VariableNames', {'metabolite', 'Percent_Error'});

% Valine, and methionine are off by more than 10%, therefore,
% these will be dropped from the constraints

mets_to_keep = percent_error_table.Percent_Error <= 10;

% Apply filter to both tables
meas_exch_fluxes = meas_exch_fluxes(mets_to_keep, :);
sim_exch_fluxes = sim_exch_fluxes(mets_to_keep, :);

clear sim_GP  percent_error_table percent_error mets_to_keep meas_GP sim_exch_fluxes root

%% Constrain metabolic models with media components + measured exchange fluxes:
% No flexibilazion factor included 

growth_rate_data = readtable('data/cell_count/growth_rate_data.csv');
% Load in MCF10A metabolic model:
load('models/tINIT_output/MCF10A_model.mat');
load('models/tINIT_output/20250329_collect_models_RobinsonCellLines.mat');

conditions = {'CORE','REGP','HS_578T','SR','NCI_H226','UO_31','MALME_3M','HT29','X786_O','MDA_MB_231_ATCC','RPMI_8226','HOP_92','HOP_62'};
    
models = {EV_0_thr1_tINIT_model, EV_0_thr1_tINIT_model, collect_models_all{:}};
% Change b to make it compatible with cobra's optimizeCbModel
for i = 1:length(models)
    models{i}.b = models{i}.b(:,1);
end

% Constrain models with components in biomass using setHamsMedium (taken
% from Robinson et al (PMID: 32209698) 
for i = 1:length(models)
    disp(conditions{i});
    models{i} = setHamsMedium_CH(models{i});
    biomass_rxn_idx = find(strcmp(models{i}.rxns,'MAR13082')==1);
    models{i}.c(biomass_rxn_idx) =1;
    sol = solveLP(models{i},'max');
    disp(sol.f);
end

% Allows no flexibility in exch. flux constraints
numeric_data = meas_exch_fluxes{:, 2:end};
lb = meas_exch_fluxes;
ub = meas_exch_fluxes;
lb{:, 2:end} = numeric_data - abs((numeric_data * 0));
ub{:, 2:end} = numeric_data + abs((numeric_data * 0));

% optimize biomass rxn
max_solutions = [];
min_solutions = [];
for i = 1:length(conditions)
    cond = conditions{i};
    disp(cond);
    % Set exchange bounds for each condition
    models{i} = setExchangeBounds(models{i}, ...
                                  meas_exch_fluxes.metabolite, ...
                                  lb.(cond), ...
                                  ub.(cond), ...
                                  false);
    sol_max = optimizeCbModel(models{i},'max');
    sol_min = optimizeCbModel(models{i},'min');
    disp(sol_max.f);
    disp(sol_min.f)
    % if the solution is infeasible
    if sol_max.stat == 0
        max_solutions = [max_solutions,sol_max.f];
    end
    if sol_max.f >= 0
        max_solutions = [max_solutions,sol_max.f];
    end
    if sol_min.stat == 0
        min_solutions = [min_solutions,sol_min.f];
        disp("Infeasible");
    end
    if sol_min.f >= 0
        min_solutions = [min_solutions,sol_min.f];
    end
    models{i}.id = cond;
end
growth_rates = table2array(growth_rate_data(1,:));

measured_fluxes_0_bounds = table(conditions', min_solutions', max_solutions', growth_rates', ...
                     'VariableNames', {'Condition', 'MinSolution', 'MaxSolution', 'GrowthRate'});

writetable(measured_fluxes_0_bounds, 'data/growth_rate_simulations/measured_fluxes_0_bounds.csv');

clear biomass_rxn_idx cond min_solutions max_solutions numeric_data i lb ub sol sol_max sol_min growth_rates measured_fluxes_0_bounds

%% 20% flexibilization included

models = {EV_0_thr1_tINIT_model, EV_0_thr1_tINIT_model, collect_models_all{:}};
% Change b to make it compatible with cobra's optimizeCbModel
for i = 1:length(models)
    models{i}.b = models{i}.b(:,1);
end

% Constrain models with components in biomass using setHamsMedium (taken
% from Robinson et al (PMID: 32209698) 
for i = 1:length(models)
    disp(conditions{i});
    models{i} = setHamsMedium_CH(models{i});
    biomass_rxn_idx = find(strcmp(models{i}.rxns,'MAR13082')==1);
    models{i}.c(biomass_rxn_idx) =1;
    sol = solveLP(models{i},'max');
    disp(sol.f);
end

% Allows no flexibility in exch. flux constraints
numeric_data = meas_exch_fluxes{:, 2:end};
lb = meas_exch_fluxes;
ub = meas_exch_fluxes;
lb{:, 2:end} = numeric_data - abs((numeric_data * 0.2));
ub{:, 2:end} = numeric_data + abs((numeric_data * 0.2));

% optimize biomass rxn
max_solutions = [];
min_solutions = [];
for i = 1:length(conditions)
    cond = conditions{i};
    disp(cond);
    % Set exchange bounds for each condition
    models{i} = setExchangeBounds(models{i}, ...
                                  meas_exch_fluxes.metabolite, ...
                                  lb.(cond), ...
                                  ub.(cond), ...
                                  false);
    sol_max = optimizeCbModel(models{i},'max');
    sol_min = optimizeCbModel(models{i},'min');
    disp(sol_max.f);
    disp(sol_min.f)
    % if the solution is infeasible
    if sol_max.stat == 0
        max_solutions = [max_solutions,sol_max.f];
    end
    if sol_max.f >= 0
        max_solutions = [max_solutions,sol_max.f];
    end
    if sol_min.stat == 0
        min_solutions = [min_solutions,sol_min.f];
        disp("Infeasible");
    end
    if sol_min.f >= 0
        min_solutions = [min_solutions,sol_min.f];
    end
    models{i}.id = cond;
end
growth_rates = table2array(growth_rate_data(1,:));

measured_fluxes_20_bounds = table(conditions', min_solutions', max_solutions', growth_rates', ...
                     'VariableNames', {'Condition', 'MinSolution', 'MaxSolution', 'GrowthRate'});

writetable(measured_fluxes_20_bounds, 'data/growth_rate_simulations/measured_fluxes_20_bounds.csv');

% Save MCF10A REGP-constrained model for FVA:
GP_model_20perc_bounds = models{2};
save('models/GP_model_20perc_bounds.mat');

clear biomass_rxn_idx cond min_solutions max_solutions numeric_data i lb ub sol sol_max sol_min growth_rates measured_fluxes_20_bounds

%% 700% flexibilization included (REGP model becomes feasible)

models = {EV_0_thr1_tINIT_model, EV_0_thr1_tINIT_model, collect_models_all{:}};
% Change b to make it compatible with cobra's optimizeCbModel
for i = 1:length(models)
    models{i}.b = models{i}.b(:,1);
end

% Constrain models with components in biomass using setHamsMedium (taken
% from Robinson et al (PMID: 32209698) 
for i = 1:length(models)
    disp(conditions{i});
    models{i} = setHamsMedium_CH(models{i});
    biomass_rxn_idx = find(strcmp(models{i}.rxns,'MAR13082')==1);
    models{i}.c(biomass_rxn_idx) =1;
    sol = solveLP(models{i},'max');
    disp(sol.f);
end

% Allows no flexibility in exch. flux constraints
numeric_data = meas_exch_fluxes{:, 2:end};
lb = meas_exch_fluxes;
ub = meas_exch_fluxes;
lb{:, 2:end} = numeric_data - abs((numeric_data * 7));
ub{:, 2:end} = numeric_data + abs((numeric_data * 7));

% optimize biomass rxn
max_solutions = [];
min_solutions = [];
for i = 1:length(conditions)
    cond = conditions{i};
    disp(cond);
    % Set exchange bounds for each condition
    models{i} = setExchangeBounds(models{i}, ...
                                  meas_exch_fluxes.metabolite, ...
                                  lb.(cond), ...
                                  ub.(cond), ...
                                  false);
    sol_max = optimizeCbModel(models{i},'max');
    sol_min = optimizeCbModel(models{i},'min');
    disp(sol_max.f);
    disp(sol_min.f)
    % if the solution is infeasible
    if sol_max.stat == 0
        max_solutions = [max_solutions,sol_max.f];
    end
    if sol_max.f >= 0
        max_solutions = [max_solutions,sol_max.f];
    end
    if sol_min.stat == 0
        min_solutions = [min_solutions,sol_min.f];
        disp("Infeasible");
    end
    if sol_min.f >= 0
        min_solutions = [min_solutions,sol_min.f];
    end
    models{i}.id = cond;
end
growth_rates = table2array(growth_rate_data(1,:));

measured_fluxes_700_bounds = table(conditions', min_solutions', max_solutions', growth_rates', ...
                     'VariableNames', {'Condition', 'MinSolution', 'MaxSolution', 'GrowthRate'});

writetable(measured_fluxes_700_bounds, 'data/growth_rate_simulations/measured_fluxes_700_bounds.csv');


