%% Make contextualized model for MCF10A RNA-seq
% Written by: Cyriel Huijer 
% Last updated: 28 March 2025

root = ''; % This should be the base directory of the github repo
cd(root);
setRavenSolver('gurobi');

%% 

load("models/Human-GEM.mat");
model = ihuman;

% Load task structure:
taskStruct = parseTaskList(strcat(root,'/data/tINIT/metabolicTasks_Essential.txt'));

% Run some preliminary steps that will allow us to skip some pre-processing
% steps in the tINIT algorithm, greatly reducing the overall run time.
[~,deletedDeadEndRxns] = simplifyModel(model,true,false,true,true,true);
cModel = removeReactions(model,deletedDeadEndRxns,false,true);
[taskReport, essentialRxnMat] = checkTasks(cModel,[],true,false,true,taskStruct);

% add pre-processing results to arrayData structure
arrayData.deletedDeadEndRxns = deletedDeadEndRxns;
arrayData.taskReport = taskReport;
arrayData.essentialRxnMat = essentialRxnMat;
% Cutoff threshold set at 1 tpm:
arrayData.threshold = 1;


rna_seq_data = readtable('data/rna_seq/nci_60_cell_lines_tpm.csv');
% Extract sample and tpm values
arrayData.tissues = rna_seq_data.Properties.VariableNames(2:end)';  % sample (tissue) names
arrayData.genes = rna_seq_data.gene;  % gene names
arrayData.levels = table2array(rna_seq_data(:, 2:end));  % gene TPM values
arrayData.threshold = 1;


% Run tINIT

n_models = numel(arrayData.tissues);
INIT_output = {};

model = ihuman; 
model = addBoundaryMets(model,false);
for i = 1:length(arrayData.tissues)
    % First try to run tINIT with shorter time limit. If it fails, then
    % try again with a longer time limit.

    try
        params.TimeLimit = 1000;
        init_model = getINITModel2(model, arrayData.tissues{i}, [], [], arrayData, [], true, [], true, true, taskStruct, params);
    catch
        params.TimeLimit = 5000;
        init_model = getINITModel2(model, arrayData.tissues{i}, [], [], arrayData, [], true, [], true, true, taskStruct, params);
    end
    % Assign results
    init_model.id = arrayData.tissues{i};
    INIT_output.id{i,1} = init_model.id;
    INIT_output.model{i,1} = init_model;
end  % End of for loop


save('models/tINIT_output/20250329_collect_models_RobinsonCellLines', 'collect_models_all');



