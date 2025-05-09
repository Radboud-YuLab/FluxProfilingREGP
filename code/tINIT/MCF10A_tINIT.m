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

% Load in RNA-seq data::

rna_seq_input = readtable(strcat(root,'/data/rna_seq/MCF10A_tpm.csv'));

rna_seq_input.Properties.VariableNames{1} = 'gene';

% Preprocessing for RNA-seq: 

% add pre-processing results to arrayData structure
arrayData.deletedDeadEndRxns = deletedDeadEndRxns;
arrayData.taskReport = taskReport;
arrayData.essentialRxnMat = essentialRxnMat;
% Cutoff threshold set at 1 tpm:
arrayData.threshold = 1;

% Formatting of RNAseq input files:

arrayData.genes = rna_seq_input.gene;
arrayData.tissues = rna_seq_input.Properties.VariableNames(2);
arrayData.levels = table2array(rna_seq_input(:, 2));

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

MCF10A_model = INIT_output.model{1,1};
save(strcat(root,'/models/tINIT_output/MCF10A_model.mat'),'MCF10A_model');


