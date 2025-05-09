%% This script runs FVA for REGP-constrained MCF10A model 
% written by: Cyriel Huijer
% Last updated: 1 April 2025


clc
clear

% Set up working directory:
root = ''; % This should be the base directory of the github repo
cd(root);

%% Sliding window + calculation of z-score:

fva_table_GP = readtable('data/fva_heatmap/fva_table_GP_simulated_exch_fluxes_20_percent_bounds_noTsp_rxns.csv');
load('models/GP_model_20perc_bounds_ATPM.mat')

% Calculate difference and sort:
fva_table_GP = fva_table_GP(fva_table_GP.min ~= fva_table_GP.max, :);
fva_table_GP.diff = fva_table_GP.max - fva_table_GP.min;
fva_table_GP = sortrows(fva_table_GP,'diff');
fva_table_GP.subsystem = GP_model.subSystems(fva_table_GP.rxn_index);
subsystems_char = cellfun(@char, fva_table_GP.subsystem, 'UniformOutput', false);
fva_table_GP = fva_table_GP(~ismember(subsystems_char, 'Transport reactions'), :);
fva_table_GP = fva_table_GP(fva_table_GP.diff > 0.1, :);

% Save supplemental table 3:

sup3 = fva_table_GP;
sup3 = removevars(sup3, 'Var6');
writetable(sup3, 'data/fva_heatmap/supplemental_table3.csv');

%% Prepare heatmap: 

window_size = 200; % windows of 200 rxns 
step_size = 10; % Sliding window with step size of 10

num_windows = floor((height(fva_table_GP) - window_size) / step_size) + 1;
reactions_in_windows = cell(num_windows, 1);

% Slide the window and extract reactions
for i = 1:num_windows
    % Determine the starting and ending index for the current window
    start_idx = (i - 1) * step_size + 1;
    end_idx = start_idx + window_size - 1;
    
    % Extract the reactions in the current window
    reactions_in_windows{i} = fva_table_GP.rxn(start_idx:end_idx);
end

% Loop over each window
for i = 1:num_windows
    % Get the reactions in the current window
    current_reactions = reactions_in_windows{i};
    
    % Find the indices of these reactions in the GP_model
    [~, idx_in_model] = ismember(current_reactions, GP_model.rxns);
    
    % Extract the subsystems for the reactions in the current window
    subsystems_in_windows{i} = GP_model.subSystems(idx_in_model(idx_in_model > 0));
end

rxn_subsystems = vertcat(GP_model.subSystems{:});
% Convert cell array to categorical array to count occurrences of each subsystem
rxn_subsystems_cat = categorical(rxn_subsystems);

% Count occurrences of each unique subsystem
subsystem_counts = countcats(rxn_subsystems_cat);

% Get the unique subsystem names
unique_subsystems = categories(rxn_subsystems_cat);

% Each subsystem should have more than 5 reactions in the model
subsystem_counts_table = table(unique_subsystems(:), subsystem_counts(:), 'VariableNames', {'Subsystem', 'Count'});
subsystem_counts_table = subsystem_counts_table(subsystem_counts_table.Count > 5, :);
subsystems_to_check = subsystem_counts_table.Subsystem;

% store the count for each window and subsystem
window_counts = zeros(length(subsystems_in_windows), height(subsystem_counts_table));

% Loop over the total number of windows, then loop over the subsystems in
% that window. In the third for loop check for the subsystem 
for i = 1:num_windows
    window = subsystems_in_windows{i};
    for j = 1:length(window)
        subsystem = char(window{j});
        for k = 1:length(subsystems_to_check)
            if strcmp(subsystem, subsystems_to_check{k})
                window_counts(i,k) = window_counts(i,k) + 1;
            end
        end
    end
end 

% Now calculate the z-score for every window.  
mu = mean(window_counts,1); % mean
sigma = std(window_counts,0,1); % standard deviation
valid_subsystems = mu ~= 0; % Check if the subsystem is actually in the dataset

% filter rows from the matrix
window_counts = window_counts(:, valid_subsystems);
subsystems_to_check = subsystems_to_check(valid_subsystems);
% Filter mean and std:
mu = mu(valid_subsystems);
sigma = sigma(valid_subsystems);

% Compute Z-score for each window and subsystem
z_scores = (window_counts - mu) ./ sigma;

% Compute the mean window index weighted by Z-score to determine enrichment trends
num_windows = size(window_counts, 1);
window_indices = (1:num_windows)';  % Column vector of window numbers

% Prepare colnames for export to R:
col_names = matlab.lang.makeValidName(subsystems_to_check);
z_scores_table = array2table(z_scores, 'VariableNames', col_names);

% Export table as csv, heatmap will be prepared in R:
writetable(z_scores_table, 'data/heatmap/z_scores_heatmap_20percent_GP.csv');




