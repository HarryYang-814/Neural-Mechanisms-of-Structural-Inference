
clear; clc; close all;
addpath 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\codes';
%% ------------------------------------------------------------------------
% SECTION 0 · DATASETA CREATION
% -------------------------------------------------------------------------
% NPVP
folder_name1 = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\8_by_conditions\oscillation\VP\';
files_info1 = dir(fullfile(folder_name1,'*.set'));
files_name1 = {files_info1.name};
folder_name2 = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\8_by_conditions\oscillation\NP\';
files_info2 = dir(fullfile(folder_name2,'*.set'));
files_name2 = {files_info2.name}; 

pac_NPVP{1,1} = 'NP';
pac_NPVP{1,2} = 'VP';
for subj = 1:length(files_name1)
    restoredefaultpath; 
    addpath(genpath('E:\Toolboxes_for_matlab\eeglab2024.2\'));
    EEG_VP = pop_loadset('filepath', folder_name1, 'filename', files_name1{1,subj});
    EEG_NP = pop_loadset('filepath', folder_name2, 'filename', files_name2{1,subj});
    restoredefaultpath;
    addpath (genpath('E:\Toolboxes_for_matlab\fieldtrip-20250915\fieldtrip-20250915'));
    data_VP = eeglab2fieldtrip( EEG_VP, 'preprocessing');
    data_NP = eeglab2fieldtrip( EEG_NP, 'preprocessing');

    pac_NPVP{subj+1,1} = data_NP;
    pac_NPVP{subj+1,2} = data_VP;
end

% VPCC
folder_name1 = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\8_by_conditions\oscillation\VP\';
files_info1 = dir(fullfile(folder_name1,'*.set'));
files_name1 = {files_info1.name};
folder_name2 = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\8_by_conditions\oscillation\control\';
files_info2 = dir(fullfile(folder_name2,'*.set'));
files_name2 = {files_info2.name}; 

pac_VPCC{1,1} = 'VP';
pac_VPCC{1,2} = 'control';
for subj = 1:length(files_name1)
    restoredefaultpath; 
    addpath(genpath('E:\Toolboxes_for_matlab\eeglab2024.2\'));
    EEG_VP = pop_loadset('filepath', folder_name1, 'filename', files_name1{1,subj});
    EEG_control = pop_loadset('filepath', folder_name2, 'filename', files_name2{1,subj});
    restoredefaultpath;
    addpath (genpath('E:\Toolboxes_for_matlab\fieldtrip-20250915\fieldtrip-20250915'));
    data_VP = eeglab2fieldtrip( EEG_VP, 'preprocessing');
    data_control = eeglab2fieldtrip( EEG_control, 'preprocessing');

    pac_VPCC{subj+1,1} = data_VP;
    pac_VPCC{subj+1,2} = data_control;
end

% NPCC
folder_name1 = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\8_by_conditions\oscillation\NP\';
files_info1 = dir(fullfile(folder_name1,'*.set'));
files_name1 = {files_info1.name};
folder_name2 = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\8_by_conditions\oscillation\control\';
files_info2 = dir(fullfile(folder_name2,'*.set'));
files_name2 = {files_info2.name}; 

pac_NPCC{1,1} = 'NP';
pac_NPCC{1,2} = 'control';
for subj = 1:length(files_name1)
    restoredefaultpath; 
    addpath(genpath('E:\Toolboxes_for_matlab\eeglab2024.2\'));
    EEG_NP = pop_loadset('filepath', folder_name1, 'filename', files_name1{1,subj});
    EEG_control = pop_loadset('filepath', folder_name2, 'filename', files_name2{1,subj});
    restoredefaultpath;
    addpath (genpath('E:\Toolboxes_for_matlab\fieldtrip-20250915\fieldtrip-20250915'));
    data_NP = eeglab2fieldtrip( EEG_NP, 'preprocessing');
    data_control = eeglab2fieldtrip( EEG_control, 'preprocessing');

    pac_NPCC{subj+1,1} = data_NP;
    pac_NPCC{subj+1,2} = data_control;
end

save pac_datasets.mat pac_NPVP pac_VPCC pac_NPCC

%% ------------------------------------------------------------------------
% SECTION 1 · USER PARAMETERS · EDIT ME!                                   
% -------------------------------------------------------------------------

cfg.srate        = 500;          % Sampling rate (Hz)

% Frequency bands of interest (Hz)
cfg.phaseBands   = struct( ...
    'delta', [1  3.5]);

cfg.ampBands     = struct( ...
    'theta', [4  7]);

% PAC pairs to compute   {phaseName, ampName}
cfg.pacPairs      = { 'delta', 'theta'};

% Latency window for mean PAC (seconds, relative to trial time axis)
cfg.latency       = [0.5 0.8];    % e.g. 200–800 ms after stimulus onset

% Channels to include (numeric indices or name cell array). Use [] for all.
cfg.chanSelect    = [];

% Permutation test settings
cfg.nPerm         = 10000;         % # shuffles
cfg.alpha         = 0.05;         % two‑tailed family‑wise alpha (Bonferroni‑Holm corrected)

% Figure aesthetics
cfg.fig.titleFont = 14;
cfg.fig.labelFont = 12;
cfg.fig.barColor  = [0.3 0.6 0.9; 0.9 0.3 0.3]; % cond 1 vs cond 2 (RGB)

% ---------------------
% ←‑‑‑ LOAD YOUR DATA ‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑▶
load('pac_datasets.mat')   % <‑‑ replace with your own .mat file

data(1,1) = pac_NPCC(2,1);
data(1,2) = pac_NPCC(2,2);


%% ------------------------------------------------------------------------
% SECTION 2 · PREP & SANITY CHECKS                                         
% -------------------------------------------------------------------------

nCond  = 2;

% Channel selection
if isempty(cfg.chanSelect)
    chanIdx = 1:numel(data{1}.label);
elseif iscell(cfg.chanSelect)
    chanIdx = find(ismember(data{1}.label, cfg.chanSelect));
else
    chanIdx = cfg.chanSelect;
end

% Latency indices
latIdx = dsearchn(data{1}.time{1}', cfg.latency');
latSmp = latIdx(1):latIdx(2);

%% ------------------------------------------------------------------------
% SECTION 3 · CORE PAC PIPELINE                                           
% -------------------------------------------------------------------------
SubjMI_NPCC = cell(length(pac_NPVP)-1,2);
for subj = 1:length(pac_NPVP)-1
fprintf('Computing PAC – Subject %d/%d …\n', subj, length(pac_NPVP)-1);
% Pre‑allocate results
allMI = cell(nCond, size(cfg.pacPairs,1));   % {cond, pair} → trials×channels MI
% Loop over conditions and PAC pairs
    for c = 1:nCond
        fprintf('Computing PAC – Condition %d/%d …\n', c, nCond);
        trials = pac_NPVP{subj+1,c}.trial;      %#ok<*USENS>
        nTrials = numel(trials);

        for t = 1:nTrials   
            ep = trials{t}(chanIdx, :);   % [chan × time]
            % Get freq bounds
            phBand = cfg.phaseBands.(cfg.pacPairs{1,1});
            amBand = cfg.ampBands.(cfg.pacPairs{1,2});
            % Compute MI per channel (helper below)
            mi = pac_MI(ep, cfg.srate, phBand, amBand, latSmp);
            allMI{c,1}(t,:) = mi;     % store [trial × chan]
        end
    end
SubjMI_NPCC(subj,:) = allMI';
end

save pac_delta_theta_statistics.mat SubjMI_NPVP SubjMI_NPCC


%% grand-averaged

Avg_MI_VP = cell(length(pac_NPVP)-1,1);
Avg_MI_NP = cell(length(pac_NPVP)-1,1);

for subj = 1:31
    Avg_MI_NP{subj} = mean(SubjMI{subj,1},1)';
    Avg_MI_VP{subj} = mean(SubjMI{subj,2},1)';
end

subj_MI_NP = zeros(31,55);
subj_MI_VP = zeros(31,55);

for subj = 1:31
    subj_MI_NP(subj,:) = Avg_MI_NP{subj,:}';
    subj_MI_VP(subj,:) = mean(SubjMI{subj,2},1)';
end

GA_MI_NP = mean(subj_MI_NP,1);
GA_MI_VP = mean(subj_MI_VP,1);


Y(1,:) = GA_MI_NP;
Y(2,:) = GA_MI_VP;
figure(1)
bar(data{1, 1}.label, Y, 'grouped');

legend({'NP','VP'}, 'Location', 'best');

xlabel('Channel');
ylabel('GA MI');
title('GA MI for NP vs VP');


[~,p,~,stats] = ttest(subj_MI_VP,subj_MI_NP);

p_fdr = mafdr(p', 'BHFDR', true);