% This is the event-related potentials analysis pipleline for EEG signals 
% used in our preprint paper:
% Yang, Q., Murphy, E., Yang, C., Liao, Y., & Hu, J. (2025). 
% Neural mechanisms of structural inference: an EEG investigation of 
% linguistic phrase structure categorization. 
% bioRxiv, 2025.2007.2001.662085. https://doi.org/10.1101/2025.07.01.662085
% 
% All analyses were carried out in MATLAB R2024b. The preprocessing
% relies on Fieldtrip:
% Oostenveld, R., Fries, P., Maris, E., & Schoffelen, J. M. (2011). 
% FieldTrip: Open source software for advanced analysis of MEG, EEG, 
% and invasive electrophysiological data. Comput Intell Neurosci, 
% 2011, 156869. https://doi.org/10.1155/2011/156869 

clear;clc

%% erp time-lock analysis

folder_name = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\8_by_conditions\erp\control\';
files_info = dir(fullfile(folder_name,'*.set'));
files_name = {files_info.name}; 
for subj = 1:length(files_name)
    restoredefaultpath; 
    addpath(genpath('E:\Toolboxes_for_matlab\eeglab2024.2\'));
    EEG = pop_loadset('filepath', folder_name, 'filename', files_name{1,subj});
    restoredefaultpath;
    addpath(genpath('E:\Toolboxes_for_matlab\fieldtrip-20250915\fieldtrip-20250915'));
    data = eeglab2fieldtrip( EEG, 'preprocessing');

    cfg = [];
    erp_CC{1,subj} = ft_timelockanalysis(cfg, data);
end

folder_name = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\8_by_conditions\erp\VP\';
files_info = dir(fullfile(folder_name,'*.set'));
files_name = {files_info.name}; 
for subj = 1:length(files_name)
    restoredefaultpath; 
    addpath(genpath('E:\Toolboxes_for_matlab\eeglab2024.2\'));
    EEG = pop_loadset('filepath', folder_name, 'filename', files_name{1,subj});
    restoredefaultpath;
    addpath(genpath('E:\Toolboxes_for_matlab\fieldtrip-20250915\fieldtrip-20250915\'));
    data = eeglab2fieldtrip( EEG, 'preprocessing');

    cfg = [];
    erp_VP{1,subj} = ft_timelockanalysis(cfg, data);
end


folder_name = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\8_by_conditions\erp\NP\';
files_info = dir(fullfile(folder_name,'*.set'));
files_name = {files_info.name}; 
for subj = 1:length(files_name)
    restoredefaultpath; 
    addpath(genpath('E:\Toolboxes_for_matlab\eeglab2024.2\'));
    EEG = pop_loadset('filepath', folder_name, 'filename', files_name{1,subj});
    restoredefaultpath;
    addpath(genpath('E:\Toolboxes_for_matlab\fieldtrip-20250915\fieldtrip-20250915\'));
    data = eeglab2fieldtrip( EEG, 'preprocessing');
  
    cfg = [];
    erp_NP{1,subj} = ft_timelockanalysis(cfg, data);
end

save erp_Averaged.mat  erp_CC erp_NP erp_VP

%% grand-average
clear;clc
load 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\2_statistical_analysis\erp_Averaged.mat\';
restoredefaultpath;
addpath 'E:\Toolboxes_for_matlab\fieldtrip-20250915\fieldtrip-20250915\';

cfg = [];   
cfg.keepindividual = 'yes';
grandavg_CC = ft_timelockgrandaverage(cfg, erp_CC{:}); 
grandavg_NP = ft_timelockgrandaverage(cfg, erp_NP{:}); 
grandavg_VP = ft_timelockgrandaverage(cfg, erp_VP{:}); 

cfg = [];
cfg.keepindividual = 'no';
grandavg_CC_plotting = ft_timelockgrandaverage(cfg, erp_CC{:}); 
grandavg_NP_plotting = ft_timelockgrandaverage(cfg, erp_NP{:}); 
grandavg_VP_plotting = ft_timelockgrandaverage(cfg, erp_VP{:});

%% basic plotting
cfg = [];
cfg.layout = 'EEG1005.lay';
cfg.interactive = 'yes';
cfg.showoutline = 'yes';
cfg.showlabels    = 'yes';
cfg.linecolor     = [0.1216, 0.4667, 0.7059; 1.0, 0.4980, 0.0549; 0.5020, 0.5020, 0.5020];
cfg.linewidth     = 1;
%Nature:
%Blue: [0.1216, 0.4667, 0.7059] Orange: [1.0, 0.4980, 0.0549] Green: [0.1725, 0.6275, 0.1725] Red: [0.8392, 0.1529, 0.1569] ...
%Purple: [0.5804, 0.4039, 0.7412]
%ft_multiplotER(cfg, grandavg_diff_plotting);
ft_multiplotER(cfg, grandavg_NP_plotting, grandavg_VP_plotting, grandavg_CC_plotting);

%% cluster-based monte-carlo permutation test 
cfg = [];
cfg.latency          = [-0.2 1.3];
cfg.avgovertime      = 'no';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT'; 
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.correcttail      = 'prob';
cfg.minnbchan        = 3;
cfg.clustertail      = 0;  %two-tailed 
cfg.tail             = 0;  %two-tailed 
% We used Bonferroni correction for three comparisons (NP, VP, control),
% therefore the critical alpha is set to be 0.05/3. We report orignial
% uncorrected p-values in the paper
cfg.alpha            = 0.05/3; 
cfg.numrandomization = 10000;
cfg_neighb = [];
cfg_neighb.method    = 'triangulation';
cfg_neighb.layout       = 'standard_1005.elc';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, grandavg_CC);

subj = length(erp_CC); 
design = zeros(2,2*subj);
for i = 1:subj
  design(1,i) = i;
end
for i = 1:subj
  design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design   = design;
cfg.uvar     = 1; 
cfg.ivar     = 2; 

%% pairwise comparisons

% NP vs control
stat.NP_vs_control = ft_timelockstatistics(cfg, grandavg_NP, grandavg_CC);
% VP vs control
stat.VP_vs_control = ft_timelockstatistics(cfg, grandavg_VP, grandavg_CC);
% VP vs NP
stat.VP_vs_NP = ft_timelockstatistics(cfg, grandavg_VP, grandavg_NP);
disp(any(stat.NP_vs_control.mask(:)))
disp(any(stat.VP_vs_control.mask(:)))
disp(any(stat.VP_vs_NP.mask(:))) 


save erp_statistics.mat stat
%%
cfg = [];
cfg.alpha  = 0.05;
cfg.parameter = 'stat';
cfg.layout = 'standard_1005.elc';
cfg.subplotsize = [4 4];
cfg.toi         = 0.1:0.1:1.3;
ft_clusterplot(cfg, stat.NP_vs_control);colorbar;clim([-5 5])   
ft_clusterplot(cfg, stat.VP_vs_control);colorbar;clim([-5 5])  
ft_clusterplot(cfg, stat.VP_vs_NP);colorbar;clim([-5 5])

%% Effect size computation -- cluster-based

cfg = [];
cfg.channel = 'EEG';
cfg.latency = [-0.2 1.3];
grandavg_a_sel = ft_selectdata(cfg,grandavg_NP);
grandavg_b_sel = ft_selectdata(cfg,grandavg_CC);

subj = length(GA_VP);

x1 = nan(subj,1);
x2 = nan(subj,1);

for i = 1:length(x1)
    sel3d = false(size(squeeze(grandavg_a_sel.individual)));
    sel3d(i,:,:) = stat_VP_vs_control.mask;
    tmp = grandavg_a_sel.individual(sel3d(:));
    x1(i) = mean(tmp);

    tmp = grandavg_b_sel.individual(sel3d(:));
    x2(i) = mean(tmp);
end
n1 = length(x1);
n2 = length(x2);

figure; plot([x1 x2]', 'o-'); xlim([0.5 2.5])
legend({'subj1', 'subj2', 'subj3', 'subj4', 'subj5', 'subj6', ...
  'subj7', 'subj8', 'subj9', 'subj10','subj11','subj12','subj13'...
  'subj14','subj15','subj16','subj17','subj18',...
  'subj19','subj20','subj21','subj22','subj23',...
  'subj24','subj25','subj26','subj27','subj28','subj29'}, 'location', 'EastOutside');
title('individual scores, averaged over cluster');

cohensd = mean(x1-x2) ./ std(x1-x2);
disp(cohensd)
