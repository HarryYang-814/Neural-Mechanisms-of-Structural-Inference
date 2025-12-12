% This is the propprocessing pipleline for EEG signals used in our 
% preprint paper:
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

%% parameter-setting 
clear;clc;

cfg                = [];                
cfg.method         = 'wavelet'; 
cfg.keeptrials     =  'no';
cfg.output         = 'pow';	
cfg.foi            = logspace(log10(2),log10(80),32);
cfg.width     = linspace(2,7,length(cfg.foi)); 
cfg.toi            = -2.6:0.01:2;
cfg.pad            = 'nextpow2';

%% VP condition
folder_name = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\8_by_conditions\oscillation\VP\';
files_info = dir(fullfile(folder_name,'*.set'));
files_name = {files_info.name}; 
for subj = 1:length(files_name)
    restoredefaultpath; 
    addpath(genpath('E:\Toolboxes_for_matlab\eeglab2024.2\'));
    EEG = pop_loadset('filepath', folder_name, 'filename', files_name{1,subj});
    restoredefaultpath;
    addpath(genpath('E:\Toolboxes_for_matlab\fieldtrip-20250915\fieldtrip-20250915'));
    data = eeglab2fieldtrip( EEG, 'preprocessing');
    
    TFRwave_VP{subj,1} = ft_freqanalysis(cfg, data);
end


%% NP condition
folder_name = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\8_by_conditions\oscillation\NP\';
files_info = dir(fullfile(folder_name,'*.set'));
files_name = {files_info.name}; 
for subj = 1:length(files_name)
    restoredefaultpath; 
    addpath(genpath('E:\Toolboxes_for_matlab\eeglab2024.2\'));
    EEG = pop_loadset('filepath', folder_name, 'filename', files_name{1,subj});
    restoredefaultpath;
    addpath (genpath('E:\Toolboxes_for_matlab\fieldtrip-20250915\fieldtrip-20250915'));
    data = eeglab2fieldtrip( EEG, 'preprocessing');

    TFRwave_NP{subj,1} = ft_freqanalysis(cfg, data);
end

%% control condition
folder_name = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\8_by_conditions\oscillation\control\';
files_info = dir(fullfile(folder_name,'*.set'));
files_name = {files_info.name}; 
for subj = 1:length(files_name)
    restoredefaultpath; 
    addpath(genpath('E:\Toolboxes_for_matlab\eeglab2024.2\'));
    EEG = pop_loadset('filepath', folder_name, 'filename', files_name{1,subj});
    restoredefaultpath;
    addpath (genpath('E:\Toolboxes_for_matlab\fieldtrip-20250915\fieldtrip-20250915'));
    data = eeglab2fieldtrip( EEG, 'preprocessing');
 
    TFRwave_CC{subj,1} = ft_freqanalysis(cfg, data);
end


save tf_power_2_80Hz.mat  TFRwave_VP TFRwave_NP TFRwave_CC


%% baseline norminalization
clc;clear 
load('G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\2_statistical_analysis\tf_power_2_80Hz.mat')
restoredefaultpath;
addpath 'E:\Toolboxes_for_matlab\fieldtrip-20250915\fieldtrip-20250915';
ft_defaults

cfg = [];

cfg.baseline = [-1.7 -1.2];
cfg.baselinetype =   'db'; % 'absolute', 'relative', 'relchange', 'normchange', 'db' or 'z-score'
cfg.parameter = 'powspctrm';
for subj = 1:length(TFRwave_VP)
    TFRwave_VP{subj,1} = ft_freqbaseline(cfg, TFRwave_VP{subj,1});
end
for subj = 1:length(TFRwave_NP)
    TFRwave_NP{subj,1} = ft_freqbaseline(cfg, TFRwave_NP{subj,1});
end
for subj = 1:length(TFRwave_CC)
    TFRwave_CC{subj,1} = ft_freqbaseline(cfg, TFRwave_CC{subj,1});
end

%% grand-average 
cfg = [];   
cfg.keepindividual = 'yes';
grandavg_VP = ft_freqgrandaverage(cfg, TFRwave_VP{:});
grandavg_NP = ft_freqgrandaverage(cfg, TFRwave_NP{:});
grandavg_control = ft_freqgrandaverage(cfg, TFRwave_CC{:});

cfg = [];   
cfg.keepindividual = 'no';
grandavg_VP_plot= ft_freqgrandaverage(cfg, TFRwave_VP{:});
grandavg_NP_plot = ft_freqgrandaverage(cfg, TFRwave_NP{:});
grandavg_control_plot = ft_freqgrandaverage(cfg, TFRwave_CC{:});

%% basic plotting
cfg = [];
cfg.showlabels   = 'yes';	        
cfg.layout       = 'EEG1005.lay';
cfg.colormap     = '*RdBu'; %The recommended colormaps include 'parula', 'cividis', 'balance', and '*RdBu'.
cfg.colorbar     = 'EastOutside';
cfg.zlim        = [-1 1];
cfg.ylim         = [30 80];
cfg.xlim         = [-1 1.3];

ft_multiplotTFR(cfg, grandavg_VP_plot);
ft_multiplotTFR(cfg, grandavg_NP_plot);
ft_multiplotTFR(cfg, grandavg_control_plot);

%% cluster-based permutation t-test ft_statistics_montecarlo
cfg = [];
cfg.latency          = [-1 1];
cfg.frequency        = [40 80]; %delta(2-3.5 Hz) theta(4-7 Hz) alpha(8-12 Hz) beta (15-30 Hz) gamma (40-80 Hz)
cfg.avgoverfreq      = 'yes';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT'; 
cfg.correctm         = 'cluster';
cfg.clusterthreshold = 'parametric';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 3;
cfg.clustertail      = 0;
cfg.tail             = 0;
cfg.correcttail      = 'prob';
cfg.alpha            = 0.05/3;
cfg.numrandomization = 10000;
% specifies with which sensors other sensors can form clusters
cfg_neighb = [];
cfg_neighb.method    = 'triangulation';
cfg_neighb.layout       = 'standard_1005.elc';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, grandavg_VP);

subj = length(TFRwave_VP);
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


%% delta
stat.sp.delta.VPvsNP = ft_freqstatistics(cfg, grandavg_NP, grandavg_VP);
stat.sp.delta.NPvsCC = ft_freqstatistics(cfg, grandavg_NP, grandavg_control);
stat.sp.delta.VPvsCC = ft_freqstatistics(cfg, grandavg_VP, grandavg_control);

disp(any(stat.sp.delta.VPvsNP.mask(:)))
disp(any(stat.sp.delta.NPvsCC.mask(:)))
disp(any(stat.sp.delta.VPvsCC.mask(:)))


%%
cfg = [];
cfg.alpha  = 0.1;
cfg.parameter = 'stat';
cfg.layout = 'standard_1005.elc';
cfg.colormap = '*RdBu';
cfg.subplotsize = [4 6];
cfg.toi         = -1:0.1:1.3;
ft_clusterplot(cfg, stat.delta.VPvsNP);colorbar;clim([-8 8])
ft_clusterplot(cfg, stat.delta.VPvsCC);colorbar;clim([-8 8])
ft_clusterplot(cfg, stat.delta.NPvsCC);colorbar;clim([-8 8])

%% theta
stat.sp.theta.VPvsNP = ft_freqstatistics(cfg, grandavg_NP, grandavg_VP);
stat.sp.theta.NPvsCC = ft_freqstatistics(cfg, grandavg_NP, grandavg_control);
stat.sp.theta.VPvsCC = ft_freqstatistics(cfg, grandavg_VP, grandavg_control);

disp(any(stat.sp.theta.VPvsNP.mask(:)))
disp(any(stat.sp.theta.NPvsCC.mask(:)))
disp(any(stat.sp.theta.VPvsCC.mask(:)))



%%
cfg = [];
cfg.alpha  = 0.1;
cfg.parameter = 'stat';
cfg.layout = 'standard_1005.elc';
cfg.colormap = '*RdBu';
cfg.subplotsize = [4 6];
cfg.toi         = -1:0.1:1.3;
ft_clusterplot(cfg, stat.theta.VPvsNP);colorbar;clim([-5 5]);
ft_clusterplot(cfg, stat.theta.VPvsCC);colorbar;clim([-8 8]);
ft_clusterplot(cfg, stat.theta.NPvsCC);colorbar;clim([-8 8]);
%% alpha
stat.sp.alpha.VPvsNP = ft_freqstatistics(cfg, grandavg_NP, grandavg_VP);
stat.sp.alpha.NPvsCC = ft_freqstatistics(cfg, grandavg_NP, grandavg_control);
stat.sp.alpha.VPvsCC = ft_freqstatistics(cfg, grandavg_VP, grandavg_control);

disp(any(stat.sp.alpha.VPvsNP.mask(:)))
disp(any(stat.sp.alpha.NPvsCC.mask(:)))
disp(any(stat.sp.alpha.VPvsCC.mask(:)))



%%
cfg = [];
cfg.alpha  = 0.05;
cfg.parameter = 'stat';
cfg.layout = 'standard_1005.elc';
cfg.colormap = '*RdBu';
cfg.subplotsize = [4 6];
cfg.toi         = -1:0.1:1.3;
ft_clusterplot(cfg, stat.alpha.VPvsNP);colorbar;clim([-3 3])
ft_clusterplot(cfg, stat.alpha.NPvsCC);colorbar;clim([-8 8])
ft_clusterplot(cfg, stat.alpha.VPvsCC);colorbar;clim([-8 8])
%% beta 
stat.sp.beta.VPvsNP = ft_freqstatistics(cfg, grandavg_NP, grandavg_VP);
stat.sp.beta.NPvsCC = ft_freqstatistics(cfg, grandavg_NP, grandavg_control);
stat.sp.beta.VPvsCC = ft_freqstatistics(cfg, grandavg_VP, grandavg_control);

disp(any(stat.sp.beta.VPvsNP.mask(:)))
disp(any(stat.sp.beta.NPvsCC.mask(:)))
disp(any(stat.sp.beta.VPvsCC.mask(:)))

%%
cfg = [];
cfg.alpha  = 0.1;
cfg.parameter = 'stat';
cfg.layout = 'standard_1005.elc';
cfg.colormap = '*RdBu';
cfg.subplotsize = [4 6];
cfg.toi         = -1:0.1:1;
ft_clusterplot(cfg, stat.beta.VPvsNP);colorbar;clim([-3 3])
ft_clusterplot(cfg, stat.beta.NPvsCC);colorbar;clim([-8 8])
ft_clusterplot(cfg, stat.beta.VPvsCC);colorbar;clim([-8 8])


%% gamma 
stat.sp.gamma.VPvsNP = ft_freqstatistics(cfg, grandavg_NP, grandavg_VP);
stat.sp.gamma.NPvsCC = ft_freqstatistics(cfg, grandavg_NP, grandavg_control);
stat.sp.gamma.VPvsCC = ft_freqstatistics(cfg, grandavg_VP, grandavg_control);

disp(any(stat.sp.gamma.VPvsNP.mask(:)))
disp(any(stat.sp.gamma.NPvsCC.mask(:)))
disp(any(stat.sp.gamma.VPvsCC.mask(:)))


cfg = [];
cfg.alpha  = 0.1;
cfg.parameter = 'stat';
cfg.layout = 'standard_1005.elc';
cfg.colormap = '*RdBu';
cfg.subplotsize = [4 6];
cfg.toi         = -1:0.1:1;
ft_clusterplot(cfg, stat.gamma.VPvsNP);colorbar;clim([-3 3])
ft_clusterplot(cfg, stat.gamma.NPvsCC);colorbar;clim([-8 8])
ft_clusterplot(cfg, stat.gamma.VPvsCC);colorbar;clim([-8 8])
%% 
stat.ex.VPvsNP = ft_freqstatistics(cfg, grandavg_NP, grandavg_VP);
stat.ex.NPvsCC = ft_freqstatistics(cfg, grandavg_NP, grandavg_control);
stat.ex.VPvsCC = ft_freqstatistics(cfg, grandavg_VP, grandavg_control);

disp(any(stat.ex.VPvsNP.mask(:)))
disp(any(stat.ex.NPvsCC.mask(:)))
disp(any(stat.ex.VPvsCC.mask(:)))

save tf_power_statistics_2_80Hz_new.mat stat
%% cluster average effect size
cfg = [];
cfg.channel = 'EEG';
cfg.latency = [-1 1];
cfg.frequency = [4 7];
cfg.avgoverfreq = 'yes';
grandavg_a_sel = ft_selectdata(cfg,grandavg_NP);
grandavg_b_sel = ft_selectdata(cfg,grandavg_VP);

subj = length(TFRwave_VP);

x1 = nan(subj,1);
x2 = nan(subj,1);


for i = 1:length(x1)
    sel3d = false(size(squeeze(grandavg_a_sel.powspctrm)));
    a = stat.ex.VPvsNP.mask;
    %a(a ~= 2) = 0;
    sel3d(i,:,:) = squeeze(a);
    tmp1 = squeeze(grandavg_a_sel.powspctrm);
    tmp2 = tmp1(sel3d(:));
    x1(i) = mean(tmp2);

    tmp1 = squeeze(grandavg_b_sel.powspctrm);
    tmp2 = tmp1(sel3d(:));
    x2(i) = mean(tmp2);
end
n1 = length(x1);
n2 = length(x2);

figure; plot([x1 x2]', 'o-'); xlim([0.5 2.5])
legend({'subj1', 'subj2', 'subj3', 'subj4', 'subj5', 'subj6', ...
  'subj7', 'subj8', 'subj9', 'subj10','subj11','subj12','subj13'...
  'subj14','subj15','subj16','subj17','subj18',...
  'subj19','subj20','subj21','subj22','subj23',...
  'subj24','subj25','subj26','subj27','subj28','subj29', 'subj30'...
  ,'subj31'}, 'location', 'EastOutside');
title('individual scores, averaged over cluster');

cohensd = mean(x1-x2) ./ std(x1-x2);
disp(cohensd)
