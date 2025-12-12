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
cfg.keeptrials     =  'yes';
cfg.output         = 'fourier';	
cfg.foi            = logspace(log10(2),log10(80),32);
cfg.width     = linspace(2,7,length(cfg.foi)); 
cfg.toi            = -2.6:0.01:2;
cfg.pad            = 'nextpow2';
cfg.channel        = {'FCz','FC1','FC2','FC3','FC4','Cz','C1','C2','C3',...
    'C4','CPz','CP1','CP2','CP3','CP4'};

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
    freq = ft_freqanalysis(cfg, data);
    
    itc = [];
    itc.label     = freq.label;
    itc.freq      = freq.freq;
    itc.time      = freq.time;
    itc.elec      = freq.elec;
    itc.dimord    = 'chan_freq_time';
    
    F = freq.fourierspctrm;   % copy the Fourier spectrum
    N = size(F,1);           % number of trials

    % compute inter-trial phase coherence (itpc)
    itc.itpc      = F./abs(F);         % divide by amplitude
    itc.itpc      = sum(itc.itpc,1);   % sum angles
    itc.itpc      = abs(itc.itpc)/N;   % take the absolute value and normalize
    itc.itpc      = squeeze(itc.itpc); % remove the first singleton dimension
    
    ITC_VP{1,subj} = itc;   
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
    freq = ft_freqanalysis(cfg, data);
    
    itc = [];
    itc.label     = freq.label;
    itc.freq      = freq.freq;
    itc.time      = freq.time;
    itc.elec      = freq.elec;
    itc.dimord    = 'chan_freq_time';
    
    F = freq.fourierspctrm;   % copy the Fourier spectrum
    N = size(F,1);           % number of trials

    % compute inter-trial phase coherence (itpc)
    itc.itpc      = F./abs(F);         % divide by amplitude
    itc.itpc      = sum(itc.itpc,1);   % sum angles
    itc.itpc      = abs(itc.itpc)/N;   % take the absolute value and normalize
    itc.itpc      = squeeze(itc.itpc); % remove the first singleton dimension
    
    ITC_NP{1,subj} = itc;   
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
 
    freq = ft_freqanalysis(cfg, data);
    
    itc = [];
    itc.label     = freq.label;
    itc.freq      = freq.freq;
    itc.time      = freq.time;
    itc.elec      = freq.elec;
    itc.dimord    = 'chan_freq_time';
    
    F = freq.fourierspctrm;   % copy the Fourier spectrum
    N = size(F,1);           % number of trials

    % compute inter-trial phase coherence (itpc)
    itc.itpc      = F./abs(F);         % divide by amplitude
    itc.itpc      = sum(itc.itpc,1);   % sum angles
    itc.itpc      = abs(itc.itpc)/N;   % take the absolute value and normalize
    itc.itpc      = squeeze(itc.itpc); % remove the first singleton dimension

    ITC_control{1,subj} = itc;   
end


save itpc_2-80Hz.mat  ITC_VP  ITC_NP  ITC_control


%% baseline 

clc;clear
load('G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\2_statistical_analysis\itpc_2-80Hz.mat')

cfg = [];
cfg.baseline = [-1.7 -1.2];
cfg.baselinetype =  'absolute';
cfg.parameter = 'itpc';
for i = 1:length(ITC_VP)
     ITC_VP{1,i} = ft_freqbaseline(cfg, ITC_VP{1,i});
end
for i = 1:length(ITC_NP)
     ITC_NP{1,i} = ft_freqbaseline(cfg, ITC_NP{1,i});
end
for i = 1:length(ITC_control)
     ITC_control{1,i} = ft_freqbaseline(cfg, ITC_control{1,i});
end

cfg = [];   
cfg.keepindividual = 'yes';
cfg.parameter = 'itpc';
grandavg_VP = ft_freqgrandaverage(cfg, ITC_VP{:});
grandavg_NP = ft_freqgrandaverage(cfg, ITC_NP{:});    
grandavg_control = ft_freqgrandaverage(cfg, ITC_control{:});    
%% create averaged channel
grandavg_NP.itpc(:,16,:,:) = mean(grandavg_NP.itpc(:,:,:,:),2);
grandavg_NP.label{16,1} = 'avg';

grandavg_VP.itpc(:,16,:,:) = mean(grandavg_VP.itpc(:,:,:,:),2);
grandavg_VP.label{16,1} = 'avg';

grandavg_control.itpc(:,16,:,:) = mean(grandavg_control.itpc(:,:,:,:),2);
grandavg_control.label{16,1} = 'avg';


%% cluster-based permutation
cfg = [];
cfg.method            = 'montecarlo';           % use the Monte Carlo Method to calculate the significance probability
cfg.statistic         = 'depsamplesT';% use the dependent samples T-statistic 
cfg.frequency        = [2 80]; %delta = 0.5 - 4 Hz, theta = 4-7 Hz, alpha = 8-12 Hz, beta = 13-30 Hz, narrowband gamma = 30-50 Hz
cfg.avgovertime       = 'no';
cfg.avgoverfreq       = 'no';
cfg.latency           = [-1 1];
cfg.correctm          = 'cluster';
cfg.clusteralpha      = 0.05;                   % alpha level of the sample-specific test statistic that will be used for thresholding
cfg.clustertail       = 0;
cfg.clusterstatistic  = 'maxsum';               % test statistic that will be evaluated under the permutation distribution.
cfg.tail              = 0;                      % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.correcttail       = 'prob';                 % the two-sided test implies that we do non-parametric two tests
cfg.alpha             = 0.05/3;                   % alpha level of the permutation test
cfg.numrandomization  = 10000;                   % number of draws from the permutation distribution
cfg.channel           = {'avg'};
cfg.neighbours        = [];                     % there are no spatial neighbours, only in time and frequency

subj = length(ITC_VP); 
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

cfg.parameter = 'itpc';

%%
stat.NP_VP = ft_freqstatistics(cfg, grandavg_NP, grandavg_VP);
stat.NP_control = ft_freqstatistics(cfg, grandavg_NP, grandavg_control);
stat.VP_control = ft_freqstatistics(cfg, grandavg_VP, grandavg_control);

save itpc_stat_2-80Hz.mat stat

cfg = [];
cfg.interactive   = 'no';
cfg.channel       = {'avg'};
cfg.renderer      = 'painters';
cfg.colorbar      = 'yes';
cfg.xlim          = [-1 1];
cfg.zlim          = [-8 8];
cfg.ylim          = [2 30];
cfg.maskparameter = 'mask';       % use significance to mask the power difference
cfg.maskstyle      = 'outline';
cfg.colormap      = '*RdBu';
cfg.maskparameter = 'mask';       % use significance to mask the power difference
cfg.parameter     = 'stat';       % display the statistical value, i.e. the t-score

ft_singleplotTFR(cfg, stat.NP_control);
title('ITPC_NP-control','FontName','Calibri','FontSize',16);
xlabel('Time (s)','FontName','Calibri','FontSize',14); 
ylabel('Frequency (Hz)','FontName','Calibri','FontSize',14);
xline(0,'--','LineWidth',1);
set(gca,'FontName','Calibri','FontSize',12)
yl = ylim(gca);
y_pos = yl(2) - 0.1*(yl(2)-yl(1));
text(0, y_pos, 'Word 2', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','top', ...
    'Rotation',90,'FontSize',12,'FontName','Calibri');


ft_singleplotTFR(cfg, stat.VP_control);
title('ITPC_VP-control','FontName','Calibri','FontSize',16);
xlabel('Time (s)','FontName','Calibri','FontSize',14); 
ylabel('Frequency (Hz)','FontName','Calibri','FontSize',14);
xline(0,'--','LineWidth',1);
set(gca,'FontName','Calibri','FontSize',12)
text(0, y_pos, 'Word 2', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','top', ...
    'Rotation',90,'FontSize',12,'FontName','Calibri');