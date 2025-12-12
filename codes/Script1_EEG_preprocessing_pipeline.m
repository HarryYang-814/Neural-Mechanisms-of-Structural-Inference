% This is the propprocessing pipleline for EEG signals used in our 
% preprint paper:
% Yang, Q., Murphy, E., Yang, C., Liao, Y., & Hu, J. (2025). 
% Neural mechanisms of structural inference: an EEG investigation of 
% linguistic phrase structure categorization. 
% bioRxiv, 2025.2007.2001.662085. https://doi.org/10.1101/2025.07.01.662085
% 
% All analyses were carried out in MATLAB R2024b. The preprocessing
% relies on the following toolboxes: EEGLAB 2024.2, ERPLAB v12.00,
% , ICLabel v1.6., CSD toolbox
%
% The current script implements the full preprocessing workflow described
% in the paper, including: 
% (1) raw data import and basic data cleaning;
% (2) ICA decomposition on 1hz high-pass filtered data; 
% (3) automatic IC components (ocular artifacts) rejection
% (4) bad channel detection;
% (5) epoching, baseline subtraction, interpolation and rereference;
% (6) bad epochs rejection.
% (7) apply surface laplacian transformation
% (8) split the datasets into conditions


clear;clc;
restoredefaultpath
addpath 'E:\Toolboxes_for_matlab\eeglab2024.2\eeglab2024.2\';
%% (1) Raw data import and basic data cleaning
%% 1Hz high-pass filtering data
data_path = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\0_raw_data\';
save_path = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\1_basic_cleaning\hp1hz\';
cd(data_path)
files = dir('*.vhdr');
file_name = {files.name};
eeglab

for subj = 1:length(file_name)
    EEG = pop_loadbv(data_path, file_name{subj});
    % we shift the event codes 13 ms earlier than original ones because
    % a photodiode test in the pilot experiment revealed a constant 13 ms
    % delay between the physical stimulus onset and the event markers
    % logged in the EEG data. We therefore correct for this fixed latency
    % at the preprocessing stage.
    EEG = pop_erplabShiftEventCodes( EEG , 'DisplayEEG', 0, 'DisplayFeedback', ...
        'summary', 'Eventcodes', { 'S 21' 'S 22' 'S 23' 'S 11' 'S 12' 'S 13' }, ...
        'Rounding', 'earlier', 'Timeshift',  13 );
    % remove segments corresponding to resting periods
    EEG  = pop_erplabDeleteTimeSegments( EEG , 'afterEventcodeBufferMS', ...
        3000, 'beforeEventcodeBufferMS', ...
        3000, 'displayEEG',  0, 'ignoreBoundary',...
  0, 'ignoreUseType', 'ignore', 'timeThresholdMS',  7000 );
    % standard_1005.elc is used for labeling channel locations
    EEG = pop_chanedit(EEG, 'lookup','E:\Toolboxes_for_matlab\eeglab2024.2\eeglab2024.2\plugins\dipfit5.5\standard_BEM\elec\standard_1005.elc');    
    % 1-80 Hz band-pass filter, non-causal, FIR
    EEG = pop_eegfiltnew(EEG, 'locutoff',1, 'hicutoff',80);
    % notch filtering to remove 50 Hz line noise
    EEG = pop_eegfiltnew(EEG, 'locutoff',48,'hicutoff',52,'revfilt',1);
    %downsampling to 500 Hz
    EEG = pop_resample(EEG, 500);
    % remove channels that are not used in later statistical analysis
    EEG = pop_select( EEG, 'rmchannel',{'I4','I5', 'I6','I3','HEOR','HEOL','VEOU','VEOL'});
    EEG.setname = [file_name{subj}(1:5),'.set'];
    EEG = pop_saveset(EEG,'filename',strcat(file_name{1,subj}(1:5),'_hp1hz.set'),'filepath',save_path);
end

%% 0.1Hz high-pass filtering data
data_path = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\0_raw_data\';
save_path = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\1_basic_cleaning\hp0.1hz\';
cd(data_path)
files = dir('*.vhdr');
file_name = {files.name};
eeglab

for subj = 1:length(file_name)
    EEG = pop_loadbv(data_path, file_name{subj});
    % we shift the event codes 13 ms earlier than original ones because
    % a photodiode test in the pilot experiment revealed a constant 13 ms
    % delay between the physical stimulus onset and the event markers
    % logged in the EEG data. We therefore correct for this fixed latency
    % at the preprocessing stage.
    EEG = pop_erplabShiftEventCodes( EEG , 'DisplayEEG', 0, 'DisplayFeedback', ...
        'summary', 'Eventcodes', { 'S 21' 'S 22' 'S 23' 'S 11' 'S 12' 'S 13' }, ...
        'Rounding', 'earlier', 'Timeshift',  13 );
    % remove segments corresponding to resting periods
    EEG  = pop_erplabDeleteTimeSegments( EEG , 'afterEventcodeBufferMS', ...
        3000, 'beforeEventcodeBufferMS', ...
        3000, 'displayEEG',  0, 'ignoreBoundary',...
  0, 'ignoreUseType', 'ignore', 'timeThresholdMS',  7000 );
    % standard_1005.elc is used for labeling channel locations
    EEG = pop_chanedit(EEG, 'lookup','E:\Toolboxes_for_matlab\eeglab2024.2\eeglab2024.2\plugins\dipfit5.5\standard_BEM\elec\standard_1005.elc');    
    % 0.1-80 Hz band-pass filter, non-causal, FIR
    EEG = pop_eegfiltnew(EEG, 'locutoff',0.2, 'hicutoff',80);
    % notch filtering to remove 50 Hz line noise
    EEG = pop_eegfiltnew(EEG, 'locutoff',48,'hicutoff',52,'revfilt',1);
    %downsampling to 500 Hz
    EEG = pop_resample(EEG, 500);
    % remove channels that are not used in later statistical analysis
    EEG = pop_select( EEG, 'rmchannel',{'I4','I5', 'I6','I3','HEOR','HEOL','VEOU','VEOL'});
    EEG.setname = [file_name{subj}(1:5),'.set'];
    EEG = pop_saveset(EEG,'filename',strcat(file_name{1,subj}(1:5),'_hp0.1hz.set'),'filepath',save_path);
end

%% (2) running ICA on 1-Hz-highpass-filtered data 
clear;clc
data_path = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\1_basic_cleaning\hp1hz';
save_path = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\2_applying_ICAdecomp\hp1Hz';
cd(data_path)
files = dir('*.set');
file_name = {files.name};
eeglab

for subj = 2:length(file_name)
    EEG = pop_loadset('filename',file_name{subj},'filepath',data_path);
    % we used PCA to reduce the data dimensionality to 50 components prior 
    % to ICA, thereby reducing computational load
    EEG = pop_runica(EEG, 'extended',1,'pca',50,'interupt','on');
    EEG = pop_saveset(EEG, 'filename',strcat(file_name{1,subj}(1:6),'hp1hz_ICA.set'),'filepath',save_path); 
end

%% transfer ICA matrix to 0.1-hz-hp-filtered data
clear;clc

data_path1 = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\1_basic_cleaning\hp0.1hz';
data_path2 = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\2_applying_ICAdecomp\hp1Hz';
save_path  = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\2_applying_ICAdecomp\hp0.1Hz';

files_info1 = dir(fullfile(data_path1,'*.set'));
files_info2 = dir(fullfile(data_path2,'*.set'));
files_name1 = {files_info1.name};
files_name2 = {files_info2.name};
eeglab


for subj = 1:length(files_name2)
    EEG_1 = pop_loadset(files_name1{subj},data_path1);
    EEG_2 = pop_loadset(files_name2{subj},data_path2);
    % copying ICA weights from EEG_2 to EEG_1
    EEG_1.icaweights = EEG_2.icaweights;
    EEG_1.icasphere = EEG_2.icasphere;
    EEG_1.icaact = EEG_2.icaact;
    EEG_1.icawinv = EEG_2.icawinv;
    EEG_1.icachansind = EEG_2.icachansind;
    EEG_1 = pop_saveset(EEG_1,'filepath',save_path,'filename',strcat(files_name1{1,subj}(1:6),'hp0.1hz_ICA.set'));
end

%% (3) automatic IC components rejection
clear;clc

data_path = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\2_applying_ICAdecomp\hp0.1Hz';
save_path = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\3_ICrejection';
cd(data_path)
files = dir('*.set');
file_name = {files.name};
num_comp_rejected{1,1} = 'num_IC_rejected';
num_comp_rejected{1,2} = 'subject';

eeglab

for subj = 1:length(file_name)
    EEG = pop_loadset('filename',file_name{subj},'filepath',data_path);
    EEG = pop_iclabel(EEG, 'default');
    % automatic rejection of ocular, muscle and channel noise artifacts with at least 90% confidence
    EEG = pop_icflag(EEG, [0 0;0.9 1; 0.9 1; 0 0; 0 0; 0.9 1; 0 0]);
    reject_comp = find(EEG.reject.gcompreject == 1);
    % save the number of removed IC components for each dataset
    num_comp_rejected{subj,1} = length(reject_comp);
    num_comp_rejected{subj,2} = file_name{1,subj}(1:5);
    EEG = pop_subcomp(EEG,reject_comp,0,0);
    EEG = pop_saveset(EEG, 'filename',strcat(file_name{1,subj}(1:6),'hp0.1hz_ICremoved.set'),'filepath',save_path); 
end
% mannually inspect other IC artifacts

save removed_ICnum.mat num_comp_rejected

%% (4) bad channels detection
clear;clc

data_path = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\3_ICrejection';
save_path = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\4_interpolation';
cd(data_path)
files = dir('*.set');
file_name = {files.name};
eeglab


% bad channels are identified by manual inspection of the raw data and
% power spectral density
for subj = 1:length(file_name)
    EEG = pop_loadset('filename',file_name{subj},'filepath',data_path);
    % bad channels are automaticlly detected with the defaul settings of
    % the function 'pop_clean_rawdata' and then mannually inspected the
    % channels
    % for double-check 
    EEG_rmchan = pop_clean_rawdata(EEG, 'FlatlineCriterion',5,...
        'ChannelCriterion',0.8,'LineNoiseCriterion',4,'Highpass','off',...
        'BurstCriterion','off','WindowCriterion','off','BurstRejection',...
        'off','Distance','Euclidian','channels_ignore',{'M2'});
    if isfield(EEG_rmchan.chaninfo, 'removedchans') && ~isempty(EEG_rmchan.chaninfo.removedchans)
    rmchan{subj,1} = EEG_rmchan.chaninfo.removedchans(1,10:length(EEG_rmchan.chaninfo.removedchans));
    end
    EEG = pop_saveset(EEG, 'filename',strcat(file_name{1,subj}(1:6),'interpolated.set'),'filepath',save_path); 
end

save rjchan_info.mat rmchan

%% (5) epoching, baseline subtraction, interpolation, rereference 
clear;clc
data_path = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\4_interpolation';
save_path = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\5_epoched';
cd(data_path)
files = dir('*.set');
file_name = {files.name};
load rjchan_info.mat

for subj = 1:length(file_name)
    EEG = pop_loadset('filename',file_name{subj},'filepath',data_path);
    EEG = pop_epoch(EEG,{'S 21'  'S 22'  'S 23'}, [-2.6  2], 'epochinfo', 'yes');
    EEG = pop_rmbase(EEG, [-2100   -2000]);
    EEG = pop_interp(EEG,[rmchan{subj,1}.urchan],'spherical');    
    EEG = pop_chanedit(EEG, 'append',56,'changefield',{57 'labels' 'M1'},...
        'lookup','E:\Toolboxes_for_matlab\eeglab2024.2\eeglab2024.2\plugins\dipfit5.5\standard_BEM\elec\standard_1005.elc');
    EEG = pop_reref( EEG, [],'refloc',struct('labels',{'M1'},'type',{''},...
        'theta',{-117.5949},'radius',{0.6944},'X',{-44.9897},'Y',{86.0761},...
        'Z',{-67.986},'sph_theta',{117.5949},'sph_phi',{-34.9916},'sph_radius',...
        {118.5549},'urchan',{57},'ref',{''},'datachan',{0}));
    EEG = pop_reref(EEG,{'M1','M2'});
    EEG = pop_saveset(EEG, 'filename',strcat(file_name{1,subj}(1:6),'epoched.set'),'filepath',save_path);     
end


%% (6) trial rejection
clear;clc
data_path = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\5_epoched';
save_path = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\6_trialrej';
cd(data_path)
files = dir('*.set');
file_name = {files.name};

rtrials{1,1} = 'remained_trials';
rtrials{1,2} = 'remained_trials(%)';
rtrials{1,3} = 'Subject';
 
for subj = 1:length(file_name)
    EEG = pop_loadset('filename',file_name{subj},'filepath',data_path);
    %EEG = pop_eegthresh(EEG,1,1:EEG.nbchan,-100,100,-2.6,1.998,0,1);
    rtrials{subj+1,1} = EEG.trials;
    rtrials{subj+1,2} = (EEG.trials/342)*100;
    rtrials{subj+1,3} = file_name{subj};    
    %EEG = pop_saveset(EEG, 'filename',strcat(file_name{1,subj}(1:6),'cleaned.set'),'filepath',save_path);     
end
% mannually inspect and reject trials 
save remained_trials.mat rtrials

%% (7) apply suface laplacian transformation (CSD)
clear;clc

data_path = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\6_trialrej';
save_path = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\7_CSD';
cd(data_path)
files = dir('*.set');
file_name = {files.name};

EEG = pop_loadset('filename',file_name{1,1},'filepath',data_path);% 
electrodes = {EEG.chanlocs.labels}'; % 

restoredefaultpath; 
addpath(genpath('E:\Toolboxes_for_matlab\csd_toolbox\'));
Montage = ExtractMontage('10-5-System_Mastoids_EGI129.csd',electrodes); 



restoredefaultpath; 
addpath(genpath('E:\Toolboxes_for_matlab\csd_toolbox\'));
[G,H] = GetGH(Montage);

addpath 'E:\Toolboxes_for_matlab\eeglab2024.2\eeglab2024.2\';
eeglab
for subj = 1:length(file_name)
    EEG = pop_loadset('filename',file_name{1,subj},'filepath',data_path); 
    EEG = eeg_checkset( EEG );     
    parfor i = 1:size(EEG.data,3)
        D = squeeze(EEG.data(:,:,i));
        X = CSD (D, G, H);
        CSDdata(:,:,i) = X;
    end
    EEG.data = CSDdata;
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset(EEG,'filename',strcat(file_name{1,subj}(1:6), 'CSD.set'), 'filepath',save_path);
    clear CSDdata
end

%% (8) By conditions
clear;clc;
restoredefaultpath
addpath 'E:\Toolboxes_for_matlab\eeglab2024.2\eeglab2024.2\';
eeglab

data_path = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\7_CSD';
cd(data_path)
files = dir('*.set');
file_name = {files.name};

save_path_VP_tf = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\8_by_conditions\oscillation\VP\';
save_path_VP_erp = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\8_by_conditions\erp\VP\';

save_path_NP_tf = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\8_by_conditions\oscillation\NP\';
save_path_NP_erp = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\8_by_conditions\erp\NP\';

save_path_control_tf = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\8_by_conditions\oscillation\control\';
save_path_control_erp = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\8_by_conditions\erp\control\';

rtrials{1,1} = 'remained_trials';
rtrials{1,2} = 'remained_trials(%)';
rtrials{1,3} = 'VP';
rtrials{1,4} = 'VP(%)';
rtrials{1,5} = 'NP';
rtrials{1,6} = 'NP(%)';
rtrials{1,7} = 'control';
rtrials{1,8} = 'control(%)';
rtrials{1,9} = 'Subject';

for subj = 1:length(file_name)
    EEG = pop_loadset('filepath', data_path, 'filename', file_name{1,subj});
    %time-frequency
    EEG_21 = pop_epoch( EEG, {  'S 21' }, [-2.6  2], 'newname', strcat(file_name{1,subj}(1:5), '_VP'), 'epochinfo', 'yes'); 
    EEG_21 = pop_rmbase( EEG_21, [-2100   -2000]);
    EEG_21.subject = file_name{1,subj}(1:5);
    EEG_21.condition = 'VP';
    EEG_21 = pop_saveset(EEG_21, 'filepath', save_path_VP_tf, 'filename', strcat(file_name{1,subj}(1:6),'VP.set'));
    %ERP
    EEG_21 = pop_epoch( EEG_21, {  'S 21'  }, [-0.2   1.3], 'newname', strcat(file_name{1,subj}(1:5), '_VP'), 'epochinfo', 'yes'); 
    EEG_21 = pop_rmbase( EEG_21, [-200  0]);
    EEG_21.subject = file_name{1,subj}(1:5);
    EEG_21.condition = 'VP';
    EEG_21 = pop_saveset(EEG_21, 'filepath', save_path_VP_erp, 'filename', strcat(file_name{1,subj}(1:6),'VP.set'));
    %time-frequency
    EEG_22 = pop_epoch( EEG, {  'S 22'  }, [-2.6  2], 'newname', strcat(file_name{1,subj}(1:5), '_NP'), 'epochinfo', 'yes'); 
    EEG_22 = pop_rmbase( EEG_22, [-2100   -2000]);
    EEG_22.subject = file_name{1,subj}(1:5);
    EEG_22.condition = 'NP';
    EEG_22 = pop_saveset(EEG_22, 'filepath', save_path_NP_tf, 'filename', strcat(file_name{1,subj}(1:6),'NP.set'));
    %ERP
    EEG_22 = pop_epoch( EEG_22, {  'S 22'  }, [-0.2   1.3], 'newname', strcat(file_name{1,subj}(1:5), '_NP'), 'epochinfo', 'yes'); 
    EEG_22 = pop_rmbase( EEG_22, [-200  0]);
    EEG_22.subject = file_name{1,subj}(1:5);
    EEG_22.condition = 'NP';
    EEG_22 = pop_saveset(EEG_22, 'filepath', save_path_NP_erp, 'filename', strcat(file_name{1,subj}(1:6),'NP.set'));
    %time-frequency
    EEG_23 = pop_epoch( EEG, {  'S 23'  }, [-2.6  2], 'newname', strcat(file_name{1,subj}(1:5), '_control'), 'epochinfo', 'yes'); 
    EEG_23 = pop_rmbase( EEG_23, [-2100   -2000]);
    EEG_23.subject = file_name{1,subj}(1:5);
    EEG_23.condition = 'CC';
    EEG_23 = pop_saveset(EEG_23, 'filepath', save_path_control_tf, 'filename', strcat(file_name{1,subj}(1:6),'control.set'));
    %ERP
    EEG_23 = pop_epoch( EEG_23, {  'S 23'  }, [-0.2   1.3], 'newname', strcat(file_name{1,subj}(1:5), '_control'), 'epochinfo', 'yes'); 
    EEG_23 = pop_rmbase( EEG_23, [-200  0]);
    EEG_23.subject = file_name{1,subj}(1:5);
    EEG_23.condition = 'CC';
    EEG_23 = pop_saveset(EEG_23, 'filepath', save_path_control_erp, 'filename', strcat(file_name{1,subj}(1:6),'control.set'));

    rtrials{subj+1,1} = EEG.trials;
    rtrials{subj+1,2} = (EEG.trials/342)*100;
    rtrials{subj+1,3} = EEG_21.trials;
    rtrials{subj+1,4} = (EEG_21.trials/114)*100;
    rtrials{subj+1,5} = EEG_22.trials;
    rtrials{subj+1,6} = (EEG_22.trials/114)*100;
    rtrials{subj+1,7} = EEG_23.trials;
    rtrials{subj+1,8} = (EEG_23.trials/114)*100;
    rtrials{subj+1,9} = file_name{1,subj}(1:5);
end

save remained_trials_for_each_subject.mat rtrials

