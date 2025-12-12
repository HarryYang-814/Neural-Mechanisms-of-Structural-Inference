% This is the propprocessing pipleline for EEG signals used in our 
% preprint paper:
% Yang, Q., Murphy, E., Yang, C., Liao, Y., & Hu, J. (2025). 
% Neural mechanisms of structural inference: an EEG investigation of 
% linguistic phrase structure categorization. 
% bioRxiv, 2025.2007.2001.662085. https://doi.org/10.1101/2025.07.01.662085
% 
% All analyses were carried out in MATLAB R2024b. The preprocessing
% relies on Fieldtrip (Oostenveld et al., 2011).
clear;clc



%% parameter-setting 
cfg                = [];                
cfg.method         = 'wavelet'; 
cfg.keeptrials     =  'yes';
cfg.output         = 'pow';	
cfg.foi            = logspace(log10(2),log10(80),32);
cfg.width     = linspace(2,7,length(cfg.foi)); 
cfg.toi            = -2.6:0.01:2;
cfg.pad            = 'nextpow2';

cfg_base = [];
cfg_base.baseline = [-1.7 -1.2];
cfg_base.baselinetype =   'db'; % 'absolute', 'relative', 'relchange', 'normchange', 'db' or 'z-score'
cfg_base.parameter = 'powspctrm';

load tf_power_2_80Hz.mat
% delta 2-3.5 Hz, theta 4-7 Hz, alpha 8-12 Hz, beta 15-30 Hz, gamma 40-80 Hz
delta_idx = find((TFRwave_VP{1, 1}.freq >= 2) & (TFRwave_VP{1, 1}.freq <= 3.5));
theta_idx = find((TFRwave_VP{1, 1}.freq >= 4) & (TFRwave_VP{1, 1}.freq <= 7));

% time of interest: -1 to 1 s
time_idx  = find((TFRwave_VP{1, 1}.time >= -1) & (TFRwave_VP{1, 1}.time <= 1));

%% VP condition
folder_name = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\8_by_conditions\oscillation\VP\';
files_info = dir(fullfile(folder_name,'*.set'));
files_name = {files_info.name}; 

VP_delta = cell(numel(files_name),1);
VP_theta = cell(numel(files_name),1);
for subj = 1:length(files_name)
    restoredefaultpath; 
    addpath(genpath('E:\Toolboxes_for_matlab\eeglab2024.2\'));
    EEG = pop_loadset('filepath', folder_name, 'filename', files_name{1,subj});
    restoredefaultpath;
    addpath(genpath('E:\Toolboxes_for_matlab\fieldtrip-20250915\fieldtrip-20250915'));
    data              = eeglab2fieldtrip( EEG, 'preprocessing');
    freq              = ft_freqanalysis(cfg, data);
    freq_base         = ft_freqbaseline(cfg_base, freq);

    delta.label       = freq_base.label;
    delta.time        = freq_base.time(time_idx);
    delta.dimord      = 'rpt_chan_time';
    delta.elec        = freq_base.elec;
    delta.powspectrm  = squeeze(mean(freq_base.powspctrm(:,:,delta_idx,time_idx),3));

    theta.label       = freq_base.label;
    theta.time        = freq_base.time(time_idx);
    theta.dimord      = 'rpt_chan_time';
    theta.elec        = freq_base.elec;
    theta.powspectrm  = squeeze(mean(freq_base.powspctrm(:,:,theta_idx,time_idx),3));
    
    VP_delta{subj,1}  = delta;
    VP_theta{subj,1}  = theta;
end


%% NP condition
folder_name = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\8_by_conditions\oscillation\NP\';
files_info = dir(fullfile(folder_name,'*.set'));
files_name = {files_info.name}; 

NP_delta = cell(numel(files_name),1);
NP_theta = cell(numel(files_name),1);
for subj = 1:length(files_name)
    restoredefaultpath; 
    addpath(genpath('E:\Toolboxes_for_matlab\eeglab2024.2\'));
    EEG = pop_loadset('filepath', folder_name, 'filename', files_name{1,subj});
    restoredefaultpath;
    addpath(genpath('E:\Toolboxes_for_matlab\fieldtrip-20250915\fieldtrip-20250915'));
    data              = eeglab2fieldtrip( EEG, 'preprocessing');
    freq              = ft_freqanalysis(cfg, data);
    freq_base         = ft_freqbaseline(cfg_base, freq);

    delta.label       = freq_base.label;
    delta.time        = freq_base.time(time_idx);
    delta.dimord      = 'rpt_chan_time';
    delta.elec        = freq_base.elec;
    delta.powspectrm  = squeeze(mean(freq_base.powspctrm(:,:,delta_idx,time_idx),3));

    theta.label       = freq_base.label;
    theta.time        = freq_base.time(time_idx);
    theta.dimord      = 'rpt_chan_time';
    theta.elec        = freq_base.elec;
    theta.powspectrm  = squeeze(mean(freq_base.powspctrm(:,:,theta_idx,time_idx),3));
    
    NP_delta{subj,1}  = delta;
    NP_theta{subj,1}  = theta;
end

%% CONTROL condition
folder_name = 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\1_preprocessed\8_by_conditions\oscillation\control\';
files_info = dir(fullfile(folder_name,'*.set'));
files_name = {files_info.name}; 

CC_delta = cell(numel(files_name),1);
CC_theta = cell(numel(files_name),1);
for subj = 1:length(files_name)
    restoredefaultpath; 
    addpath(genpath('E:\Toolboxes_for_matlab\eeglab2024.2\'));
    EEG = pop_loadset('filepath', folder_name, 'filename', files_name{1,subj});
    restoredefaultpath;
    addpath(genpath('E:\Toolboxes_for_matlab\fieldtrip-20250915\fieldtrip-20250915'));
    data              = eeglab2fieldtrip( EEG, 'preprocessing');
    freq              = ft_freqanalysis(cfg, data);
    freq_base         = ft_freqbaseline(cfg_base, freq);

    delta.label       = freq_base.label;
    delta.time        = freq_base.time(time_idx);
    delta.dimord      = 'rpt_chan_time';
    delta.elec        = freq_base.elec;
    delta.powspectrm  = squeeze(mean(freq_base.powspctrm(:,:,delta_idx,time_idx),3));

    theta.label       = freq_base.label;
    theta.time        = freq_base.time(time_idx);
    theta.dimord      = 'rpt_chan_time';
    theta.elec        = freq_base.elec;
    theta.powspectrm  = squeeze(mean(freq_base.powspctrm(:,:,theta_idx,time_idx),3));
    
    CC_delta{subj,1}  = delta;
    CC_theta{subj,1}  = theta;
end

save tf_power2power_delta2theta.mat CC_theta CC_delta NP_theta NP_delta VP_theta VP_delta


%% delta-theta power2power coupling analysis
clear;clc
load tf_power2power_delta2theta.mat

LP = [20:22 27:30 36:39];

d_t_coup_z = zeros(numel(VP_delta),3);

% VP
for subj = 1:numel(VP_delta)
    d = mean(VP_delta{subj, 1}.powspectrm(:,LP),2);
    t = mean(VP_theta{subj, 1}.powspectrm(:,LP),2);

    d = zscore(d);
    t = zscore(t);

    r = corr(d, t, 'type','Pearson');
    d_t_coup_z(subj,1) = atanh(r);
end

% NP
for subj = 1:numel(NP_delta)
    d = mean(NP_delta{subj, 1}.powspectrm(:,LP),2);
    t = mean(NP_theta{subj, 1}.powspectrm(:,LP),2);

    d = zscore(d);
    t = zscore(t);

    r = corr(d, t, 'type','Pearson');
    d_t_coup_z(subj,2) = atanh(r);
end

% control
for subj = 1:numel(CC_delta)
    d = mean(CC_delta{subj, 1}.powspectrm(:,LP),2);
    t = mean(CC_theta{subj, 1}.powspectrm(:,LP),2);

    d = zscore(d);
    t = zscore(t);

    r = corr(d, t, 'type','Pearson');
    d_t_coup_z(subj,3) = atanh(r);
end

%% statistical analysis 
% NP vs VP
[~, p_NP_VP, ~, stats_NP_VP] = ttest(d_t_coup_z(:,1), d_t_coup_z(:,2));

% NP vs CC
[~, p_NP_CC, ~, stats_NP_CC] = ttest(d_t_coup_z(:,2), d_t_coup_z(:,3));

% VP vs CC
[~, p_VP_CC, ~, stats_VP_CC] = ttest(d_t_coup_z(:,1), d_t_coup_z(:,3));


%% moving-window analysis

nSub  = numel(NP_delta);   % 被试数量，应该是 31
condNames = {'NP','VP','CC'};
nCond = numel(condNames);

% 把三种条件的 delta/theta 放在一起，方便循环
ALL_delta = {NP_delta, VP_delta, CC_delta};
ALL_theta = {NP_theta, VP_theta, CC_theta};

% 以第一个被试的 NP_delta 当模板
tmp = NP_delta{1};
time = tmp.time;               % 1 x nTime
dt   = time(2) - time(1);      % 采样间隔（秒）
nChan = numel(tmp.label);

roi_idx = [20:22 27:30 36:39];

% 分析时间范围（单位：秒）
tmin = -0.5;
tmax =  0.8;

% 找出对应的时间 index 范围
idx_min = find(time >= tmin, 1, 'first');
idx_max = find(time <= tmax, 1, 'last');

% 滑动时间窗参数
win_len   = 0.100;   % 每个窗口长度 = 200 ms
win_step  = 0.020;   % 步长 = 20 ms

win_samp  = round(win_len  / dt);  % 窗口长度对应的采样点数
step_samp = round(win_step / dt);  % 步长对应的采样点数

% 生成所有时间窗的 index
win_idx   = {};
win_center = [];

start_idx = idx_min;
cnt = 0;
while (start_idx + win_samp - 1) <= idx_max
    cnt = cnt + 1;
    this_idx = start_idx : (start_idx + win_samp - 1);
    win_idx{cnt} = this_idx;
    win_center(cnt) = mean(time(this_idx));  % 记录窗口中心时刻
    start_idx = start_idx + step_samp;
end

nWin = numel(win_idx);
fprintf('总共得到 %d 个时间窗\n', nWin);

% coup_z: nSub x nCond x nWin
coup_z = nan(nSub, nCond, nWin);

for s = 1:nSub
    fprintf('Subject %d / %d\n', s, nSub);
    
    for c = 1:nCond
        
        % 取出该被试该条件的 delta / theta 数据
        S_delta = ALL_delta{c}{s};   % struct, dimord: rpt_chan_time
        S_theta = ALL_theta{c}{s};
        
        % 确认时间向量一致（一般是一致的）
        % assert(isequal(S_delta.time, S_theta.time));
        
        % powspctrm: nTrial x nChan x nTime
        PowD = S_delta.powspectrm;
        PowT = S_theta.powspectrm;
        
        % 在 ROI 内平均 → nTrial x nTime
        PowD_roi = squeeze(mean(PowD(:, roi_idx, :), 2));
        PowT_roi = squeeze(mean(PowT(:, roi_idx, :), 2));
        % 如果 nTrial=1 时 squeeze 会变成 1 x nTime，这个时候可以再包一层处理
        
        nTrial = size(PowD_roi, 1);
        
        for w = 1:nWin
            idx = win_idx{w};
            
            % 在时间窗内对 power 再平均 → 每个 trial 一个值
            d = mean(PowD_roi(:, idx), 2);   % nTrial x 1
            t = mean(PowT_roi(:, idx), 2);   % nTrial x 1
            
            % 去掉可能的 NaN trial
            good = ~isnan(d) & ~isnan(t);
            d = d(good);
            t = t(good);
            
            if numel(d) > 2   % 至少要有 3 个 trial 才能算相关
                % 可选：在条件内 z-score，消除均值/方差的影响
                d = zscore(d);
                t = zscore(t);
                
                r = corr(d, t, 'type', 'Pearson');
                coup_z(s, c, w) = atanh(r);  % Fisher z 变换
            end
        end
    end
end



p_t = nan(1, nWin);
tval = nan(1, nWin);

for w = 1:nWin
    datNP = coup_z(:, 1, w);   % nSub x 1
    datVP = coup_z(:, 2, w);
    
    % 去掉这一窗里是 NaN 的被试
    okSub = ~isnan(datNP) & ~isnan(datVP);
    datNP = datNP(okSub);
    datVP = datVP(okSub);
    
    if numel(datNP) > 2
        [~, p_t(w), ~, stats] = ttest(datNP, datVP); % 配对 t 检验
        tval(w) = stats.tstat;
    end
end

% 对所有时间窗的 p 值做 FDR 校正
q = 0.05;   % FDR 水平
[h_fdr, crit_p] = fdr_bh_simple(p_t, q);

sig_win = find(h_fdr);  % FDR 显著的时间窗 index
fprintf('FDR q=%.2f 下 NP vs VP 显著时间窗数量: %d\n', q, numel(sig_win));