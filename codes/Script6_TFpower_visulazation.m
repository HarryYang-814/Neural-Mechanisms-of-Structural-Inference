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

clear;clc
    restoredefaultpath; 
load 'G:\Research Data\1_Projects\2_Merge and composition\Data and analysis\data_analysis_revised\data\2_statistical_analysis\tf_power_statistics_2_80Hz.mat\';
addpath E:\Toolboxes_for_matlab\fieldtrip-20250915\fieldtrip-20250915
ft_defaults

%% band-specific information
% delta 2-3.5 Hz, theta 4-7 Hz, alpha 8-12 Hz, beta 15-30 Hz, gamma 40-80 Hz
delta    = find((stat.ex.VPvsNP.freq >= 2) & (stat.ex.VPvsNP.freq <= 3.5));
theta    = find((stat.ex.VPvsNP.freq >= 4) & (stat.ex.VPvsNP.freq <= 7));
theta_lo = find((stat.ex.VPvsNP.freq >= 4) & (stat.ex.VPvsNP.freq <= 5));
theta_hg = find((stat.ex.VPvsNP.freq >= 6) & (stat.ex.VPvsNP.freq <= 7));
alpha    = find((stat.ex.VPvsNP.freq >= 8) & (stat.ex.VPvsNP.freq <= 12));
alpha_lo = find((stat.ex.VPvsNP.freq >= 8) & (stat.ex.VPvsNP.freq <= 9.5));
alpha_hg = find((stat.ex.VPvsNP.freq >= 10) & (stat.ex.VPvsNP.freq <= 12));
beta     = find((stat.ex.VPvsNP.freq >= 15) & (stat.ex.VPvsNP.freq <= 30));
gamma    = find((stat.ex.VPvsNP.freq >= 40) & (stat.ex.VPvsNP.freq <= 80));

%% delta
% NPVP
time.delta.NPVP      = double(squeeze(any(any(stat.ex.VPvsNP.mask(:,delta,:) == 1,1),2)));
time.delta.NPVP(:,2) = stat.ex.VPvsNP.time';
chan.delta.NPVP      = num2cell(double(squeeze(any(any(stat.ex.VPvsNP.mask(:,delta,:) == 1,3),2))));
chan.delta.NPVP(:,2) = stat.ex.VPvsNP.label;
% NPCC
time.delta.NPCC         = double(squeeze(any(any(stat.ex.NPvsCC.mask(:,delta,:) == 1,1),2)));
time.delta.NPCC(:,2)    = stat.ex.NPvsCC.time';
chan.delta.NPCC_c1      = num2cell(double(squeeze(any(any(stat.ex.NPvsCC.mask(:,delta,14:56) == 1,3),2))));
chan.delta.NPCC_c1(:,2) = stat.ex.NPvsCC.label;
chan.delta.NPCC_c2      = num2cell(double(squeeze(any(any(stat.ex.NPvsCC.mask(:,delta,143:186) == 1,3),2))));
chan.delta.NPCC_c2(:,2) = stat.ex.NPvsCC.label;
% VPCC
time.delta.VPCC         = double(squeeze(any(any(stat.ex.VPvsCC.mask(:,delta,:) == 1,1),2)));
time.delta.VPCC(:,2)    = stat.ex.VPvsCC.time';
chan.delta.VPCC_c1      = num2cell(double(squeeze(any(any(stat.ex.VPvsCC.mask(:,delta,10:57) == 1,3),2))));
chan.delta.VPCC_c1(:,2) = stat.ex.VPvsCC.label;
chan.delta.VPCC_c2      = num2cell(double(squeeze(any(any(stat.ex.VPvsCC.mask(:,delta,145:177) == 1,3),2))));
chan.delta.VPCC_c2(:,2) = stat.ex.VPvsCC.label;

%% theta (4-7 Hz)
% NPVP
time.theta.NPVP      = double(squeeze(any(any(stat.ex.VPvsNP.mask(:,theta,:) == 1,1),2)));
time.theta.NPVP(:,2) = stat.ex.VPvsNP.time';
chan.theta.NPVP      = num2cell(double(squeeze(any(any(stat.ex.VPvsNP.mask(:,theta,:) == 1,3),2))));
chan.theta.NPVP(:,2) = stat.ex.VPvsNP.label;
% NPCC
time.theta.NPCC         = double(squeeze(any(any(stat.ex.NPvsCC.mask(:,theta,:) == 1,1),2)));
time.theta.NPCC(:,2)    = stat.ex.NPvsCC.time';
chan.theta.NPCC_c1      = num2cell(double(squeeze(any(any(stat.ex.NPvsCC.mask(:,theta,14:72) == 1,3),2))));
chan.theta.NPCC_c1(:,2) = stat.ex.NPvsCC.label;
chan.theta.NPCC_c2      = num2cell(double(squeeze(any(any(stat.ex.NPvsCC.mask(:,theta,141:191) == 1,3),2))));
chan.theta.NPCC_c2(:,2) = stat.ex.NPvsCC.label;
% VPCC
time.theta.VPCC         = double(squeeze(any(any(stat.ex.VPvsCC.mask(:,theta,:) == 1,1),2)));
time.theta.VPCC(:,2)    = stat.ex.VPvsCC.time';
chan.theta.VPCC_c1      = num2cell(double(squeeze(any(any(stat.ex.VPvsCC.mask(:,theta,10:69) == 1,3),2))));
chan.theta.VPCC_c1(:,2) = stat.ex.VPvsCC.label;
chan.theta.VPCC_c2      = num2cell(double(squeeze(any(any(stat.ex.VPvsCC.mask(:,theta,83:118) == 1,3),2))));
chan.theta.VPCC_c2(:,2) = stat.ex.VPvsCC.label;
chan.theta.VPCC_c3      = num2cell(double(squeeze(any(any(stat.ex.VPvsCC.mask(:,theta,147:184) == 1,3),2))));
chan.theta.VPCC_c3(:,2) = stat.ex.VPvsCC.label;
%% low-theta (4-5 Hz)
% NPVP
time.theta_lo.NPVP      = double(squeeze(any(any(stat.ex.VPvsNP.mask(:,theta_lo,:) == 1,1),2)));
time.theta_lo.NPVP(:,2) = stat.ex.VPvsNP.time';
chan.theta_lo.NPVP      = num2cell(double(squeeze(any(any(stat.ex.VPvsNP.mask(:,theta_lo,:) == 1,3),2))));
chan.theta_lo.NPVP(:,2) = stat.ex.VPvsNP.label;
% NPCC
time.theta_lo.NPCC         = double(squeeze(any(any(stat.ex.NPvsCC.mask(:,theta_lo,:) == 1,1),2)));
time.theta_lo.NPCC(:,2)    = stat.ex.NPvsCC.time';
chan.theta_lo.NPCC_c1      = num2cell(double(squeeze(any(any(stat.ex.NPvsCC.mask(:,theta_lo,14:43) == 1,3),2))));
chan.theta_lo.NPCC_c1(:,2) = stat.ex.NPvsCC.label;
chan.theta_lo.NPCC_c2      = num2cell(double(squeeze(any(any(stat.ex.NPvsCC.mask(:,theta_lo,141:187) == 1,3),2))));
chan.theta_lo.NPCC_c2(:,2) = stat.ex.NPvsCC.label;
% VPCC
time.theta_lo.VPCC         = double(squeeze(any(any(stat.ex.VPvsCC.mask(:,theta_lo,:) == 1,1),2)));
time.theta_lo.VPCC(:,2)    = stat.ex.VPvsCC.time';
chan.theta_lo.VPCC_c1      = num2cell(double(squeeze(any(any(stat.ex.VPvsCC.mask(:,theta_lo,10:69) == 1,3),2))));
chan.theta_lo.VPCC_c1(:,2) = stat.ex.VPvsCC.label;
chan.theta_lo.VPCC_c2      = num2cell(double(squeeze(any(any(stat.ex.VPvsCC.mask(:,theta_lo,147:177) == 1,3),2))));
chan.theta_lo.VPCC_c2(:,2) = stat.ex.VPvsCC.label;

%% high-theta (6-7 Hz)
% NPVP
time.theta_hg.NPVP      = double(squeeze(any(any(stat.ex.VPvsNP.mask(:,theta_lo,:) == 1,1),2)));
time.theta_hg.NPVP(:,2) = stat.ex.VPvsNP.time';
chan.theta_hg.NPVP      = num2cell(double(squeeze(any(any(stat.ex.VPvsNP.mask(:,theta_hg,:) == 1,3),2))));
chan.theta_hg.NPVP(:,2) = stat.ex.VPvsNP.label;
% NPCC
time.theta_hg.NPCC         = double(squeeze(any(any(stat.ex.NPvsCC.mask(:,theta_hg,:) == 1,1),2)));
time.theta_hg.NPCC(:,2)    = stat.ex.NPvsCC.time';

% VPCC
time.theta_hg.VPCC         = double(squeeze(any(any(stat.ex.VPvsCC.mask(:,theta_hg,:) == 1,1),2)));
time.theta_hg.VPCC(:,2)    = stat.ex.VPvsCC.time';
chan.theta_hg.VPCC         = num2cell(double(squeeze(any(any(stat.ex.VPvsCC.mask(:,theta_hg,83:118) == 1,3),2))));
chan.theta_hg.VPCC(:,2)    = stat.ex.VPvsCC.label;

%% alpha
% NPVP
time.alpha.NPVP      = double(squeeze(any(any(stat.ex.VPvsNP.mask(:,alpha,:) == 1,1),2)));
time.alpha.NPVP(:,2) = stat.ex.VPvsNP.time';
chan.alpha.NPVP      = num2cell(double(squeeze(any(any(stat.ex.VPvsNP.mask(:,alpha,:) == 1,3),2))));
chan.alpha.NPVP(:,2) = stat.ex.VPvsNP.label;
% NPCC
time.alpha.NPCC = double(squeeze(any(any(stat.ex.NPvsCC.mask(:,alpha,:) == 1,1),2)));
time.alpha.NPCC(:,2) = stat.ex.NPvsCC.time';
chan.alpha.NPCC_c1   = num2cell(double(squeeze(any(any(stat.ex.NPvsCC.mask(:,alpha,11:75) == 1,3),2))));
chan.alpha.NPCC_c1(:,2) = stat.ex.NPvsCC.label;
chan.alpha.NPCC_c2   = num2cell(double(squeeze(any(any(stat.ex.NPvsCC.mask(:,alpha,140:192) == 1,3),2))));
chan.alpha.NPCC_c2(:,2) = stat.ex.NPvsCC.label;
% VPCC
time.alpha.VPCC = double(squeeze(any(any(stat.ex.VPvsCC.mask(:,alpha,:) == 1,1),2)));
time.alpha.VPCC(:,2) = stat.ex.VPvsCC.time';
chan.alpha.VPCC_c1   = num2cell(double(squeeze(any(any(stat.ex.VPvsCC.mask(:,alpha,19:69) == 1,3),2))));
chan.alpha.VPCC_c1(:,2) = stat.ex.VPvsCC.label;
chan.alpha.VPCC_c2   = num2cell(double(squeeze(any(any(stat.ex.VPvsCC.mask(:,alpha,73:122) == 1,3),2))));
chan.alpha.VPCC_c2(:,2) = stat.ex.VPvsCC.label;
chan.alpha.VPCC_c3   = num2cell(double(squeeze(any(any(stat.ex.VPvsCC.mask(:,alpha,146:185) == 1,3),2))));
chan.alpha.VPCC_c3(:,2) = stat.ex.VPvsCC.label;

%% beta
% NPVP
time.beta.NPVP = double(squeeze(any(any(stat.ex.VPvsNP.mask(:,beta,:) == 1,1),2)));
time.beta.NPVP(:,2) = stat.ex.VPvsNP.time';
chan.beta.NPVP      = num2cell(double(squeeze(any(any(stat.ex.VPvsNP.mask(:,beta,:) == 1,3),2))));
chan.beta.NPVP(:,2) = stat.ex.VPvsNP.label;
% NPCC
time.beta.NPCC = double(squeeze(any(any(stat.ex.NPvsCC.mask(:,beta,:) == 1,1),2)));
time.beta.NPCC(:,2) = stat.ex.NPvsCC.time';
chan.beta.NPCC_c1   = num2cell(double(squeeze(any(any(stat.ex.NPvsCC.mask(:,beta,12:61) == 1,3),2))));
chan.beta.NPCC_c1(:,2) = stat.ex.NPvsCC.label;
chan.beta.NPCC_c2   = num2cell(double(squeeze(any(any(stat.ex.NPvsCC.mask(:,beta,138:189) == 1,3),2))));
chan.beta.NPCC_c2(:,2) = stat.ex.NPvsCC.label;
% VPCC
time.beta.VPCC = double(squeeze(any(any(stat.ex.VPvsCC.mask(:,beta,:) == 1,1),2)));
time.beta.VPCC(:,2) = stat.ex.VPvsCC.time';
chan.beta.VPCC_c1   = num2cell(double(squeeze(any(any(stat.ex.VPvsCC.mask(:,beta,20:62) == 1,3),2))));
chan.beta.VPCC_c1(:,2) = stat.ex.VPvsCC.label;
chan.beta.VPCC_c2   = num2cell(double(squeeze(any(any(stat.ex.VPvsCC.mask(:,beta,73:118) == 1,3),2))));
chan.beta.VPCC_c2(:,2) = stat.ex.VPvsCC.label;
chan.beta.VPCC_c3   = num2cell(double(squeeze(any(any(stat.ex.VPvsCC.mask(:,beta,141:176) == 1,3),2))));
chan.beta.VPCC_c3(:,2) = stat.ex.VPvsCC.label;
%% gamma
% NPVP
time.gamma.NPVP = double(squeeze(any(any(stat.ex.VPvsNP.mask(:,gamma,:) == 1,1),2)));
time.gamma.NPVP(:,2) = stat.ex.VPvsNP.time';
% NPCC
time.gamma.NPCC = double(squeeze(any(any(stat.ex.NPvsCC.mask(:,gamma,:) == 1,1),2)));
time.gamma.NPCC(:,2) = stat.ex.NPvsCC.time';
% VPCC
time.gamma.VPCC = double(squeeze(any(any(stat.ex.VPvsCC.mask(:,gamma,:) == 1,1),2)));
time.gamma.VPCC(:,2) = stat.ex.VPvsCC.time';

%%
save tf_power_statistics_2_80Hz.mat stat time chan
%%  ft_singleplotTFR
% display raw effect
cfg = [];
%cfg.channel       = NPvsVP_sig_elecs_cluster1;
cfg.renderer      = 'painters';
cfg.colorbar      = 'yes';
cfg.zlim          = [-0.7 0.7];
cfg.ylim          = [2 30];
cfg.xlim          = [-0.3 1.2];
cfg.colormap      = '*RdBu';
cfg.parameter     = 'powspctrm'; % display the statistical value, i.e. the t-score
ft_singleplotTFR(cfg, grandavg_NP);ft_singleplotTFR(cfg, grandavg_control);



% display t-values
cfg = [];
cfg.channel       = {'TP7'};
cfg.renderer      = 'painters';
cfg.colorbar      = 'yes';
cfg.zlim          = [-5 5];
cfg.ylim          = [2 30];
cfg.maskparameter = 'mask';       % use significance to mask the power difference
cfg.maskstyle      = 'outline';
cfg.colormap      = '*RdBu';
cfg.maskparameter = 'mask';       % use significance to mask the power difference
cfg.parameter     = 'stat'; % display the statistical value, i.e. the t-score
ft_singleplotTFR(cfg, stat.ex.VPvsNP);
title('t-score (not corrected)')
% multiplot display t-values

cfg              = [];
cfg.renderer     = 'painters';
cfg.parameter   = 'stat';
cfg.maskparameter= 'mask';
cfg.maskstyle    = 'outline';
cfg.showlabels   = 'yes';	        
cfg.layout       = 'EEG1005.lay';
cfg.colormap     = '*RdBu'; %The recommended colormaps include 'parula', 'cividis', 'balance', and '*RdBu'.
cfg.colorbar     = 'EastOutside';
cfg.zlim         = [-5 5];
cfg.ylim         = [2 30];
cfg.xlim         = [-1 1];
ft_multiplotTFR(cfg,stat.ex.VPvsNP)
ft_multiplotTFR(cfg,stat.ex.NPvsCC)
ft_multiplotTFR(cfg,stat.ex.VPvsCC)

% display net power effect
%chansel  = find(strcmp(grandavg_NP.label, 'CPz'));
timesel  = find((grandavg_NP.time>=-0.3)&(grandavg_NP.time<=1.2));
power_NP = mean(grandavg_NP.powspctrm(:,:,:,timesel),1);    
power_VP = mean(grandavg_VP.powspctrm(:,:,:,timesel),1);
power_CC = mean(grandavg_control.powspctrm(:,:,:,timesel),1);
power_net = power_NP - power_VP;
siz    = size(power_net);
power_net = reshape(power_net, siz(2:end));
all_stat_VPvsNP.effect = power_net;

power_net_VPCC = power_VP - power_CC;
siz    = size(power_net_VPCC);
power_net_VPCC = reshape(power_net_VPCC, siz(2:end));
all_stat_VPvscontrol.effect = power_net_VPCC;

power_net_NPCC = power_NP - power_CC;
siz    = size(power_net_NPCC);
power_net_NPCC = reshape(power_net_NPCC, siz(2:end));
all_stat_NPvscontrol.effect = power_net_NPCC;

cfg = [];
cfg.channel       = {'PO4'};
cfg.ylim          = [2 30];
cfg.xlim          = [-0.3 1.2];
cfg.colormap      = '*RdBu';
cfg.renderer      = 'openGL';     % painters does not support opacity, openGL does
cfg.colorbar      = 'no';
cfg.parameter     = 'effect';     % display the power
cfg.maskparameter = 'mask';       % use significance to mask the power difference
cfg.maskstyle      = 'outline';
cfg.maskalpha     = 1;          % make non-significant regions 0% visible
cfg.zlim          = [-0.5 0.5];
ft_singleplotTFR(cfg, all_stat_VPvscontrol);
title('PO4')
cb = colorbar;
cb.YTick = [-0.5, 0, 0.5];set(gca, 'FontSize', 20);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
%% number of channels in freq, time points

Num_channels_NPVP = squeeze(sum(stat.ex.VPvsNP.mask, 1));
Num_channels_NPCC = squeeze(sum(stat.ex.NPvsCC.mask, 1));
Num_channels_VPCC = squeeze(sum(stat.ex.VPvsCC.mask, 1));

%% number of significant channels on frequency bands
% NPVP
figure(1);
imagesc(stat.ex.VPvsNP.time, 1:28, Num_channels_NPVP); 
axis xy; 
xlabel('Time (s)'); 
ylabel('Frequency Bands');
title('NP-VP');
cmap = ft_colormap('*RdBu');
colormap(cmap);
colorbar; 
clim([0, 10]);

% NPCC
figure(2);
imagesc(stat.ex.NPvsCC.time, 1:28, Num_channels_NPCC); 
axis xy; 
xlabel('Time (s)'); 
ylabel('Frequency Bands');
title('NP-CC');
cmap = ft_colormap('*RdBu');
colormap(cmap);
colorbar; 
clim([0, 30]);

% VPCC
figure(3);
imagesc(stat.ex.NPvsCC.time, 1:28, Num_channels_VPCC); 
axis xy; 
xlabel('Time (s)'); 
ylabel('Frequency Bands');
title('VP-CC');
cmap = ft_colormap('*RdBu');
colormap(cmap);
colorbar; 
clim([0, 30]);
%% plotting: topographies of the effects size
%% Step 1: 创建新的变量 stat_masked
% 将 all_stat_VPvsNP.stat 中 mask 为0的点置为0
stat_masked = all_stat_VPvsNP.stat;
stat_masked(all_stat_VPvsNP.mask == 0) = 0;

% 假设 stat_masked 的维度为 [channels x freq x time]
[nchan, nfreq, ntime] = size(stat_masked);

%% Step 2: 对每个电极计算 t 值平均值（只考虑非0数据点）
%% 假设 stat_masked 已经由以下语句生成：
% stat_masked = all_stat_VPvsNP.stat;
% stat_masked(all_stat_VPvsNP.mask == 0) = 0;

% 初始化存储每个电极正值和负值之和的变量
sum_t_pos = nan(nchan, 1);
sum_t_neg = nan(nchan, 1);

toi = find((all_stat_VPvsNP.time >= -0.5)& (all_stat_VPvsNP.time <= 1.3));
foi = find((all_stat_VPvsNP.freq >= 8) & (all_stat_VPvsNP.freq <= 12));
for ch = 1:nchan
    % 取出该电极所有频率和时间点的 t 值
    t_vals = squeeze(stat_masked(ch,foi,toi));
    % 筛选出所有非0 t 值
    nonzero_t = t_vals(t_vals ~= 0);
    
    if ~isempty(nonzero_t)
        % 筛选正值，并计算总和
        pos_vals = nonzero_t(nonzero_t > 0);
        if ~isempty(pos_vals)
            sum_t_pos(ch) = sum(pos_vals);
        else
            sum_t_pos(ch) = 0;
        end
        
        % 筛选负值，并计算总和
        neg_vals = nonzero_t(nonzero_t < 0);
        if ~isempty(neg_vals)
            sum_t_neg(ch) = sum(neg_vals);
        else
            sum_t_neg(ch) = 0;
        end
    else
        sum_t_pos(ch) = 0;
        sum_t_neg(ch) = 0;
    end
end

%% topoplot
    restoredefaultpath; 
    addpath(genpath('E:\Toolboxes_for_matlab\eeglab2024.2\'));
%定义数据所在的路径
data_path = 'E:\SUDA\论文撰写\Merge and composition\Data and analysis\data_analysis_new\data\TFR\NP';
%将数据所在的路径定义为工作路径
cd(data_path)
%筛选当前路径下所有的vhdr结尾的文件
files = dir('*.set');
%提取文件名
fn = {files.name};
EEG = pop_loadset('filename',fn{2},'filepath',data_path);

chanlocs = EEG.chanlocs;

topoplot(sum_t_pos, chanlocs, 'electrodes', 'off','style','map');
cmap = ft_colormap('*RdBu');
colormap(cmap);
cb = colorbar;
cb.YTick = [-500, 0, 500];clim([-500, 500]);set(gca, 'FontSize', 20);
