% compute modulation index WN responses. last editted 01.21.2020 ira
% variables gethered with getWNresponses.m

% three states: sit + small pupil, sit + large pupil, running  + large
%

clear; close all; dbstop if error
variables_dir = 'C:\Users\lab\Resilio Sync\Paper1Figures\code\variables';

cd(variables_dir);
load('WNdataLaserOFF.mat'); % 'WNdataLaserOFF.mat (OFF1 -  longer responses), dynamic pupil threshold, but varified.
data = WNdataLaserOFF;
load('WNdataLaserON.mat');
data1 = WNdataLaserON;
load('CellsQualityStats.mat')

maxFRall =[];
SP = nan(length(data),2);
SPL = nan(length(data),2);
meanON = nan(length(data),2); % On response only (0 - 100 ms)
meanONL = nan(length(data),2);
meanWN = nan(length(data),2); % Full response (0 - 600 ms)
meanWNL = nan(length(data),2);
meanPreStim = nan(length(data),2);
meanPreStimL = nan(length(data),2);
depths = nan(length(data),1);
fs = zeros(length(data),1);
Rs = zeros(length(data),1);
WNdirsM =[]; WNcellsM = 0; SSdirsM =[]; SScellsM = 0;
recs = []; cells=[];

evoked = zeros(length(data),4); zstats = zeros(length(data),4);
color1 = [0.7 0.7 0.7]; color2 = [0.2 0.2 0.2];
nreps_check_running = 6; %number of repetitions in each condition for comparison
nreps_check_pupil = 8;
CL = {[0 301], [300 401], [400 601], [600 2000]};

cdVIP; load('Silence_DistanceCorr_dirs.mat')

for cc =1:length(data)
    if data(cc).dir<47
        try
            meanSpikeCount = nanmean([data(cc).SpikeCountWN data1(cc).SpikeCountWN ]); %data(cc).SpikeCountSS data1(cc).SpikeCountSS
        catch
            meanSpikeCount = NaN;
        end
        % exclude cells with very low spikecount, they usually have very
        % large effects; count spikes to WN only
        if meanSpikeCount > 1 && CellsQualityStats.SNR(cc)>.5 && CellsQualityStats.uQ(cc)>10
            
            Spont = [data(cc).mSSon; data(cc).mSSoff]; % spont trials in all states
            PreStim = [data(cc).mNson; data(cc).mNsoff]; %pre stimulus fr
            
            ON = [data(cc).mNON_on; data(cc).mNON_off]; % data(cc).mNON_on; data(cc).mNON_off]; %ON = nanmean(ON);
            [rs, h, stats] = ranksum( ON(:), Spont(:));
            evoked(cc,1) = h;
            zstats(cc,1) = stats.zval;
            
            Sust = [data(cc).mNSustained_on; data(cc).mNSustained_off];
            [rs, h, stats] = ranksum( Sust(:), Spont(:));
            evoked(cc,2) = h;
            zstats(cc,2) = stats.zval;
            
            OFF = [data(cc).mNOFF_on; data(cc).mNOFF_off];
            [rs, h, stats] = ranksum(OFF(:), Spont(:));
            evoked(cc,3) = h;
            zstats(cc,3) = stats.zval;
            
            [rs, h, stats] = ranksum( data(cc).WNresponse, data(cc).SSresponse); %WN response is 0-600ms window
            evoked(cc,4) = h;
            zstats(cc,4) = stats.zval;
            
            % WN laser off trials
            if data(cc).nrepsWNMon > nreps_check_running && data(cc).nrepsWNMoff > nreps_check_running %&& data1(cc).nrepsWNMon > nreps_check_running && data1(cc).nrepsWNMoff > nreps_check_running
                meanWN(cc,1)= nanmean([nanmean(nanmean(data(cc).mWNon))]);% nanmean(nanmean(data(cc).mNSustained_on)) nanmean(nanmean(data(cc).mNOFF_on))]);
                meanWN(cc,2)= nanmean([nanmean(nanmean(data(cc).mWNoff))]);% nanmean(nanmean(data(cc).mNSustained_off)) nanmean(nanmean(data(cc).mNOFF_off))]);
                meanON(cc,1) = nanmean([nanmean(nanmean(data(cc).mNON_on))]);%-nanmean(nanmean(data(cc).mNson));
                meanON(cc,2) = nanmean([nanmean(nanmean(data(cc).mNON_off))]);%-nanmean(nanmean(data(cc).mNsoff));
                meanPreStim(cc,1) = nanmean([nanmean(nanmean(data(cc).mNson))]);
                meanPreStim(cc,2) = nanmean([nanmean(nanmean(data(cc).mNsoff))]);
                WNdirsM = [WNdirsM data(cc).dir]; % collect directory (rec#) for all recordings included in the analysis to make sure N is ok.
                WNcellsM = WNcellsM + 1; % count number of cells included in the analysis
            end
            
            % WN laser on trials
            if data1(cc).nrepsWNMon > nreps_check_running && data1(cc).nrepsWNMoff > nreps_check_running && data(cc).nrepsWNMon > nreps_check_running && data(cc).nrepsWNMoff > nreps_check_running
                meanWNL(cc,1)= nanmean(nanmean(data1(cc).mWNon));%nanmean([nanmean(nanmean(data1(cc).mNON_on)) nanmean(nanmean(data1(cc).mNSustained_on)) nanmean(nanmean(data1(cc).mNOFF_on))]);
                meanWNL(cc,2)= nanmean(nanmean(data1(cc).mWNoff));%nanmean([nanmean(nanmean(data1(cc).mNON_off)) nanmean(nanmean(data1(cc).mNSustained_off)) nanmean(nanmean(data1(cc).mNOFF_off))]);
                meanONL(cc,1) =nanmean([nanmean(nanmean(data1(cc).mNON_on))]);%-nanmean(nanmean(data1(cc).mNson));
                meanONL(cc,2) =nanmean([nanmean(nanmean(data1(cc).mNON_off))]);%-nanmean(nanmean(data1(cc).mNsoff));
                meanPreStimL(cc,1) = nanmean([nanmean(nanmean(data1(cc).mNson))]);
                meanPreStimL(cc,2) = nanmean([nanmean(nanmean(data1(cc).mNsoff))]);
            end
            
            % silent sound laser off
            if data(cc).nrepsSSMon > nreps_check_running && data(cc).nrepsSSMoff > nreps_check_running %&& data1(cc).nrepsSSMon > nreps_check_running && data1(cc).nrepsSSMoff > nreps_check_running
                SP(cc,1) = nanmean(nanmean(data(cc).mSSon)); SP(cc,2) = nanmean(nanmean(data(cc).mSSoff));
                
                WNdirsM = [WNdirsM data(cc).dir]; % collect directory (rec#) for all recordings included in the analysis to make sure N is ok.
                WNcellsM = WNcellsM + 1; % count number of cells included in the analysis
            end
            
            % silent sound laser on
            if data1(cc).nrepsSSMon > nreps_check_running && data1(cc).nrepsSSMoff > nreps_check_running && data(cc).nrepsSSMon > nreps_check_running && data(cc).nrepsSSMoff > nreps_check_running
                SPL(cc,1) = nanmean(nanmean(data1(cc).mSSon)); SPL(cc,2) = nanmean(nanmean(data1(cc).mSSoff));
            end
        else
            evoked(cc,:,:,:) = [0 0 0 0];
            zstats(cc, :,:,:) = [NaN NaN NaN NaN];
        end % inclusion criteria
        
        if isempty(data(cc).maxFR) || data(cc).maxFR == 0
            data(cc).maxFR= NaN;
        end
        
        maxFRall = [maxFRall data(cc).maxFR];
        
        try
            depths(cc) = data(cc).depth;
            if data(cc).width < 0.71 && data(cc).endslope<0.01 || data(cc).width<0.31
                fs(cc) = 1;
            else
                Rs(cc) = 1;
            end
        end
    end
    recs = [recs data(cc).dir];
    cells = [cells data(cc).cell];
end

allDirsM = unique(WNdirsM);

% which cells to include? 1- On response, 2 - Sustained, 3-  Off response,
% 4 - full 600 ms. zstats > 0 activated, <0 suppressed
% 
evoked1 = logical(evoked(:,1) & zstats(:,4)>0);
%evoked1 = logical(ones(length(evoked),1));
WN1 = meanON(evoked1,:);
WN1L = meanONL(evoked1,:);
%WN1 = meanWN(evoked1,:);
%WN1L = meanWNL(evoked1,:);
SP1 = SP(evoked1,:);
SP1L = SPL(evoked1,:);
meanPreStim1 = meanPreStim(evoked1,:);
meanPreStim1L = meanPreStimL(evoked1,:);
depths1 = depths(evoked1);
fs1 = fs(evoked1);
rs1 = Rs(evoked1);
rs1 = logical(rs1); fs1 = logical(fs1);

EvokedWN = WN1;
modulation_indx1 = (EvokedWN(:,1) - EvokedWN(:,2))./(EvokedWN(:,1) + EvokedWN(:,2));
modulation_indx1_sp = (SP1(:,1) - SP1(:,2))./(SP1(:,1) + SP1(:,2));

%modulation index by cell type
modulation_indx1_rs = (WN1(rs1,1) -WN1(rs1,2))./(WN1(rs1,1) + WN1(rs1,2)); %running
modulation_indx1_fs = (WN1(fs1,1) -WN1(fs1,2))./(WN1(fs1,1) + WN1(fs1,2)); %running

modulation_indx1_sp_rs = (SP1(rs1,1) -SP1(rs1,2))./(SP1(rs1,1) + SP1(rs1,2)); %running
modulation_indx1_sp_fs = (SP1(fs1,1) -SP1(fs1,2))./(SP1(fs1,1) + SP1(fs1,2)); %running

modulation_indx1_rsL = (WN1L(rs1,1) -WN1L(rs1,2))./(WN1L(rs1,1) + WN1L(rs1,2)); %running
modulation_indx1_fsL = (WN1L(fs1,1) -WN1L(fs1,2))./(WN1L(fs1,1) + WN1L(fs1,2)); %running

modulation_indx1_sp_rsL = (SP1L(rs1,1) -SP1L(rs1,2))./(SP1L(rs1,1) + SP1L(rs1,2)); %running
modulation_indx1_sp_fsL = (SP1L(fs1,1) -SP1L(fs1,2))./(SP1L(fs1,1) + SP1L(fs1,2)); %running

% subtract spontaneous activity from evoked
WN2 = WN1 - meanPreStim1; WN2L = WN1L - meanPreStim1L;
modulation_indx2 = (WN2(:,1) -WN2(:,2))./(WN2(:,1) + WN2(:,2));
modulation_indx2L = (WN2L(:,1) -WN2L(:,2))./(WN2L(:,1) + WN2L(:,2));

%% plot the results

figure; hold on
bar([1 2], [nanmean(modulation_indx1_sp) nanmean(modulation_indx1)], 'BarWidth', .4)
xticks([1:2])
xticklabels({'spont', 'response'})
ylabel('Modulation Index (mean/SEM)')
errorbar( [1 2], [nanmean(modulation_indx1_sp) nanmean(modulation_indx1)],  ...
    [sem(modulation_indx1_sp) sem(modulation_indx1)])
title('Effects of running on spont and evoked activity')

% plot distribution of modulation index for spont and evoked
figure; subplot(1,2,1); hold on
hist(modulation_indx1_sp,[-1: .05:1], 'Color', [0.8 0.8 0.8])
xlim([-1 1]); xlabel('Modulation Index'); ylabel('Number of cells')
plot([nanmean(modulation_indx1_sp) nanmean(modulation_indx1_sp)], [0 max(ylim)], 'r-', 'Linewidth', 2)
plot([nanmedian(modulation_indx1_sp) nanmedian(modulation_indx1_sp)], [0 max(ylim)], 'c-', 'Linewidth', 2)
title('spont'); legend('hist','mean','median')
subplot(1,2,2); hold on
hist(modulation_indx1, [-1: .05:1], 'Color', [0.8 0.8 0.8])
xlim([-1 1]); xlabel('Modulation Index'); ylabel('Number of cells')
plot([nanmean(modulation_indx1) nanmean(modulation_indx1)], [0 max(ylim)], 'r-', 'Linewidth', 2)
plot([nanmedian(modulation_indx1) nanmedian(modulation_indx1)], [0 max(ylim)], 'c-', 'Linewidth', 2)
title('response'); legend('hist','mean','median')
set(gcf, 'Position',  [17 558 1223 420])

% Is the effect of running the same across recordings? Cell types?
recs1= (recs(evoked1));
recs_rs= (recs1(rs1));
recs_fs=(recs1(fs1));

figure; plot( recs_rs,modulation_indx1_sp_rs, 'ko');
nanindx = find(isnan(modulation_indx1_sp_rs)==1);
x=modulation_indx1_sp_rs;
[p,tbl1,stats] = kruskalwallis(x, recs_rs);
c = multcompare(stats);
title('Spont RS, running')
x=modulation_indx1_rs;
[p,tbl1,stats] = kruskalwallis(x, recs_rs);
c = multcompare(stats);
title('Evoked RS, running')

x=modulation_indx1_sp_fs;
[p,tbl1,stats] = kruskalwallis(x, recs_fs);
c = multcompare(stats);
title('Spont FS, running')
x=modulation_indx1_fs;
[p,tbl1,stats] = kruskalwallis(x, recs_fs);
c = multcompare(stats);
title('Evoked FS, running')

x=modulation_indx1_sp;
[p,tbl1,stats] = kruskalwallis(x, recs1);
c = multcompare(stats);
title('Spont all cells, running')

x=modulation_indx1;
[p,tbl1,stats] = kruskalwallis(x, recs1);
c = multcompare(stats);
title('Evoked all cells, running')
[m,s]=grpstats(x,recs1,{'mean','sem'})

%[m,s]=grpstats(x,recs_rs,{'mean','sem'})
r =unique(recs_rs);
%%
MODULATION_INDX1 = []; MODULATION_INDX1L = [];  MODULATION_INDX1Laser = []; %  for stats across cortical layers, evoked
MODULATION_sp_INDX1 = []; MODULATION_sp_INDX1L = []; MODULATION_sp_INDX1Laser = []; %  for stats across cortical layers, spont
cells1 = cells(evoked1); recs1 = recs(evoked1); depths1 = depths(evoked1);

indx1 = find(modulation_indx2<-1); modulation_indx2(indx1) = -1;
indx2 = find(modulation_indx2> 1);modulation_indx2(indx2) = 1;

for cl = 1:length(CL)
    layer = CL{cl}; % layer depth limits
    indx = find(depths1 >layer(1) & depths1 <layer(2)); % indices of cells within this layer
    recs2 = recs1(indx); % recordings in this layer
    cells2 = cells1(indx); % cells in this layer
    fs2 = fs1(indx); fs2 = logical(fs2); % fast spiking cells in this layer
    rs2 = rs1(indx); rs2 = logical(rs2); % regular spiking cells in this layer
    
    % % evoked (WN respone)
    modulation_indx1 = (WN1(indx,1) -WN1(indx,2))./(WN1(indx,1) + WN1(indx,2)); %running, laser off trials
    modulation_indx1L = (WN1L(indx,1) -WN1L(indx,2))./(WN1L(indx,1) + WN1L(indx,2)); %running, laser on trials
    modulation_indx1Laser = (nanmean(WN1L(indx,:),2) - nanmean(WN1(indx,:),2))./(nanmean(WN1L(indx,:),2) + nanmean(WN1(indx,:),2)); % effect of layer in this layer across states
    
    % means
    meanMI1(cl) = nanmean(modulation_indx1); % mean, laser off
    meanMI1L(cl) = nanmean(modulation_indx1L); % mean, laser on
    meanMI1_Laser(cl) = nanmean(modulation_indx1Laser); % mean laser effect
    
    meanMI1_rs(cl) = nanmean(modulation_indx1(rs2));
    meanMI1_fs(cl) = nanmean(modulation_indx1(fs2));
    meanMI1_rsL(cl) = nanmean(modulation_indx1L(rs2));
    meanMI1_fsL(cl) = nanmean(modulation_indx1L(fs2));
    
    semMI1(cl) = sem(modulation_indx1);
    semMI1L(cl) = sem(modulation_indx1L);
    semMI1_Laser(cl) = sem(modulation_indx1Laser);
    
    semMI1_rs(cl) = sem(modulation_indx1(rs2));
    semMI1_fs(cl) = sem(modulation_indx1(fs2));
    semMI1_rsL(cl) = sem(modulation_indx1L(rs2));
    semMI1_fsL(cl) = sem(modulation_indx1L(fs2));
    
    % collect modulation indices across layers for further stats [MI cl]
    MODULATION_INDX1 = [MODULATION_INDX1; modulation_indx1 ones(length(indx),1)*cl];
    MODULATION_INDX1L = [MODULATION_INDX1L; modulation_indx1L ones(length(indx),1)*cl];
    MODULATION_INDX1Laser = [MODULATION_INDX1Laser; modulation_indx1Laser ones(length(indx),1)*cl];
    
    n_layer_evoked(cl) = length(find(~isnan(modulation_indx1)==1));
    n_layer_evokedL(cl) = length(find(~isnan(modulation_indx1L)==1));
    n_layer_evokedLaser(cl) = length(find(~isnan(modulation_indx1Laser)==1));
    
    rs_cl(cl) = sum(~isnan(modulation_indx1(rs2))); % number of regular spiking cells in this layer
    fs_cl(cl) = sum(~isnan(modulation_indx1(fs2))); % number of fast spiking cells in this layer
    
    % addative effect
    MI_evoked_plus = (WN1L(indx,1) - WN1(indx,2))./(WN1L(indx,1) + WN1(indx,2));
    meanMI_evoked_plus(cl) = nanmean(MI_evoked_plus);
    semMI_evoked_plus(cl) = sem(MI_evoked_plus);
    n_layer_evoked_plus(cl) = length(find(~isnan(MI_evoked_plus)==1)); %number of cell
    
    % evoked - spont
    
    meanMI2(cl) = nanmean(modulation_indx2(indx)); % mean, laser off
    meanMI2L(cl) = nanmean(modulation_indx2L(indx)); % mean, laser on
    
    meanMI2(cl) = nanmean(modulation_indx2(indx));
    meanMI2L(cl) = nanmean(modulation_indx2L(indx));
    
    semMI2(cl) = sem(modulation_indx2(indx));
    semMI2L(cl) = sem(modulation_indx2L(indx));
    
    % % spont activity
    modulation_indx1_sp = (SP1(indx,1) -SP1(indx,2))./(SP1(indx,1) + SP1(indx,2)); %running + laser
    modulation_indx1_spL = (SP1L(indx,1) -SP1L(indx,2))./(SP1L(indx,1) + SP1L(indx,2)); %running + laser
    modulation_indx1_spLaser = (nanmean(SP1L(indx,:),2) - nanmean(SP1(indx,:),2))./(nanmean(SP1L(indx,:),2) + nanmean(SP1(indx,:),2)); %effects of laser across condition
    
    meanMI1_sp(cl) = nanmean(modulation_indx1_sp);
    meanMI1_spL(cl) = nanmean(modulation_indx1_spL);
    meanMI1_spLaser(cl) = nanmean(modulation_indx1_spLaser);
    meanMI1_sp_rs(cl) = nanmean(modulation_indx1_sp(rs2));
    meanMI1_spL_rs(cl) = nanmean(modulation_indx1_spL(rs2));
    meanMI1_spLaser_rs(cl) = nanmean(modulation_indx1_spLaser(rs2));
    meanMI1_sp_fs(cl) = nanmean(modulation_indx1_sp(fs2));
    meanMI1_spL_fs(cl) = nanmean(modulation_indx1_spL(fs2));
    meanMI1_spLaser_fs(cl) = nanmean(modulation_indx1_spLaser(fs2));
    
    semMI1_sp(cl) = sem(modulation_indx1_sp);
    semMI1_spL(cl) = sem(modulation_indx1_spL);
    semMI1_spLaser(cl) = sem(modulation_indx1_spLaser);
    semMI1_sp_rs(cl) = sem(modulation_indx1_sp(rs2));
    semMI1_spL_rs(cl) = sem(modulation_indx1_spL(rs2));
    semMI1_spLaser_rs(cl) = sem(modulation_indx1_spLaser(rs2));
    semMI1_sp_fs(cl) = sem(modulation_indx1_sp(fs2));
    semMI1_spL_fs(cl) = sem(modulation_indx1_spL(fs2));
    semMI1_spLaser_fs(cl) = sem(modulation_indx1_spLaser(fs2));
    
    n_layer_sp(cl) = length(find(~isnan(modulation_indx1_sp)==1)); %number of cells in each layer
    n_layer_spL(cl) = length(find(~isnan(modulation_indx1_spL)==1));
    n_layer_spLaser(cl) = length(find(~isnan(modulation_indx1_spLaser)==1));
    
    rs_sp_cl(cl) = sum(~isnan(modulation_indx1_sp(rs2)));
    fs_sp_cl(cl) = sum(~isnan(modulation_indx1_sp(fs2)));
    
    % additive effect
    MI_sp_plus = (SP1L(indx,1) - SP1(indx,2))./(SP1L(indx,1) + SP1(indx,2));
    meanMI_sp_plus(cl) = nanmean(MI_sp_plus);
    semMI_sp_plus(cl) = sem(MI_sp_plus);
    n_layer_sp_plus(cl) = length(find(~isnan(MI_sp_plus)==1));
    
    MODULATION_sp_INDX1 = [MODULATION_sp_INDX1; modulation_indx1_sp ones(length(indx),1)*cl];
    MODULATION_sp_INDX1L = [MODULATION_sp_INDX1L; modulation_indx1_spL ones(length(indx),1)*cl];
    MODULATION_sp_INDX1Laser = [MODULATION_sp_INDX1Laser; modulation_indx1_spLaser ones(length(indx),1)*cl];
end

% state by layer without layer
figure; hold on; subplot(2,1,1); hold on
errorbar([1:4], meanMI1, semMI1, 'ko-');
xlabel('Cortical Layer')
plot([0 5], [0 0], '--')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Modulation Index');
title_string = sprintf('WN On response, n = %d, %d, %d, %d',  n_layer_evoked);
title(title_string)
% [p,tbl1,stats] = kruskalwallis(mi, layers);
% c = multcompare(stats);
subplot(2,1,2); hold on
errorbar([1:4], meanMI1_sp, semMI1_sp, 'ko-');
xlabel('Cortical Layer')
plot([0 5], [0 0], '--')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Modulation Index');
title_string = sprintf('Spont On response, n = %d, %d, %d, %d',  n_layer_sp);
title(title_string)

%State by layer with laser
figure; hold on; subplot(2,1,1); hold on
errorbar([1:4], meanMI1L, semMI1L, 'bo-');
xlabel('Cortical Layer')
plot([0 5], [0 0], 'k--')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Modulation Index');
title_string = sprintf('WN On response, laser on, n = %d, %d, %d, %d',  n_layer_evokedL);
title(title_string)
% [p,tbl1,stats] = kruskalwallis(mi, layers);
% c = multcompare(stats);
subplot(2,1,2); hold on
errorbar([1:4], meanMI1_spL, semMI1_spL, 'bo-');
xlabel('Cortical Layer')
plot([0 5], [0 0], 'k--')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Modulation Index');
title_string = sprintf('Spont On response, laser on, n = %d, %d, %d, %d',  n_layer_spL);
title(title_string)

figure; hold on; subplot(2,1,1); hold on
errorbar([1:4], meanMI1_Laser, semMI1_Laser, 'co-');
xlabel('Cortical Layer')
plot([0 5], [0 0], 'k--')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Modulation Index - Laser effect');
title_string = sprintf('WN On response, n = %d, %d, %d, %d',  n_layer_evokedLaser);
title(title_string)
% [p,tbl1,stats] = kruskalwallis(mi, layers);
% c = multcompare(stats);
subplot(2,1,2); hold on
errorbar([1:4], meanMI1_spLaser, semMI1_spLaser, 'co-');
xlabel('Cortical Layer')
plot([0 5], [0 0], 'k--')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Modulation Index - Laser effect');
title_string = sprintf('Spont On response, n = %d, %d, %d, %d',  n_layer_spLaser);
title(title_string)

% RS FS WN SP1 by layer
figure; hold on; subplot(2,1,1); hold on
errorbar([1:4], meanMI1_rs, semMI1_rs, 'ko');
errorbar([1.1:4.1], meanMI1_fs, semMI1_fs, 'k*');
plot([0 5], [0 0], 'k--')
legend({ 'run rs', 'run fs'})
xlabel('Cortical Layer')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Modulation Index'); title('WN responses')

hold on; subplot(2,1,2); hold on
errorbar([1:4], meanMI1_sp_rs, semMI1_sp_rs, 'ko');
errorbar([1.1:4.1], meanMI1_sp_fs, semMI1_sp_fs, 'k*');
plot([0 5], [0 0], 'k--')
legend({ 'run rs', 'run fs'})
xlabel('Cortical Layer')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Modulation Index'); title('Spont Activity')

% compare running, laser, and running + layer by layer
figure; hold on; subplot(2,1,1); hold on
errorbar([1:4], meanMI1, semMI1, 'ko');
errorbar([1.1:4.1], meanMI1_Laser, semMI1_Laser, 'bo');
errorbar([1.2:4.2], meanMI_evoked_plus, semMI_evoked_plus, 'go');
plot([0 5], [0 0], 'k--')
legend({ 'run', 'laser' 'run + laser' })
xlabel('Cortical Layer')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Modulation Index - Run'); title('Evoked responses')
subplot(2,1,2); hold on
errorbar([1:4], meanMI1_sp, semMI1_sp, 'ko');
errorbar([1.1:4.1], meanMI1_spLaser, semMI1_spLaser, 'bo');
errorbar([1.2:4.2], meanMI_sp_plus, semMI_sp_plus, 'go');
plot([0 5], [0 0], 'k--')
legend({ 'run', 'laser' 'run + laser' })
xlabel('Cortical Layer')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Modulation Index - Run'); title('Spont responses')

% subtract spont activity from evoked.
figure; hold on;  hold on
errorbar([1:4], meanMI2, semMI2, 'ko-');
xlabel('Cortical Layer')
plot([0 5], [0 0], 'k--')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Modulation Index');
title_string = sprintf('WN On response (w/o spont), n = %d, %d, %d, %d',  n_layer_evoked);
title(title_string)


% plot firing rate distributions
figure; 
h1 = histogram(WN1(:,1));
h1.Normalization = 'probability';
h1.BinWidth = 3;
hold on;
h2 = histogram(WN1(:,2));
h2.Normalization = 'probability';
h2.BinWidth = 3;

figure; 
h1 = histogram(SP1(:,1));
h1.BinWidth = .5;
hold on;
h2 = histogram(SP1(:,2));
h2.BinWidth = .5;
ylabel('Number of cells')
xlabel('Firing Rate (Hz)')
legend('running', 'sitting')

%%
% sound modulation index
MI_sound_run = (WN1(:,1) - SP1(:,1))./ (WN1(:,1) +SP1(:,1));
MI_sound_sit = (WN1(:,2) - SP1(:,2))./ (WN1(:,2) +SP1(:,2));
MI_sound_run_laser = (WN1L(:,1) - SP1L(:,1))./ (WN1L(:,1) +SP1L(:,1));
MI_sound_sit_laser = (WN1L(:,2) - SP1L(:,2))./ (WN1L(:,2) +SP1L(:,2));
[h1,x] = hist(MI_sound_sit, [-1:.1:1]);
h = smooth(h1,3);
figure; hold on
plot(x,h);
[h2,x] = hist(MI_sound_run, [-1:.1:1]);
h = smooth(h2,3);
plot(x,h, '--');
title(' all On responses, laser off')
legend('sitting - laser off', 'running - laser off')
xlabel('Modulation Index - Sound')
ylabel('Number of cells')
figure; hold on;
[h3,x] = hist(MI_sound_sit_laser, [-1:.1:1]);
h = smooth(h3,3);
plot(x,h);
[h4,x] = hist(MI_sound_run_laser, [-1:.1:1]);
h = smooth(h4,3);
plot(x,h, '--');

legend( 'sitting - laser on', 'running - laser on')
xlabel('Modulation Index - Sound')
ylabel('Number of cells')
title(' all On responses, laser off')

fprintf('running vs sitting, laser off')
[p,h stats] =ranksum(MI_sound_run, MI_sound_sit)

fprintf('running vs sitting, laser on')
[p,h stats] =ranksum(MI_sound_run_laser, MI_sound_sit_laser)

for cl = 1:length(CL)
    layer = CL{cl}; % layer depth limits
    indx = find(depths1 >layer(1) & depths1 <layer(2)); % indices of cells within this layer
    meanMI_sound_run(cl) = nanmean(MI_sound_run(indx));
    meanMI_sound_sit(cl) = nanmean(MI_sound_sit(indx));
    semMI_sound_run(cl) = sem(MI_sound_run(indx));
    semMI_sound_sit(cl) = sem(MI_sound_sit(indx));
    n_layer(cl) = sum(isnan(MI_sound_run(indx)));
    
    meanMI_sound_runL(cl) = nanmean(MI_sound_run_laser(indx));
    meanMI_sound_sitL(cl) = nanmean(MI_sound_sit_laser(indx));
    semMI_sound_runL(cl) = sem(MI_sound_run_laser(indx));
    semMI_sound_sitL(cl) = sem(MI_sound_sit_laser(indx));
    n_layerL(cl) = sum(isnan(MI_sound_run_laser(indx)));
end

figure; hold on
errorbar([1:4], meanMI_sound_run, semMI_sound_run, 'ko--');
errorbar([1:4], meanMI_sound_sit, semMI_sound_sit, 'ko-');
xlabel('Cortical Layer')
plot([0 5], [0 0], 'k--')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Modulation Index (mean/SEM)- sound effect');
title_string = sprintf( 'On response MI, n = %d, %d, %d, %d',  n_layer);
title(title_string)
legend('running', 'sitting')

figure; hold on
errorbar([1:4], meanMI_sound_runL, semMI_sound_runL, 'co--');
errorbar([1:4], meanMI_sound_sitL, semMI_sound_sitL, 'co-');
xlabel('Cortical Layer')
plot([0 5], [0 0], 'k--')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Modulation Index (mean/SEM)- sound effect');
title_string = sprintf( 'On response MI, n = %d, %d, %d, %d',  n_layer);
title(title_string)
legend('running', 'sitting')

figure; hold on
errorbar([1:4], mean([meanMI_sound_sitL; meanMI_sound_runL]), mean([semMI_sound_sitL; semMI_sound_runL]), 'co-');
errorbar([1.1:4.1], mean([meanMI_sound_sit; meanMI_sound_run]), mean([semMI_sound_sit; semMI_sound_run]), 'ko-');
xlabel('Cortical Layer')
plot([0 5], [0 0], 'k--')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Modulation Index (mean/SEM)- sound effect');
title_string = sprintf( 'On response MI, n = %d, %d, %d, %d',  n_layer);
title(title_string)
legend('laser on', 'laser off')

figure; hold on
errorbar([1:4], meanMI_sound_sitL, semMI_sound_sitL, 'co-');
errorbar([1.1:4.1], meanMI_sound_sit, semMI_sound_sit, 'ko-');
xlabel('Cortical Layer')
plot([0 5], [0 0], 'k--')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Modulation Index (mean/SEM)- sound effect');
title_string = sprintf( 'On response MI, n = %d, %d, %d, %d',  n_layer);
title(title_string)
legend('laser on', 'laser off')








fs1 = logical(fs1); rs1 = logical(rs1);
MI1_rs = (WN1(rs1,1) - WN1(rs1,2))./ (WN1(rs1,1) +WN1(rs1,2));
MI1_fs = (WN1(fs1,1) - WN1(fs1,2))./ (WN1(fs1,1) +WN1(fs1,2));

SP1_rs = (SP1(rs1,1) - SP1(rs1,4))./ (SP1(rs1,1) +SP1(rs1,4));
SP2_rs = (SP1(rs1,2) - SP1(rs1,3))./ (SP1(rs1,2) +SP1(rs1,3));
SP1_fs = (SP1(fs1,1) - SP1(fs1,4))./ (SP1(fs1,1) +SP1(fs1,4));
SP2_fs = (SP1(fs1,2) - SP1(fs1,3))./ (SP1(fs1,2) +SP1(fs1,3));

%laser VIP
MI1_rsL = (WN1L(rs1,1) - WN1L(rs1,4))./ (WN1L(rs1,1) +WN1L(rs1,4));
MI2_rsL = (WN1L(rs1,2) - WN1L(rs1,3))./ (WN1L(rs1,2) +WN1L(rs1,3));
MI1_fsL = (WN1L(fs1,1) - WN1L(fs1,4))./ (WN1L(fs1,1) +WN1L(fs1,4));
MI2_fsL = (WN1L(fs1,2) - WN1L(fs1,3))./ (WN1L(fs1,2) +WN1L(fs1,3));
SP1_rsL = (SP1L(rs1,1) - SP1L(rs1,4))./ (SP1L(rs1,1) +SP1L(rs1,4));
SP2_rsL = (SP1L(rs1,2) - SP1L(rs1,3))./ (SP1L(rs1,2) +SP1L(rs1,3));
SP1_fsL = (SP1L(fs1,1) - SP1L(fs1,4))./ (SP1L(fs1,1) +SP1L(fs1,4));
SP2_fsL = (SP1L(fs1,2) - SP1L(fs1,3))./ (SP1L(fs1,2) +SP1L(fs1,3));

figure(106); hold on
bar([0.8 1.8], [nanmean(SP1_rs) nanmean(SP2_rs)], 'BarWidth', .1)
bar([0.9 1.9], [nanmean(MI1_rs) nanmean(MI2_rs)], 'BarWidth', .1)

bar([1  2 ], [ nanmean(SP1_fs) nanmean(SP2_fs) ], 'BarWidth', .1)
bar([1.1 2.1], [nanmean(MI1_fs) nanmean(MI2_fs)], 'BarWidth', .1)

bar([2.8 3.8], [nanmean(SP1_rsL) nanmean(SP2_rsL)], 'BarWidth', .1)
bar([2.9 3.9], [nanmean(MI1_rsL) nanmean(MI2_rsL)], 'BarWidth', .1)

bar([3  4 ], [ nanmean(SP1_fsL) nanmean(SP2_fsL) ], 'BarWidth', .1)
bar([3.1 4.1], [nanmean(MI1_fsL) nanmean(MI2_fsL)], 'BarWidth', .1)

% bar([2.8 2.9 3 3.1 3.8 3.9 4 4.1] , [nanmean(SP1_rsL) nanmean(MI1_rsL) nanmean(SP1_fsL) nanmean(MI1_fsL) ...
%     nanmean(SP2_rsL) nanmean(MI2_rsL) nanmean(SP2_fsL) nanmean(MI2_fsL)])

errorbar( [0.8 0.9 1 1.1 1.8 1.9 2 2.1],[nanmean(SP1_rs) nanmean(MI1_rs) nanmean(SP1_fs) nanmean(MI1_fs)  ...
    nanmean(SP2_rs) nanmean(MI2_rs) nanmean(SP2_fs) nanmean(MI2_fs)], ...
    [sem(SP1_rs) sem(MI1_rs) sem(SP1_fs) sem(MI1_fs)  ...
    sem(SP2_rs) sem(MI2_rs) sem(SP2_fs) sem(MI2_fs)])

errorbar( [2.8 2.9 3 3.1 3.8 3.9 4 4.1],[nanmean(SP1_rsL) nanmean(MI1_rsL) nanmean(SP1_fsL) nanmean(MI1_fsL)  ...
    nanmean(SP2_rsL) nanmean(MI2_rsL) nanmean(SP2_fsL) nanmean(MI2_fsL)], ...
    [sem(SP1_rsL) sem(MI1_rsL) sem(SP1_fsL) sem(MI1_fsL)  ...
    sem(SP2_rsL) sem(MI2_rsL) sem(SP2_fsL) sem(MI2_fsL)])
legend({'spont rs', 'evoked rs', 'spont fs', 'evoked fs'})
xticks([1:4])
xticklabels({'run + laser off', 'sit + large pupil + laser off', 'run + laser on', 'sit + large pupil + laser on'})
ylabel('Modulation Index')

modulation_indx1_rs = (WN1L(rs1,1) -WN1(rs1,1))./(WN1L(rs1,1) + WN1(rs1,1)); %running
modulation_indx2_rs = (WN1L(rs1,2) -WN1(rs1,2))./(WN1L(rs1,2) + WN1(rs1,2)); %large pupil + sit
modulation_indx3_rs = (WN1L(rs1,4) -WN1(rs1,4))./(WN1L(rs1,4) + WN1(rs1,4)); %small pupil + sit
modulation_indx1_fs = (WN1L(fs1,1) -WN1(fs1,1))./(WN1L(fs1,1) + WN1(fs1,1)); %running
modulation_indx2_fs = (WN1L(fs1,2) -WN1(fs1,2))./(WN1L(fs1,2) + WN1(fs1,2)); %large pupil + sit
modulation_indx3_fs = (WN1L(fs1,4) -WN1(fs1,4))./(WN1L(fs1,4) + WN1(fs1,4)); %small pupil + sit

modulation_indx1_sp_rs = (SP1L(rs1,1) -SP1(rs1,1))./(SP1L(rs1,1) + SP1(rs1,1)); %running
modulation_indx2_sp_rs = (SP1L(rs1,2) -SP1(rs1,2))./(SP1L(rs1,2) + SP1(rs1,2)); %large pupil + sit
modulation_indx3_sp_rs = (SP1L(rs1,4) -SP1(rs1,4))./(SP1L(rs1,4) + SP1(rs1,4)); %small pupil + sit
modulation_indx1_sp_fs = (SP1L(fs1,1) -SP1(fs1,1))./(SP1L(fs1,1) + SP1(fs1,1)); %running
modulation_indx2_sp_fs = (SP1L(fs1,2) -SP1(fs1,2))./(SP1L(fs1,2) + SP1(fs1,2)); %large pupil + sit
modulation_indx3_sp_fs = (SP1L(fs1,4) -SP1(fs1,3))./(SP1L(fs1,4) + SP1(fs1,4)); %small pupil + sit

data = [modulation_indx1_rs modulation_indx2_rs modulation_indx2_rs];
cate = [ones(length(modulation_indx1_rs),1)*1 ones(length(modulation_indx2_rs),1)*2 ones(length(modulation_indx3_rs),1)*3];

data = [];
cate = [];
data= [nanmean(modulation_indx3_sp_rs), nanmean(modulation_indx2_sp_rs), nanmean(modulation_indx1_sp_rs); nanmean(modulation_indx3_rs) nanmean(modulation_indx2_rs) nanmean(modulation_indx1_rs);...
    nanmean(modulation_indx3_sp_fs) nanmean(modulation_indx2_sp_fs) nanmean(modulation_indx1_sp_fs);nanmean(modulation_indx3_fs) nanmean(modulation_indx2_fs) nanmean(modulation_indx1_fs)];
cate = 2;
[p,tbl] = anova2(data,cate);
% VIP effect in ech state
figure(210); hold on
bar([.8  1.8  2.8 ], [nanmean(modulation_indx3_sp_rs) nanmean(modulation_indx2_sp_rs) nanmean(modulation_indx1_sp_rs)], 'BarWidth', .1)
bar([.9  1.9  2.9 ], [nanmean(modulation_indx3_sp_fs) nanmean(modulation_indx2_sp_fs) nanmean(modulation_indx1_sp_fs)], 'BarWidth', .1)
bar([1  2  3 ], [nanmean(modulation_indx3_rs) nanmean(modulation_indx2_rs) nanmean(modulation_indx1_rs)], 'BarWidth', .1)
bar([1.1  2.1 3.1 ], [nanmean(modulation_indx3_fs) nanmean(modulation_indx2_fs) nanmean(modulation_indx1_fs)], 'BarWidth', .1)

errorbar([.8 .9 1. 1.1 1.8 1.9 2 2.1 2.8 2.9 3 3.1], [nanmean(modulation_indx3_sp_rs) nanmean(modulation_indx3_sp_fs) nanmean(modulation_indx3_rs) nanmean(modulation_indx3_fs) ...
    nanmean(modulation_indx2_sp_rs) nanmean(modulation_indx2_sp_fs)  nanmean(modulation_indx2_rs)  nanmean(modulation_indx2_fs) ...
    nanmean(modulation_indx1_sp_rs) nanmean(modulation_indx1_sp_fs) nanmean(modulation_indx1_rs) nanmean(modulation_indx1_fs)],...
    [sem(modulation_indx3_sp_rs) sem(modulation_indx3_sp_fs) sem(modulation_indx3_rs) sem(modulation_indx3_fs) ...
    sem(modulation_indx2_sp_rs) sem(modulation_indx2_sp_fs)  sem(modulation_indx2_rs)  sem(modulation_indx2_fs) ...
    sem(modulation_indx1_sp_rs) sem(modulation_indx1_sp_fs) sem(modulation_indx1_rs) sem(modulation_indx1_fs)])
xticks([1:3])
xticklabels({'sit + small pupil', 'sit + large pupil', 'run'})
ylabel('Modulation Index'); title('VIP activation')
legend({'spont rs', 'spont fs', 'evoked rs', 'evoked fs'})

figure(213); all_MI=[];
mi_cl = []; mi_sp_cl = [];
for cl = 1:length(CL)
    layer = CL{cl};
    indx = find(depths2 >layer(1) & depths2 <layer(2));
    fs1 = fs1(indx); fs1 = logical(fs1);
    rs1 = rs1(indx); rs1 = logical(rs1);
    layer_indx(cl).indx = indx;
    mi = (nanmean(WN1L(indx,:),2) - nanmean(WN1(indx,:),2))./(nanmean(WN1L(indx,:),2) + nanmean(WN1(indx,:),2));
    mi_run = (WN1(indx,1) - WN1(indx,4))./(WN1(indx,1) + WN1(indx,4));
    
    meanMI(cl) = nanmean(mi);
    semMI(cl) = sem(mi);
    nonnan_indx = find(~isnan(mi)==1);
    fs_layer(cl) = sum(fs1(nonnan_indx));
    rs_layer(cl) = sum(rs1(nonnan_indx));
    n_layer_VIP(cl,1) = sum(~isnan(mi));
    mi_cl = [ mi_cl; mi ones(length(mi),1)*cl fs1 rs1];
    
    
    figure(300); hold on
    plot(ones(length(mi),1)*cl, mi, 'k.');
    plot(cl, nanmean(mi), 'ro')
    
    meanMI_rs(cl) = nanmean(mi(rs1));
    semMI_rs(cl) = sem(mi(rs1));
    meanMI_fs(cl) = nanmean(mi(fs1));
    semMI_fs(cl) = sem(mi(fs1));
    
    modulation_indx1 = (WN1L(indx,1) -WN1(indx,1))./(WN1L(indx,1) + WN1(indx,1)); %running
    modulation_indx2 = (WN1L(indx,2) -WN1(indx,2))./(WN1L(indx,2) + WN1(indx,2)); %large pupil + sit
    modulation_indx3 = (WN1L(indx,4) -WN1(indx,4))./(WN1L(indx,4) + WN1(indx,4));
    
    
    meanMI1(cl) = nanmean(modulation_indx1);
    semMI1(cl) = sem(modulation_indx1);
    meanMI2(cl) = nanmean(modulation_indx2);
    semMI2(cl) = sem(modulation_indx2);
    meanMI3(cl) = nanmean(modulation_indx3);
    semMI3(cl) = sem(modulation_indx3);
    
    figure(213);
    subplot(1,4,cl)
    plot(mi, mi_run, 'b.')
    lsline; title(sprintf('Layer %d', cl))
    xlabel('VIP MI'); ylabel('Running MI')
    
    %spont
    mi_sp_run = (SP1(indx,1) - SP1(indx,4))./(SP1(indx,1) + SP1(indx,4));
    mi_sp = (nanmean(SP1L(indx,:),2) - nanmean(SP1(indx,:),2))./(nanmean(SP1L(indx,:),2) + nanmean(SP1(indx,:),2));
    all_MI = [all_MI; mi mi_run mi_sp mi_sp_run recs1(indx)' cells1(indx)' ones(length(indx),1)*cl fs1 rs1];
    
    n_layer_VIP(cl,2) = sum(~isnan(mi_sp));
    mi_sp_cl = [ mi_sp_cl; mi_sp ones(length(mi_sp),1)*cl fs1 rs1];
    nonnan_indx = find(~isnan(mi_sp)==1);
    fs_sp_layer(cl) = sum(fs1(nonnan_indx));
    rs_sp_layer(cl) = sum(rs1(nonnan_indx));
    
    figure(213); hold on
    subplot(1,4,cl)
    plot(mi_sp, mi_sp_run, 'k.')
    lsline; title(sprintf('Layer %d', cl))
    [rhoe,pe]=corr(mi, mi_run, 'Type', 'Spearman', 'rows', 'pairwise');
    [rhos,ps]=corr(mi_sp, mi_sp_run, 'Type', 'Spearman', 'rows', 'pairwise');
    title(sprintf('e:%.2f %.4f; sp:%.2f %.4f',rhoe,pe, rhos,ps));
    
    meanMI_sp(cl) = nanmean(mi_sp);
    semMI_sp(cl) = sem(mi_sp);
    
    meanMI_sp_rs(cl) = nanmean(mi_sp(rs1));
    semMI_sp_rs(cl) = sem(mi_sp(rs1));
    meanMI_sp_fs(cl) = nanmean(mi_sp(fs1));
    semMI_sp_fs(cl) = sem(mi_sp(fs1));
    
    modulation_indx1_sp = (SP1L(indx,1) -SP1(indx,1))./(SP1L(indx,1) + SP1(indx,1)); %running
    modulation_indx2_sp = (SP1L(indx,2) -SP1(indx,2))./(SP1L(indx,2) + SP1(indx,2)); %large pupil + sit
    modulation_indx3_sp = (SP1L(indx,4) -SP1(indx,4))./(SP1L(indx,4) + SP1(indx,4)); %small pupil pupil + sit
    
    meanMI1_sp(cl) = nanmean(modulation_indx1_sp);
    semMI1_sp(cl) = sem(modulation_indx1_sp);
    meanMI2_sp(cl) = nanmean(modulation_indx2_sp);
    semMI2_sp(cl) = sem(modulation_indx2_sp);
    meanMI3_sp(cl) = nanmean(modulation_indx3_sp);
    semMI3_sp(cl) = sem(modulation_indx3_sp);
end
figure; hold on
WN1 = nanmean(WN1,2);
errorbar([1:4], [nanmean(WN1(layer_indx(1).indx,:)) nanmean(WN1(layer_indx(2).indx)) nanmean(WN1(layer_indx(3).indx)) nanmean(WN1(layer_indx(4).indx))], ...
    [sem(WN1(layer_indx(1).indx)) sem(WN1(layer_indx(2).indx)) sem(WN1(layer_indx(3).indx)) sem(WN1(layer_indx(4).indx))], 'ro')
xticks([1:4]); xticklabels({'2/3', '4' '5' '6'})
xlabel('Cortical Layers')
ylabel('mean FR (Hz), evoked and spont')

WN1 = nanmean(SP1,2);
errorbar([1:4], [nanmean(WN1(layer_indx(1).indx,:)) nanmean(WN1(layer_indx(2).indx)) nanmean(WN1(layer_indx(3).indx)) nanmean(WN1(layer_indx(4).indx))], ...
    [sem(WN1(layer_indx(1).indx)) sem(WN1(layer_indx(2).indx)) sem(WN1(layer_indx(3).indx)) sem(WN1(layer_indx(4).indx))], 'ko')
xticks([1:4]); xticklabels({'2/3', '4' '5' '6'})
xlabel('Cortical Layers')



mi = (nanmean(WN1L(:,1:4),2) - nanmean(WN1(:,1:4),2))./(nanmean(WN1L(:,1:4),2) + nanmean(WN1(:,1:4),2));
mi_sp = (nanmean(SP1L(:,1:4),2) - nanmean(SP1(:,1:4),2))./(nanmean(SP1L(:,1:4),2) + nanmean(SP1(:,1:4),2));

%mean Laser Effect
figure; hold on
bar([1,2], [nanmean(mi_sp), nanmean(mi)])
errorbar([1,2], [nanmean(mi_sp), nanmean(mi)], [sem(mi_sp), sem(mi)])
xticks([1 2])
xticklabels({'spont','evoked'})
ylabel('VIP MI (mean/sem)')

%mean laser effect by cell type
figure;hold on
bar([1 2], [nanmean(mi_sp(rs1)),  nanmean(mi(rs1))], 'BarWidth', .3)
bar([1.3, 2.3], [ nanmean(mi_sp(fs1)) nanmean(mi(fs1))], 'BarWidth', .3)
errorbar([1,1.3,2,2.3], [nanmean(mi_sp(rs1)), nanmean(mi_sp(fs1)), nanmean(mi(rs1)), nanmean(mi(fs1))], [sem(mi_sp(rs1)),sem(mi_sp(fs1)),sem(mi(rs1)) sem(mi(fs1))])
xticks([1.15 2.15])
xticklabels({'spont','evoked'})
ylabel('VIP MI (mean/sem)')
legend({'rs', 'fs'})

%VIP effect by cell type and layer WN + SP1
figure(211); hold on
errorbar([0.9 1.9 2.9 3.9], [meanMI_sp_rs], [semMI_sp_rs], 'ko')
errorbar([1 2 3 4], [meanMI_sp_fs], [semMI_sp_fs], 'k*')
errorbar([1.1 2.1 3.1 4.1], [meanMI_rs], [semMI_rs], 'ro')
errorbar([1.2 2.2 3.2 4.2], [meanMI_fs], [semMI_fs], 'r*')
xticks([1:4]); xticklabels({'2/3', '4' '5' '6'})
xlabel('Cortical Layers')
ylabel('Modulation Index - VIP'); title('WN responses')
plot([0 5], [0 0], 'k--')
legend({'spont rs', 'spont fs', 'evoked rs', 'evoked fs'})

%VIP effect by layer across cells
figure(212); hold on
errorbar([0.9 1.9 2.9 3.9], [meanMI_sp], [semMI_sp], 'ko')
errorbar([1 2 3 4], [meanMI], [semMI], 'bo')
xticks([1:4]); xticklabels({'2/3', '4' '5' '6'})
xlabel('Cortical Layers')
ylabel('Modulation Index - VIP');
plot([0 5], [0 0], 'k--')
legend({'spont', 'evoked'})

% indx1 = 2:2:4;
% indx2 = 3:4;
indx1 = 4;

VIP_MI_rs_evoked = (nanmean(WN1L(rs1,indx1),2)- nanmean(WN1(rs1,indx1),2))./(nanmean(WN1L(rs1,indx1),2) + nanmean(WN1(rs1,indx1),2));
VIP_MI_fs_evoked = (nanmean(WN1L(fs1,indx1),2)- nanmean(WN1(fs1,indx1),2))./(nanmean(WN1L(fs1,indx1),2) + nanmean(WN1(fs1,indx1),2));
VIP_MI_rs_sp = (nanmean(SP1L(rs1,indx1),2)- nanmean(SP1(rs1,indx1),2))./(nanmean(SP1L(rs1,indx1),2) + nanmean(SP1(rs1,indx1),2));
VIP_MI_fs_sp = (nanmean(SP1L(fs1,indx1),2)- nanmean(SP1(fs1,indx1),2))./(nanmean(SP1L(fs1,indx1),2) + nanmean(SP1(fs1,indx1),2));

RUN_MI_rs_evoked = (WN1(rs1,1) - WN1(rs1,4))./(WN1(rs1,1) + WN1(rs1,4));
RUN_MI_fs_evoked = (WN1(fs1,1) - WN1(fs1,4))./(WN1(fs1,1) + WN1(fs1,4));
RUN_MI_rs_sp = (SP1(rs1,1) - SP1(rs1,4))./(SP1(rs1,1) + SP1(rs1,4));
RUN_MI_fs_sp = (SP1(fs1,1) - SP1(fs1,4))./(SP1(fs1,1) + SP1(fs1,4));
figure; hold on
layer23_indx= find(mi_cl(:,2,:,:)==1);
layer4_indx= find(mi_cl(:,2,:,:)==2);
layer5_indx= find(mi_cl(:,2,:,:)==3);
layer6_indx= find(mi_cl(:,2,:,:)==4);
plot(ones(length(mi_cl(layer23_indx)),1)*1,mi_cl(layer23_indx), 'k.')
plot(ones(length(mi_cl(layer4_indx)),1)*2,mi_cl(layer4_indx), 'k.')
plot(ones(length(mi_cl(layer5_indx)),1)*3,mi_cl(layer5_indx), 'k.')
plot(ones(length(mi_cl(layer6_indx)),1)*4,mi_cl(layer6_indx), 'k.')

figure(216); hold on
subplot(2,1,1); hold on;
plot(VIP_MI_rs_evoked, RUN_MI_rs_evoked, 'k.');
plot(VIP_MI_fs_evoked, RUN_MI_fs_evoked, 'b.');
lsline;
xlabel('Modulation Index - VIP')
ylabel('Modulation Index - Run')
title('Evoked Activity')
legend({'rs', 'fs'})

subplot(2,1,2); hold on;
plot(VIP_MI_rs_sp, RUN_MI_rs_sp, 'k.');
plot(VIP_MI_fs_sp, RUN_MI_fs_sp, 'b.');
lsline;
xlabel('Modulation Index - VIP')
ylabel('Modulation Index - Run')
title('Spont Activity')
legend({'rs', 'fs'})

figure(217); hold on
subplot(2,1,1); hold on;
plot([VIP_MI_rs_evoked; VIP_MI_fs_evoked], [RUN_MI_rs_evoked; RUN_MI_fs_evoked], 'k.');
lsline;
xlabel('Modulation Index - VIP')
ylabel('Modulation Index - Run')
title('Evoked Activity')

subplot(2,1,2); hold on;
plot([VIP_MI_rs_sp; VIP_MI_fs_sp], [RUN_MI_rs_sp; RUN_MI_fs_sp], 'k.');
lsline;
xlabel('Modulation Index - VIP')
ylabel('Modulation Index - Run')
title('Spont Activity')

fprintf('\nCorr between VIP and run evoked and spont all\n')
[rho,p]=corr([VIP_MI_rs_evoked; VIP_MI_fs_evoked], [RUN_MI_rs_evoked; RUN_MI_fs_evoked], 'Type', 'Spearman', 'rows', 'pairwise')
[rho,p]=corr([VIP_MI_rs_sp; VIP_MI_fs_sp], [RUN_MI_rs_sp; RUN_MI_fs_sp], 'Type', 'Spearman', 'rows', 'pairwise')

fprintf('\nCorr between VIP and run evoked and spont fs\n')
[rho,p]=corr([VIP_MI_fs_evoked], [RUN_MI_fs_evoked], 'Type', 'Spearman', 'rows', 'pairwise')
[rho,p]=corr([VIP_MI_fs_sp],  [RUN_MI_fs_sp], 'Type', 'Spearman', 'rows', 'pairwise')

fprintf('\nCorr between VIP and run evoked and spont rs\n')
[rho,p]=corr([VIP_MI_rs_evoked], [RUN_MI_rs_evoked], 'Type', 'Spearman', 'rows', 'pairwise')
[rho,p]=corr([VIP_MI_rs_sp],  [RUN_MI_rs_sp], 'Type', 'Spearman', 'rows', 'pairwise')

% evoked and spont modulation are correlated for running
figure; subplot(2,1,1);
a = (WN1(:,1)-WN1(:,4))./(WN1(:,1)+WN1(:,4));
b = (SP1(:,1)-SP1(:,4))./(SP1(:,1)+SP1(:,4));
[rho,p]=corr(b, a, 'Type', 'Spearman', 'rows', 'pairwise')
plot(b(rs1),a(rs1),'k.')
hold on; plot(b(fs1),a(fs1),'b.')
legend('rs', 'fs')
xlabel('running MI spont'); ylabel('running MI evoked')

% evoked and spont modulation are correlated for pupil
subplot(2,1,2); hold on
a = (WN1(:,2)-WN1(:,3))./(WN1(:,2)+WN1(:,3));
b = (SP1(:,2)-SP1(:,3))./(SP1(:,2)+SP1(:,3));
[rho,p]=corr(b, a, 'Type', 'Spearman', 'rows', 'pairwise')
plot(b(rs1),a(rs1),'k.')
hold on; plot(b(fs1),a(fs1),'b.')
legend('rs', 'fs')
xlabel('pupil MI spont'); ylabel('pupil MI evoked')


% running
WN_rs = (WN1(rs1,1) - WN1(rs1,4))./(WN1(rs1,1) + WN1(rs1,4));
WN_fs = (WN1(fs1,1) - WN1(fs1,4))./(WN1(fs1,1) + WN1(fs1,4));
WN_rsL = (WN1L(rs1,4) - WN1(rs1,4))./(WN1L(rs1,4) + WN1(rs1,4));
WN_fsL = (WN1L(fs1,4) - WN1(fs1,4))./(WN1L(fs1,4) + WN1(fs1,4));
WN_rs_plus = (WN1L(rs1,1) - WN1(rs1,4))./(WN1L(rs1,1) + WN1(rs1,4));
WN_fs_plus = (WN1L(fs1,1) - WN1(fs1,4))./(WN1L(fs1,1) + WN1(fs1,4));
WN_exp_rs = sum([WN_rs WN_rsL],2);
WN_exp_fs = sum([WN_fs WN_fsL],2);
[p,h,zstat] = ranksum(WN_rs_plus, WN_exp_rs);
[p,h,zstat] = ranksum(WN_fs_plus, WN_exp_fs);

SP_rs = (SP1(rs1,1) - SP1(rs1,4))./(SP1(rs1,1) + SP1(rs1,4));
SP_fs = (SP1(fs1,1) - SP1(fs1,4))./(SP1(fs1,1) + SP1(fs1,4));
SP_rsL = (nanmean(SP1L(rs1,2:3),2) - nanmean(SP1(rs1,2:3),2))./(nanmean(SP1L(rs1,2:3),2) + nanmean(SP1(rs1,2:3),2));
SP_fsL = (nanmean(SP1L(fs1,2:3),2) - nanmean(SP1(fs1,2:3),2))./(nanmean(SP1L(fs1,2:3),2) + nanmean(SP1(fs1,2:3),2));
SP_rs_plus = (SP1L(rs1,1) - SP1(rs1,4))./(SP1L(rs1,1) + SP1(rs1,4));
SP_fs_plus = (SP1L(fs1,1) - SP1(fs1,4))./(SP1L(fs1,1) + SP1(fs1,4));
SP_exp_rs = sum([SP_rs SP_rsL],2);
SP_exp_fs = sum([SP_fs SP_fsL],2);

[p,h,zstat] = ranksum(SP_rs_plus, SP_exp_rs);
[p,h,zstat] = ranksum(SP_fs_plus, SP_exp_fs);

figure; subplot(2,1,1);hold on
plot(WN_rs, WN_rsL, 'k.')
plot(WN_fs, WN_fsL, 'b.')
[rho_rs,p_rs]=corr([WN_rs], [WN_rsL], 'Type', 'Spearman', 'rows', 'pairwise')
[rho_fs,p_fs]=corr([WN_fs], [WN_fsL], 'Type', 'Spearman', 'rows', 'pairwise')
title(sprintf('Evoked, rs rho=%.4f,p=%.4f fs rho=%.4f,p=%.4f', rho_rs, p_rs, rho_fs, p_fs))
xlabel('run - VIP'); ylabel('VIP - run'); lsline
subplot(2,1,2);hold on
plot(SP_rs, SP_rsL, 'k.')
plot(SP_fs, SP_fsL, 'b.')
xlabel('run - VIP'); ylabel('VIP - run'); lsline
[rho_rs,p_rs]=corr([SP_rs], [SP_rsL], 'Type', 'Spearman', 'rows', 'pairwise')
[rho_fs,p_fs]=corr([SP_fs], [SP_fsL], 'Type', 'Spearman', 'rows', 'pairwise')
title(sprintf('Spont, rs rho=%.4f,p=%.4f fs rho=%.4f,p=%.4f', rho_rs, p_rs, rho_fs, p_fs))


figure; subplot(2,1,1);hold on
plot(WN_rs_plus, WN_exp_rs, 'k.');
plot(WN_fs_plus, WN_exp_fs, 'b.');
xlabel('VIP + run'); ylabel('Expected'); lsline
[rho_rs,p_rs]=corr(WN_rs_plus, WN_exp_rs, 'Type', 'Spearman', 'rows', 'pairwise')
[rho_fs,p_fs]=corr(WN_fs_plus, WN_exp_fs, 'Type', 'Spearman', 'rows', 'pairwise')
title(sprintf('Evoked, rs rho=%.4f,p=%.4f fs rho=%.4f,p=%.4f', rho_rs, p_rs, rho_fs, p_fs))
subplot(2,1,2);hold on
plot(SP_rs_plus, SP_exp_rs, 'k.');
plot(SP_fs_plus, SP_exp_fs, 'b.');
xlabel('VIP + run'); ylabel('Expected'); lsline
[rho_rs,p_rs]=corr(SP_rs_plus, SP_exp_rs, 'Type', 'Spearman', 'rows', 'pairwise')
[rho_fs,p_fs]=corr(SP_fs_plus, SP_exp_fs, 'Type', 'Spearman', 'rows', 'pairwise')
title(sprintf('Spont, rs rho=%.4f,p=%.4f fs rho=%.4f,p=%.4f', rho_rs, p_rs, rho_fs, p_fs))


numCells = [sum(~isnan(SP_rs)) sum(~isnan(SP_rsL)) sum(~isnan(SP_rs_plus)) sum(~isnan(SP_exp_rs));sum(~isnan(SP_fs)) sum(~isnan(SP_fsL)) sum(~isnan(SP_fs_plus)) sum(~isnan(SP_exp_fs))];
numCells_evoked = [sum(~isnan(WN_rs)) sum(~isnan(WN_rsL)) sum(~isnan(WN_rs_plus)) sum(~isnan(WN_exp_rs));sum(~isnan(WN_fs)) sum(~isnan(WN_fsL)) sum(~isnan(WN_fs_plus)) sum(~isnan(WN_exp_fs))];

figure; hold on
bar([1 2 3], [nanmean(WN_rs), nanmean(WN_rsL), nanmean(WN_rs_plus)],'BarWidth', .3)
bar([1.3 2.3 3.3], [nanmean(WN_fs), nanmean(WN_fsL), nanmean(WN_fs_plus)], 'BarWidth', .3)
bar([4 5 6], [nanmean(SP_rs), nanmean(SP_rsL), nanmean(SP_rs_plus)],'BarWidth', .3)
bar([4.3 5.3 6.3], [nanmean(SP_fs), nanmean(SP_fsL), nanmean(SP_fs_plus)], 'BarWidth', .3)
xticks(1:6);
xticklabels({'run MI - VIP', 'VIP MI - run', 'run + VIP MI', 'run MI - VIP', 'VIP MI - run', 'run + VIP MI'})
errorbar([1 2 3 1.3 2.3 3.3 4 5 6 4.3 5.3 6.3], [nanmean(WN_rs), nanmean(WN_rsL), nanmean(WN_rs_plus), nanmean(WN_fs), nanmean(WN_fsL), nanmean(WN_fs_plus),...
    nanmean(SP_rs), nanmean(SP_rsL), nanmean(SP_rs_plus), nanmean(SP_fs), nanmean(SP_fsL), nanmean(SP_fs_plus)], ...
    [sem(WN_rs), sem(WN_rsL), sem(WN_rs_plus), sem(WN_fs), sem(WN_fsL), sem(WN_fs_plus),...
    sem(SP_rs), sem(SP_rsL), sem(SP_rs_plus), sem(SP_fs), sem(SP_fsL), sem(SP_fs_plus)])
plot([3 6], [nanmean(WN_exp_rs)  nanmean(SP_exp_rs)], 'ko', 'MarkerSize', 10)
plot([3.3 6.3], [nanmean(WN_exp_fs)  nanmean(SP_exp_fs)], 'ko', 'MarkerSize', 10)
legend({'evoked rs', 'evoked fs', 'spont rs', 'spont fs'})
ylabel('Modulation Index')

[p,h,zstat] = ranksum(SP_rs_plus, SP_exp_rs);
[p,h,zstat] = ranksum(SP_fs_plus, SP_exp_fs);

n_boots =100;
m_exp = bootstrp(n_boots, @nanmean, WN_exp_rs);
m_plus = bootstrp(n_boots, @nanmean, WN_rs_plus);
[p,h,zstat] = ranksum(m_plus, m_exp);
Stats(1,:,:) = [p, zstat.zval];
m_exp = bootstrp(n_boots, @nanmean, WN_exp_fs);
m_plus = bootstrp(n_boots, @nanmean, WN_fs_plus);
[p,h,zstat] = ranksum(m_plus, m_exp);
Stats(2,:,:) = [p, zstat.zval];

m_exp = bootstrp(n_boots, @nanmean, SP_exp_rs);
m_plus = bootstrp(n_boots, @nanmean, SP_rs_plus);
[p,h,zstat] = ranksum(m_plus, m_exp);
Stats(3,:,:) = [p, zstat.zval];
m_exp = bootstrp(n_boots, @nanmean, SP_exp_fs);
m_plus = bootstrp(n_boots, @nanmean, SP_fs_plus);
[p,h,zstat] = ranksum(m_plus, m_exp);
Stats(4,:,:) = [p, zstat.zval];

% pupil

WN_rs = (WN1(rs1,2) - WN1(rs1,3))./(WN1(rs1,2) + WN1(rs1,3));
WN_fs = (WN1(fs1,2) - WN1(fs1,3))./(WN1(fs1,2) + WN1(fs1,3));
WN_rsL = (WN1L(rs1,3) - WN1(rs1,3))./(WN1L(rs1,3) + WN1(rs1,3));
WN_fsL = (WN1L(fs1,3) - WN1(fs1,3))./(WN1L(fs1,3) + WN1(fs1,3));
WN_rs_plus = (WN1L(rs1,2) - WN1(rs1,3))./(WN1L(rs1,2) + WN1(rs1,3));
WN_fs_plus = (WN1L(fs1,2) - WN1(fs1,3))./(WN1L(fs1,2) + WN1(fs1,3));
% WN_exp_rs = [nanmean(WN_rs) nanmean(WN_rsL)];
% WN_exp_fs = [nanmean(WN_fs) nanmean(WN_fsL)];
WN_exp_rs = sum([WN_rs WN_rsL],2);
WN_exp_fs = sum([WN_fs WN_fsL],2);

[p,h,zstat] = ranksum(WN_rs_plus, WN_exp_rs);
[p,h,zstat] = ranksum(WN_fs_plus, WN_exp_fs);

SP_rs = (SP1(rs1,2) - SP1(rs1,3))./(SP1(rs1,2) + SP1(rs1,3));
SP_fs = (SP1(fs1,2) - SP1(fs1,3))./(SP1(fs1,2) + SP1(fs1,3));
SP_rsL = (SP1L(rs1,3) - SP1(rs1,3))./(SP1L(rs1,3) + SP1(rs1,3));
SP_fsL = (SP1L(fs1,3) - SP1(fs1,3))./(SP1L(fs1,3) + SP1(fs1,3));
SP_rs_plus = (SP1L(rs1,2) - SP1(rs1,3))./(SP1L(rs1,2) + SP1(rs1,3));
SP_fs_plus = (SP1L(fs1,2) - SP1(fs1,3))./(SP1L(fs1,2) + SP1(fs1,3));
%SP_exp_rs = [nanmean(SP_rs) nanmean(SP_rsL)];
%SP_exp_fs = [nanmean(SP_fs) nanmean(SP_fsL)];
SP_exp_rs = sum([SP_rs SP_rsL],2);
SP_exp_fs = sum([SP_fs SP_fsL],2);

figure; hold on
bar([1 2 3], [nanmean(WN_rs), nanmean(WN_rsL), nanmean(WN_rs_plus)],'BarWidth', .3)
bar([1.3 2.3 3.3], [nanmean(WN_fs), nanmean(WN_fsL), nanmean(WN_fs_plus)], 'BarWidth', .3)
bar([4 5 6], [nanmean(SP_rs), nanmean(SP_rsL), nanmean(SP_rs_plus)],'BarWidth', .3)
bar([4.3 5.3 6.3], [nanmean(SP_fs), nanmean(SP_fsL), nanmean(SP_fs_plus)], 'BarWidth', .3)
xticks(1:6);
xticklabels({'pupil - VIP', 'VIP - pupil', 'pupil + VIP', 'pupil - VIP', 'VIP - pupil', 'pupil + VIP'})
errorbar([1 2 3 1.3 2.3 3.3 4 5 6 4.3 5.3 6.3], [nanmean(WN_rs), nanmean(WN_rsL), nanmean(WN_rs_plus), nanmean(WN_fs), nanmean(WN_fsL), nanmean(WN_fs_plus),...
    nanmean(SP_rs), nanmean(SP_rsL), nanmean(SP_rs_plus), nanmean(SP_fs), nanmean(SP_fsL), nanmean(SP_fs_plus)], ...
    [sem(WN_rs), sem(WN_rsL), sem(WN_rs_plus), sem(WN_fs), sem(WN_fsL), sem(WN_fs_plus),...
    sem(SP_rs), sem(SP_rsL), sem(SP_rs_plus), sem(SP_fs), sem(SP_fsL), sem(SP_fs_plus)])
plot([3 6], [nanmean(WN_exp_rs)  nanmean(SP_exp_rs)], 'ko', 'MarkerSize', 10)
plot([3.3 6.3], [nanmean(WN_exp_fs)  nanmean(SP_exp_fs)], 'ko', 'MarkerSize', 10)
legend({'evoked rs', 'evoked fs', 'spont rs', 'spont fs'})
ylabel('Modulation Index')

[p,h,zstat] = ranksum(SP_rs_plus, SP_exp_rs);
[p,h,zstat] = ranksum(SP_fs_plus, SP_exp_fs);

figure; subplot(2,1,1);hold on
plot(WN_rs_plus, WN_exp_rs, 'k.');
plot(WN_fs_plus, WN_exp_fs, 'b.');
xlabel('VIP + pupil'); ylabel('Expected'); lsline
[rho_rs,p_rs]=corr(WN_rs_plus, WN_exp_rs, 'Type', 'Spearman', 'rows', 'pairwise');
[rho_fs,p_fs]=corr(WN_fs_plus, WN_exp_fs, 'Type', 'Spearman', 'rows', 'pairwise');
title(sprintf('Evoked, rs rho=%.4f,p=%.4f fs rho=%.4f,p=%.4f', rho_rs, p_rs, rho_fs, p_fs))
subplot(2,1,2);hold on
plot(SP_rs_plus, SP_exp_rs, 'k.');
plot(SP_fs_plus, SP_exp_fs, 'b.');
xlabel('VIP + pupil'); ylabel('Expected'); lsline
[rho_rs,p_rs]=corr(SP_rs_plus, SP_exp_rs, 'Type', 'Spearman', 'rows', 'pairwise');
[rho_fs,p_fs]=corr(SP_fs_plus, SP_exp_fs, 'Type', 'Spearman', 'rows', 'pairwise');
title(sprintf('Spont, rs rho=%.4f,p=%.4f fs rho=%.4f,p=%.4f', rho_rs, p_rs, rho_fs, p_fs))

figure; subplot(2,1,1);hold on
plot(WN_rs, WN_rsL, 'k.')
plot(WN_fs, WN_fsL, 'b.')
[rho_rs,p_rs]=corr([WN_rs], [WN_rsL], 'Type', 'Spearman', 'rows', 'pairwise');
[rho_fs,p_fs]=corr([WN_fs], [WN_fsL], 'Type', 'Spearman', 'rows', 'pairwise');
title(sprintf('Evoked, rs rho=%.4f,p=%.4f fs rho=%.4f,p=%.4f', rho_rs, p_rs, rho_fs, p_fs))
xlabel('pupil - VIP'); ylabel('VIP - pupil'); lsline
subplot(2,1,2);hold on
plot(SP_rs, SP_rsL, 'k.')
plot(SP_fs, SP_fsL, 'b.')
xlabel('pupil - VIP'); ylabel('VIP - pupil'); lsline
[rho_rs,p_rs]=corr([SP_rs], [SP_rsL], 'Type', 'Spearman', 'rows', 'pairwise');
[rho_fs,p_fs]=corr([SP_fs], [SP_fsL], 'Type', 'Spearman', 'rows', 'pairwise');
title(sprintf('Spont, rs rho=%.4f,p=%.4f fs rho=%.4f,p=%.4f', rho_rs, p_rs, rho_fs, p_fs))

m_exp = bootstrp(n_boots, @nanmean, WN_exp_rs);
m_plus = bootstrp(n_boots, @nanmean, WN_rs_plus);
[p,h,zstat] = ranksum(m_plus, m_exp);
Stats(5,:,:) = [p, zstat.zval];
m_exp = bootstrp(n_boots, @nanmean, WN_exp_fs);
m_plus = bootstrp(n_boots, @nanmean, WN_fs_plus);
[p,h,zstat] = ranksum(m_plus, m_exp);
Stats(6,:,:) = [p, zstat.zval];

m_exp = bootstrp(n_boots, @nanmean, SP_exp_rs);
m_plus = bootstrp(n_boots, @nanmean, SP_rs_plus);
[p,h,zstat] = ranksum(m_plus, m_exp);
Stats(7,:,:) = [p, zstat.zval];
m_exp = bootstrp(n_boots, @nanmean, SP_exp_fs);
m_plus = bootstrp(n_boots, @nanmean, SP_fs_plus);
[p,h,zstat] = ranksum(m_plus, m_exp);
Stats(8,:,:) = [p, zstat.zval];



% cells
figure; hold on
bar([1 2 3 4], [sum(zstats(:,1)>0 & evoked(:,1)) sum(zstats(:,2)>0 & evoked(:,2)) sum(zstats(:,3)>0 & evoked(:,3)) sum(zstats(:,4)>0 & evoked(:,4))], 'BarWidth', .3)
bar([1.3 2.3 3.3 4.3], [sum(zstats(:,1)<0 & evoked(:,1)) sum(zstats(:,2)<0 & evoked(:,2)) sum(zstats(:,3)<0 & evoked(:,3)) sum(zstats(:,4)<0 & evoked(:,4))], 'BarWidth', .3)
bar(5, sum(evoked(:,4)), 'BarWidth', .3)
bar(5.3, sum(evoked(:,4)==0), 'BarWidth', .3)
xticks([1:5])
xticklabels({'ON', 'Sust', 'OFF', 'Combined', 'All'})
legend({'activated', 'suppressed', 'responsive', 'not responsive'})
ylabel('number of cells')
%%
% bootstramp results

x = SP_rs_plus;
SP_rs_plus1 = bootstrp(100, @(x) nanmean(x), x);
x = SP_exp_rs;
SP_exp_rs1 = bootstrp(100, @(x) nanmean(x), x);
x = SP_fs_plus;
SP_fs_plus1 = bootstrp(100, @(x) nanmean(x), x);
x = SP_exp_fs;
SP_exp_fs1 = bootstrp(100, @(x) nanmean(x), x);
x = WN_rs_plus;
WN_rs_plus1 = bootstrp(100, @(x) nanmean(x), x);
x = WN_exp_rs;
WN_exp_rs1 = bootstrp(100, @(x) nanmean(x), x);
x = WN_fs_plus;
WN_fs_plus1 = bootstrp(100, @(x) nanmean(x), x);
x = SP_exp_fs;
WN_exp_fs1 = bootstrp(100, @(x) nanmean(x), x);

figure; hold on
plot(SP_rs_plus1, SP_exp_rs1, 'k.')
plot(SP_fs_plus1, SP_exp_fs1, 'b.')


figure; hold on
bar([1 2 3], [nanmean(WN_rs), nanmean(WN_rsL), nanmean(WN_rs_plus)],'BarWidth', .3)
bar([1.3 2.3 3.3], [nanmean(WN_fs), nanmean(WN_fsL), nanmean(WN_fs_plus)], 'BarWidth', .3)
bar([4 5 6], [nanmean(SP_rs), nanmean(SP_rsL), nanmean(SP_rs_plus)],'BarWidth', .3)
bar([4.3 5.3 6.3], [nanmean(SP_fs), nanmean(SP_fsL), nanmean(SP_fs_plus)], 'BarWidth', .3)
xticks(1:6);
xticklabels({'pupil - VIP', 'VIP - pupil', 'pupil + VIP', 'pupil - VIP', 'VIP - pupil', 'pupil + VIP'})
errorbar([1 2 3 1.3 2.3 3.3 4 5 6 4.3 5.3 6.3], [nanmean(WN_rs), nanmean(WN_rsL), nanmean(WN_rs_plus), nanmean(WN_fs), nanmean(WN_fsL), nanmean(WN_fs_plus1),...
    nanmean(SP_rs), nanmean(SP_rsL), nanmean(SP_rs_plus), nanmean(SP_fs), nanmean(SP_fsL), nanmean(SP_fs_plus1)], ...
    [sem(WN_rs), sem(WN_rsL), sem(WN_rs_plus1), sem(WN_fs), sem(WN_fsL), sem(WN_fs_plus1),...
    sem(SP_rs), sem(SP_rsL), sem(SP_rs_plus1), sem(SP_fs), sem(SP_fsL), sem(SP_fs_plus1)])
plot([3 6], [nanmean(WN_exp_rs1)  nanmean(SP_exp_rs1)], 'ko', 'MarkerSize', 10)
plot([3.3 6.3], [nanmean(WN_exp_fs1)  nanmean(SP_exp_fs1)], 'ko', 'MarkerSize', 10)
legend({'evoked rs', 'evoked fs', 'spont rs', 'spont fs'})
ylabel('Modulation Index')


