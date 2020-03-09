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
load ('C:\Users\lab\Resilio Sync\Paper1Figures\code\colorblind_colormap.mat'); % use color blind friendly colors

maxFRall =[];
SP = nan(length(data),2);
SPL = nan(length(data),2);
meanON = nan(length(data),2); % On response only (0 - 100 ms)
meanONL = nan(length(data),2);
meanWN = nan(length(data),2); % Full response (0 - 600 ms)
meanWNL = nan(length(data),2);
meanPreStim = nan(length(data),2);
meanPreStimL = nan(length(data),2);
meanSilentSound = nan(length(data),2);
depths = nan(length(data),1);
fs = zeros(length(data),1);
Rs = zeros(length(data),1);
WNdirsM =[]; WNcellsM = 0; SSdirsM =[]; SScellsM = 0;
recs = []; cells=[]; cell_number= nan(length(data),1);

evoked = zeros(length(data),4); zstats = zeros(length(data),4);
color1 = [0.7 0.7 0.7]; color2 = [0.2 0.2 0.2];
nreps_check_running = 6; %number of repetitions in each condition for comparison
CL = {[0 301], [300 401], [400 601], [600 2000]};

cdVIP; load('Silence_DistanceCorr_dirs.mat')

for cc =1:length(data)
    if data(cc).dir < 30
        if data(cc).dir~=0 && data(cc).dir~=0 % ecluding putlier recordings
        try
            meanSpikeCount = nanmean([data(cc).SpikeCountWN data1(cc).SpikeCountWN ]); %data(cc).SpikeCountSS data1(cc).SpikeCountSS
        catch
            meanSpikeCount = NaN;
        end
        % exclude cells with very low spikecount, they usually have very
        % large effects; count spikes to WN only
        if meanSpikeCount > 2  && CellsQualityStats.SNR(cc)>.5 && CellsQualityStats.uQ(cc)>10
            
            Spont = [data(cc).mSSon; data(cc).mSSoff]; % spont trials in all states
            PreStim = [data(cc).mNson; data(cc).mNsoff]; %pre stimulus fr
            
            ON = [data(cc).mNON_on; data(cc).mNON_off]; % data(cc).mNON_on; data(cc).mNON_off]; %ON = nanmean(ON);
            [rs, h, stats] = ranksum( ON(:), Spont(:));
            evoked(cc,1) = h;
            zstats(cc,1) = stats.zval;
            %zstats(cc,1) = stats.tstat;
            
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
                meanWN(cc,1)= nanmean(nanmean(data(cc).mWNon));% nanmean(nanmean(data(cc).mNSustained_on)) nanmean(nanmean(data(cc).mNOFF_on))]);
                meanWN(cc,2)= nanmean(nanmean(data(cc).mWNoff));% nanmean(nanmean(data(cc).mNSustained_off)) nanmean(nanmean(data(cc).mNOFF_off))]);
                meanON(cc,1) = nanmean(nanmean(data(cc).mNON_on));%-nanmean(nanmean(data(cc).mNson));
                meanON(cc,2) = nanmean(nanmean(data(cc).mNON_off));%-nanmean(nanmean(data(cc).mNsoff));
                meanPreStim(cc,1) = nanmean(nanmean(data(cc).mNson));
                meanPreStim(cc,2) = nanmean(nanmean(data(cc).mNsoff));
                WNdirsM = [WNdirsM data(cc).dir]; % collect directory (rec#) for all recordings included in the analysis to make sure N is ok.
                WNcellsM = WNcellsM + 1; % count number of cells included in the analysis
            end
            
            % WN laser on trials
            if data1(cc).nrepsWNMon > nreps_check_running && data1(cc).nrepsWNMoff > nreps_check_running % && data(cc).nrepsWNMon > nreps_check_running && data(cc).nrepsWNMoff > nreps_check_running
                meanWNL(cc,1)= nanmean(nanmean(data1(cc).mWNon));%nanmean([nanmean(nanmean(data1(cc).mNON_on)) nanmean(nanmean(data1(cc).mNSustained_on)) nanmean(nanmean(data1(cc).mNOFF_on))]);
                meanWNL(cc,2)= nanmean(nanmean(data1(cc).mWNoff));%nanmean([nanmean(nanmean(data1(cc).mNON_off)) nanmean(nanmean(data1(cc).mNSustained_off)) nanmean(nanmean(data1(cc).mNOFF_off))]);
                meanONL(cc,1) =nanmean([nanmean(nanmean(data1(cc).mNON_on))]);%-nanmean(nanmean(data1(cc).mNson));
                meanONL(cc,2) =nanmean([nanmean(nanmean(data1(cc).mNON_off))]);%-nanmean(nanmean(data1(cc).mNsoff));
                meanPreStimL(cc,1) = nanmean([nanmean(nanmean(data1(cc).mNson))]);
                meanPreStimL(cc,2) = nanmean([nanmean(nanmean(data1(cc).mNsoff))]);
                cell_number(cc) = cc;
            end
            
            % silent sound laser off
            if data(cc).nrepsSSMon > nreps_check_running && data(cc).nrepsSSMoff > nreps_check_running %&& data1(cc).nrepsSSMon > nreps_check_running && data1(cc).nrepsSSMoff > nreps_check_running
                SP(cc,1) = nanmean(nanmean(data(cc).mSSon)); SP(cc,2) = nanmean(nanmean(data(cc).mSSoff));
                meanSilentSound(cc,2) = nanmean([nanmean(nanmean(data(cc).mSSoff))]);
                meanSilentSound(cc,1) = nanmean([nanmean(nanmean(data(cc).mSSon))]);
                WNdirsM = [WNdirsM data(cc).dir]; % collect directory (rec#) for all recordings included in the analysis to make sure N is ok.
                WNcellsM = WNcellsM + 1; % count number of cells included in the analysis
            end
            
            % silent sound laser on
            if data1(cc).nrepsSSMon > nreps_check_running && data1(cc).nrepsSSMoff > nreps_check_running % && data(cc).nrepsSSMon > nreps_check_running && data(cc).nrepsSSMoff > nreps_check_running
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
    end
    recs = [recs data(cc).dir];
    cells = [cells data(cc).cell];
    
end

allDirsM = unique(WNdirsM);

% which cells to include? 1- On response, 2 - Sustained, 3-  Off response,
% 4 - full 600 ms. zstats > 0 activated, <0 suppressed
% 
evoked1 = logical(evoked(:,1) & zstats(:,1)>0); 
cd(variables_dir)
save('evoked_indx_epistatic.mat', 'evoked1')

%evoked1 = logical(ones(length(evoked),1));

responses =  (meanON - SP)./(meanON+SP);
figure; hist(responses)
legend('running', 'sitting')
[p, h, stats]= signrank(responses(:,1), responses(:,2));
title_string = sprintf( 'Sound Modulation Index all cells On responses z = %.2f, p = %d', stats.zval, p);
title(title_string)
xlabel('Sound MI')
ylabel('Number of cells')

responses =  (meanWN - SP)./(meanWN+SP);
figure; hist(responses)
legend('running', 'sitting')
[p, h, stats]= signrank(responses(:,1), responses(:,2));
title_string = sprintf( 'Sound Modulation Index all cells Full responses z = %.2f, p = %d', stats.zval, p);
title(title_string)
xlabel('Sound MI')
ylabel('Number of cells')


WN1 = meanON(evoked1,:);
WN1L = meanONL(evoked1,:);
%WN1 = meanWN(evoked1,:);
%WN1L = meanWNL(evoked1,:);
SP1 = SP(evoked1,:);
SP1L = SPL(evoked1,:);
meanPreStim1 = meanPreStim(evoked1,:);
meanPreStim1L = meanPreStimL(evoked1,:);
meanSilentSound1 = meanSilentSound(evoked1,:);
depths1 = depths(evoked1);
fs1 = fs(evoked1);
rs1 = Rs(evoked1);
rs1 = logical(rs1); fs1 = logical(fs1);
cell_number1 = cell_number(evoked1);

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
indx = find(WN1(:,2)< meanPreStim1(:,2));
indx1 = find(WN1(:,1) < meanPreStim1(:,1));
WN2 = WN1 - meanPreStim1; WN2L = WN1L - meanPreStim1L;
WN2(indx,:) = NaN; WN2(indx1,:) = NaN;
%WN2 = WN1 - meanSilentSound1; WN2L = WN1L - meanPreStim1L;
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
cells1 = cells(evoked1); recs1 = recs(evoked1); 


modulation_indx1 = (WN1(:,1) -WN1(:,2))./(WN1(:,1) + WN1(:,2)); %running, laser off trials
modulation_indx1L = (WN1L(:,1) -WN1L(:,2))./(WN1L(:,1) + WN1L(:,2)); %running, laser on trials
modulation_indx1Laser = (WN1L(:,2) - WN1(:,2))./(WN1L(:,2) + WN1(:,2)); % effect of layer in this layer across states

modulation_indx1_sp = (SP1(:,1) -SP1(:,2))./(SP1(:,1) + SP1(:,2)); %spont running effect laser off
modulation_indx1_spL = (SP1L(:,1) -SP1L(:,2))./(SP1L(:,1) + SP1L(:,2)); %spont running effect laser on
modulation_indx1_spLaser = (SP1L(:,2) - SP1(:,2))./(SP1L(:,2) + SP1(:,2)); %effects of laser in sitting

MI_sp_plus = (SP1L(:,1) - SP1(:,2))./(SP1L(:,1) + SP1(:,2));
MI_evoked_plus = (WN1L(:,1) - WN1(:,2))./(WN1L(:,1) + WN1(:,2));


for cl = 1:length(CL)
    layer = CL{cl}; % layer depth limits
    indx = find(depths1 >layer(1) & depths1 <layer(2)); % indices of cells within this layer
    recs2 = recs1(indx); % recordings in this layer
    cells2 = cells1(indx); % cells in this layer
    fs2 = fs1(indx); fs2 = logical(fs2); % fast spiking cells in this layer
    rs2 = rs1(indx); rs2 = logical(rs2); % regular spiking cells in this layer
    
    % % evoked (WN respone)
    % means
    meanMI1(cl) = nanmean(modulation_indx1(indx)); % mean, laser off
    meanMI1L(cl) = nanmean(modulation_indx1L(indx)); % mean, laser on
    meanMI1_Laser(cl) = nanmean(modulation_indx1Laser(indx)); % mean laser effect
    meanMI2(cl) = nanmean(modulation_indx2(indx));
    
    meanMI1_rs(cl) = nanmean(modulation_indx1(indx(rs2)));
    meanMI1_fs(cl) = nanmean(modulation_indx1(indx(fs2)));
    meanMI1_rsL(cl) = nanmean(modulation_indx1L(indx(rs2)));
    meanMI1_fsL(cl) = nanmean(modulation_indx1L(indx(fs2)));
    
    semMI1(cl) = sem(modulation_indx1(indx));
    semMI1L(cl) = sem(modulation_indx1L(indx));
    semMI1_Laser(cl) = sem(modulation_indx1Laser(indx));
    semMI2(cl) = sem(modulation_indx2(indx));
    
    semMI1_rs(cl) = sem(modulation_indx1(indx(rs2)));
    semMI1_fs(cl) = sem(modulation_indx1(indx(fs2)));
    semMI1_rsL(cl) = sem(modulation_indx1L(indx(rs2)));
    semMI1_fsL(cl) = sem(modulation_indx1L(indx(fs2)));
    
    % collect modulation indices across layers for further stats [MI cl]
    MODULATION_INDX1 = [MODULATION_INDX1; modulation_indx1(indx) ones(length(indx),1)*cl];
    MODULATION_INDX1L = [MODULATION_INDX1L; modulation_indx1L(indx) ones(length(indx),1)*cl];
    MODULATION_INDX1Laser = [MODULATION_INDX1Laser; modulation_indx1Laser(indx) ones(length(indx),1)*cl];
    
    n_layer_evoked(cl) = length(find(~isnan(modulation_indx1(indx))==1));
    n_layer_evokedL(cl) = length(find(~isnan(modulation_indx1L(indx))==1));
    n_layer_evokedLaser(cl) = length(find(~isnan(modulation_indx1Laser(indx))==1));
    
    rs_cl(cl) = sum(~isnan(modulation_indx1(indx(rs2)))); % number of regular spiking cells in this layer
    fs_cl(cl) = sum(~isnan(modulation_indx1(indx(fs2)))); % number of fast spiking cells in this layer
    
    % addative effect
    meanMI_evoked_plus(cl) = nanmean(MI_evoked_plus(indx));
    semMI_evoked_plus(cl) = sem(MI_evoked_plus(indx));
    n_layer_evoked_plus(cl) = length(find(~isnan(MI_evoked_plus(indx))==1)); %number of cell
    
    % evoked - spont
    meanMI2(cl) = nanmean(modulation_indx2(indx)); % mean, laser off
    meanMI2L(cl) = nanmean(modulation_indx2L(indx)); % mean, laser on
    
    meanMI2(cl) = nanmean(modulation_indx2(indx));
    meanMI2L(cl) = nanmean(modulation_indx2L(indx));
    
    semMI2(cl) = sem(modulation_indx2(indx));
    semMI2L(cl) = sem(modulation_indx2L(indx));
    
    % % spont activity
    meanMI1_sp(cl) = nanmean(modulation_indx1_sp(indx));
    meanMI1_spL(cl) = nanmean(modulation_indx1_spL(indx));
    meanMI1_spLaser(cl) = nanmean(modulation_indx1_spLaser(indx));
    meanMI1_sp_rs(cl) = nanmean(modulation_indx1_sp(indx(rs2)));
    meanMI1_spL_rs(cl) = nanmean(modulation_indx1_spL(indx(rs2)));
    meanMI1_spLaser_rs(cl) = nanmean(modulation_indx1_spLaser(indx(rs2)));
    meanMI1_sp_fs(cl) = nanmean(modulation_indx1_sp(indx(fs2)));
    meanMI1_spL_fs(cl) = nanmean(modulation_indx1_spL(indx(fs2)));
    meanMI1_spLaser_fs(cl) = nanmean(modulation_indx1_spLaser(indx(fs2)));
    
    semMI1_sp(cl) = sem(modulation_indx1_sp(indx));
    semMI1_spL(cl) = sem(modulation_indx1_spL(indx));
    semMI1_spLaser(cl) = sem(modulation_indx1_spLaser(indx));
    semMI1_sp_rs(cl) = sem(modulation_indx1_sp(indx(rs2)));
    semMI1_spL_rs(cl) = sem(modulation_indx1_spL(indx(rs2)));
    semMI1_spLaser_rs(cl) = sem(modulation_indx1_spLaser(indx(rs2)));
    semMI1_sp_fs(cl) = sem(modulation_indx1_sp(indx(fs2)));
    semMI1_spL_fs(cl) = sem(modulation_indx1_spL(indx(fs2)));
    semMI1_spLaser_fs(cl) = sem(modulation_indx1_spLaser(indx(fs2)));
    
    n_layer_sp(cl) = length(find(~isnan(modulation_indx1_sp(indx))==1)); %number of cells in each layer
    n_layer_spL(cl) = length(find(~isnan(modulation_indx1_spL(indx))==1));
    n_layer_spLaser(cl) = length(find(~isnan(modulation_indx1_spLaser(indx))==1));
    
    rs_sp_cl(cl) = sum(~isnan(modulation_indx1_sp(indx(rs2))));
    fs_sp_cl(cl) = sum(~isnan(modulation_indx1_sp(indx(fs2))));
    
    % additive effect
    meanMI_sp_plus(cl) = nanmean(MI_sp_plus(indx));
    semMI_sp_plus(cl) = sem(MI_sp_plus(indx));
    n_layer_sp_plus(cl) = length(find(~isnan(MI_sp_plus(indx))==1));
    
    MODULATION_sp_INDX1 = [MODULATION_sp_INDX1; modulation_indx1_sp(indx) ones(length(indx),1)*cl];
    MODULATION_sp_INDX1L = [MODULATION_sp_INDX1L; modulation_indx1_spL(indx) ones(length(indx),1)*cl];
    MODULATION_sp_INDX1Laser = [MODULATION_sp_INDX1Laser; modulation_indx1_spLaser(indx) ones(length(indx),1)*cl];
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

figure; hold on
plot(modulation_indx1(rs1), modulation_indx1_sp(rs1), 'ko')
plot(modulation_indx1(fs1), modulation_indx1_sp(fs1), 'go')
lsline; xlabel('response MI'); ylabel('spont MI')
plot([0 0], [-1 1], 'k--')
plot([-1 1], [0 0], 'k--')
[r1, p1] = corr(modulation_indx1(rs1), modulation_indx1_sp(rs1), 'Type','Spearman','Rows', 'complete');
[r2, p2] = corr(modulation_indx1(fs1), modulation_indx1_sp(fs1), 'Type','Spearman','Rows', 'complete');
title_str = sprintf('rs rho = %.2f, p=%.4f; fs rho = %.2f, p = %.4f', r1, p1, r2, p2);
title(title_str); pbaspect([1 1 1]);

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
 set(gcf, 'PaperPositionMode', 'auto');

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
errorbar([1:4], meanMI1_rs, semMI1_rs, 'ko-');
errorbar([1.1:4.1], meanMI1_fs, semMI1_fs, 'k*-');
plot([0 5], [0 0], 'k--')
legend({ 'run rs', 'run fs'})
xlabel('Cortical Layer')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Modulation Index'); title('WN responses')
hold on; subplot(2,1,2); hold on
errorbar([1:4], meanMI1_sp_rs, semMI1_sp_rs, 'ko-');
errorbar([1.1:4.1], meanMI1_sp_fs, semMI1_sp_fs, 'k*-');
plot([0 5], [0 0], 'k--')
legend({ 'run rs', 'run fs'})
xlabel('Cortical Layer')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Modulation Index'); title('Spont Activity')

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

%%
% sound modulation index
MI_sound_run = (WN1(:,1) - SP1(:,1))./ (WN1(:,1) +SP1(:,1));
MI_sound_sit = (WN1(:,2) - SP1(:,2))./ (WN1(:,2) +SP1(:,2));
MI_sound_run_laser = (WN1L(:,1) - SP1L(:,1))./ (WN1L(:,1) +SP1L(:,1));
MI_sound_sit_laser = (WN1L(:,2) - SP1L(:,2))./ (WN1L(:,2) +SP1L(:,2));
MI_sound_predicted = MI_sound_sit_laser + MI_sound_run;

% plot MI
figure; hold on
plot(MI_sound_sit(rs1) ,MI_sound_run(rs1), 'ko')
plot(MI_sound_sit(fs1) ,MI_sound_run(fs1), 'go')
ylabel('sound MI running'); xlabel('sound MI sitting')
[r, p] = corr(MI_sound_run,MI_sound_sit, 'Type','Spearman','Rows', 'complete')
plot([-1 1], [-1 1], 'r-')
plot([-1 1], [0 0], 'k--')
plot([0 0], [-1 1], 'k--')
pbaspect([1 1 1]);  set(gcf, 'PaperPositionMode', 'auto');

% plot difference in MI
run_diff =  MI_sound_run - MI_sound_sit;
laser_diff = MI_sound_sit_laser - MI_sound_sit;
predicted_diff = run_diff + laser_diff;
actual_diff =  MI_sound_run_laser - MI_sound_sit;

figure; hold
plot(run_diff(rs1), laser_diff(rs1), 'o', 'Color', [.8 .80 .80], 'MarkerSize', 8)
plot(run_diff(fs1), laser_diff(fs1), 'o', 'Color', [0 .255 0], 'MarkerSize', 8)
xlabel('Running Effect (MI diff)')
ylabel('Laser Effect (MI diff)')
plot([-2 2], [0 0], 'k--')
plot([0 0], [-2 2], 'k--')
plot([-2 2], [-2 2], 'k--')
[r, p] = corr(run_diff,laser_diff, 'Type','Spearman','Rows', 'complete')
title_string = sprintf( 'rho = %.4f, p = %.4f', r,p);
title(title_string); legend('Regular spiking', 'Narrow spiking')
pbaspect([1 1 1]);  set(gcf, 'PaperPositionMode', 'auto');
ylim([-2 2])
xlim([-2 2])

figure; hold on
plot(actual_diff(rs1), predicted_diff(rs1), 'ko')
plot(actual_diff(fs1), predicted_diff(fs1), 'go')
xlabel('Combined Effect (running laser on trials)')
ylabel('Predicted combined Effect (running effect + laser effect)')
plot([-4 2], [0 0], 'k--')
plot([0 0], [-4 2], 'k--')
plot([-4 4], [-4 4], 'k--')
[r, p] = corr(actual_diff,predicted_diff, 'Type','Spearman','Rows', 'complete')
title_string = sprintf( 'rho = %.4f, p = %.4f', r,p);
ylim([-4 2])
xlim([-4 2])
title(title_string); legend('Regular spiking', 'Narrow spiking')
pbaspect([1 1 1]);  set(gcf, 'PaperPositionMode', 'auto');

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
title(' all On responses, laser on')

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
    n_layer(cl) = sum(~isnan(MI_sound_run(indx)));
    
    meanMI_sound_runL(cl) = nanmean(MI_sound_run_laser(indx));
    meanMI_sound_sitL(cl) = nanmean(MI_sound_sit_laser(indx));
    semMI_sound_runL(cl) = sem(MI_sound_run_laser(indx));
    semMI_sound_sitL(cl) = sem(MI_sound_sit_laser(indx));
    n_layerL(cl) = sum(~isnan(MI_sound_run_laser(indx)));
end

figure; hold on
errorbar([1:4], meanMI_sound_sit, semMI_sound_sit, 'ko-');
errorbar([1:4], meanMI_sound_run, semMI_sound_run, 'ko--');
xlabel('Cortical Layer')
plot([0 5], [0 0], 'k--')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Modulation Index (mean/SEM)- sound effect');
title_string = sprintf( 'On response MI, n = %d, %d, %d, %d',  n_layer);
title(title_string)
legend( 'sitting', 'running')

figure; hold on
bar([1 2], [nanmean(MI_sound_sit) nanmean(MI_sound_run)])
errorbar([1 2], [nanmean(MI_sound_sit) nanmean(MI_sound_run)], [sem(MI_sound_sit) sem(MI_sound_run)])
xticks([1 2])
xlim([0 3])
xticklabels({'sitting', 'running'})
ylabel('Mean (SEM) Sound Modulation Index')

figure; hold on
errorbar([1:4], meanMI_sound_sitL, semMI_sound_sitL, 'co-');
errorbar([1:4], meanMI_sound_runL, semMI_sound_runL, 'co--');
xlabel('Cortical Layer')
plot([0 5], [0 0], 'k--')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Modulation Index (mean/SEM)- sound effect');
title_string = sprintf( 'On response MI, n = %d, %d, %d, %d',  n_layer);
title(title_string)
legend('sitting','running')

% distribution scatter plot
figure; hold on
plot(WN1(rs1,2), WN1(rs1,1), 'ko')
plot(WN1(fs1,2), WN1(fs1,1), 'go')
maxFR = max(max(WN1));
plot([0 maxFR], [0 maxFR], 'r-')
plot(nanmean(WN1(:,2)), nanmean(WN1(:,1)), 'ro')
[p,h,stats] = signrank(WN1(:,1),  WN1(:,2));
title_str = sprintf('Evoked Activity, p = %d', p);
title(title_str)
xlabel('FR sitting (Hz)'); ylabel('FR running (Hz)')
pbaspect([1 1 1]); set(gcf, 'PaperPositionMode', 'auto');

figure; hold on;
plot(SP1(rs1,2), SP1(rs1,1), 'ko')
plot(SP1(fs1,2), SP1(fs1,1), 'go')
maxFR = max(max([SP1(:,2), SP1(:,1)]));
plot([0 maxFR], [0 maxFR], 'r-')
plot(nanmean(SP1(:,2)), nanmean(SP1(:,1)), 'ro')
[p,h,stats] = signrank(SP1(:,1),  SP1(:,2));
title_str = sprintf('Spont Activity, p = %d', p);
title(title_str)
xlabel('FR running (Hz)'); ylabel('FR running (Hz)')
pbaspect([1 1 1]); set(gcf, 'PaperPositionMode', 'auto');

fs1 = logical(fs1); rs1 = logical(rs1);
MI1_rs = (WN1(rs1,1) - WN1(rs1,2))./ (WN1(rs1,1) +WN1(rs1,2));
MI1_fs = (WN1(fs1,1) - WN1(fs1,2))./ (WN1(fs1,1) +WN1(fs1,2));

SP1_rs = (SP1(rs1,1) - SP1(rs1,2))./ (SP1(rs1,1) +SP1(rs1,2));
SP1_fs = (SP1(fs1,1) - SP1(fs1,2))./ (SP1(fs1,1) +SP1(fs1,2));

figure(106); hold on
bar([0.8 1.8], [nanmean(SP1_rs) nanmean(MI1_rs)], 'BarWidth', .1)
bar([1  2 ], [ nanmean(SP1_fs) nanmean(MI1_fs) ], 'BarWidth', .1)

errorbar( [0.8 1 1.8 2],[nanmean(SP1_rs) nanmean(SP1_fs) nanmean(MI1_rs)  nanmean(MI1_fs)], ...
    [sem(SP1_rs) sem(SP1_fs) sem(MI1_rs) sem(MI1_fs)])
legend({'RS', 'FS'})
xticks([1:2])
xticklabels({'Spont run + laser off', 'Evoked run + laser off'})
ylabel('Modulation Index')


%% laser modulation
figure; hold on
errorbar([1.1:4.1], meanMI_sound_sit, semMI_sound_sit, 'ko-');
errorbar([1:4], meanMI_sound_sitL, semMI_sound_sitL, 'co-');
xlabel('Cortical Layer')
plot([0 5], [0 0], 'k--')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Modulation Index (mean/SEM)- sound effect');
legend('laser off', 'laser on')
fprintf('laer effect')
[p, h, stats] = signrank(MI_sound_sit_laser, MI_sound_sit);
title_string = sprintf( 'On response MI, n = %d, %d, %d, %d, laser effect, p = %d',  n_layer, p);
title(title_string)
set(gcf, 'PaperPositionMode', 'auto');

figure; hold on
plot(MI_sound_sit(rs1), MI_sound_sit_laser(rs1), 'ko')
plot(MI_sound_sit(fs1), MI_sound_sit_laser(fs1), 'go')
plot([-1 1], [-1 1], 'r--')
plot([-1 1], [0 0], 'k--')
plot( [0 0], [-1 1],'k--')
xlabel(' Sound MI laser off')
ylabel('Sound MI laser on')
pbaspect([1 1 1]); set(gcf, 'PaperPositionMode', 'auto');

figure; hold on;
plot(WN1(rs1,2), WN1L(rs1,2), 'ko')
plot(WN1(fs1,2), WN1L(fs1,2), 'go')
maxFR = max(max([WN1(:,2), WN1L(:,2)]));
plot(nanmean(WN1(:,2)), nanmean(WN1L(:,2)), 'ro')
plot([0 maxFR], [0 maxFR], 'r-')
title(title_str)
xlabel('FR laser off'); ylabel('FR laser on')
pbaspect([1 1 1]); set(gcf, 'PaperPositionMode', 'auto');

figure; hold on;
plot(SP1(rs1,2), SP1L(rs1,2), 'ko')
plot(SP1(fs1,2), SP1L(fs1,2), 'go')
maxFR = max(max([SP1(:,2), SP1L(:,2)]));
plot(nanmean(SP1(:,2)), nanmean(SP1L(:,2)), 'ro')
plot([0 maxFR], [0 maxFR], 'r-')
[p,h,stats] = signrank(SP1L(:,2),  SP1(:,2));
title_str = sprintf('Spont Activity, p = %d', p);
title(title_str)
xlabel('FR laser off'); ylabel('FR laser on')
pbaspect([1 1 1]); set(gcf, 'PaperPositionMode', 'auto');

figure; hold on;
plot(MI_sound_sit(rs1), MI_sound_sit_laser(rs1), 'ko')
plot(MI_sound_sit(fs1), MI_sound_sit_laser(fs1), 'go')
maxFR = max(max([MI_sound_sit, MI_sound_sit_laser]));
plot([-1 1], [-1 1], 'r-')
plot([-1 1], [0 0], 'k--')
plot([0 0], [-1 1], 'k--')
[p,h,stats] = signrank(MI_sound_sit_laser,MI_sound_sit)
title_str = sprintf('Evoked Activity, p = %d', p);
title(title_str)
xlabel('sound modulation index - sit laser off');
ylabel('sound modulation index - sit laser on')
legend('Regular spiking', 'Narrow spiking')
pbaspect([1 1 1]); set(gcf, 'PaperPositionMode', 'auto');

meanMI_sound_laser = nanmean(MI_sound_sit_laser);
semMI_sound_laser = nansem(MI_sound_sit_laser);

meanMI_sound = nanmean(MI_sound_sit);
semMI_sound = nansem(MI_sound_sit);

figure;hold on;
bar([1 2], [meanMI_sound_laser meanMI_sound])
errorbar(1, meanMI_sound_laser, semMI_sound_laser , 'k-' );
errorbar(2, meanMI_sound, semMI_sound , 'k-' )
xlim([0 3]); ylabel('Modulation Index - sound')
xticks([1 2]); xticklabels({'laser on', 'laser off'})
title('effect off sound in laser on and off conditions')

%% Linear interaction test with running MI and VIP MI
MI1_rsL = (WN1L(rs1,2) -WN1(rs1,2))./(WN1L(rs1,2) + WN1(rs1,2)); % laser rs
MI1_fsL = (WN1L(fs1,2) -WN1(fs1,2))./(WN1L(fs1,2) + WN1(fs1,2)); % laser fs

MI1_sp_rsL = (SP1L(rs1,2) -SP1(rs1,2))./(SP1L(rs1,2) + SP1(rs1,2)); % laser
MI1_sp_fsL = (SP1L(fs1,2) -SP1(fs1,2))./(SP1L(fs1,2) + SP1(fs1,2)); % laser

% VIP effect in ech state
figure; hold on
bar([.8  1.8 ], [nanmean(MI1_sp_rsL) nanmean(MI1_rsL)], 'BarWidth', .1)
bar([.9  1.9], [nanmean(MI1_sp_fsL) nanmean(MI1_fsL)], 'BarWidth', .1)
errorbar([.8 .9 1.8 1.9], [nanmean(MI1_sp_rsL) nanmean(MI1_sp_fsL) nanmean(MI1_rsL) nanmean(MI1_fsL)],...
    [sem(MI1_sp_rsL) sem(MI1_sp_fsL) sem(MI1_rsL) sem(MI1_fsL)])
xticks([1:2])
xticklabels({'spont', 'evoked'})
ylabel('Modulation Index'); title('VIP activation')
legend({'RS', 'FS'})

% running
MI1 = (WN1(:,1) - WN1(:,2))./(WN1(:,1) + WN1(:,2));
MI1_sp = (SP1(:,1) - SP1(:,2))./(SP1(:,1) + SP1(:,2));

% VIP
MI1_Laser = (WN1L(:,2) - WN1(:,2))./(WN1L(:,2) + WN1(:,2)); % laser all cells
MI1_sp_Laser = (SP1L(:,2) - SP1(:,2))./(SP1L(:,2) + SP1(:,2)); % laser all cells

% combined effect of VIP and running
MI1_runLaser = (WN1L(:,1) - WN1(:,2))./(WN1L(:,1) + WN1(:,2));
MI1_sp_runLaser = (SP1L(:,1) - SP1(:,2))./(SP1L(:,1) + SP1(:,2));

% predicted effect
MI1_predicted = MI1 + MI1_Laser;
MI1_sp_predicted = MI1_sp + MI1_sp_Laser;

figure; subplot(2,1,1); hold on
plot(MI1(rs1), MI1_Laser(rs1), 'ko')
plot(MI1(fs1), MI1_Laser(fs1), 'go')
plot([-1 1], [-1 1], 'r--')
xlabel('running MI'); ylabel('VIP MI')
[r, p] = corr(MI1,MI1_Laser, 'Type','Spearman','Rows', 'complete');
title_string = sprintf( 'Evoked, rho = %.4f, p = %.4f', r, p);
title(title_string)
subplot(2,1,2); hold on
plot(MI1_sp(rs1), MI1_sp_Laser(rs1), 'ko')
plot(MI1_sp(fs1), MI1_sp_Laser(fs1), 'go')
plot([-1 1], [-1 1], 'r--')
xlabel('running MI'); ylabel('VIP MI')
[r, p] = corr(MI1_sp,MI1_sp_Laser, 'Type','Spearman','Rows', 'complete');
title_string = sprintf( 'Spont, rho = %.4f, p = %.4f', r, p);
title(title_string)
set(gcf, 'PaperPositionMode', 'auto');

% linearity analysis
figure; subplot(2,1,1); hold on
plot(MI1_predicted(rs1), MI1_runLaser(rs1), 'ko')
plot(MI1_predicted(fs1), MI1_runLaser(fs1), 'go')
plot([0 0], [-2 2], 'k--'); plot( [-2 2], [0 0], 'k--')
lsline
xlabel('predicted MI'); ylabel('VIP + running MI')
[r, p] = corr(MI1_predicted, MI1_runLaser, 'Type','Spearman','Rows', 'complete');
title_string = sprintf( 'Evoked, rho = %.4f, p = %d', r, p);
title(title_string)
subplot(2,1,2); hold on
plot(MI1_sp_predicted(rs1), MI1_sp_runLaser(rs1), 'ko')
plot(MI1_sp_predicted(fs1), MI1_sp_runLaser(fs1), 'go')
plot([0 0], [-2 2], 'k--'); plot( [-2 2], [0 0], 'k--')
lsline
xlabel('predicted MI'); ylabel('VIP + running MI')
[r, p] = corr(MI1_sp_predicted,MI1_sp_runLaser, 'Type','Spearman','Rows', 'complete');
title_string = sprintf( 'Spont, rho = %.4f, p = %d', r, p);
title(title_string)
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Position',  [260 124 902 864])

%% plot examples, FR data
% for cc = 1:length(cell_number1)
%     c =  cell_number1(cc);
%     if ~isnan(c)
%         figure;  hold on
%         sitting_trace = smooth(nanmean(data(c).mNON_off),3);
%         running_trace = smooth(nanmean(data(c).mNON_on),3);
%         plot(sitting_trace, 'k')
%         plot(running_trace, 'r--')
%         laser_on_sitting_trace = smooth(nanmean(data1(c).mNON_off),3);
%         plot( laser_on_sitting_trace, 'c')
%         predicted_trace = nanmean([running_trace laser_on_sitting_trace],2);
%         laser_on_running_trace = smooth(nanmean(data1(c).mNON_on),3);
%         plot(predicted_trace, 'g')
%         plot(laser_on_running_trace, 'c--')
%         legend('sitting', 'running','laser on','predicted', 'runnning +laser on')
%         title_string = sprintf
%     end
% end
% plot examples, spiking data data



%%  STOP HERE, THE REST IS OLD


