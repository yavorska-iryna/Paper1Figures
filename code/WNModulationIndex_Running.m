 % compute sound modulation index WN responses, 
 % effect of running only
 
% variables gethered with getWNresponses.m

% three states: sit + small pupil, sit + large pupil, running  + large


clear; close all; dbstop if error
variables_dir = 'C:\Users\lab\Documents\GitHub\Paper1Figures\code\variables';

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
meanSilentSound = nan(length(data),2);
depths = nan(length(data),1);
fs = zeros(length(data),1);
Rs = zeros(length(data),1);
WNdirsM =[]; WNcellsM = 0; SSdirsM =[]; SScellsM = 0;
recs = []; cells=[]; mouse_ID = []; mouse = [];

evoked = zeros(length(data),4); zstats = zeros(length(data),4);
color1 = [0.7 0.7 0.7]; color2 = [0.2 0.2 0.2];
nreps_check_running = 6; %number of repetitions in each condition for comparison
CL = { [129 380], [379 525], [524 806], [805 2000]};  id = 0;

cdVIP; load('Silence_DistanceCorr_dirs.mat')
cdPV;  load('allPINPdirs.mat') % load all 45 dirs
PINPdirs = PINPDIRS;

% collect waveform, SNR, and uQ from cells. Must match length of data.
try
    cd(variables_dir)
    load('CellsInfo.mat')
catch
    WFs = []; SNRs = []; uQs = [];
    for d = 1:length(PINPdirs)
        if ~isempty (PINPdirs{d})
            [WF, SNR, uQ] = getCellsInfo(PINPdirs{d});
            WFs = [WFs; WF];
            SNRs = [SNRs; SNR];
            uQs = [uQs; uQ];
        end
    end
    cd(variables_dir)
    save('CellsInfo.mat', 'WFs','SNRs','uQs')
end

for cc =1:length(data)
    if data(cc).dir < 47
        if data(cc).dir~=33 && data(cc).dir~=0 % excluding putlier recordings
        try
            meanSpikeCount = nanmean([data(cc).SpikeCountWN data1(cc).SpikeCountWN ]); 
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
        d = data(cc).dir;
        try
            runsit_reps(d,:,:,:) = [data(cc).nrepsWNMon data(cc).nrepsWNMoff data(cc).nrepsSSMon data(cc).nrepsSSMoff];
        catch
            runsit_reps(d,:,:,:) = [NaN NaN NaN NaN];
        end
        end
    end
    
    recs = [recs data(cc).dir];
    cells = [cells data(cc).cell];
    try
        cd(PINPdirs{data(cc).dir})
        load('notebook.mat')
        mouse_ID(cc) = str2num(nb.mouseID);
    catch
        mouse_ID(cc) = NaN;
    end
end
%% ANALYSIS %%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculate fraction of running trials
movement_probability = nansum(runsit_reps(:,1:2:end),2)./nansum(runsit_reps,2);
grand_movement_probability = nanmean(movement_probability);
[CI,BOOTSTAT] = bootci(100,{@nanmean, movement_probability}, 'type', 'per');
allDirsM = unique(WNdirsM);

% which cells to include? 1- On response, 2 - Sustained, 3-  Off response,
% 4 - full 600 ms. zstats > 0 activated, <0 suppressed

evoked1 = logical(evoked(:,1) & zstats(:,1)>0); 

responses_sit =  meanON(:,2) - SP(:,2);
responses_run =  meanON(:,1) - SP(:,1);
figure; hold on;
hist([responses_sit(evoked1) responses_run(evoked1)], 20)
legend('running', 'sitting')
[p, h, stats]= signrank(responses_run, responses_sit);
title_string = sprintf( 'Sound responses (evoked - spont) z = %.2f, p = %d', stats.zval, p);
title(title_string)
xlabel('Sound responses, FR (Hz)')
legend('sit', 'run')
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
mouse_ID1 = mouse_ID(evoked1);

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

% check running effects on raw FR in RS and FS cells, evoked and spont;
evokedChange = WN1(:,1) - WN1(:,2);
spontChange = SP1(:,1) - SP1(:,2);
[p, h, stats]= ranksum(evokedChange(rs1), evokedChange(fs1));
N = sum(~isnan(evokedChange));
r = stats.zval/sqrt(N);
fprintf('\n Mean+/-SEM RS = %.4f +/- %.4f, NS = %.4f +/- %.4f', nanmean(evokedChange(rs1)), sem(evokedChange(rs1)), ...
    nanmean(evokedChange(fs1)), sem(evokedChange(fs1)))
fprintf('\n Evoked Change FR between RS and FS: z= %.4f, p = %d, r = %.4f', stats.zval, p, r)

[p, h, stats]= ranksum(spontChange(rs1), spontChange(fs1));
N = sum(~isnan(spontChange));
r = stats.zval/sqrt(N);
fprintf('\n Mean+/-SEM RS = %.4f +/- %.4f, NS = %.4f +/- %.4f', nanmean(spontChange(rs1)),sem(spontChange(rs1)), ...
    nanmean(spontChange(fs1)), sem(spontChange(fs1)))
fprintf('\n Spont Change FR between RS and FS: z= %.4f, p = %d, r = %.4f', stats.zval, p, r)

depth1 = depths(evoked1);
layer1 = nan(length(depth1),1);
for l = 1:length(CL)
    layer = CL{l};
    for d = 1:length(depth1)
        if depth1(d) > layer(1) & depth1(d) < layer(2)
            layer1(d) = l;
        end
    end
end

% change in FR all cells across layer
[p,tbl1,stats] = kruskalwallis(spontChange, layer1);
c = multcompare(stats);
title('Spont all cells, running FR change')
x=modulation_indx1_rs;
[p,tbl1,stats] = kruskalwallis(evokedChange, layer1);
c = multcompare(stats);
title('Evoked all cells, running FR change')

% change in FR RS cells
[p,tbl1,stats] = kruskalwallis(spontChange(rs1), layer1(rs1));
c = multcompare(stats);
title('Spont RS, running FR change')
x=modulation_indx1_rs;
[p,tbl1,stats] = kruskalwallis(evokedChange(rs1), layer1(rs1));
c = multcompare(stats);
title('Evoked RS, running FR change')

% change in FR FS cells
[p,tbl1,stats] = kruskalwallis(spontChange(fs1), layer1(fs1));
c = multcompare(stats);
title('Spont FS, running FR change')
x=modulation_indx1_rs;
[p,tbl1,stats] = kruskalwallis(evokedChange(fs1), layer1(fs1));
c = multcompare(stats);
title('Evoked FS, running FR change')

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

layers = nan(sum(evoked1),1);
for cl = 1:length(CL)
    layer = CL{cl}; % layer depth limits
    indx = find(depths1 >layer(1) & depths1 <layer(2)); % indices of cells within this layer
    recs2 = recs1(indx); % recordings in this layer
    cells2 = cells1(indx); % cells in this layer
    fs2 = fs1(indx); fs2 = logical(fs2); % fast spiking cells in this layer
    rs2 = rs1(indx); rs2 = logical(rs2); % regular spiking cells in this layer
    layers(indx) = 1*cl;
    
    % % evoked (WN respone)
    % means
    meanMI1(cl) = nanmean(modulation_indx1(indx)); % mean, laser off
    meanMI1L(cl) = nanmean(modulation_indx1L(indx)); % mean, laser on
    meanMI1_Laser(cl) = nanmean(modulation_indx1Laser(indx)); % mean laser effect
    
    meanMI1_rs(cl) = nanmean(modulation_indx1(indx(rs2)));
    meanMI1_fs(cl) = nanmean(modulation_indx1(indx(fs2)));
    meanMI1_rsL(cl) = nanmean(modulation_indx1L(indx(rs2)));
    meanMI1_fsL(cl) = nanmean(modulation_indx1L(indx(fs2)));
    
    semMI1(cl) = sem(modulation_indx1(indx));
    semMI1L(cl) = sem(modulation_indx1L(indx));
    semMI1_Laser(cl) = sem(modulation_indx1Laser(indx));
    
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

% Plot firing rates by cortical layer to identify what drives running effect
H = [];
for l = 1:4
    FR_evoked_means_sitting(l) = nanmean(WN1(layers == l,2));
    FR_evoked_means_running(l) = nanmean(WN1(layers == l,1));
    FR_spont_means_sitting(l) = nanmean(SP1(layers == l,2));
    FR_spont_means_running(l) = nanmean(SP1(layers == l,1));
    
    
    FR_evoked_medians_sitting(l) = nanmedian(WN1(layers == l,2));
    FR_evoked_medians_running(l) = nanmedian(WN1(layers == l,1));
    FR_spont_medians_sitting(l) = nanmedian(SP1(layers == l,2));
    FR_spont_medians_running(l) = nanmedian(SP1(layers == l,1));
    
    FR_evoked_sems_sitting(l) = sem(WN1(layers == l,2));
    FR_evoked_sems_running(l) = sem(WN1(layers == l,1));
    FR_spont_sems_sitting(l) = sem(SP1(layers == l,2));
    FR_spont_sems_running(l) = sem(SP1(layers == l,1));
    
    % stats
    [p,h,STATS] = signrank(WN1(layers==l,2), WN1(layers==l,1));
    if p < 0.0125
        H(1,l) = 1.15;
    else
         H(1,l) = NaN;
    end
    
    [p,h,STATS] = signrank(SP1(layers==l,2), SP1(layers==l,1));
    if p < 0.0125
        H(2,l) = 1.15;
    else
        H(2,l) = NaN;
    end
end

figure; subplot(2,1,1); hold on
errorbar([1:4], FR_evoked_means_sitting, FR_evoked_sems_sitting, 'ko-')
errorbar([1.2:4.2], FR_evoked_means_running, FR_evoked_sems_running, 'ko--')
plot([1.1:4.1], max(FR_evoked_means_running).*H(1,:), '*r')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Mean FR /SEM')
title('Evoked')
subplot(2,1,2); hold on;
errorbar([1:4], FR_spont_means_sitting, FR_spont_sems_sitting, 'ko-')
errorbar([1.2:4.2], FR_spont_means_running, FR_spont_sems_running, 'ko--')
plot([1.1:4.1], max(FR_spont_means_running).*H(2,:), '*r')
xticks([1:4]); xlim([0 5])
title('Spont')
xticklabels({'2/3', '4', '5', '6'});
ylabel('Mean/SEM FR')

figure; subplot(2,1,1); hold on
errorbar([1:4], FR_evoked_medians_sitting, FR_evoked_sems_sitting, 'ko-')
errorbar([1.2:4.2], FR_evoked_medians_running, FR_evoked_sems_running, 'ko--')
plot([1.1:4.1], max(FR_evoked_medians_sitting).*H(1,:), '*r')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
title('Evoked')
ylabel('Median/SEM FR')

subplot(2,1,2); hold on;
errorbar([1:4], FR_spont_medians_sitting, FR_spont_sems_sitting, 'ko-')
errorbar([1.2:4.2], FR_spont_medians_running, FR_spont_sems_running, 'ko--')
plot([1.1:4.1], max(FR_spont_medians_running).*H(2,:), '*r')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
title('Spont')
ylabel('Median/SEM FR')
xlabel('Cortical layers')

% stats
SP_change = SP1(:,1) - SP1(:,2);
WN_change = WN1(:,1) - WN1(:,2);
[p,tbl1,stats] = kruskalwallis(SP_change, layers);
c = multcompare(stats);
title('Evoked all cells, running')
[m,s]=grpstats(SP_change, layers,{'mean','sem'})
title('SP FR change across layers')

[p,tbl1,stats] = kruskalwallis(WN_change, layers);
c = multcompare(stats);
title('Evoked all cells, running')
 [m,s]=grpstats(WN_change, layers,{'mean','sem'})
title('WN FR change across layers')

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


% plot firing rate distributions
figure; 
h1 = histogram(SP1(:,1));
h1.BinWidth = .5;
hold on;
h2 = histogram(SP1(:,2));
h2.BinWidth = .5;
ylabel('Number of cells')
xlabel('Firing Rate (Hz)')
legend('running', 'sitting')
[p,h,stats] = signrank(SP1(:,1), SP1(:,2))
title_str = sprintf('Spontaneous Activity, p = %d', p)
title(title_str)

figure; 
h1 = histogram(WN1(:,1));
h1.BinWidth = 3;
hold on;
h2 = histogram(WN1(:,2));
h2.BinWidth = 3;
ylabel('Number of cells')
xlabel('Firing Rate (Hz)')
legend('running', 'sitting')
[p,h,stats] = signrank(WN1(:,1), WN1(:,2))
title_str = sprintf('Evoked Activity, p = %d', p)
title(title_str)

x=WN1(:,1) - WN1(:,2);
[p,tbl1,stats] = kruskalwallis(x, layers);
c = multcompare(stats);
title('Evoked FR change all cells, layers')
[m,s]=grpstats(x,layers,{'mean','sem'});

x=SP1(:,1) - SP1(:,2);
[p,tbl1,stats] = kruskalwallis(x, layers);
c = multcompare(stats);
title('Spont FR change all cells, layers')
[m,s]=grpstats(x,layers,{'mean','sem'});

%%
% sound modulation index
MI_sound_run = (WN1(:,1) - SP1(:,1))./ (WN1(:,1) +SP1(:,1));
MI_sound_sit = (WN1(:,2) - SP1(:,2))./ (WN1(:,2) +SP1(:,2));
MI_sound_run_laser = (WN1L(:,1) - SP1L(:,1))./ (WN1L(:,1) +SP1L(:,1));
MI_sound_sit_laser = (WN1L(:,2) - SP1L(:,2))./ (WN1L(:,2) +SP1L(:,2));
MI_sound_predicted = MI_sound_sit_laser + MI_sound_run;

% does running effect differ across recordings?
x=MI_sound_run - MI_sound_sit;
[p,tbl1,stats] = kruskalwallis(x, recs1);
c = multcompare(stats);
title('Running Effect (Sound MI run - Sound MI sit), recordings')
[m,s]=grpstats(x, recs1,{'mean','sem'});

% does sound responsiveness differ across recordings
x= MI_sound_sit;
[p,tbl1,stats] = kruskalwallis(x, recs1);
c = multcompare(stats);
title('Sound MI across , recordings')
[m,s]=grpstats(x,recs1,{'mean','sem'});

x=MI_sound_run - MI_sound_sit;
[p,tbl1,stats] = kruskalwallis(x, mouse_ID1);
c = multcompare(stats);
title('Running Effect (Sound MI run - Sound MI sit), mice')
[m,s]=grpstats(x,mouse_ID1,{'mean','sem'});


% plot MI
figure; hold on
plot(MI_sound_sit(rs1) ,MI_sound_run(rs1), 'ko')
plot(MI_sound_sit(fs1) ,MI_sound_run(fs1), 'go')
plot(nanmean(MI_sound_sit), nanmean(MI_sound_run), 'ro', 'MarkerSize', 10)
plot(nanmean(MI_sound_sit), nanmean(MI_sound_run), 'ro', 'MarkerSize', 10)
ylabel('sound MI running'); xlabel('sound MI sitting')
[r, p] = corr(MI_sound_run,MI_sound_sit, 'Type','Spearman','Rows', 'complete')
plot([-1 1], [-1 1], 'r-')
plot([-1 1], [0 0], 'k--')
plot([0 0], [-1 1], 'k--')
pbaspect([1 1 1]);  set(gcf, 'PaperPositionMode', 'auto');


[h1,x] = hist(MI_sound_sit, [-1:.1:1]);
h = smooth(h1,3);
figure; hold on
plot(x,h);
[h2,x] = hist(MI_sound_run, [-1:.1:1]);
h = smooth(h2,3);
plot(x,h, '--');
title(' On responses, laser off')
legend('sitting - laser off', 'running - laser off')
xlabel('Sound Modulation Index ')
ylabel('Number of cells')


fprintf('running vs sitting, laser off')
[p,h stats] =ranksum(MI_sound_run, MI_sound_sit)

fprintf('running vs sitting, laser on')
[p,h stats] =ranksum(MI_sound_run_laser, MI_sound_sit_laser)


for cl = 1:length(CL)
    layer = CL{cl}; % layer depth limits
    indx = find(depths1 >layer(1) & depths1 <layer(2)); % indices of cells within this layer
    rs2 = rs1(indx); fs2 = fs1(indx);
    
    
    meanMI_sound_run(cl) = nanmean(MI_sound_run(indx));
    meanMI_sound_sit(cl) = nanmean(MI_sound_sit(indx));
    m1=MI_sound_run(indx);
    m2 = MI_sound_sit(indx);
    meanMI_sound_run_rs(cl) = nanmean(m1(rs2));
    meanMI_sound_sit_rs(cl) = nanmean(m2(rs2));
    meanMI_sound_run_fs(cl) = nanmean(m1(fs2));
    meanMI_sound_sit_fs(cl) = nanmean(m2(fs2));
    
    
    semMI_sound_run(cl) = sem(MI_sound_run(indx));
    semMI_sound_sit(cl) = sem(MI_sound_sit(indx));
    semMI_sound_run_rs(cl) = sem(m1(rs2));
    semMI_sound_sit_rs(cl) = sem(m2(rs2));
    semMI_sound_run_fs(cl) = sem(m1(fs2));
    semMI_sound_sit_fs(cl) = sem(m2(fs2));
    
    n_layer(cl) = sum(~isnan(MI_sound_run(indx)));
    
    [p,h,STATS] = signrank(MI_sound_run(indx), MI_sound_sit(indx));
    fprintf('\n Layer %d , signedrank = %d, p = %d\n', cl, STATS.signedrank, p)
    
    meanMI_sound_runL(cl) = nanmean(MI_sound_run_laser(indx));
    meanMI_sound_sitL(cl) = nanmean(MI_sound_sit_laser(indx));
    semMI_sound_runL(cl) = sem(MI_sound_run_laser(indx));
    semMI_sound_sitL(cl) = sem(MI_sound_sit_laser(indx));
    n_layerL(cl) = sum(~isnan(MI_sound_run_laser(indx)));
end

figure; hold on
errorbar([1:4], meanMI_sound_sit, semMI_sound_sit, 'ko-');
errorbar([1.1:4.1], meanMI_sound_run, semMI_sound_run, 'ko--');
xlabel('Cortical Layer')
plot([0 5], [0 0], 'k--')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Sound Modulation Index (mean/SEM)');
title_string = sprintf( ' Sound MI, n = %d, %d, %d, %d',  n_layer);
title(title_string)
legend( 'sitting', 'running')

figure; hold on
errorbar([1:4], meanMI_sound_sit_rs, semMI_sound_sit_rs, 'ko-');
errorbar([1:4], meanMI_sound_run_rs, semMI_sound_run_rs, 'ko--');
errorbar([1:4], meanMI_sound_sit_fs, semMI_sound_sit_fs, 'go-');
errorbar([1:4], meanMI_sound_run_fs, semMI_sound_run_fs, 'go--');
xlabel('Cortical Layer')
plot([0 5], [0 0], 'k--')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Sound Modulation Index (mean/SEM)');
title_string = sprintf( 'On response MI, n = %d, %d, %d, %d',  n_layer);
title(title_string)
legend( 'sitting', 'running')

% distribution scatter plot
figure; hold on
plot(WN1(rs1,2), WN1(rs1,1), 'ko')
plot(WN1(fs1,2), WN1(fs1,1), 'go')
plot(nanmean(WN1(:,2)), nanmean(WN1(:,1)), 'ro', 'MarkerSize', 10)
plot(nanmedian(WN1(:,2)), nanmedian(WN1(:,1)), 'mo', 'MarkerSize', 10)
maxFR = max(max(WN1));
plot([0 maxFR], [0 maxFR], 'r-')
plot(nanmean(WN1(:,2)), nanmean(WN1(:,1)), 'ro')
[p,h,stats] = signrank(WN1(:,1),  WN1(:,2));
title_str = sprintf('Evoked Activity, p = %d', p);
title(title_str)
xlabel('FR sitting (Hz)'); ylabel('FR running (Hz)')
legend('Regular spiking', 'Narrow spiking', 'Mean', 'Median')
pbaspect([1 1 1]); set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Position',  [260 124 902 864])

figure; hold on;
plot(SP1(rs1,2), SP1(rs1,1), 'ko')
plot(SP1(fs1,2), SP1(fs1,1), 'go')
plot(nanmean(SP1(:,2)), nanmean(SP1(:,1)), 'ro', 'MarkerSize', 10)
plot(nanmedian(SP1(:,2)), nanmedian(SP1(:,1)), 'mo', 'MarkerSize', 10)
maxFR = max(max([SP1(:,2), SP1(:,1)]));
plot([0 maxFR], [0 maxFR], 'r-')
plot(nanmean(SP1(:,2)), nanmean(SP1(:,1)), 'ro')
[p,h,stats] = signrank(SP1(:,1),  SP1(:,2));
title_str = sprintf('Spont Activity, p = %d', p);
title(title_str)
xlabel('FR sitting (Hz)'); ylabel('FR running (Hz)')
legend('Regular spiking', 'Narrow spiking', 'Mean', 'Median')
pbaspect([1 1 1]); set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Position',  [260 124 902 864])

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
legend({'Regular spiking', 'Narrow spiking'})
xticks([1:2])
xticklabels({'Spont run + laser off', 'Evoked run + laser off'})
ylabel('Modulation Index')

meanMI_sound = nanmean(MI_sound_sit);
semMI_sound = sem(MI_sound_sit);

%% Test results on different repetitions

%%  STOP HERE, THE REST IS OLD
for i = 1:100
    cc = 0;
    for c = 1 : length(data)
        if evoked1(c) ==1
            FR_sit = data(c).mNON_off;
            FR_run = data(c).mNON_on;
            
            SP_sit = data(c).mSSoff;
            SP_run = data(c).mSSon;
            
            num_reps_sit = min([size(FR_sit,1) size(SP_sit,1)]);
            num_reps_run = min([size(FR_run,1) size(SP_run,1)]);
            
            if num_reps_sit > 6 && num_reps_run > 6
                indx = randi(num_reps_sit, 1, num_reps_run); % generate indices to match running trials
                cc = cc+1;
                MI_sit(cc) = (nanmean(nanmean(FR_sit(indx,:))) - nanmean(nanmean(SP_sit(indx,:)))) / (nanmean(nanmean(FR_sit(indx,:))) + nanmean(nanmean(SP_sit(indx,:))));
                MI_run(cc) = (nanmean(nanmean(FR_run)) - nanmean(nanmean(SP_run)))/ (nanmean(nanmean(FR_run)) + nanmean(nanmean(SP_run)));
                
            end
        end
    end
    [P,H,STATS] = signrank(MI_run, MI_sit);
    MI_sit_all(i, 1) = nanmean(MI_sit);
    MI_sit_all(i, 2) = sem(MI_sit);
    p_values(i) = P;
    z_values(i) = STATS.zval;
end
effect_size = abs(z_values)/sqrt(length(MI_run))

figure; hold on
bar([1 2], [nanmean(MI_sit_all(:,1)) nanmean(MI_run)], 'BarWidth', .4);
errorbar([1 2], [nanmean(MI_sit_all(:,1)) nanmean(MI_run)], [nanmean(MI_sit_all(:,2)) sem(MI_run)])
xticks([1 2])
xticklabels({'sitting', 'running'})
title('matched sample size')
ylabel('Sound Modulation Index (Mean/Sem)')

% data stats
fprintf('\nN cells Sound MI sitting =%d, N running = %d \n',  sum(~isnan(MI_sound_sit)),  sum(~isnan(MI_sound_run)))
fprintf('\nN cells WN =%d, N running = %d \n',  sum(~isnan(WN1(:,2))),  sum(~isnan(WN1(:,1))))
fprintf('\nN cells Spont sitting =%d, N running = %d \n',  sum(~isnan(SP1(:,2))),  sum(~isnan(SP1(:,1))))

