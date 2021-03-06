% compute modulation index WN responses. last editted 01.21.2020 ira
% variables gethered with getWNresponses.m

% three states: sit + small pupil, sit + large pupil, running  + large
%

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
recs = []; cells=[]; cell_number= 1:length(data);

evoked = zeros(length(data),4); zstats = zeros(length(data),4);
color1 = [0.7 0.7 0.7]; color2 = [0.2 0.2 0.2];
nreps_check_running = 6; %number of repetitions in each condition for comparison
CL = {[0 301], [300 401], [400 601], [600 2000]};

cdVIP; load('Silence_DistanceCorr_dirs.mat')
% collect waveform, SNR, and uQ from cells. Must match length of data.
try
    cd(variables_dir)
    load('CellsInfo.mat') % run WNModularionIndex_VIP or _runnning to get this variable saved
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

evoked1 = logical(evoked(:,1) & zstats(:,1)>0); 
cd(variables_dir)
save('evoked_indx_epistatic.mat', 'evoked1')

% On responses, suppressed and activated distributions
responses =  (meanON - SP)./(meanON+SP);
figure; hist(responses)
legend('running', 'sitting')
[p, h, stats]= signrank(responses(:,1), responses(:,2));
title_string = sprintf( 'Sound Modulation Index all cells On responses z = %.2f, p = %d', stats.zval, p);
title(title_string)
xlabel('Sound MI')
ylabel('Number of cells')

% Full responses, suppressed and activated distributions
responses =  (meanWN - SP)./(meanWN+SP);
figure; hist(responses)
legend('running', 'sitting')
[p, h, stats]= signrank(responses(:,1), responses(:,2));
title_string = sprintf( 'Sound Modulation Index all cells Full responses z = %.2f, p = %d', stats.zval, p);
title(title_string)
xlabel('Sound MI')
ylabel('Number of cells')

% subselect data with evoked responses only
WN1 = meanON(evoked1,:);
WN1L = meanONL(evoked1,:);
SP1 = SP(evoked1,:);
SP1L = SPL(evoked1,:);
meanPreStim1 = meanPreStim(evoked1,:);
meanPreStim1L = meanPreStimL(evoked1,:);
meanSilentSound1 = meanSilentSound(evoked1,:);
depths1 = depths(evoked1);
fs1 = fs(evoked1);
rs1 = Rs(evoked1);
rs1 = logical(rs1); fs1 = logical(fs1);
cell_number1 = cell_number(:,evoked1);
cell_number1 = cell_number1';

EvokedWN = WN1;
modulation_indx1 = (EvokedWN(:,1) - EvokedWN(:,2))./(EvokedWN(:,1) + EvokedWN(:,2));
modulation_indx1_sp = (SP1(:,1) - SP1(:,2))./(SP1(:,1) + SP1(:,2));

%modulation index by cell type
modulation_indx1 = (WN1(:,1) -WN1(:,2))./(WN1(:,1) + WN1(:,2)); %running, laser off trials
modulation_indx1L = (WN1L(:,2) -WN1(:,2))./(WN1L(:,2) + WN1(:,2)); %running, laser on trials

modulation_indx1_sp = (SP1(:,1) -SP1(:,2))./(SP1(:,1) + SP1(:,2)); %spont running effect laser off
modulation_indx1_spL = (SP1L(:,2) -SP1(:,2))./(SP1(:,2) + SP1(:,2)); %spont running effect laser on

MI_sp_plus = (SP1L(:,1) - SP1(:,2))./(SP1L(:,1) + SP1(:,2));
MI_evoked_plus = (WN1L(:,1) - WN1(:,2))./(WN1L(:,1) + WN1(:,2));

modulation_indx1_rs = modulation_indx1(rs1); %running
modulation_indx1_fs = modulation_indx1(fs1); %running

modulation_indx1_sp_rs = modulation_indx1_sp(rs1); %running
modulation_indx1_sp_fs = modulation_indx1_sp(fs1); %running

modulation_indx1_rsL = modulation_indx1L(rs1); %running
modulation_indx1_fsL = modulation_indx1L(fs1); %running

modulation_indx1_sp_rsL = modulation_indx1_spL(rs1); %running
modulation_indx1_sp_fsL = modulation_indx1_spL(fs1); %running

% Is the effect of running the same across recordings? Cell types?
recs1= (recs(evoked1)); cells1 = cells(evoked1);

x=modulation_indx1_sp_rs;
[p,tbl1,stats] = kruskalwallis(x, recs1(rs1));
c = multcompare(stats);
title('Spont RS, running')
x=modulation_indx1_rs;
[p,tbl1,stats] = kruskalwallis(x, recs1(rs1));
c = multcompare(stats);
title('Evoked RS, running')

x=modulation_indx1_sp_fs;
[p,tbl1,stats] = kruskalwallis(x, recs1(fs1));
c = multcompare(stats);
title('Spont FS, running')
x=modulation_indx1_fs;
[p,tbl1,stats] = kruskalwallis(x, recs1(fs1));
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
[m,s]=grpstats(x,recs1,{'mean','sem'});

%% Layer analysis
MODULATION_INDX1 = []; MODULATION_INDX1L = []; %  for stats across cortical layers, evoked
MODULATION_sp_INDX1 = []; MODULATION_sp_INDX1L = []; %  for stats across cortical layers, spont
rs_layers = []; fs_layers = [];
for cl = 1:length(CL)
    layer = CL{cl}; % layer depth limits
    indx = find(depths1 >layer(1) & depths1 <layer(2)); % indices of cells within this layer
    recs2 = recs1(indx); % recordings in this layer
    cells2 = cells1(indx); % cells in this layer
    fs2 = fs1(indx); fs2 = logical(fs2); % fast spiking cells in this layer
    rs2 = rs1(indx); rs2 = logical(rs2); % regular spiking cells in this layer
    rs_layers = [rs_layers; rs2];
    fs_layers = [fs_layers; fs2];
    
    % % evoked (WN respone)
    % means
    meanMI1(cl) = nanmean(modulation_indx1(indx)); % mean, laser off
    meanMI1L(cl) = nanmean(modulation_indx1L(indx)); % mean, laser on
    
    meanMI1_rs(cl) = nanmean(modulation_indx1(indx(rs2)));
    meanMI1_fs(cl) = nanmean(modulation_indx1(indx(fs2)));
    meanMI1_rsL(cl) = nanmean(modulation_indx1L(indx(rs2)));
    meanMI1_fsL(cl) = nanmean(modulation_indx1L(indx(fs2)));
    
    semMI1(cl) = sem(modulation_indx1(indx));
    semMI1L(cl) = sem(modulation_indx1L(indx));
    
    semMI1_rs(cl) = sem(modulation_indx1(indx(rs2)));
    semMI1_fs(cl) = sem(modulation_indx1(indx(fs2)));
    semMI1_rsL(cl) = sem(modulation_indx1L(indx(rs2)));
    semMI1_fsL(cl) = sem(modulation_indx1L(indx(fs2)));
    
    % collect modulation indices across layers for further stats [MI cl]
    MODULATION_INDX1 = [MODULATION_INDX1; modulation_indx1(indx) ones(length(indx),1)*cl];
    MODULATION_INDX1L = [MODULATION_INDX1L; modulation_indx1L(indx) ones(length(indx),1)*cl];
    
    n_layer_evoked(cl) = length(find(~isnan(modulation_indx1(indx))==1));
    n_layer_evokedL(cl) = length(find(~isnan(modulation_indx1L(indx))==1));
    
    rs_cl(cl) = sum(~isnan(modulation_indx1(indx(rs2)))); % number of regular spiking cells in this layer
    fs_cl(cl) = sum(~isnan(modulation_indx1(indx(fs2)))); % number of fast spiking cells in this layer
    
    % addative effect
    meanMI_evoked_plus(cl) = nanmean(MI_evoked_plus(indx));
    semMI_evoked_plus(cl) = sem(MI_evoked_plus(indx));
    n_layer_evoked_plus(cl) = length(find(~isnan(MI_evoked_plus(indx))==1)); %number of cell
    
    % % spont activity
    meanMI1_sp(cl) = nanmean(modulation_indx1_sp(indx));
    meanMI1_spL(cl) = nanmean(modulation_indx1_spL(indx));
    meanMI1_sp_rs(cl) = nanmean(modulation_indx1_sp(indx(rs2)));
    meanMI1_spL_rs(cl) = nanmean(modulation_indx1_spL(indx(rs2)));
    meanMI1_sp_fs(cl) = nanmean(modulation_indx1_sp(indx(fs2)));
    meanMI1_spL_fs(cl) = nanmean(modulation_indx1_spL(indx(fs2)));
    
    semMI1_sp(cl) = sem(modulation_indx1_sp(indx));
    semMI1_spL(cl) = sem(modulation_indx1_spL(indx));
    semMI1_sp_rs(cl) = sem(modulation_indx1_sp(indx(rs2)));
    semMI1_spL_rs(cl) = sem(modulation_indx1_spL(indx(rs2)));
    semMI1_sp_fs(cl) = sem(modulation_indx1_sp(indx(fs2)));
    semMI1_spL_fs(cl) = sem(modulation_indx1_spL(indx(fs2)));
    
    n_layer_sp(cl) = length(find(~isnan(modulation_indx1_sp(indx))==1)); %number of cells in each layer
    n_layer_spL(cl) = length(find(~isnan(modulation_indx1_spL(indx))==1));
    
    rs_sp_cl(cl) = sum(~isnan(modulation_indx1_sp(indx(rs2))));
    fs_sp_cl(cl) = sum(~isnan(modulation_indx1_sp(indx(fs2))));
    
    % additive effect
    meanMI_sp_plus(cl) = nanmean(MI_sp_plus(indx));
    semMI_sp_plus(cl) = sem(MI_sp_plus(indx));
    n_layer_sp_plus(cl) = length(find(~isnan(MI_sp_plus(indx))==1));
    
    MODULATION_sp_INDX1 = [MODULATION_sp_INDX1; modulation_indx1_sp(indx) ones(length(indx),1)*cl];
    MODULATION_sp_INDX1L = [MODULATION_sp_INDX1L; modulation_indx1_spL(indx) ones(length(indx),1)*cl];
end
rs_layers = logical(rs_layers);
fs_layers = logical(fs_layers);

% stats
[p,tbl1,stats] = kruskalwallis(MODULATION_INDX1(:,1), MODULATION_INDX1(:,2));
 c = multcompare(stats);
 title('Running MI, evoked, all cells')
 
 [p,tbl1,stats] = kruskalwallis(MODULATION_sp_INDX1(:,1), MODULATION_sp_INDX1(:,2));
 c = multcompare(stats);
 title('Running MI, spont, all cells')
 
 [p,tbl1,stats] = kruskalwallis(MODULATION_INDX1(rs_layers,1), MODULATION_INDX1(rs_layers,2));
 c = multcompare(stats);
 title('Running MI, evoked, RS cells')
 
 [p,tbl1,stats] = kruskalwallis(MODULATION_sp_INDX1(rs_layers,1), MODULATION_sp_INDX1(rs_layers,2));
 c = multcompare(stats);
 title('Running MI, spont, RS cells')
 
  [p,tbl1,stats] = kruskalwallis(MODULATION_INDX1(fs_layers,1), MODULATION_INDX1(fs_layers,2));
 c = multcompare(stats);
 title('Running MI, evoked, FS cells')
 
 [p,tbl1,stats] = kruskalwallis(MODULATION_sp_INDX1(fs_layers,1), MODULATION_sp_INDX1(fs_layers,2));
 c = multcompare(stats);
 title('Running MI, spont, FS cells')

% state by layer without layer
figure; hold on; subplot(2,1,1); hold on
errorbar([1:4], meanMI1, semMI1, 'ko-');
xlabel('Cortical Layer')
plot([0 5], [0 0], '--')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Modulation Index - Run');
title_string = sprintf('WN On response, n = %d, %d, %d, %d',  n_layer_evoked);
title(title_string)

subplot(2,1,2); hold on
errorbar([1:4], meanMI1_sp, semMI1_sp, 'ko-');
xlabel('Cortical Layer')
plot([0 5], [0 0], '--')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Modulation Index - Run');
title_string = sprintf('Spont On response, n = %d, %d, %d, %d',  n_layer_sp);
title(title_string)

%run modulation by cell type
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
ylabel('Modulation Index - Laser effect');
title_string = sprintf('WN On response, laser on, n = %d, %d, %d, %d',  n_layer_evokedL);
title(title_string)
subplot(2,1,2); hold on
errorbar([1:4], meanMI1_spL, semMI1_spL, 'bo-');
xlabel('Cortical Layer')
plot([0 5], [0 0], 'k--')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Modulation Index - Laser effect');
title_string = sprintf('Spont On response, laser on, n = %d, %d, %d, %d',  n_layer_spL);
title(title_string)
 set(gcf, 'PaperPositionMode', 'auto');

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


%%
% sound modulation index
MI_sound_run = (WN1(:,1) - SP1(:,1))./ (WN1(:,1) +SP1(:,1));
MI_sound_sit = (WN1(:,2) - SP1(:,2))./ (WN1(:,2) +SP1(:,2));
MI_sound_run_laser = (WN1L(:,1) - SP1L(:,1))./ (WN1L(:,1) +SP1L(:,1));
MI_sound_sit_laser = (WN1L(:,2) - SP1L(:,2))./ (WN1L(:,2) +SP1L(:,2));
MI_sound_predicted = MI_sound_sit_laser + MI_sound_run;

% for regression
y = MI_sound_predicted;
indx = isnan(MI_sound_predicted);
y(indx)=[];
x = MI_sound_run_laser;
x(indx) = [];
x=[x ones(length(x),1)];
[B,BINT,R,RINT,STATS] = regress(y,x) ;

%find 99 cells used in this analysis, save their indices, 
indx1 = find(~isnan(MI_sound_run_laser)==1); 
indx2 = find(~isnan(MI_sound_predicted)==1);
indx = intersect(indx1, indx2); 
cell_number1 = cell_number1(indx);
cd(variables_dir)
save('ModulationIndices.mat', 'MI_sound_run', 'MI_sound_sit', 'MI_sound_run_laser', 'MI_sound_sit_laser', 'MI_sound_predicted')
save('example_cells_epistatic.mat', 'cell_number1')

%%%%% plot sound modulation indices %%%%%%%%
% These plots are not very meaningful because they will both show high
% correlations: responsive cells will continue to be responsive in all
% conditions
figure; hold on
plot(MI_sound_predicted(rs1) ,MI_sound_run_laser(rs1), 'ko')
plot(MI_sound_predicted(fs1) ,MI_sound_run_laser(fs1), 'go')
xlabel('sound MI combined (runnning + laser)'); ylabel('sound MI running + laser on trials')
[r, p] = corr(MI_sound_predicted,MI_sound_run_laser, 'Type','Spearman','Rows', 'complete')
plot([-2 2], [-2 2], 'r-')
plot([-2 2], [0 0], 'k--')
plot([0 0], [-2 2], 'k--')
pbaspect([1 1 1]);  set(gcf, 'PaperPositionMode', 'auto');
title_string = sprintf( 'rho = %.4f, p = %.4f', r,p);
title(title_string); legend('Regular spiking', 'Narrow spiking')

figure; hold on
plot(MI_sound_run(rs1) ,MI_sound_sit_laser(rs1), 'ko')
plot(MI_sound_run(fs1) ,MI_sound_sit_laser(fs1), 'go')
xlabel('sound MI run'); ylabel('sound MI laser on')
[r, p] = corr(MI_sound_run,MI_sound_sit_laser, 'Type','Spearman','Rows', 'complete')
plot([-2 2], [-2 2], 'r-')
plot([-2 2], [0 0], 'k--')
plot([0 0], [-2 2], 'k--')
pbaspect([1 1 1]);  set(gcf, 'PaperPositionMode', 'auto');
title_string = sprintf( 'rho = %.4f, p = %.4f', r,p);
title(title_string); legend('Regular spiking', 'Narrow spiking')


%%%% plot difference in MI %%%%
% To look more closely at differences in MI in different conditions 
run_diff =  MI_sound_run - MI_sound_sit;
laser_diff = MI_sound_sit_laser - MI_sound_sit;
predicted_diff = run_diff + laser_diff;
actual_diff =  MI_sound_run_laser - MI_sound_sit;

% for regression
y = actual_diff;
indx = isnan(actual_diff);
y(indx)=[];
x = run_diff;
x(indx) = [];
X = [ones(length(x),1) x];
[B,BINT,R,RINT,STATS] = regress(y,X) ;

y = actual_diff;
indx = isnan(actual_diff);
y(indx)=[];
x = laser_diff;
x(indx) = [];
X = [ones(length(x),1) x];
[B,BINT,R,RINT,STATS] = regress(y,X);

% for regression
y = laser_diff;
indx = isnan(laser_diff);
y(indx)=[];
x = run_diff;
x(indx) = [];
X = [ones(length(x),1) x];
[B,BINT,R,RINT,STATS] = regress(y,X) ;
regressline = laser_diff*B(2) + B(1);

figure; hold
plot(run_diff(rs1), laser_diff(rs1), 'ko', 'MarkerSize', 8)
plot(run_diff(fs1), laser_diff(fs1), 'go', 'MarkerSize', 8)
plot(run_diff, regressline, 'r-')
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


% for regression
y = actual_diff;
indx = isnan(actual_diff);
y(indx)=[];
x = predicted_diff;
x(indx) = [];
X=[ones(length(x),1) x];
[B,BINT,R,RINT,STATS] = regress(y,X) ;
regressline = predicted_diff*B(2) + B(1);

figure; hold on
plot(predicted_diff(rs1), actual_diff(rs1), 'ko', 'MarkerSize', 8)
plot(predicted_diff(fs1), actual_diff(fs1), 'go', 'MarkerSize', 8)
xlabel({'Computed Combined Effect';'(running effect +  laser effect)'})
ylabel({'Recorded Combined Effect';'(running + laser on trials)'})
plot(predicted_diff, regressline)
plot([-4 2], [0 0], 'k--')
plot([0 0], [-4 2], 'k--')
plot([-4 4], [-4 4], 'k--')
[r, p] = corr(actual_diff,predicted_diff, 'Type','Spearman','Rows', 'complete')
title_string = sprintf( 'rho = %.4f, p = %.4f', r,p);
ylim([-4 2])
xlim([-4 2])
title(title_string); legend('Regular spiking', 'Narrow spiking', 'Regression Line')
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

%% Firing Rate Analysis %%%
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

laser_effect_evokedFR = WN1L(:,2) - WN1(:,2);
laser_effect_spontFR = SP1L(:,2) - SP1(:,2);

run_effect_evokedFR = WN1(:,1) - WN1(:,2);
run_effect_spontFR = SP1L(:,1) - SP1(:,2);

run_laser_effect_evokedFR = WN1L(:,1) - WN1(:,2);
run_laser_effect_spontFR = SP1L(:,1) - SP1(:,2);

predicted_effect_evokedFR = laser_effect_evokedFR  + run_effect_evokedFR;
predicted_effect_spontFR = laser_effect_spontFR + run_effect_spontFR ;

[r, p] = corr(run_laser_effect_evokedFR, predicted_effect_evokedFR, 'Type','Spearman','Rows', 'complete');
[r, p] = corr(run_laser_effect_spontFR, predicted_effect_spontFR, 'Type','Spearman','Rows', 'complete');

% laser and runnign effect
figure; 
subplot(2,1,1); hold on
plot(run_effect_evokedFR(rs1), laser_effect_evokedFR(rs1),  'ko', 'MarkerSize', 8)
plot(run_effect_evokedFR(fs1), laser_effect_evokedFR(fs1), 'go', 'MarkerSize', 8)
xlabel('Change in FR (running)')
ylabel('Change in FR (laser on)')
plot([min(run_effect_evokedFR) max(run_effect_evokedFR)], [0 0], 'k--')
plot([0 0], [min(laser_effect_evokedFR) max(run_effect_evokedFR)], 'k--')
[r, p] = corr(run_effect_evokedFR, laser_effect_evokedFR, 'Type','Spearman','Rows', 'complete')
title_string = sprintf( 'rho = %.4f, p = %.4f', r,p);
title(title_string); legend('Regular spiking', 'Narrow spiking')
pbaspect([1 1 1]);  set(gcf, 'PaperPositionMode', 'auto');
xlim([min(run_effect_evokedFR) max(run_effect_evokedFR)]);
ylim([min(laser_effect_evokedFR) max(run_effect_evokedFR)])

subplot(2,1,2); hold on
plot(run_effect_spontFR(rs1), laser_effect_spontFR(rs1),  'ko', 'MarkerSize', 8)
plot(run_effect_spontFR(fs1), laser_effect_spontFR(fs1), 'go', 'MarkerSize', 8)
xlabel('Change in FR (running)')
ylabel('Change in FR (laser on)')
plot([min(run_effect_spontFR) max(run_effect_spontFR)], [0 0], 'k--')
plot([0 0], [min(laser_effect_spontFR) max(run_effect_spontFR)], 'k--')
[r, p] = corr(run_effect_spontFR, laser_effect_spontFR, 'Type','Spearman','Rows', 'complete')
title_string = sprintf( 'rho = %.4f, p = %.4f', r,p);
title(title_string); legend('Regular spiking', 'Narrow spiking')
pbaspect([1 1 1]);  set(gcf, 'PaperPositionMode', 'auto');
xlim([min(run_effect_spontFR) max(run_effect_spontFR)])
ylim([min(laser_effect_spontFR) max(run_effect_spontFR)])

% predicted and actual

figure; 
subplot(2,1,1); hold on
plot(predicted_effect_evokedFR(rs1), run_laser_effect_evokedFR(rs1),  'ko', 'MarkerSize', 8)
plot(predicted_effect_evokedFR(fs1), run_laser_effect_evokedFR(fs1), 'go', 'MarkerSize', 8)
xlabel('Change in FR (running)')
ylabel('Change in FR (laser on)')
plot([min(predicted_effect_evokedFR) max(run_laser_effect_evokedFR)], [0 0], 'k--')
plot([0 0], [min(predicted_effect_evokedFR) max(run_laser_effect_evokedFR)], 'k--')
[r, p] = corr(predicted_effect_evokedFR, run_laser_effect_evokedFR, 'Type','Spearman','Rows', 'complete');
title_string = sprintf( 'rho = %.4f, p = %.4f', r,p);
title(title_string); legend('Regular spiking', 'Narrow spiking')
pbaspect([1 1 1]);  set(gcf, 'PaperPositionMode', 'auto');
xlim([min(predicted_effect_evokedFR) max(run_laser_effect_evokedFR)]);
ylim([min(predicted_effect_evokedFR) max(run_laser_effect_evokedFR)])

subplot(2,1,2); hold on
plot(predicted_effect_spontFR(rs1), run_laser_effect_spontFR(rs1),  'ko', 'MarkerSize', 8)
plot(predicted_effect_spontFR(fs1), run_laser_effect_spontFR(fs1), 'go', 'MarkerSize', 8)
xlabel('Change in FR (running)')
ylabel('Change in FR (laser on)')
plot([min(predicted_effect_spontFR) max(run_laser_effect_spontFR)], [0 0], 'k--')
plot([0 0], [min(predicted_effect_spontFR) max(run_laser_effect_spontFR)], 'k--')
[r, p] = corr(predicted_effect_spontFR, run_laser_effect_spontFR, 'Type','Spearman','Rows', 'complete');
title_string = sprintf( 'rho = %.4f, p = %.4f', r,p);
title(title_string); legend('Regular spiking', 'Narrow spiking')
pbaspect([1 1 1]);  set(gcf, 'PaperPositionMode', 'auto');
xlim([min(predicted_effect_spontFR) max(run_laser_effect_spontFR)])
ylim([min(predicted_effect_spontFR) max(run_laser_effect_spontFR)])

[r, p] = corr(predicted_effect_evokedFR, run_laser_effect_evokedFR, 'Type','Spearman','Rows', 'complete')
[r, p] = corr(predicted_effect_evokedFR(rs1), run_laser_effect_evokedFR(rs1), 'Type','Spearman','Rows', 'complete')
[r, p] = corr(predicted_effect_evokedFR(fs1), run_laser_effect_evokedFR(fs1), 'Type','Spearman','Rows', 'complete')
[r, p] = corr(predicted_effect_spontFR, run_laser_effect_spontFR, 'Type','Spearman','Rows', 'complete')
[r, p] = corr(predicted_effect_spontFR(rs1), run_laser_effect_spontFR(rs1), 'Type','Spearman','Rows', 'complete')
[r, p] = corr(predicted_effect_spontFR(fs1), run_laser_effect_spontFR(fs1), 'Type','Spearman','Rows', 'complete')
%% modulation index separately for evoked and spont
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
semMI_sound_laser = sem(MI_sound_sit_laser);

meanMI_sound = nanmean(MI_sound_sit);
semMI_sound = sem(MI_sound_sit);

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
plot( MI1_runLaser(rs1), MI1_predicted(rs1), 'ko')
plot( MI1_runLaser(fs1), MI1_predicted(fs1), 'go')
plot([0 0], [-2 2], 'k--'); plot( [-2 2], [0 0], 'k--')
lsline
xlabel('predicted MI'); ylabel('VIP + running MI')
[r, p] = corr( MI1_runLaser, MI1_predicted, 'Type','Spearman','Rows', 'complete');
title_string = sprintf( 'Evoked, rho = %.4f, p = %d', r, p);
title(title_string)
subplot(2,1,2); hold on
plot( MI1_sp_runLaser(rs1),MI1_sp_predicted(rs1), 'ko')
plot( MI1_sp_runLaser(fs1), MI1_sp_predicted(fs1), 'go')
plot([0 0], [-2 2], 'k--'); plot( [-2 2], [0 0], 'k--')
lsline
xlabel('predicted MI'); ylabel('VIP + running MI')
[r, p] = corr(MI1_sp_runLaser, MI1_sp_predicted, 'Type','Spearman','Rows', 'complete');
title_string = sprintf( 'Spont, rho = %.4f, p = %d', r, p);
title(title_string)
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Position',  [260 124 902 864])

% regression
y = MI_sound_predicted;
indx = isnan(MI_sound_predicted);
y(indx)=[];
x = MI_sound_run_laser;
x(indx) = [];
x=[x ones(length(x),1)];
 [B,BINT,R,RINT,STATS] = regress(y,x) ;
 

 



