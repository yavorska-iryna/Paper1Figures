% compute sound modulation index of WN responses
% VIP effect only (dirs = 1:30 )

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
sit_color = 4;
run_color = 2;
laser_color = 11;

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
recs = []; cells=[]; mouse_ID = {};

evoked = zeros(length(data),4); zstats = zeros(length(data),4);
color1 = [0.7 0.7 0.7]; color2 = [0.2 0.2 0.2];
nreps_check_running = -1; %number of repetitions in each condition for comparison
CL = {[0 301], [300 401], [400 601], [600 2000]}; id =0;

cdVIP; load('Silence_DistanceCorr_dirs.mat'); 
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
    if data(cc).dir < 30
        if data(cc).dir~=0  && ~isempty(data(cc).nrepsWNMoff)
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
            
            recs = [recs data(cc).dir]; % unique recordings, all
            cells = [cells data(cc).cell]; % cell numbers
            
            if cc == 1 || recs(cc) - recs(cc-1) ~= 0
                id = id +1;
                cd(PINPdirs{recs(cc)})
                load('notebook.mat')
                if strcmp(nb.mouseID, '') == 0
                    mouse_ID{id} = nb.mouseID;
                else
                    load('dirs.mat')
                    for d = 1: length(dirs)
                        cd(dirs{d})
                        load('notebook.mat')
                        if strcmp(nb.mouseID,'') == 0
                            mouse_ID{id} = nb.mouseID;
                        end
                    end
                end   
            end
            
        else
            recs = [recs NaN]; % unique recordings, all
            cells = [cells data(cc).cell]; % cell numbers
            
        end
    end

end

allDirsM = unique(WNdirsM); % recording that passed min running trials
all_mice = unique(mouse_ID);

%% ANALYSIS %%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% evoked responses indices

% which cells to include? 1- On response, 2 - Sustained, 3-  Off response,
% 4 - full 600 ms. zstats > 0 activated, <0 suppressed
evoked1 = logical(evoked(:,1) & zstats(:,1)>0);

% Looked at the evoked responses without spontaneous activity
% figure 1 - On responses
responses1 =  meanON(:,2) - SP(:,2);
responses2 =  meanONL(:,2) - SPL(:,2);
figure; hold on
hist([responses1(evoked1) responses2(evoked1)], 100)
legend('laser off', 'laser on')
[p, h, stats]= signrank(responses1(evoked1), responses2(evoked1));
title_string = sprintf( 'On responses without spont z = %.2f, p = %d', stats.zval, p);
title(title_string)
xlabel('Firing Rate (Hz)')
ylabel('Number of cells')
xlim([0 max(responses1(evoked1))])

% figure 2 - Full WN responses
evoked2 = logical(evoked(:,4) & zstats(:,4)>0);

responses1 =  meanWN(:,2) - SP(:,2);
responses2 =  meanWNL(:,2) - SPL(:,2);
figure; hist([responses1(evoked2) responses2(evoked2)], 100)
legend('laser off', 'laser on')
[p, h, stats]= signrank(responses1(evoked2), responses2(evoked2));
title_string = sprintf( 'Full responses without spont z = %.2f, p = %d', stats.zval, p);
title(title_string)
xlabel('Firing Rate (Hz)')
ylabel('Number of cells')
xlim([0 max(responses1(evoked2))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subset evoked responses
WN1 = meanON(evoked1,:);
WN1L = meanONL(evoked1,:);
SP1 = SP(evoked1,:);
SP1L = SPL(evoked1,:);

depths1 = depths(evoked1);
fs1 = fs(evoked1); rs1 = Rs(evoked1);
rs1 = logical(rs1); fs1 = logical(fs1);

recs1= (recs(evoked1));
cells1 = cells(evoked1); 

% plot firing rate distributions
figure;
h1 = histogram(SP1L(:,2));
h1.BinWidth = .5;
hold on;
h2 = histogram(SP1(:,2));
h2.BinWidth = .5;
ylabel('Number of cells')
xlabel('Firing Rate (Hz)')
legend('laser on', 'laser off')
[p,h,stats] = signrank(SP1L(:,2), SP1(:,2));
r = stats.zval/sqrt(sum(~isnan(SP1L(:,2)))+sum(~isnan(SP1(:,2)))); %z / sqrt(N)
title_str = sprintf('Spontaneous Activity,  p = %d, r = %.4f', p, r)
title(title_str)

figure;
h1 = histogram(WN1L(:,2));
h1.BinWidth = 3;
hold on;
h2 = histogram(WN1(:,2));
h2.BinWidth = 3;
ylabel('Number of cells')
xlabel('Firing Rate (Hz)')
legend('laser on', 'laser off')
[p,h,stats] = signrank(WN1L(:,2), WN1(:,2))
r = stats.zval/sqrt(sum(~isnan(WN1L(:,2)))+sum(~isnan(WN1(:,2)))); %z / sqrt(N)
title_str = sprintf('Evoked Activity, p = %d, r = %.4f', p, r)
title(title_str)

figure; hold on;
plot(WN1(rs1,2), WN1L(rs1,2), 'ko')
plot(WN1(fs1,2), WN1L(fs1,2), 'go')
maxFR = max(max([WN1(:,2), WN1L(:,2)]));
plot(nanmean(WN1(:,2)), nanmean(WN1L(:,2)), 'ro','MarkerSize', 10)
plot(nanmedian(WN1(:,2)), nanmedian(WN1L(:,2)), 'ro','MarkerSize', 10)
plot([0 maxFR], [0 maxFR], 'r-')
title(title_str)
xlabel('FR laser off'); ylabel('FR laser on')
legend('Regular spiking', 'Narrow spiking', 'Mean', 'Median')
pbaspect([1 1 1]); set(gcf, 'PaperPositionMode', 'auto');

figure; hold on;
plot(SP1(rs1,2), SP1L(rs1,2), 'ko')
plot(SP1(fs1,2), SP1L(fs1,2), 'go')
plot(nanmean(SP1), nanmean(SP1L), 'ro', 'MarkerSize', 10)
plot(nanmedian(SP1), nanmedian(SP1L), 'mo', 'MarkerSize', 10)
maxFR = max(max([SP1(:,2), SP1L(:,2)]));
plot(nanmean(SP1(:,2)), nanmean(SP1L(:,2)), 'ro')
plot([0 maxFR], [0 maxFR], 'r-')
[p,h,stats] = signrank(SP1L(:,2),  SP1(:,2));
title_str = sprintf('Spont Activity, p = %d', p);
title(title_str)
legend('Regular spiking', 'Narrow spiking', 'Mean', 'Median')
xlabel('FR laser off'); ylabel('FR laser on')
pbaspect([1 1 1]); set(gcf, 'PaperPositionMode', 'auto');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot waveforms of cells
WFs1 = WFs(evoked1,:); %subset evoeked cells only that were used in the analysis
SNRs1 = SNRs(evoked1,:);
uQs1 = uQs(evoked1,:);
figure; hold on

x = 1:size(WFs1,2);
y = WFs1(rs1,:);
N = size(y,1);                                      % Number of ‘Experiments’ In Data Set
yMean = mean(y);                                    % Mean Of All Experiments At Each Value Of ‘x’
ySEM = std(y)/sqrt(N);                              % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
yCI95 = bsxfun(@times, ySEM, CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’

plot(x, yCI95+yMean, 'Color', [0.7 0.7 0.7]) 
plot(x, yMean, 'k', 'LineWidth', 4)                                      % Plot Mean WF Of RS cells

y = WFs1(fs1,:);
N = size(y,1);                                      
yMean = mean(y);                                    
ySEM = std(y)/sqrt(N);                              
CI95 = tinv([0.025 0.975], N-1);                    
yCI95 = bsxfun(@times, ySEM, CI95(:));              
           
plot(x, yCI95+yMean, 'Color', 'g') 
plot(x, yMean, 'g', 'LineWidth', 4)    % Plot Mean WF Of RS cells

xticks([0 30 60 90 120]);
xticklabels({'0', '1', '2', '3', '4'})
xlim([0 100])
xlabel('time (ms)')
ylabel('normalized spike amplitude')
 legend('','','mean RS waveform', '','','mean NS waveform')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test if VIP activation is the same across all recordings and across regular and fast spiking cells

VIP_MI_evoked = (WN1L(:,2) - WN1(:,2))./ (WN1L(:,2) + WN1(:,2));
VIP_MI_spont = (SP1L(:,2) - SP1(:,2))./ (SP1L(:,2) + SP1(:,2));

% STATS / ANOVA across recordings
x = VIP_MI_evoked;
[p,tbl1,stats] = kruskalwallis(x, recs1);
c = multcompare(stats);
title('VIP effect Evoked all cells')
[m, s] = grpstats(x,recs1,{'mean','sem'});

figure;
errorbar([1:length(unique(recs1))], m, s, 'ko-');
xlabel('Recording number')
ylabel('Mean MI /SEM')
title('VIP effect Evoked all cells')

x = VIP_MI_spont;
[p,tbl1,stats] = kruskalwallis(x, recs1);
c = multcompare(stats);
title('VIP effect Spont all cells')
[m, s] = grpstats(x,recs1,{'mean','sem'});

figure;
errorbar([1:length(unique(recs1))], m, s, 'ko-');
xlabel('Recording number')
ylabel('Mean MI /SEM')
title('VIP effect Spont all cells')

% split by cell type, RS
x = VIP_MI_evoked(rs1);
[p,tbl1,stats] = kruskalwallis(x, recs1(rs1));
c = multcompare(stats);
title('VIP effect Evoked RS cells')

x = VIP_MI_spont(rs1);
[p,tbl1,stats] = kruskalwallis(x, recs1(rs1));
c = multcompare(stats);
title('VIP effect Spont RS cells')

% split by cell type, FS
x = VIP_MI_evoked(fs1);
[p,tbl1,stats] = kruskalwallis(x, recs1(fs1));
c = multcompare(stats);
title('VIP effect Evoked FS cells')

x = VIP_MI_spont(fs1);
[p,tbl1,stats] = kruskalwallis(x, recs1(fs1));
c = multcompare(stats);
title('VIP effect Spont FS cells')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Effect of VIP activation on evoked and spontaneous activity. 
% 
LASER_EFFECT_EVOKED = []; %  all cells, for stats across cortical layers, evoked
LASER_EFFECT_SPONT = []; %  all cells, for stats across cortical layers, spont
layers1 = nan(sum(evoked1),1);

% old way to calculate laser effect, because it corr with FR we switched to
% modulation index
% LaserEffect_evoked = (WN1L(:,2) - WN1(:,2))./ max([WN1L(:,2)'; WN1(:,2)'])'; % Laser effect on evoked FR
% LaserEffect_spont = (SP1L(:,2) - SP1(:,2))./ max([SP1L(:,2)'; SP1(:,2)'])'; % Laser effect on spont FR

LaserEffect_evoked = (WN1L(:,2) - WN1(:,2))./(WN1L(:,2) + WN1(:,2)); % Laser effect on evoked FR
LaserEffect_spont = (SP1L(:,2) - SP1(:,2))./(SP1L(:,2) + SP1(:,2)); % Laser effect on spont FR


for cl = 1:length(CL)
    layer = CL{cl}; % layer depth limits
    indx = find(depths1 > layer(1) & depths1 <layer(2)); % indices of cells within this layer
    recs2 = recs1(indx); % recordings in this layer
    cells2 = cells1(indx); % cells in this layer
    fs2 = fs1(indx); fs2 = logical(fs2); % fast spiking cells in this layer
    rs2 = rs1(indx); rs2 = logical(rs2); % regular spiking cells in this layer
    layers1(indx) = 1*cl;
    
    % % evoked (WN response)
    meanLaserEffect_evoked(cl) = nanmean(LaserEffect_evoked(indx)); % mean laser effect
    meanLaserEffect_evoked_rs(cl) = nanmean(LaserEffect_evoked(indx(rs2)));
    meanLaserEffect_evoked_fs(cl) = nanmean(LaserEffect_evoked(indx(fs2)));
    
    semLaserEffect_evoked(cl) = sem(LaserEffect_evoked(indx)); % sem of laser effect
    semLaserEffect_evoked_rs(cl) = sem(LaserEffect_evoked(indx(rs2)));
    semLaserEffect_evoked_fs(cl) = sem(LaserEffect_evoked(indx(fs2)));

    % % spont activity
    meanLaserEffect_spont(cl) = nanmean(LaserEffect_spont(indx));
    meanLaserEffect_spont_rs(cl) = nanmean(LaserEffect_spont(indx(rs2)));
    meanLaserEffect_spont_fs(cl) = nanmean(LaserEffect_spont(indx(fs2)));
    
    semLaserEffect_spont(cl) = sem(LaserEffect_spont(indx));
    semLaserEffect_spont_rs(cl) = sem(LaserEffect_spont(indx(rs2)));
    semLaserEffect_spont_fs(cl) = sem(LaserEffect_spont(indx(fs2)));
    
    % collect modulation indices across layers for further stats [MI cl]
    LASER_EFFECT_EVOKED = [LASER_EFFECT_EVOKED; LaserEffect_evoked(indx) ones(length(indx),1)*cl];
    n_layer_evokedLaser(cl) = length(find(~isnan(LaserEffect_evoked(indx))==1));
    LASER_EFFECT_SPONT = [LASER_EFFECT_SPONT; LaserEffect_spont(indx) ones(length(indx),1)*cl];
    n_layer_spLaser(cl) = length(find(~isnan(LaserEffect_spont(indx))==1));
       
    n_rs_layer(cl) = sum(~isnan(LaserEffect_evoked(indx(rs2)))); % number of regular spiking cells in this layer
    n_fs_layer(cl) = sum(~isnan(LaserEffect_evoked(indx(fs2)))); % number of fast spiking cells in this layer   
end

% state by layer without layer
figure; hold on; subplot(2,1,1); hold on
errorbar([1:4], meanLaserEffect_evoked, semLaserEffect_evoked, 'ko-');
xlabel('Cortical Layer')
plot([0 5], [0 0], '--')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Laser Effect');
title_string = sprintf('WN On response, n = %d, %d, %d, %d',  n_layer_evokedLaser);
title(title_string)

subplot(2,1,2); hold on
errorbar([1:4], meanLaserEffect_spont, semLaserEffect_spont, 'ko-');
xlabel('Cortical Layer')
plot([0 5], [0 0], '--')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Laser Effect');
title_string = sprintf('Spont On response, n = %d, %d, %d, %d',  n_layer_spLaser);
title(title_string)

x = LASER_EFFECT_EVOKED(:,1); % effect
y = LASER_EFFECT_EVOKED(:,2); % layer
[p,tbl1,stats] = kruskalwallis(x, y);
title('Laser Effect by layer, Evoked')
c = multcompare(stats);
title('Laser Effect by layer, Evoked')

x = LASER_EFFECT_SPONT(:,1); % effect 
y = LASER_EFFECT_SPONT(:,2); % layer
[p,tbl1,stats] = kruskalwallis(x, y);
title('Laser Effect by layer, Spont')
c = multcompare(stats);
title('Laser Effect by layer, Spont')


x = LASER_EFFECT_SPONT(rs1,1); % effect 
y = LASER_EFFECT_SPONT(rs1,2); % layer
[p,tbl1,stats] = kruskalwallis(x, y);
title('Laser Effect by layer RS, Spont')
c = multcompare(stats);
title('Laser Effect by layer, Spont')

% % 
% PROBABILITY PLOTS
%
suppressed_evoked_indx = logical(LaserEffect_evoked < 0);
disinhibited_evoked_indx = logical(LaserEffect_evoked > 0);
suppressed_spont_indx = logical(LaserEffect_spont < 0);
disinhibited_spont_indx = logical(LaserEffect_spont > 0);

figure; 
subplot(2,1,1); hold on % evoked
[Ns, x] = hist(depths1(suppressed_evoked_indx), [0:50:800]);
[Nd, x] = hist(depths1(disinhibited_evoked_indx), [0:50:800]);
bar(x,[Ns' Nd'], 1,'grouped')
title('Distribution of VIP effects on evoked activity')
ylabel('Number of cells')
xlabel('Cortical Depth (\mum)')
legend('suppressed', 'disinihbited')

subplot(2,1,2); hold on % spont
[Ns, x] = hist(depths1(suppressed_spont_indx), [0:50:800]);
[Nd, x] = hist(depths1(disinhibited_spont_indx), [0:50:800]);
bar(x,[Ns' Nd'], 1,'grouped')
title('Distribution of VIP effects on spont activity')
ylabel('Number of cells')
xlabel('Cortical Depth (\mum)')
legend('suppressed', 'disinihbited')


figure; hold on
x = depths1(suppressed_evoked_indx);
x = x(~isnan(x));
h1 = histogram(x, [0:25:850], 'Normalization','probability');
y = depths1(disinhibited_evoked_indx);
y = y(~isnan(y));
h2 = histogram(y, [0:25:850],  'Normalization', 'probability');

x = depths1(suppressed_spont_indx);
x = x(~isnan(x));
h3 = histogram(x, [0:25:850], 'Normalization','probability');
y = depths1(disinhibited_spont_indx);
y = y(~isnan(y));
h4 = histogram(y, [0:25:850],  'Normalization', 'probability');


% % Probability Evoked and Spont 
figure; subplot(2,1,1); hold on
plot(h1.BinEdges(1:end-1), smooth(h1.Values,5), 'b')
plot(h2.BinEdges(1:end-1), smooth(h2.Values,5), 'r')
ylabel('Probability')
legend('suppressed', 'disinihbited')
title('Distribution of probability of VIP activation on evoked activity')
xlim([0 850])

subplot(2,1,2); hold on
plot(h3.BinEdges(1:end-1), smooth(h3.Values,5), 'b')
plot(h4.BinEdges(1:end-1), smooth(h4.Values,5), 'r')
ylabel('Probability')
xlabel('Cortical Depth (\mum)')
legend('suppressed', 'disinihbited')
title('Distribution of probability of VIP activation on spont activity')
xlim([0 850])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % look at the difference between RS and FS cells

figure; hold on
plot(LaserEffect_evoked(rs1), LaserEffect_spont(rs1), 'ko')
plot(LaserEffect_evoked(fs1), LaserEffect_spont(fs1), 'go')
lsline; xlabel('Laser Effect Evoked'); ylabel('Laser Effect Spont')
plot(nanmean(LaserEffect_evoked(rs1)), nanmean(LaserEffect_spont(rs1)), 'r*', 'MarkerSize', 15)
plot(nanmean(LaserEffect_evoked(fs1)), nanmean(LaserEffect_spont(fs1)), 'm*', 'MarkerSize', 15)
plot([0 0], [-1 1], 'k--')
plot([-1 1], [0 0], 'k--')
legend('RS', 'FS', '', '', 'mean RS', 'mean FS')
[r1, p1] = corr(LaserEffect_evoked(rs1), LaserEffect_spont(rs1), 'Type','Spearman','Rows', 'complete');
[r2, p2] = corr(LaserEffect_evoked(fs1), LaserEffect_spont(fs1), 'Type','Spearman','Rows', 'complete');
title_str = sprintf('rs rho = %.2f, p=%.4f; fs rho = %.2f, p = %.4f', r1, p1, r2, p2);
title(title_str); pbaspect([1 1 1]);

% RS and FS evoked and spont by layers
figure; % evoked
subplot(2,1,1); hold on
errorbar([0.9:3.9], meanLaserEffect_evoked_rs, semLaserEffect_evoked_rs, 'ko-');
errorbar([1.1:4.1], meanLaserEffect_evoked_fs, semLaserEffect_evoked_fs, 'k*-');
xlabel('Cortical Layer')
plot([0 5], [0 0], 'k--')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Laser Effect'); 
legend('RS', 'FS')
title_string = sprintf('Laser Effect, Evoked, RS n = %d, %d, %d, %d, n = %d, %d, %d, %d',  n_rs_layer, n_fs_layer);
title(title_string)
% spont
subplot(2,1,2); hold on
errorbar([0.9:3.9], meanLaserEffect_spont_rs, semLaserEffect_spont_rs, 'ko-');
errorbar([1.1:4.1], meanLaserEffect_spont_fs, semLaserEffect_spont_fs, 'k*-');
xlabel('Cortical Layer')
plot([0 5], [0 0], 'k--')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Laser Effect');
legend('RS', 'FS')
title_string = sprintf('Laser Effectt, Spont activity');
title(title_string)
set(gcf, 'PaperPositionMode', 'auto');


% plot FR vs Laser Effect
figure;
subplot(2,1,1); hold on
plot(WN1(rs1,2), LaserEffect_evoked(rs1), 'ko')
plot(WN1(fs1,2), LaserEffect_evoked(fs1), 'go')
plot([0 max(WN1(:,2))], [0 0 ], 'k--')
xlabel('Evoked FR'); ylabel('Laser Effect - evoked')
legend('rs', 'fs')
[r1, p1] = corr(abs(LaserEffect_evoked(rs1)), WN1(rs1,2), 'Type','Spearman','Rows', 'complete');
[r2, p2] = corr(abs(LaserEffect_evoked(fs1)), WN1(fs1,2), 'Type','Spearman','Rows', 'complete');
title(sprintf('Evoked, RS r = %.4f, p = %d, FS = %.4f, p = %d', r1, p1, r2, p2))
subplot(2,1,2); hold on
plot(SP1(rs1,2), LaserEffect_spont(rs1), 'ko')
plot(SP1(fs1,2), LaserEffect_spont(fs1), 'go')
plot([0 max(SP1(:,2))], [0 0 ], 'k--')
xlabel('Spont FR'); ylabel('Laser Effect - evoked')
legend('rs', 'fs')
[r1, p1] = corr(abs(LaserEffect_spont(rs1)), SP1(rs1,2), 'Type','Spearman','Rows', 'complete');
[r2, p2] = corr(abs(LaserEffect_spont(fs1)), SP1(fs1,2), 'Type','Spearman','Rows', 'complete');
title(sprintf('Spont, RS r = %.4f, p = %d, FS = %.4f, p = %d', r1, p1, r2, p2))

% STATS
% not sure what stats to do here yet

% Plot firing rates by cortical layer to identify what drives laser effect
H = []; FRs_evoked = []; FRs_spont =[]; FSindx = []; RSindx = [];
FRs_evoked_laser = []; FRs_spont_laser =[];
for cl = 1:4
    layer = CL{cl}; % layer depth limits
    indx = find(depths1 >layer(1) & depths1 < layer(2)); % indices of cells within this layer
    fs2 = fs1(indx); fs2 = logical(fs2); % fast spiking cells in this layer
    rs2 = rs1(indx); rs2 = logical(rs2); % regular spiking cells in this layer
    FSindx = [FSindx; fs2]; RSindx = [RSindx; rs2];
    
    if cl == 2 % save index of cells in layer 4. to look at them closer
        layer4_index = indx;
    end
    
    FR_evoked_means_laser_off(cl) = nanmean(WN1(indx,2));
    FR_evoked_means_laser_on(cl) = nanmean(WN1L(indx,2));
    FR_spont_means_laser_off(cl) = nanmean(SP1(indx,2));
    FR_spont_means_laser_on(cl) = nanmean(SP1L(indx,2));
    
    FRs_evoked = [FRs_evoked; WN1(indx,2) ones(length(WN1(indx,2)),1)*cl];
    FRs_spont = [FRs_spont; SP1(indx,2) ones(length(SP1(indx,2)),1)*cl];
    FRs_evoked_laser = [FRs_evoked_laser; WN1L(indx,2) ones(length(WN1L(indx,2)),1)*cl];
    FRs_spont_laser = [FRs_spont_laser; SP1L(indx,2) ones(length(SP1L(indx,2)),1)*cl];
    
    FR_evoked_medians_laser_off(cl) = nanmedian(WN1(indx,2));
    FR_evoked_medians_laser_on(cl) = nanmedian(WN1L(indx,2));
    FR_spont_medians_laser_off(cl) = nanmedian(SP1(indx,2));
    FR_spont_medians_laser_on(cl) = nanmedian(SP1L(indx,2));
    
    FR_evoked_means_laser_off_rs(cl) = nanmean(WN1(indx(rs2),2));
    FR_evoked_means_laser_off_fs(cl) = nanmean(WN1(indx(fs2),2));
    FR_evoked_means_laser_on_rs(cl) = nanmean(WN1L(indx(rs2),2));
    FR_evoked_means_laser_on_fs(cl) = nanmean(WN1L(indx(fs2),2));
    
    FR_spont_means_laser_off_rs(cl) = nanmean(SP1(indx(rs2),2));
    FR_spont_means_laser_on_rs(cl) = nanmean(SP1L(indx(rs2),2));
    FR_spont_means_laser_off_fs(cl) = nanmean(SP1(indx(fs2),2));
    FR_spont_means_laser_on_fs(cl) = nanmean(SP1L(indx(fs2),2));
    
    FR_evoked_sems_laser_off(cl) = sem(WN1(indx,2));
    FR_evoked_sems_laser_on(cl) = sem(WN1L(indx,2));
    FR_spont_sems_laser_off(cl) = sem(SP1(indx,2));
    FR_spont_sems_laser_on(cl) = sem(SP1L(indx,2));
    
    FR_evoked_sems_laser_off_rs(cl) = sem(WN1(indx(rs2),2));
    FR_evoked_sems_laser_on_rs(cl) = sem(WN1L(indx(rs2),2));
    FR_spont_sems_laser_off_rs(cl) = sem(SP1(indx(rs2),2));
    FR_spont_sems_laser_on_rs(cl) = sem(SP1L(indx(rs2),2));
    FR_evoked_sems_laser_off_fs(cl) = sem(WN1(indx(fs2),2));
    FR_evoked_sems_laser_on_fs(cl) = sem(WN1L(indx(fs2),2));
    FR_spont_sems_laser_off_fs(cl) = sem(SP1(indx(fs2),2));
    FR_spont_sems_laser_on_fs(cl) = sem(SP1L(indx(fs2),2));
    
    % stats
    [p,h,STATS] = signrank(WN1(indx,2), WN1L(indx,2));
    H(1,cl) = p;
    
    [p,h,STATS] = signrank(SP1(indx,2), SP1L(indx,2));
    H(2,cl) = p;
end

% Overall
[p,tbl1,stats] = kruskalwallis(FRs_evoked(:,1), FRs_evoked(:,2));
c = multcompare(stats);
title('Evoked FR by layer')

[p,tbl1,stats] = kruskalwallis(FRs_spont(:,1), FRs_spont(:,2));
c = multcompare(stats);
title('Spont FR by layer')

RSindx = logical(RSindx);
% Regular spiking cells
[p,tbl1,stats] = kruskalwallis(FRs_evoked(RSindx,1), FRs_evoked(RSindx,2));
c = multcompare(stats);
title('Evoked FR RS by layer')

[p,tbl1,stats] = kruskalwallis(FRs_spont(RSindx,1), FRs_spont(RSindx,2));
c = multcompare(stats);
title('Spont FR RS by layer')

FSindx = logical(FSindx);
% Narrow spiking cells
[p,tbl1,stats] = kruskalwallis(FRs_evoked(FSindx,1), FRs_evoked(FSindx,2));
c = multcompare(stats);
title('Evoked FR FS by layer')

[p,tbl1,stats] = kruskalwallis(FRs_spont(FSindx,1), FRs_spont(FSindx,2));
c = multcompare(stats);
title('Spont FR FS by layer')

% % FR change across layers % %
% Regular spiking cells
[p,tbl1,stats] = kruskalwallis(FRs_evoked_laser(RSindx,1) - FRs_evoked(RSindx,1), FRs_evoked(RSindx,2));
c = multcompare(stats);
title('Evoked FR change RS by layer')

[p,tbl1,stats] = kruskalwallis(FRs_spont_laser(RSindx,1) - FRs_spont(RSindx,1), FRs_spont(RSindx,2));
c = multcompare(stats);
title('Spont FR change RS by layer')

% Narrow spiking cells
[p,tbl1,stats] = kruskalwallis(FRs_evoked_laser(FSindx,1) - FRs_evoked(FSindx,1), FRs_evoked(FSindx,2));
c = multcompare(stats);
title('Evoked FR change FS by layer')

[p,tbl1,stats] = kruskalwallis(FRs_spont_laser(FSindx,1) - FRs_spont(FSindx,1), FRs_spont(FSindx,2));
c = multcompare(stats);
title('Spont FR change FS by layer')

% look at layer 4 fs spiking cells closer
fs_layer4_index =  layer4_index(fs1(layer4_index));
recs_layer4_index = recs1(fs_layer4_index);
cells_layer4_index = cells1(fs_layer4_index);

figure; hold on
for i = 1:length(recs_layer4_index)
    plot(WFs1(fs_layer4_index(i),:)')
    title(sprintf('uQ = %.4f, SNR = %.4f', uQs1(fs_layer4_index(i)), SNRs1(fs_layer4_index(i))))
end

% plot FR across layers
figure; subplot(2,1,1); hold on
errorbar([1:4], FR_evoked_means_laser_off, FR_evoked_sems_laser_off, 'ko-')
plot([1:4], FR_evoked_medians_laser_off, 'k*-')
errorbar([1.2:4.2], FR_evoked_means_laser_on, FR_evoked_sems_laser_on, 'co-')
plot([1.1:4.1], FR_evoked_medians_laser_on, 'c*-')
plot([1.1:4.1], max(FR_evoked_means_laser_on).*(H(1,:) <0.0125), '*r')
legend('mean laser off', 'median laser off', 'mean laser on', 'median laser on')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Mean (median) FR /SEM')
title('Evoked')
subplot(2,1,2); hold on;
errorbar([1:4], FR_spont_means_laser_off, FR_spont_sems_laser_off, 'ko-')
plot([1:4], FR_spont_medians_laser_off, 'k*-')
errorbar([1.2:4.2], FR_spont_means_laser_on, FR_spont_sems_laser_on, 'co-')
plot([1.2:4.2], FR_spont_medians_laser_on, 'c*-')
plot([1.1:4.1], max(FR_spont_means_laser_on).*(H(2,:) <0.0125), '*r')
legend('mean laser off', 'median laser off', 'mean laser on', 'median laser on')
xticks([1:4]); xlim([0 5])
title('Spont')
xticklabels({'2/3', '4', '5', '6'});
ylabel('Mean (median)/SEM FR')

% plot changes in FR in RS and FS by layer
figure; subplot(2,1,1); hold on
errorbar([1:4], FR_evoked_means_laser_off_rs, FR_evoked_sems_laser_off_rs, 'ko-')
errorbar([1.1:4.1], FR_evoked_means_laser_off_fs, FR_evoked_sems_laser_off_fs, 'k*-')
errorbar([1.2:4.2], FR_evoked_means_laser_on_rs, FR_evoked_sems_laser_on_rs, 'co-')
errorbar([1.3:4.3], FR_evoked_means_laser_on_fs, FR_evoked_sems_laser_on_fs, 'c*-')
xticks([1:4]); xlim([0 5])
legend('rs laser off', 'fs laser off', 'rs laser on', 'fs laser on')
xticklabels({'2/3', '4', '5', '6'});
ylabel('Mean FR /SEM')
title('Evoked')
subplot(2,1,2); hold on;
errorbar([1:4], FR_spont_means_laser_off_rs, FR_spont_sems_laser_off_rs, 'ko-')
errorbar([1.1:4.1], FR_spont_means_laser_off_fs, FR_spont_sems_laser_off_fs, 'k*-')
errorbar([1.2:4.2], FR_spont_means_laser_on_rs, FR_spont_sems_laser_on_rs, 'co-')
errorbar([1.3:4.3], FR_spont_means_laser_on_fs, FR_spont_sems_laser_on_fs, 'c*-')
legend('rs laser off', 'fs laser off', 'rs laser on', 'fs laser on')
xticks([1:4]); xlim([0 5])
title('Spont')
xticklabels({'2/3', '4', '5', '6'});
ylabel('Mean /SEM FR')

%% Sound Modulation Index, Figure 3 plots
% compute sound modulation index
MI_sound_run = (WN1(:,1) - SP1(:,1))./ (WN1(:,1) +SP1(:,1)); % running laser off
MI_sound_sit = (WN1(:,2) - SP1(:,2))./ (WN1(:,2) +SP1(:,2)); % sitting laser off (control)
MI_sound_run_laser = (WN1L(:,1) - SP1L(:,1))./ (WN1L(:,1) +SP1L(:,1)); % running laser on
MI_sound_sit_laser = (WN1L(:,2) - SP1L(:,2))./ (WN1L(:,2) +SP1L(:,2)); % running laser on
MI_sound_predicted = MI_sound_sit_laser + MI_sound_run; % linear sum of running and laser effects

% plot distributions laser off vs laser on
figure; hold on;
[h3,x] = hist(MI_sound_sit, [-1:.1:1]);
h = smooth(h3,3);
plot(x,h, 'k-');
[h4,x] = hist(MI_sound_sit_laser, [-1:.1:1]);
h = smooth(h4,3);
plot(x,h, 'c-');
legend( 'sitting - laser off', 'sitting - laser on')
xlabel('Modulation Index - Sound')
ylabel('Number of cells')
[p,h stats] = signrank(MI_sound_sit , MI_sound_sit_laser);
title(sprintf(' Sound MI, laser off/on (smoothed), signrank = %.4f, p = %.4f', stats.zval, p))

% print results of regression
[B,BINT,R,RINT,STATS] = regress(MI_sound_predicted, [MI_sound_run_laser ones(length(MI_sound_run_laser),1)]);
fprintf('\n %%%%%%%%%%%%%%%%%%%%%')
fprintf('\n Regressing Running and VIP effect on Predicted Sound MI values')
fprintf('\n B = %.4f, Rsqrd = %.4f, F = %.4f p = %d', B, STATS(1), STATS(2), STATS(3))
fprintf('\n %%%%%%%%%%%%%%%%%%%%%')

% plot individual cells sound MI, RS and FS separately
figure; hold on
plot(MI_sound_sit(rs1) ,MI_sound_sit_laser(rs1), 'ko')
plot(MI_sound_sit(fs1) ,MI_sound_sit_laser(fs1), 'go')
xlabel('sound MI sitting lasr off'); ylabel('sound MI sitting laser on')
[r1, p1] = corr(MI_sound_run(rs1),MI_sound_sit(rs1), 'Type','Spearman','Rows', 'complete');
[r2, p2] = corr(MI_sound_run(fs1),MI_sound_sit(fs1), 'Type','Spearman','Rows', 'complete');
title_str = sprintf('RS r = %4.f, p = %d, FS r = %.4f, p = %d', r1, p1, r2, p2);
title(title_str)
plot([-1 1], [-1 1], 'r-'); plot([-1 1], [0 0], 'k--'); plot([0 0], [-1 1], 'k--')
pbaspect([1 1 1]);  set(gcf, 'PaperPositionMode', 'auto');

%look at the effect of laser on sound modulation by layer
soundMI_laserOFF = []; soundMI_laserON = [];
for cl = 1:length(CL)
    layer = CL{cl}; % layer depth limits
    indx = find(depths1 >layer(1) & depths1 <layer(2)); % indices of cells within this layer
    
    meanMI_sound_sit_laserOFF(cl) = nanmean(MI_sound_sit(indx)); % mean sound MI sitting in this layer
    semMI_sound_sit_laserOFF(cl) = sem(MI_sound_sit(indx)); % SEM sound MI sitting
    n_layer_laserOFF(cl) = sum(~isnan(MI_sound_sit(indx))); % number of non NaN sound MI in this layer
    soundMI_laserOFF = [soundMI_laserOFF; MI_sound_sit(indx) ones(length(MI_sound_sit(indx)),1)*cl]; % collect individual sound MI with layer assignment
    
    meanMI_sound_sit_laserON(cl) = nanmean(MI_sound_sit_laser(indx));
    semMI_sound_sit_laserON(cl) = sem(MI_sound_sit_laser(indx));
    n_layer_laserON(cl) = sum(~isnan(MI_sound_sit_laser(indx)));
    soundMI_laserON = [soundMI_laserON; MI_sound_sit_laser(indx) ones(length(MI_sound_sit_laser(indx)),1)*cl]; % collect individual sound MI with layer assignment
end

% ANOVA to check if sound MI is different across the layers
[p,tbl1,stats] = kruskalwallis(soundMI_laserOFF(:,1)', soundMI_laserOFF(:,2)');
c = multcompare(stats);
title('Sound modulation index by layer, laser off')
[p,tbl1,stats] = kruskalwallis(soundMI_laserON(:,1)', soundMI_laserON(:,2)');
c = multcompare(stats);
title('Sound modulation index by layer, laser on')
[p,tbl1,stats] = kruskalwallis(soundMI_laserON(:,1)-soundMI_laserOFF(:,1), soundMI_laserOFF(:,2));
c = multcompare(stats);
title('Sound MI diff (laser on - laser off) , VIP')

% plot FR vs Sound MI
figure;
subplot(2,1,1); hold on
plot(WN1(rs1,2), MI_sound_sit(rs1), 'ko')
plot(WN1(fs1,2), MI_sound_sit(fs1), 'go')
plot([0 max(WN1(:,2))], [0 0 ], 'k--')
xlabel('Evoked FR'); ylabel('Sound MI sit + laser off')
legend('rs', 'fs')
[r1, p1] = corr(abs(MI_sound_sit(rs1)), WN1(rs1,2), 'Type','Spearman','Rows', 'complete');
[r2, p2] = corr(abs(MI_sound_sit(fs1)), WN1(fs1,2), 'Type','Spearman','Rows', 'complete');
title(sprintf('Evoked, RS r = %.4f, p = %d, FS = %.4f, p = %d', r1, p1, r2, p2))
subplot(2,1,2); hold on
plot(SP1(rs1,2), MI_sound_sit(rs1), 'ko')
plot(SP1(fs1,2), MI_sound_sit(fs1), 'go')
plot([0 max(SP1(:,2))], [0 0 ], 'k--')
xlabel('Spont FR'); ylabel('Sound MI sit + laser off')
legend('rs', 'fs')
[r1, p1] = corr(abs(MI_sound_sit(rs1)), SP1(rs1,2), 'Type','Spearman','Rows', 'complete');
[r2, p2] = corr(abs(MI_sound_sit(fs1)), SP1(fs1,2), 'Type','Spearman','Rows', 'complete');
title(sprintf('Spont, RS r = %.4f, p = %d, FS = %.4f, p = %d', r1, p1, r2, p2))

figure; hold on
errorbar([1.1:4.1], meanMI_sound_sit_laserOFF, semMI_sound_sit_laserOFF, 'ko-');
errorbar([1:4], meanMI_sound_sit_laserON, semMI_sound_sit_laserON, 'co-');
xlabel('Cortical Layer')
plot([0 5], [0 0], 'k--')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Modulation Index (mean/SEM)- sound effect');
legend('laser off', 'laser on')
fprintf('laer effect')
[p, h, stats] = signrank(soundMI_laserON(:,1), soundMI_laserOFF(:,1));
title_string = sprintf( 'On response MI, n = %d, %d, %d, %d, laser effect, p = %d',  n_layer_laserOFF, p);
title(title_string)
set(gcf, 'PaperPositionMode', 'auto');

figure; hold on;
plot(MI_sound_sit(rs1), MI_sound_sit_laser(rs1), 'ko')
plot(MI_sound_sit(fs1), MI_sound_sit_laser(fs1), 'go')
plot(nanmean(MI_sound_sit), nanmean(MI_sound_sit_laser), 'ro', 'MarkerSize', 10)
plot(nanmedian(MI_sound_sit), nanmedian(MI_sound_sit_laser), 'mo', 'MarkerSize', 10)
maxFR = max(max([MI_sound_sit, MI_sound_sit_laser]));
plot([-1 1], [-1 1], 'r-')
plot([-1 1], [0 0], 'k--')
plot([0 0], [-1 1], 'k--')
[p,h,stats] = signrank(MI_sound_sit_laser, MI_sound_sit)
title_str = sprintf('Evoked Activity, p = %d', p);
title(title_str)
xlabel('sound modulation index - sit laser off');
ylabel('sound modulation index - sit laser on')
legend('Regular spiking', 'Narrow spiking', 'Mean', 'Median')
pbaspect([1 1 1]); set(gcf, 'PaperPositionMode', 'auto');

meanMI_sound_laser = nanmean(soundMI_laserON(:,1));
semMI_sound_laser = sem(soundMI_laserON(:,1));

meanMI_sound = nanmean(soundMI_laserOFF(:,1));
semMI_sound = sem(soundMI_laserOFF(:,1));

figure;hold on;
bar([1 2], [meanMI_sound_laser meanMI_sound])
errorbar(1, meanMI_sound_laser, semMI_sound_laser , 'k-' );
errorbar(2, meanMI_sound, semMI_sound , 'k-' )
xlim([0 3]); ylabel('Modulation Index - sound')
xticks([1 2]); xticklabels({'laser on', 'laser off'})
title('effect off sound in laser on and off conditions')

%% Distance Correlation analysis to test if Laser Effect is dependent on FR

% Test with laser effect first
% calculate distance corr for RS and FS separately
DistCorr_evoked_rs = distcorr(LaserEffect_evoked(rs1), WN1(rs1,2));
DistCorr_evoked_fs = distcorr(LaserEffect_evoked(fs1), WN1(fs1,2));
DistCorr_spont_rs = distcorr(LaserEffect_spont(rs1), SP1(rs1,2));
DistCorr_spont_fs = distcorr(LaserEffect_spont(fs1), SP1(fs1,2));

total = [DistCorr_evoked_rs DistCorr_evoked_fs DistCorr_spont_rs DistCorr_spont_fs];
% bootstrap shuffle to assess noise
n = 500;
for i = 1:n
    DistCorr_evoked_rs_noise(i) = distcorr(LaserEffect_evoked(rs1), Shuffle(WN1(rs1,2)));
    DistCorr_evoked_fs_noise(i) = distcorr(LaserEffect_evoked(fs1), Shuffle(WN1(fs1,2)));
    DistCorr_spont_rs_noise(i) = distcorr(LaserEffect_spont(rs1), Shuffle(SP1(rs1,2)));
    DistCorr_spont_fs_noise(i) = distcorr(LaserEffect_spont(fs1), Shuffle(SP1(fs1,2)));
end

noise = [median(DistCorr_evoked_rs_noise) median(DistCorr_evoked_fs_noise) median(DistCorr_spont_rs_noise) median(DistCorr_spont_fs_noise)];

figure;
bar(1:4, total - noise)
xticklabels({'Evoked RS', 'Evoked FS', 'Spont RS', 'Spont FS'})
ylabel('distance corr')
title('Dependence of Laser Effect on FR')

% Test with sound Modulation Index
% calculate distance corr for RS and FS separately

run_diff =  MI_sound_run - MI_sound_sit;
laser_diff = MI_sound_sit_laser - MI_sound_sit;
predicted_diff = run_diff + laser_diff;
actual_diff =  MI_sound_run_laser - MI_sound_sit;

DistCorr_evoked_rs = distcorr(laser_diff(rs1), WN1(rs1,2));
DistCorr_evoked_fs = distcorr(laser_diff(fs1), WN1(fs1,2));
DistCorr_spont_rs = distcorr(laser_diff(rs1), SP1(rs1,2));
DistCorr_spont_fs = distcorr(laser_diff(fs1), SP1(fs1,2));

total = [DistCorr_evoked_rs DistCorr_evoked_fs DistCorr_spont_rs DistCorr_spont_fs];
% bootstrap shuffle to assess noise
n = 500;
for i = 1:n
    DistCorr_evoked_rs_noise(i) = distcorr(laser_diff(rs1), Shuffle(WN1(rs1,2)));
    DistCorr_evoked_fs_noise(i) = distcorr(laser_diff(fs1), Shuffle(WN1(fs1,2)));
    DistCorr_spont_rs_noise(i) = distcorr(laser_diff(rs1), Shuffle(SP1(rs1,2)));
    DistCorr_spont_fs_noise(i) = distcorr(laser_diff(fs1), Shuffle(SP1(fs1,2)));
end

noise = [median(DistCorr_evoked_rs_noise) median(DistCorr_evoked_fs_noise) median(DistCorr_spont_rs_noise) median(DistCorr_spont_fs_noise)];

figure;
bar(1:4, total - noise)
xticklabels({'Evoked RS', 'Evoked FS', 'Spont RS', 'Spont FS'})
ylabel('distance corr')
title('Dependence of Laser Effect on Sound MI with FR')

% It appears that changes in sound MI are less dependent on the FR of
% neurons and since different layers have different mean FR, it is a better
% measure to check for layer specific effects
LaserDiffs = [];
for cl = 1:4
    layer = CL{cl}; % layer depth limits
    indx = find(depths1 > layer(1) & depths1 < layer(2)); % indices of cells within this layer
     fs2 = fs1(indx); fs2 = logical(fs2); % fast spiking cells in this layer
    rs2 = rs1(indx); rs2 = logical(rs2); % regular spiking cells in this layer
    
    mean_diff_rs(cl) = nanmean(laser_diff(rs2));
    mean_diff_fs(cl) = nanmean(laser_diff(fs2));
    
    sem_diff_rs(cl) = sem(laser_diff(rs2));
    sem_diff_fs(cl) = sem(laser_diff(fs2));
    
    LaserDiffs = [ LaserDiffs; laser_diff(indx), ones(length(indx),1)*cl ];
end
figure; hold on
errorbar([1:4], mean_diff_rs, sem_diff_rs, 'ko-')
errorbar([1:4], mean_diff_fs, sem_diff_fs, 'k*-')
xlabel('Cortical Layer')
plot([0 5], [0 0], 'k--')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('difference in sound MI (laser on - laser off)')

[p,tbl1,stats] = kruskalwallis(LaserDiffs(:,1), LaserDiffs(:,2));
c = multcompare(stats);
title('Difference in sound MI by layer')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Are any effects correlated with cluster quality?
figure; hold on
plot(MI_sound_sit, uQs1, 'ko')
[rho, p]=corr(MI_sound_sit, uQs1, 'Type', 'Spearman', 'rows', 'pairwise');
title_str = sprintf('rs rho = %.2f, p=%.4f', rho, p);
title(title_str)
ylabel('uQ'); xlabel('Sound MI')

figure; hold on
plot(MI_sound_sit_laser - MI_sound_sit, uQs1, 'ko')
[rho, p]=corr(MI_sound_sit_laser - MI_sound_sit, uQs1, 'Type', 'Spearman', 'rows', 'pairwise');
title_str = sprintf('rs rho = %.2f, p=%.4f', rho, p);
title(title_str)
ylabel('uQ'); xlabel('Sound MI laser - Sound MI')

figure; hold on
plot(MI_sound_sit_laser - MI_sound_sit, SNRs1, 'ko')
[rho, p]=corr(MI_sound_sit_laser - MI_sound_sit, SNRs1, 'Type', 'Spearman', 'rows', 'pairwise');
title_str = sprintf('rs rho = %.2f, p=%.4f', rho, p);
title(title_str)
ylabel('SNR'); xlabel('Sound MI laser - Sound MI')


%% Linear interaction test with running MI and VIP MI
% This is not going to be very accurate because we do not control number of
% running trials. this is just for curiosity For more accurate analysis see
% WNModulationIndex_epistatic file
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

figure; hold
plot(run_diff(rs1), laser_diff(rs1), 'o', 'Color', [.8 .80 .80], 'MarkerSize', 8)
plot(run_diff(fs1), laser_diff(fs1), 'o', 'Color', [0 .255 0], 'MarkerSize', 8)
xlabel('Running Effect (MI diff)')
ylabel('Laser Effect (MI diff)')
plot([-2 1], [0 0], 'k--')
plot([0 0], [-2 1], 'k--')
[r, p] = corr(run_diff,laser_diff, 'Type','Spearman','Rows', 'complete')
title_string = sprintf( 'rho = %.4f, p = %.4f', r,p);
title(title_string)
pbaspect([1 1 1]);  set(gcf, 'PaperPositionMode', 'auto');

figure; hold on
plot(actual_diff(rs1), predicted_diff(rs1), 'ko')
plot(actual_diff(fs1), predicted_diff(fs1), 'ro')
xlabel('Combined Effect (running laser on trials)')
ylabel('Predicted combined Effect (running effect + laser effect)')
plot([-2 2], [0 0], 'k--')
plot([0 0], [-2 2], 'k--')
[r, p] = corr(actual_diff,predicted_diff, 'Type','Spearman','Rows', 'complete')
title_string = sprintf( 'rho = %.4f, p = %.4f', r,p);
title(title_string)
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

%%  print
fprintf('\nN cells Sound MI sitting laser on =%d, N laser off = %d \n',  sum(~isnan(MI_sound_sit_laser)),  sum(~isnan(MI_sound_sit)))
fprintf('\nN cells WN laser on =%d, N sitting laser off = %d \n',  sum(~isnan(WN1L(:,2))),  sum(~isnan(WN1(:,2))))
fprintf('\nN cells Spont sitting laser on =%d, N sitting laser off = %d \n',  sum(~isnan(SP1L(:,2))),  sum(~isnan(SP1(:,2))))
fprintf('\nSound MI laser off by layer =%.2f +/- %.2f, %.2f +/- %.2f, %.2f +/- %.2f, %.2f +/- %.2f,  \n',  nanmean(MI_sound_sit(layers1==1)), sem(MI_sound_sit(layers1==1)), ...
    nanmean(MI_sound_sit(layers1==2)), sem(MI_sound_sit(layers1==2)), nanmean(MI_sound_sit(layers1==3)), sem(MI_sound_sit(layers1==3)), nanmean(MI_sound_sit(layers1==4)), sem(MI_sound_sit(layers1==4)));

fprintf('\nSound MI laser on by layer =%.2f +/- %.2f, %.2f +/- %.2f, %.2f +/- %.2f, %.2f +/- %.2f,  \n',  nanmean(MI_sound_sit_laser(layers1==1)), sem(MI_sound_sit_laser(layers1==1)), ...
    nanmean(MI_sound_sit_laser(layers1==2)), sem(MI_sound_sit_laser(layers1==2)), nanmean(MI_sound_sit_laser(layers1==3)), sem(MI_sound_sit_laser(layers1==3)), nanmean(MI_sound_sit_laser(layers1==4)), sem(MI_sound_sit_laser(layers1==4)));

y = [MI_sound_sit_laser; MI_sound_sit];
g1 = [ones(length(MI_sound_sit_laser),1); ones(length(MI_sound_sit_laser),1) * 2]; % laser on / off grouping varibale
g2 = [layers1; layers1]; %cortical layers grouping variable
[p,tbl,stats] = anovan(y,{g1,g2}, 'model', 2, 'varnames', {'laser', 'layer'});

