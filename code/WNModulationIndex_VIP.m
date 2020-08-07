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

cdVIP; load('Silence_DistanceCorr_dirs.mat'); load('WNdirsVIP.mat')


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
                cd(WNdirs{recs(cc)})
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
responses =  meanON - SP;
figure; hist(responses(evoked1,:), 100)
legend('running', 'sitting')
[p, h, stats]= signrank(responses(evoked1,1), responses(evoked1,2));
title_string = sprintf( 'On responses without spont z = %.2f, p = %d', stats.zval, p);
title(title_string)
xlabel('Firing Rate (Hz)')
ylabel('Number of cells')
xlim([0 max(responses(evoked1))])

% figure 2 - Full WN responses
evoked2 = logical(evoked(:,4) & zstats(:,4)>0);

responses =  meanWN - SP;
figure; hist(responses(evoked2,:), 100)
legend('running', 'sitting')
[p, h, stats]= signrank(responses(evoked2,1), responses(evoked2,2));
title_string = sprintf( 'Full responses without spont z = %.2f, p = %d', stats.zval, p);
title(title_string)
xlabel('Firing Rate (Hz)')
ylabel('Number of cells')
xlim([0 max(responses(evoked2))])

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

LASER_EFFECT_EVOKED = []; %  all cells, for stats across cortical layers, evoked
LASER_EFFECT_SPONT = []; %  all cells, for stats across cortical layers, spont
layers1 = nan(sum(evoked1),1);

LaserEffect_evoked = (WN1L(:,2) - WN1(:,2))./ max([WN1L(:,2)'; WN1(:,2)'])'; % Laser effect on evoked FR
LaserEffect_spont = (SP1L(:,2) - SP1(:,2))./ max([SP1L(:,2)'; SP1(:,2)'])'; % Laser effect on spont FR

for cl = 1:length(CL)
    layer = CL{cl}; % layer depth limits
    indx = find(depths1 >layer(1) & depths1 <layer(2)); % indices of cells within this layer
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

x = LASER_EFFECT_EVOKED(:,1);
y = LASER_EFFECT_EVOKED(:,2);
[p,tbl1,stats] = kruskalwallis(x, y);
title('Laser Effect by layer, Evoked')
c = multcompare(stats);
title('Laser Effect by layer, Evoked')

x = LASER_EFFECT_SPONT(:,1);
y = LASER_EFFECT_SPONT(:,2);
[p,tbl1,stats] = kruskalwallis(x, y);
title('Laser Effect by layer, Spont')
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



figure; hold on
plot(LaserEffect_evoked(rs1), LaserEffect_spont(rs1), 'ko')
plot(LaserEffect_evoked(fs1), LaserEffect_spont(fs1), 'go')
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
errorbar([1:4], meanLaserEffect_evoked, semLaserEffect_evoked, 'co-');
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
errorbar([1:4], meanLaserEffect_spont, semLaserEffect_spont, 'co-');
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

%%
% sound modulation index
MI_sound_run = (WN1(:,1) - SP1(:,1))./ (WN1(:,1) +SP1(:,1));
MI_sound_sit = (WN1(:,2) - SP1(:,2))./ (WN1(:,2) +SP1(:,2));
MI_sound_run_laser = (WN1L(:,1) - SP1L(:,1))./ (WN1L(:,1) +SP1L(:,1));
MI_sound_sit_laser = (WN1L(:,2) - SP1L(:,2))./ (WN1L(:,2) +SP1L(:,2));
MI_sound_predicted = MI_sound_sit_laser + MI_sound_run;

[B,BINT,R,RINT,STATS] = regress(MI_sound_run_laser,MI_sound_predicted)


% indx = find(MI_sound_sit<0);

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
title(' On responses, laser on')

fprintf('running vs sitting, laser off')
[p,h stats] =ranksum(MI_sound_run, MI_sound_sit)

fprintf('running vs sitting, laser on')
[p,h stats] =ranksum(MI_sound_run_laser, MI_sound_sit_laser)

mi = []; miL=[]; layers = []; layersL=[];
for cl = 1:length(CL)
    layer = CL{cl}; % layer depth limits
    indx = find(depths1 >layer(1) & depths1 <layer(2)); % indices of cells within this layer
    meanMI_sound_run(cl) = nanmean(MI_sound_run(indx));
    meanMI_sound_sit(cl) = nanmean(MI_sound_sit(indx));
    semMI_sound_run(cl) = sem(MI_sound_run(indx));
    semMI_sound_sit(cl) = sem(MI_sound_sit(indx));
    n_layer(cl) = sum(~isnan(MI_sound_run(indx)));
    layers = [layers; ones(length(MI_sound_sit(indx)),1)*cl];
    mi = [mi; MI_sound_sit(indx)];
    
    meanMI_sound_runL(cl) = nanmean(MI_sound_run_laser(indx));
    meanMI_sound_sitL(cl) = nanmean(MI_sound_sit_laser(indx));
    semMI_sound_runL(cl) = sem(MI_sound_run_laser(indx));
    semMI_sound_sitL(cl) = sem(MI_sound_sit_laser(indx));
    n_layerL(cl) = sum(~isnan(MI_sound_run_laser(indx)));
    layersL=[layersL; ones(length(MI_sound_run_laser(indx)),1)*cl];
    miL = [miL; MI_sound_sit_laser(indx)];
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

[p,tbl1,stats] = kruskalwallis(mi', layers');
c = multcompare(stats);
title('Sound modulation index by layer, laser off')
[p,tbl1,stats] = kruskalwallis(miL', layersL');
c = multcompare(stats);
title('Sound modulation index by layer, laser on')
[p,tbl1,stats] = kruskalwallis(miL-mi, layers');
c = multcompare(stats);
title('Sound modulation index change, VIP')


figure; hold on
errorbar([1:4], meanMI_sound_sitL, semMI_sound_sitL, 'co-');
errorbar([1:4], meanMI_sound_runL, semMI_sound_runL, 'co--');
xlabel('Cortical Layer')
plot([0 5], [0 0], 'k--')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Sound Modulation Index (mean/SEM)');
title_string = sprintf( 'On response MI laser on trials, n = %d, %d, %d, %d',  n_layer);
title(title_string)
legend('sitting','running')

fs1 = logical(fs1); rs1 = logical(rs1);
MI1_rs = (WN1(rs1,1) - WN1(rs1,2))./ (WN1(rs1,1) +WN1(rs1,2));
MI1_fs = (WN1(fs1,1) - WN1(fs1,2))./ (WN1(fs1,1) +WN1(fs1,2));

SP1_rs = (SP1(rs1,1) - SP1(rs1,2))./ (SP1(rs1,1) +SP1(rs1,2));
SP1_fs = (SP1(fs1,1) - SP1(fs1,2))./ (SP1(fs1,1) +SP1(fs1,2));

figure(106); hold on
bar([0.8 1.8], [nanmean(SP1_rs) nanmean(MI1_rs)], 'BarWidth', .2)
bar([1  2 ], [ nanmean(SP1_fs) nanmean(MI1_fs) ], 'BarWidth', .2)

errorbar( [0.8 1 1.8 2],[nanmean(SP1_rs) nanmean(SP1_fs) nanmean(MI1_rs)  nanmean(MI1_fs)], ...
    [sem(SP1_rs) sem(SP1_fs) sem(MI1_rs) sem(MI1_fs)])
legend({'RS', 'FS'})
xticks([1:2])
xticklabels({'Spont run + laser off', 'Evoked run + laser off'})
ylabel('Modulation Index')


% Plot firing rates by cortical layer to identify what drives laser effect
H = [];
for l = 1:4
    FR_evoked_means_laser_off(l) = nanmean(WN1(layers1 == l,2));
    FR_evoked_means_laser_on(l) = nanmean(WN1L(layers1 == l,2));
    FR_spont_means_laser_off(l) = nanmean(SP1(layers1 == l,2));
    FR_spont_means_laser_on(l) = nanmean(SP1L(layers1 == l,2));
    
    
    FR_evoked_medians_laser_off(l) = nanmedian(WN1(layers1 == l,2));
    FR_evoked_medians_laser_on(l) = nanmedian(WN1L(layers1 == l,2));
    FR_spont_medians_laser_off(l) = nanmedian(SP1(layers1 == l,2));
    FR_spont_medians_laser_on(l) = nanmedian(SP1L(layers1 == l,2));
    
    FR_evoked_sems_laser_off(l) = sem(WN1(layers1 == l,2));
    FR_evoked_sems_laser_on(l) = sem(WN1L(layers1 == l,2));
    FR_spont_sems_laser_off(l) = sem(SP1(layers1 == l,2));
    FR_spont_sems_laser_on(l) = sem(SP1L(layers1 == l,2));
    
    % stats
    [p,h,STATS] = signrank(WN1(layers1==l,2), WN1L(layers1==l,2));
    if p < 0.0125
        H(1,l) = 1.15;
    else
        H(1,l) = NaN;
    end
    
    
    [p,h,STATS] = signrank(SP1(layers1==l,2), SP1L(layers1==l,2));
    if p < 0.0125
        H(2,l) = 1.15;
    else
        H(2,l) = NaN;
    end
end

figure; subplot(2,1,1); hold on
errorbar([1:4], FR_evoked_means_laser_off, FR_evoked_sems_laser_off, 'ko-')
errorbar([1.2:4.2], FR_evoked_means_laser_on, FR_evoked_sems_laser_on, 'co-')
plot([1.1:4.1], max(FR_evoked_means_laser_on).*H(1,:), '*r')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Mean FR /SEM')
title('Evoked')
subplot(2,1,2); hold on;
errorbar([1:4], FR_spont_means_laser_off, FR_spont_sems_laser_off, 'ko-')
errorbar([1.2:4.2], FR_spont_means_laser_on, FR_spont_sems_laser_on, 'co-')
plot([1.1:4.1], max(FR_spont_means_laser_on).*H(2,:), '*r')
xticks([1:4]); xlim([0 5])
title('Spont')
xticklabels({'2/3', '4', '5', '6'});
ylabel('Mean/SEM FR')

figure; subplot(2,1,1); hold on
errorbar([1:4], FR_evoked_medians_laser_off, FR_evoked_sems_laser_off, 'ko-')
errorbar([1.2:4.2], FR_evoked_medians_laser_on, FR_evoked_sems_laser_on, 'co-')
plot([1.1:4.1], max(FR_evoked_medians_laser_on).*H(1,:), '*r')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
title('Evoked')
ylabel('Mean/SEM FR')

subplot(2,1,2); hold on;
errorbar([1:4], FR_spont_medians_laser_off, FR_spont_sems_laser_off, 'ko-')
errorbar([1.2:4.2], FR_spont_medians_laser_on, FR_spont_sems_laser_on, 'co-')
plot([1.1:4.1], max(FR_spont_medians_laser_on).*H(2,:), '*r')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
title('Spont')
ylabel('Mean/SEM FR')
xlabel('Cortical layers')


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

figure; hold on;
plot(MI_sound_sit(rs1), MI_sound_sit_laser(rs1), 'ko')
plot(MI_sound_sit(fs1), MI_sound_sit_laser(fs1), 'go')
plot(nanmean(MI_sound_sit), nanmean(MI_sound_sit_laser), 'ro', 'MarkerSize', 10)
plot(nanmedian(MI_sound_sit), nanmedian(MI_sound_sit_laser), 'mo', 'MarkerSize', 10)
maxFR = max(max([MI_sound_sit, MI_sound_sit_laser]));
plot([-1 1], [-1 1], 'r-')
plot([-1 1], [0 0], 'k--')
plot([0 0], [-1 1], 'k--')
[p,h,stats] = signrank(MI_sound_sit_laser,MI_sound_sit)
title_str = sprintf('Evoked Activity, p = %d', p);
title(title_str)
xlabel('sound modulation index - sit laser off');
ylabel('sound modulation index - sit laser on')
legend('Regular spiking', 'Narrow spiking', 'Mean', 'Median')
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

% VIP
MI1_Laser = (WN1L(:,2) - WN1(:,2))./(WN1L(:,2) + WN1(:,2)); % laser all cells
MI1_sp_Laser = (SP1L(:,2) - SP1(:,2))./(SP1L(:,2) + SP1(:,2)); % laser all cells



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

x=[MI_sound_sit_laser MI_sound_sit];
[p,tbl,stats] = kruskalwallis(x, [g1,g2])
