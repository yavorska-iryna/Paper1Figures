% compute modulation index WN responses. last editted 01.21.2020 ira
% variables gethered with getWNresponses.m

% three states: sit + small pupil, sit + large pupil, running  + large
% 

clear; close all
variables_dir = 'C:\Users\lab\Resilio Sync\Paper1Figures\code\variables';

cd(variables_dir);
load('WNdataLaserOFF.mat'); % 'WNdataLaserOFF.mat (OFF1 -  longer responses), dynamic pupil threshold, but varified.
data = WNdataLaserOFF;
load('WNdataLaserON.mat');
data1 = WNdataLaserON;
load('CellsQualityStats.mat')

 maxFRall =[];
SP = nan(length(data,4);
SP = nan(length(data,4));
meanON = nan(length(data),1); % On response only (0 - 100 ms)
meanONL = nan(length(data),1); 
meanWN = nan(length(data),4); % Full response (0 - 600 ms)
meanWNL = nan(length(data),4);
depths = nan(length(data),1); 
fs = zeros(length(data),1); 
Rs = zeros(length(data),1);
dirsM =[]; cellsM = 0;  dirsP =[]; cellsP = 0; 
recs = []; cells=[];

evoked = nan(length(data),4); 
color1 = [0.7 0.7 0.7]; color2 = [0.2 0.2 0.2];
nreps_check_running = 7; %number of repetitions in each condition for comparison
nreps_check_pupil = 8;
CL = {[0 301], [300 401], [400 601], [600 2000]};


cdVIP; load('Silence_DistanceCorr_dirs.mat')

for cc =1:length(data)
    if data(cc).dir<46
        try
            meanSpikeCount = nanmean([data(cc).SpikeCountWN data1(cc).SpikeCountWN ]); %data(cc).SpikeCountSS data1(cc).SpikeCountSS
        catch
            meanSpikeCount = NaN;
        end
        % exclude cells with very low spikecount, they usually have very
        % large effects; count spikes to WN only
        if meanSpikeCount > 0 && CellsQualityStats.SNR(cc)>.5 && CellsQualityStats.uQ(cc)>10 
            
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
            
            if data(cc).nrepsWNMon > nreps_check_running && data(cc).nrepsWNMoff > nreps_check_running %&& data1(cc).nrepsWNMon > nreps_check_running && data1(cc).nrepsWNMoff > nreps_check_running
                meanWN(cc,1)= nanmean([nanmean(nanmean(data(cc).mWNon))]);% nanmean(nanmean(data(cc).mNSustained_on)) nanmean(nanmean(data(cc).mNOFF_on))]);
                meanWN(cc,4)= nanmean([nanmean(nanmean(data(cc).mWNoff))]);% nanmean(nanmean(data(cc).mNSustained_off)) nanmean(nanmean(data(cc).mNOFF_off))]);
                meanON(cc,1) =nanmean([nanmean(nanmean(data(cc).mNON_on))])-nanmean(nanmean(data(cc).mNson));
                meanON(cc,4) =nanmean([nanmean(nanmean(data(cc).mNON_off))])-nanmean(nanmean(data(cc).mNsoff));
                
                dirsM = [dirsM data(cc).dir];
                cellsM = cellsM+1;
            end
            if data(cc).dir~= 0 
                if data(cc).nrepsWNPon > nreps_check_pupil && data(cc).nrepsWNPoff > nreps_check_pupil %&& data1(cc).nrepsWNPon > nreps_check_pupil && data1(cc).nrepsWNPoff > nreps_check_pupil
                    meanWN(cc,2)= nanmean(nanmean(data(cc).pWNon));%nanmean([nanmean(nanmean(data(cc).pNON_on)) nanmean(nanmean(data(cc).pNSustained_on)) nanmean(nanmean(data(cc).pNOFF_on))]);
                    meanWN(cc,3)= nanmean(nanmean(data(cc).pWNoff));%nanmean([nanmean(nanmean(data(cc).pNON_off))] nanmean(nanmean(data(cc).pNSustained_off)) nanmean(nanmean(data(cc).pNOFF_off))]);
                    meanON(cc,2)= nanmean(nanmean(data(cc).pNON_on))-nanmean(nanmean(data(cc).pNson));%nanmean([nanmean(nanmean(data(cc).pNON_on)) nanmean(nanmean(data(cc).pNSustained_on)) nanmean(nanmean(data(cc).pNOFF_on))]);
                    meanON(cc,3)= nanmean(nanmean(data(cc).pNON_off))-nanmean(nanmean(data(cc).pNsoff));
                end
            end
            
            % laser on trials
            if data1(cc).nrepsWNMon > nreps_check_running && data1(cc).nrepsWNMoff > nreps_check_running && data(cc).nrepsWNMon > nreps_check_running && data(cc).nrepsWNMoff > nreps_check_running
                meanWNL(cc,1)= nanmean(nanmean(data1(cc).mWNon));%nanmean([nanmean(nanmean(data1(cc).mNON_on)) nanmean(nanmean(data1(cc).mNSustained_on)) nanmean(nanmean(data1(cc).mNOFF_on))]);
                meanWNL(cc,4)= nanmean(nanmean(data1(cc).mWNoff));%nanmean([nanmean(nanmean(data1(cc).mNON_off)) nanmean(nanmean(data1(cc).mNSustained_off)) nanmean(nanmean(data1(cc).mNOFF_off))]);
                meanONL(cc,1) =nanmean([nanmean(nanmean(data1(cc).mNON_on))])-nanmean(nanmean(data1(cc).mNson));
                meanONL(cc,4) =nanmean([nanmean(nanmean(data1(cc).mNON_off))])-nanmean(nanmean(data1(cc).mNsoff));
            
            end
            if data(cc).dir~= 0 
                if data1(cc).nrepsWNPon > nreps_check_pupil && data1(cc).nrepsWNPoff > nreps_check_pupil && data(cc).nrepsWNPon > nreps_check_pupil && data(cc).nrepsWNPoff > nreps_check_pupil
                    meanWNL(cc,2)= nanmean(nanmean(data1(cc).pWNon));%nanmean([nanmean(nanmean(data1(cc).pNON_on)) nanmean(nanmean(data1(cc).pNSustained_on)) nanmean(nanmean(data1(cc).pNOFF_on))]);
                    meanWNL(cc,3)= nanmean(nanmean(data1(cc).pWNoff));%nanmean([nanmean(nanmean(data1(cc).pNON_off)) nanmean(nanmean(data1(cc).pNSustained_off)) nanmean(nanmean(data1(cc).pNOFF_off))]);
                    meanONL(cc,2) =nanmean([nanmean(nanmean(data1(cc).pNON_on))])-nanmean(nanmean(data1(cc).pNson));
                    meanONL(cc,3) =nanmean([nanmean(nanmean(data1(cc).pNON_off))])-nanmean(nanmean(data1(cc).pNsoff));
                end
            end
            
            %silent sound
            if data(cc).nrepsSSMon > nreps_check_running && data(cc).nrepsSSMoff > nreps_check_running %&& data1(cc).nrepsSSMon > nreps_check_running && data1(cc).nrepsSSMoff > nreps_check_running
                SP(cc,1) = nanmean(nanmean(data(cc).mSSon)); SP(cc,4) = nanmean(nanmean(data(cc).mSSoff));
            end
             if data1(cc).nrepsSSMon > nreps_check_running && data1(cc).nrepsSSMoff > nreps_check_running && data(cc).nrepsSSMon > nreps_check_running && data(cc).nrepsSSMoff > nreps_check_running
                SPL(cc,1) = nanmean(nanmean(data1(cc).mSSon)); SPL(cc,4) = nanmean(nanmean(data1(cc).mSSoff));
             end
            
            if data(cc).dir~= 0 
                if data(cc).nrepsSSPon > nreps_check_pupil && data(cc).nrepsSSPoff > nreps_check_pupil %&& data1(cc).nrepsSSPon > nreps_check_pupil && data1(cc).nrepsSSPoff > nreps_check_pupil
                    SP(cc,2) =  nanmean(nanmean(data(cc).pSSon)); SP(cc,3) =  nanmean(nanmean(data(cc).pSSoff));
                    dirsP = [dirsP data(cc).dir];
                    cellsP = cellsP+1;
                end
                if data1(cc).nrepsSSPon > nreps_check_pupil && data1(cc).nrepsSSPoff > nreps_check_pupil && data(cc).nrepsSSPon > nreps_check_pupil && data(cc).nrepsSSPoff > nreps_check_pupil
                    SPL(cc, 2) = nanmean(nanmean(data1(cc).pSSon)); SPL(cc, 3) = nanmean(nanmean(data1(cc).pSSoff));
                end
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

allDirsM = unique(dirsM);
allDirsP = unique(dirsP);

% full time window response
evoked2 = evoked(:,4); % evoked is 1 and 0; 1) ON, 2) Sustained, 3) OFF, 4) Full Response comparison 
notresp = find(evoked2==0);
evoked2 = logical(evoked2);

fs2 = fs(evoked2); % only evoked fs and rs
rs2 = Rs(evoked2);
recs1 = recs(evoked2); cells1 = cells(evoked2);
WN1 = meanWN(evoked2,1:4);
WN1L = meanWNL(evoked2,1:4);
SP1 = SP(evoked2,:);
SP1L = SPL(evoked2,:);
depths2 = depths(evoked2);
rs2 = logical(rs2); fs2 = logical(fs2);

%modulation index by cell type
modulation_indx1_rs = (WN1(rs2,1) -WN1(rs2,4))./(WN1(rs2,1) + WN1(rs2,4)); %running
modulation_indx2_rs = (WN1(rs2,2) -WN1(rs2,3))./(WN1(rs2,2) + WN1(rs2,3)); %large pupil + sit
modulation_indx1_fs = (WN1(fs2,1) -WN1(fs2,4))./(WN1(fs2,1) + WN1(fs2,4)); %running
modulation_indx2_fs = (WN1(fs2,2) -WN1(fs2,3))./(WN1(fs2,2) + WN1(fs2,3)); %large pupil + sit

modulation_indx1_sp_rs = (SP1(rs2,1) -SP1(rs2,4))./(SP1(rs2,1) + SP1(rs2,4)); %running
modulation_indx2_sp_rs = (SP1(rs2,2) -SP1(rs2,3))./(SP1(rs2,2) + SP1(rs2,3)); %large pupil + sit
modulation_indx1_sp_fs = (SP1(fs2,1) -SP1(fs2,4))./(SP1(fs2,1) + SP1(fs2,4)); %running
modulation_indx2_sp_fs = (SP1(fs2,2) -SP1(fs2,3))./(SP1(fs2,2) + SP1(fs2,3)); %large pupil + sit

modulation_indx1_rsL = (WN1L(rs2,1) -WN1L(rs2,4))./(WN1L(rs2,1) + WN1L(rs2,4)); %running
modulation_indx2_rsL = (WN1L(rs2,2) -WN1L(rs2,3))./(WN1L(rs2,2) + WN1L(rs2,3)); %large pupil + sit
modulation_indx1_fsL = (WN1L(fs2,1) -WN1L(fs2,4))./(WN1L(fs2,1) + WN1L(fs2,4)); %running
modulation_indx2_fsL = (WN1L(fs2,2) -WN1L(fs2,3))./(WN1L(fs2,2) + WN1L(fs2,3)); %large pupil + sit

modulation_indx1_sp_rsL = (SP1L(rs2,1) -SP1L(rs2,4))./(SP1L(rs2,1) + SP1L(rs2,4)); %running
modulation_indx2_sp_rsL = (SP1L(rs2,2) -SP1L(rs2,3))./(SP1L(rs2,2) + SP1L(rs2,3)); %large pupil + sit
modulation_indx1_sp_fsL = (SP1L(fs2,1) -SP1L(fs2,4))./(SP1L(fs2,1) + SP1L(fs2,4)); %running
modulation_indx2_sp_fsL = (SP1L(fs2,2) -SP1L(fs2,3))./(SP1L(fs2,2) + SP1L(fs2,3)); %large pupil + sit

% modulation index total
modulation_indx1 = (WN1(:,1) -WN1(:,4))./(WN1(:,1) + WN1(:,4)); %running
modulation_indx1_sp= (SP1(:,1) -SP1(:,4))./(SP1(:,1) + SP1(:,4));

% activated cell only
evoked1 = evoked(:,1) & zstats(:,4)>0;
evoked1 = ones(length(evoked),1);
evoked1 = logical(evoked1);
actWN1 = meanON(evoked1,1:4);
actWN1L = meanONL(evoked1,1:4);
actSP1 = SP(evoked1,1:4);
actSP1L = SPL(evoked1,1:4);
depth1 = depths(evoked1);
fs1 = fs(evoked1); % only evoked fs and rs
rs1 = Rs(evoked1);
rs1 = logical(rs1); fs1 = logical(fs1);

EvokedWN = actWN1;
modulation_indx1_evoked = (EvokedWN(:,1) - EvokedWN(:,4))./(EvokedWN(:,1) + EvokedWN(:,4));
modulation_indx1_sp_act = (actSP1(:,1) - actSP1(:,4))./(actSP1(:,1) + actSP1(:,4));

%% plot the results
% state means by RS FS WN SP
figure(1); hold on
bar([0.8 1.8], [nanmean(modulation_indx1_sp_rs) nanmean(modulation_indx2_sp_rs)], 'BarWidth', .1)
bar([0.9 1.9], [nanmean(modulation_indx1_rs) nanmean(modulation_indx2_rs)], 'BarWidth', .1)
bar([1  2 ], [ nanmean(modulation_indx1_sp_fs) nanmean(modulation_indx2_sp_fs ) ], 'BarWidth', .1)
bar([1.1 2.1], [nanmean(modulation_indx1_fs) nanmean(modulation_indx2_fs)], 'BarWidth', .1)
xticks([1:2])
xticklabels({'run', 'sit + large pupil'})
ylabel('Modulation Index')
errorbar( [0.8 0.9 1 1.1 1.8 1.9 2 2.1],[nanmean(modulation_indx1_sp_rs) nanmean(modulation_indx1_rs) nanmean(modulation_indx1_sp_fs) nanmean(modulation_indx1_fs)  ...
    nanmean(modulation_indx2_sp_rs) nanmean(modulation_indx2_rs) nanmean(modulation_indx2_sp_fs) nanmean(modulation_indx2_fs)], ...
    [sem(modulation_indx1_sp_rs) sem(modulation_indx1_rs) sem(modulation_indx1_sp_fs) sem(modulation_indx1_fs)  ...
    sem(modulation_indx2_sp_rs) sem(modulation_indx2_rs) sem(modulation_indx2_sp_fs) sem(modulation_indx2_fs)])
legend({'spont rs', 'evoked rs', 'spont fs', 'evoked fs'})

figure(2); hold on
bar([1 2], [nanmean(modulation_indx1_sp) nanmean(modulation_indx1)], 'BarWidth', .4)
xticks([1:2])
xticklabels({'spont', 'response'})
ylabel('Modulation Index (mean/SEM)')
errorbar( [1 2], [nanmean(modulation_indx1_sp) nanmean(modulation_indx1)],  ...
    [sem(modulation_indx1_sp) sem(modulation_indx1)])
title('Effects of running on spont and evoked activity')

figure(3); hold on
bar([1 2], [nanmean(modulation_indx1_sp_act) nanmean(modulation_indx1_evoked)], 'BarWidth', .4)
xticks([1:2])
xticklabels({'spont (act resp)', 'evoked (response-sp)'})
ylabel('Modulation Index (mean/SEM)')
errorbar( [1 2], [nanmean(modulation_indx1_sp_act) nanmean(modulation_indx1_evoked)],  ...
    [sem(modulation_indx1_sp_act) sem(modulation_indx1_evoked)])
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

figure; subplot(1,2,1); hold on
hist(modulation_indx1_sp,[-1: .05:1], 'Color', [0.8 0.8 0.8])
xlim([-1 1]); xlabel('Modulation Index'); ylabel('Number of cells')
plot([nanmean(modulation_indx1_sp) nanmean(modulation_indx1_sp)], [0 max(ylim)], 'r-', 'Linewidth', 2)
plot([nanmedian(modulation_indx1_sp) nanmedian(modulation_indx1_sp)], [0 max(ylim)], 'c-', 'Linewidth', 2)
title('spont'); legend('hist','mean','median')
subplot(1,2,2); hold on
hist(modulation_indx1_evoked, [-1: .05:1], 'Color', [0.8 0.8 0.8])
xlim([-1 1]); xlabel('Modulation Index'); ylabel('Number of cells')
plot([nanmean(modulation_indx1_evoked) nanmean(modulation_indx1_evoked)], [0 max(ylim)], 'r-', 'Linewidth', 2)
plot([nanmedian(modulation_indx1_evoked) nanmedian(modulation_indx1_evoked)], [0 max(ylim)], 'c-', 'Linewidth', 2)
title('evoked'); legend('hist','mean','median')
set(gcf, 'Position',  [17 558 1223 420])
%%

recs1= (recs(evoked2));
recs_rs= (recs1(rs2));
recs_fs=(recs1(fs2));

figure; plot( recs_rs,modulation_indx1_sp_rs, 'ko');
nanindx = find(isnan(modulation_indx1_sp_rs)==1);
x=modulation_indx1_sp_rs;
[p,tbl1,stats] = kruskalwallis(x, recs_rs);
c = multcompare(stats);
title('Spont RS, running')
x=modulation_indx1_fs;
[p,tbl1,stats] = kruskalwallis(x, recs_fs);
c = multcompare(stats);
title('Evoked RS, running')
x=modulation_indx1_sp;
[p,tbl1,stats] = kruskalwallis(x, recs1);
c = multcompare(stats);
title('Spont RS, pupil')
[m,s]=grpstats(x,recs_rs,{'mean','sem'})

x=modulation_indx2_sp_fs;
[p,tbl1,stats] = kruskalwallis(x, recs_fs);
c = multcompare(stats);
title('Spont FS, pupil')
x=modulation_indx2_sp_rs;
[p,tbl1,stats] = kruskalwallis(x, recs_rs);
c = multcompare(stats);
title('Spont RS, pupil')

%[m,s]=grpstats(x,recs_rs,{'mean','sem'})
r =unique(recs_rs);
%%
MODULATION_INDX1 = []; MODULATION_INDX2 =[]; MODULATION_INDX1L = []; MODULATION_INDX2L = [];
MODULATION_sp_INDX1 = []; MODULATION_sp_INDX2 =[]; MODULATION_sp_INDX1L = []; MODULATION_sp_INDX2L = [];
cells1 = cells(evoked2); recs1 = recs(evoked2);
for cl = 1:length(CL)
    layer = CL{cl};
    indx = find(depths2 >layer(1) & depths2 <layer(2));
    modulation_indx1 = (WN1(indx,1) -WN1(indx,4))./(WN1(indx,1) + WN1(indx,4)); %running
    modulation_indx2 = (WN1(indx,2) -WN1(indx,3))./(WN1(indx,2) + WN1(indx,3)); %large pupil + sit
    
    best_examples = find(modulation_indx1>0.2);
    recs2 = recs1(indx);
    cells2 = cells1(indx);
    recs2(best_examples)
    cells2(best_examples)
    
    modulation_indx1L = (WN1L(indx,1) -WN1L(indx,3))./(WN1L(indx,1) + WN1L(indx,3)); %running + laser
    modulation_indx2L = (WN1L(indx,2) -WN1L(indx,3))./(WN1L(indx,2) + WN1L(indx,3)); %large pupil + sit +laser
    modulation_indx1Laser = (nanmean(WN1L(indx,2:2:4),2) - nanmean(WN1(indx,2:2:4),2))./(nanmean(WN1L(indx,2:2:4),2) + nanmean(WN1(indx,2:2:4),2)); %running + laser
    modulation_indx2Laser = (nanmean(WN1L(indx,2:3),2) - nanmean(WN1(indx,2:3),2))./(nanmean(WN1L(indx,2:3),2) + nanmean(WN1(indx,2:3),2)); %pupil
    
    effect(cl).mi = modulation_indx1;
    effect(cl).cells = cells2;
    effect(cl).recs = recs2;
    
    
    fs2 = fs2(indx); fs2 = logical(fs2);
    rs2 = rs2(indx); rs2 = logical(rs2);
    
    recs2 = recs1(indx);
    cells2 = cells1(indx);
    %stats
    MODULATION_INDX1 = [MODULATION_INDX1; modulation_indx1 ones(length(indx),1)*cl];
    MODULATION_INDX2 = [MODULATION_INDX2; modulation_indx2 ones(length(indx),1)*cl];
    MODULATION_INDX1L = [MODULATION_INDX1L; modulation_indx1L ones(length(indx),1)*cl];
    MODULATION_INDX2L = [MODULATION_INDX2L; modulation_indx2L ones(length(indx),1)*cl];
    %
    meanMI1_Laser(cl) = nanmean(modulation_indx1Laser);
    meanMI2_Laser(cl) = nanmean(modulation_indx2Laser);
    semMI1_Laser(cl) = sem(modulation_indx1Laser);
    semMI2_Laser(cl) = sem(modulation_indx2Laser);
    
    meanMI1(cl) = nanmean(modulation_indx1);
    n_layer_evoked(cl,1) = length(find(~isnan(modulation_indx1)==1));
    n_layer_evokedL(cl) = length(find(~isnan(modulation_indx1L)==1));
    meanMI2(cl) = nanmean(modulation_indx2);
    n_layer_evoked(cl,2) = length(find(~isnan(modulation_indx2)==1));
    semMI1(cl) = sem(modulation_indx1);
    semMI2(cl) = sem(modulation_indx2);
    meanMI1_fs(cl) = nanmean(modulation_indx1(fs2));
    meanMI2_fs(cl) = nanmean(modulation_indx2(fs2));
    semMI1_fs(cl) = sem(modulation_indx1(fs2));
    semMI2_fs(cl) = sem(modulation_indx2(fs2));
    meanMI1_rs(cl) = nanmean(modulation_indx1(rs2));
    meanMI2_rs(cl) = nanmean(modulation_indx2(rs2));
    semMI1_rs(cl) = sem(modulation_indx1(rs2));
    semMI2_rs(cl) = sem(modulation_indx2(rs2));
    
    rs_cl(cl,1) = sum(~isnan(modulation_indx1(rs2)));
    fs_cl(cl,1) = sum(~isnan(modulation_indx1(fs2)));
    rs_cl(cl,2) = sum(~isnan(modulation_indx2(rs2)));
    fs_cl(cl,2) = sum(~isnan(modulation_indx2(fs2)));
    
    % addative effect
    MI_evoked_plus = (WN1L(indx,1) - WN1(indx,4))./(WN1L(indx,1) + WN1(indx,4));
    meanMI_evoked_plus(cl) = nanmean(MI_evoked_plus);
    semMI_evoked_plus(cl) = sem(MI_evoked_plus);
    n_layer_evoked_plus(cl) = length(find(~isnan(MI_evoked_plus)==1));
    MI_evoked_plus = (WN1L(indx,2) - WN1(indx,3))./(WN1L(indx,2) + WN1(indx,3));
    meanMI2_evoked_plus(cl) = nanmean(MI_evoked_plus);
    semMI2_evoked_plus(cl) = sem(MI_evoked_plus);
    
    %laser on
    meanMI1L(cl) = nanmean(modulation_indx1L);
    meanMI2L(cl) = nanmean(modulation_indx2L);
    semMI1L(cl) = sem(modulation_indx1L);
    semMI2L(cl) = sem(modulation_indx2L);
    meanMI1_fsL(cl) = nanmean(modulation_indx1L(fs2));
    meanMI2_fsL(cl) = nanmean(modulation_indx2L(fs2));
    semMI1_fsL(cl) = sem(modulation_indx1L(fs2));
    semMI2_fsL(cl) = sem(modulation_indx2L(fs2));
    meanMI1_rsL(cl) = nanmean(modulation_indx1L(rs2));
    meanMI2_rsL(cl) = nanmean(modulation_indx2L(rs2));
    semMI1_rsL(cl) = sem(modulation_indx1L(rs2));
    semMI2_rsL(cl) = sem(modulation_indx2L(rs2));
    
    % spont act
    indx = find(depths2 >layer(1) & depths2 <layer(2));
    fs2 = fs2(indx); fs2 = logical(fs2);
    rs2 = Rs(indx); rs2 = logical(rs2);
    modulation_indx1_sp = (SP1(indx,1) -SP1(indx,4))./(SP1(indx,1) + SP1(indx,4)); %running + laser
    modulation_indx2_sp = (SP1(indx,2) - SP1(indx,3))./(SP1(indx,2) +SP1(indx,3)); %large pupil + sit +laser
    modulation_indx1_spL = (SP1L(indx,1) -SP1L(indx,4))./(SP1L(indx,1) + SP1L(indx,4)); %running + laser
    modulation_indx2_spL = (SP1L(indx,2) - SP1L(indx,3))./(SP1L(indx,2) +SP1L(indx,3)); %large pupil + sit +laser
    modulation_indx1_spLaser = (nanmean(SP1L(indx,2:2:4),2) - nanmean(SP1(indx,2:2:4),2))./(nanmean(SP1(indx,2:2:4),2) + nanmean(SP1(indx,2:2:4),2)); %running + laser
    modulation_indx2_spLaser = (nanmean(SP1L(indx,2:3),2) - nanmean(SP1(indx,2:3),2))./(nanmean(SP1(indx,2:3),2) + nanmean(SP1(indx,2:3),2)); %running + laser
    
    MI_sp_plus = (SP1L(indx,1) - SP1(indx,4))./(SP1L(indx,1) + SP1(indx,4));
    meanMI_sp_plus(cl) = nanmean(MI_sp_plus);
    semMI_sp_plus(cl) = sem(MI_sp_plus);
    n_layer_sp_plus(cl) = length(find(~isnan(MI_sp_plus)==1));
    MI2_sp_plus = (SP1L(indx,2) - SP1(indx,3))./(SP1L(indx,2) + SP1(indx,3));
    meanMI2_sp_plus(cl) = nanmean(MI_sp_plus);
    semMI2_sp_plus(cl) = sem(MI_sp_plus);
    
    meanMI1_sp(cl) = nanmean(modulation_indx1_sp);
    n_layer_sp(cl,1) = length(find(~isnan(modulation_indx1_sp)==1)); %number of cells in each layer
    n_layer_spL(cl) = length(find(~isnan(modulation_indx1_spL)==1));
    meanMI1_spL(cl) = nanmean(modulation_indx1_spL);
    meanMI1_spLaser(cl) = nanmean(modulation_indx1_spLaser);
    semMI1_spLaser(cl) = sem(modulation_indx1_spLaser);
    meanMI2_spLaser(cl) = nanmean(modulation_indx2_spLaser);
    semMI2_spLaser(cl) = sem(modulation_indx2_spLaser);
    
    meanMI2_sp(cl) = nanmean(modulation_indx2_sp);
    n_layer_sp(cl,2) = length(find(~isnan(modulation_indx2_sp)==1)); %number of cells in each layer
    meanMI2_spL(cl) = nanmean(modulation_indx2_spL);
    semMI1_sp(cl) = sem(modulation_indx1_sp);
    semMI2_sp(cl) = sem(modulation_indx2_sp);
    semMI1_spL(cl) = sem(modulation_indx1_spL);
    semMI2_spL(cl) = sem(modulation_indx2_spL);
    
    meanMI1_sp_fs(cl) = nanmean(modulation_indx1_sp(fs2));
    rs_sp_cl(cl,1) = sum(~isnan(modulation_indx1_sp(rs2)));
    fs_sp_cl(cl,1) = sum(~isnan(modulation_indx1_sp(fs2)));
    rs_sp_cl(cl,2) = sum(~isnan(modulation_indx2_sp(rs2)));
    fs_sp_cl(cl,2) = sum(~isnan(modulation_indx2_sp(fs2)));
    
    meanMI1_spL_fs(cl) = nanmean(modulation_indx1_spL(fs2));
    meanMI2_sp_fs(cl) = nanmean(modulation_indx2_sp(fs2));
    meanMI2_spL_fs(cl) = nanmean(modulation_indx2_spL(fs2));
    semMI1_sp_fs(cl) = sem(modulation_indx1_sp(fs2));
    semMI2_sp_fs(cl) = sem(modulation_indx2_sp(fs2));
    semMI1_spL_fs(cl) = sem(modulation_indx1_spL(fs2));
    semMI2_spL_fs(cl) = sem(modulation_indx2_spL(fs2));
    
    meanMI1_sp_rs(cl) = nanmean(modulation_indx1_sp(rs2));
    meanMI1_spL_rs(cl) = nanmean(modulation_indx1_spL(rs2));
    meanMI2_sp_rs(cl) = nanmean(modulation_indx2_sp(rs2));
    meanMI2_spL_rs(cl) = nanmean(modulation_indx2_spL(rs2));
    semMI1_sp_rs(cl) = sem(modulation_indx1_sp(rs2));
    semMI2_sp_rs(cl) = sem(modulation_indx2_sp(rs2));
    semMI1_spL_rs(cl) = sem(modulation_indx1_spL(rs2));
    semMI2_spL_rs(cl) = sem(modulation_indx2_spL(rs2));
    
    MODULATION_sp_INDX1 = [MODULATION_sp_INDX1; modulation_indx1_sp ones(length(indx),1)*cl];
    MODULATION_sp_INDX2 = [MODULATION_sp_INDX2; modulation_indx2_sp ones(length(indx),1)*cl];
    MODULATION_sp_INDX1L = [MODULATION_sp_INDX1L; modulation_indx1_spL ones(length(indx),1)*cl];
    MODULATION_sp_INDX2L = [MODULATION_sp_INDX2L; modulation_indx2_spL ones(length(indx),1)*cl];
    
end

% state by layer without layer
figure(101); subplot(2,1,1); hold on
errorbar([1:4], meanMI1, semMI1, 'ko');
errorbar([1.2:4.2], meanMI2, semMI2, 'ro');
xlabel('Cortical Layer')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Modulation Index'); title('WN response')
% [p,tbl1,stats] = kruskalwallis(mi, layers);
% c = multcompare(stats);
subplot(2,1,2); hold on
errorbar([1:4], meanMI1_sp, semMI1_sp, 'ko');
errorbar([1.2:4.2], meanMI2_sp, semMI2_sp, 'ro');

xlabel('Cortical Layer')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Modulation Index'); title('Spont activity');

%State by layer with laser
figure(101); subplot(2,1,2); hold on
errorbar([1.1:4.1], meanMI1_spL, semMI1_spL, 'k>');
errorbar([1.3:4.3], meanMI2_spL, semMI2_spL, 'b>');
plot([0 5], [0 0], 'k--')
legend({ 'run laser off', 'sit + large pupil laser off', 'run laser on', 'sit + large pupil laser on'})
xlabel('Cortical Layer')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Modulation Index'); title('Spont activity - VIP laser on')
subplot(2,1,1); hold on
errorbar([1.1:4.1], meanMI1L, semMI1L, 'k>');
errorbar([1.3:4.3], meanMI2L, semMI2L, 'b>');
plot([0 5], [0 0], 'k--')
legend({ 'run laser off', 'sit + large pupil laser off', 'run laser on', 'sit + large pupil laser on'})
xlabel('Cortical Layer')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Modulation Index'); title('Evoked activity - VIP laser on')


% RS FS WN SP1 by layer
figure(202); subplot(2,1,2); hold on
errorbar([1:4], meanMI1_sp_rs, semMI1_sp_rs, 'ko');
errorbar([1.2:4.2], meanMI2_sp_rs, semMI2_sp_rs, 'ro');
errorbar([1.1:4.1], meanMI1_sp_fs, semMI1_sp_fs, 'k*');
errorbar([1.3:4.3], meanMI2_sp_fs, semMI2_sp_fs, 'b*');
plot([0 5], [0 0], 'k--')
legend({ 'run rs', 'sit + large pupil rs', 'run fs', 'sit + large pupil fs'})
xlabel('Cortical Layer')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Modulation Index'); title('Spont Activity')
subplot(2,1,1); hold on
errorbar([1:4], meanMI1_rs, semMI1_rs, 'ko');
errorbar([1.2:4.2], meanMI2_rs, semMI2_rs, 'ro');
errorbar([1.1:4.1], meanMI1_fs, semMI1_fs, 'k*');
errorbar([1.3:4.3], meanMI2_fs, semMI2_fs, 'b*');
plot([0 5], [0 0], 'k--')
legend({ 'run rs', 'sit + large pupil rs', 'run fs', 'sit + large pupil fs'})
xlabel('Cortical Layer')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Modulation Index'); title('WN responses')


% compare laser
%running
figure(103); subplot(2,1,1); hold on
errorbar([1:4], meanMI1, semMI1, 'ko');
errorbar([1.1:4.1], meanMI1_Laser, semMI1_Laser, 'bo');
errorbar([1.2:4.2], meanMI_evoked_plus, semMI_evoked_plus, 'go');
plot([0 5], [0 0], 'k--')
legend({ 'run', 'laser' 'run + laser' })
xlabel('Cortical Layer')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Modulation Index - Run'); title('Evoked responses')
figure(103); subplot(2,1,2); hold on
errorbar([1:4], meanMI1_sp, semMI1_sp, 'ko');
errorbar([1.1:4.1], meanMI1_spLaser, semMI1_spLaser, 'bo');
errorbar([1.2:4.2], meanMI_sp_plus, semMI_sp_plus, 'go');
plot([0 5], [0 0], 'k--')
legend({ 'run', 'laser' 'run + laser' })
xlabel('Cortical Layer')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Modulation Index - Run'); title('Spont responses')
%pupils
figure(104); subplot(2,1,1); hold on
errorbar([1:4], meanMI2, semMI2, 'ko');
errorbar([1.1:4.1], meanMI2_Laser, semMI2_Laser, 'bo');
errorbar([1.2:4.2], meanMI2_evoked_plus, semMI2_evoked_plus, 'go');
plot([0 5], [0 0], 'k--')
legend({ 'pupil', 'laser' 'pupil + laser' })
xlabel('Cortical Layer')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Modulation Index - Pupil'); title('Evoked responses')
subplot(2,1,2); hold on
errorbar([1:4], meanMI2_sp, semMI2_sp, 'ko');
errorbar([1.1:4.1], meanMI2_spLaser, semMI2_spLaser, 'bo');
errorbar([1.2:4.2], meanMI2_sp_plus, semMI2_sp_plus, 'go');
plot([0 5], [0 0], 'k--')
legend({ 'pupil', 'laser' 'pupil + laser' }); xlabel('Cortical Layer')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Modulation Index - Pupil'); title('Spont responses')

fs2 = logical(fs2); rs2 = logical(rs2);
MI1_rs = (WN1(rs2,1) - WN1(rs2,4))./ (WN1(rs2,1) +WN1(rs2,4));
MI2_rs = (WN1(rs2,2) - WN1(rs2,3))./ (WN1(rs2,2) +WN1(rs2,3));
MI1_fs = (WN1(fs2,1) - WN1(fs2,4))./ (WN1(fs2,1) +WN1(fs2,4));
MI2_fs = (WN1(fs2,2) - WN1(fs2,3))./ (WN1(fs2,2) +WN1(fs2,3));

SP1_rs = (SP1(rs2,1) - SP1(rs2,4))./ (SP1(rs2,1) +SP1(rs2,4));
SP2_rs = (SP1(rs2,2) - SP1(rs2,3))./ (SP1(rs2,2) +SP1(rs2,3));
SP1_fs = (SP1(fs2,1) - SP1(fs2,4))./ (SP1(fs2,1) +SP1(fs2,4));
SP2_fs = (SP1(fs2,2) - SP1(fs2,3))./ (SP1(fs2,2) +SP1(fs2,3));

%laser VIP
MI1_rsL = (WN1L(rs2,1) - WN1L(rs2,4))./ (WN1L(rs2,1) +WN1L(rs2,4));
MI2_rsL = (WN1L(rs2,2) - WN1L(rs2,3))./ (WN1L(rs2,2) +WN1L(rs2,3));
MI1_fsL = (WN1L(fs2,1) - WN1L(fs2,4))./ (WN1L(fs2,1) +WN1L(fs2,4));
MI2_fsL = (WN1L(fs2,2) - WN1L(fs2,3))./ (WN1L(fs2,2) +WN1L(fs2,3));
SP1_rsL = (SP1L(rs2,1) - SP1L(rs2,4))./ (SP1L(rs2,1) +SP1L(rs2,4));
SP2_rsL = (SP1L(rs2,2) - SP1L(rs2,3))./ (SP1L(rs2,2) +SP1L(rs2,3));
SP1_fsL = (SP1L(fs2,1) - SP1L(fs2,4))./ (SP1L(fs2,1) +SP1L(fs2,4));
SP2_fsL = (SP1L(fs2,2) - SP1L(fs2,3))./ (SP1L(fs2,2) +SP1L(fs2,3));


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

modulation_indx1_rs = (WN1L(rs2,1) -WN1(rs2,1))./(WN1L(rs2,1) + WN1(rs2,1)); %running
modulation_indx2_rs = (WN1L(rs2,2) -WN1(rs2,2))./(WN1L(rs2,2) + WN1(rs2,2)); %large pupil + sit
modulation_indx3_rs = (WN1L(rs2,4) -WN1(rs2,4))./(WN1L(rs2,4) + WN1(rs2,4)); %small pupil + sit
modulation_indx1_fs = (WN1L(fs2,1) -WN1(fs2,1))./(WN1L(fs2,1) + WN1(fs2,1)); %running
modulation_indx2_fs = (WN1L(fs2,2) -WN1(fs2,2))./(WN1L(fs2,2) + WN1(fs2,2)); %large pupil + sit
modulation_indx3_fs = (WN1L(fs2,4) -WN1(fs2,4))./(WN1L(fs2,4) + WN1(fs2,4)); %small pupil + sit

modulation_indx1_sp_rs = (SP1L(rs2,1) -SP1(rs2,1))./(SP1L(rs2,1) + SP1(rs2,1)); %running
modulation_indx2_sp_rs = (SP1L(rs2,2) -SP1(rs2,2))./(SP1L(rs2,2) + SP1(rs2,2)); %large pupil + sit
modulation_indx3_sp_rs = (SP1L(rs2,4) -SP1(rs2,4))./(SP1L(rs2,4) + SP1(rs2,4)); %small pupil + sit
modulation_indx1_sp_fs = (SP1L(fs2,1) -SP1(fs2,1))./(SP1L(fs2,1) + SP1(fs2,1)); %running
modulation_indx2_sp_fs = (SP1L(fs2,2) -SP1(fs2,2))./(SP1L(fs2,2) + SP1(fs2,2)); %large pupil + sit
modulation_indx3_sp_fs = (SP1L(fs2,4) -SP1(fs2,3))./(SP1L(fs2,4) + SP1(fs2,4)); %small pupil + sit

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
    fs2 = fs2(indx); fs2 = logical(fs2);
    rs2 = rs2(indx); rs2 = logical(rs2);
    layer_indx(cl).indx = indx;
    mi = (nanmean(WN1L(indx,:),2) - nanmean(WN1(indx,:),2))./(nanmean(WN1L(indx,:),2) + nanmean(WN1(indx,:),2));
    mi_run = (WN1(indx,1) - WN1(indx,4))./(WN1(indx,1) + WN1(indx,4));
    
    meanMI(cl) = nanmean(mi);
    semMI(cl) = sem(mi);
    nonnan_indx = find(~isnan(mi)==1);
    fs_layer(cl) = sum(fs2(nonnan_indx));
    rs_layer(cl) = sum(rs2(nonnan_indx));
    n_layer_VIP(cl,1) = sum(~isnan(mi));
    mi_cl = [ mi_cl; mi ones(length(mi),1)*cl fs2 rs2];
    
    
    figure(300); hold on
    plot(ones(length(mi),1)*cl, mi, 'k.');
    plot(cl, nanmean(mi), 'ro')
    
    meanMI_rs(cl) = nanmean(mi(rs2));
    semMI_rs(cl) = sem(mi(rs2));
    meanMI_fs(cl) = nanmean(mi(fs2));
    semMI_fs(cl) = sem(mi(fs2));
    
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
    all_MI = [all_MI; mi mi_run mi_sp mi_sp_run recs1(indx)' cells1(indx)' ones(length(indx),1)*cl fs2 rs2];
    
    n_layer_VIP(cl,2) = sum(~isnan(mi_sp));
    mi_sp_cl = [ mi_sp_cl; mi_sp ones(length(mi_sp),1)*cl fs2 rs2];
    nonnan_indx = find(~isnan(mi_sp)==1);
    fs_sp_layer(cl) = sum(fs2(nonnan_indx));
    rs_sp_layer(cl) = sum(rs2(nonnan_indx));
    
    figure(213); hold on
    subplot(1,4,cl)
    plot(mi_sp, mi_sp_run, 'k.')
    lsline; title(sprintf('Layer %d', cl))
    [rhoe,pe]=corr(mi, mi_run, 'Type', 'Spearman', 'rows', 'pairwise');
    [rhos,ps]=corr(mi_sp, mi_sp_run, 'Type', 'Spearman', 'rows', 'pairwise');
    title(sprintf('e:%.2f %.4f; sp:%.2f %.4f',rhoe,pe, rhos,ps));
    
    meanMI_sp(cl) = nanmean(mi_sp);
    semMI_sp(cl) = sem(mi_sp);
    
    meanMI_sp_rs(cl) = nanmean(mi_sp(rs2));
    semMI_sp_rs(cl) = sem(mi_sp(rs2));
    meanMI_sp_fs(cl) = nanmean(mi_sp(fs2));
    semMI_sp_fs(cl) = sem(mi_sp(fs2));
    
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
WN2 = nanmean(WN1,2);
errorbar([1:4], [nanmean(WN2(layer_indx(1).indx,:)) nanmean(WN2(layer_indx(2).indx)) nanmean(WN2(layer_indx(3).indx)) nanmean(WN2(layer_indx(4).indx))], ...
    [sem(WN2(layer_indx(1).indx)) sem(WN2(layer_indx(2).indx)) sem(WN2(layer_indx(3).indx)) sem(WN2(layer_indx(4).indx))], 'ro')
xticks([1:4]); xticklabels({'2/3', '4' '5' '6'})
xlabel('Cortical Layers')
ylabel('mean FR (Hz), evoked and spont')

WN2 = nanmean(SP1,2);
errorbar([1:4], [nanmean(WN2(layer_indx(1).indx,:)) nanmean(WN2(layer_indx(2).indx)) nanmean(WN2(layer_indx(3).indx)) nanmean(WN2(layer_indx(4).indx))], ...
    [sem(WN2(layer_indx(1).indx)) sem(WN2(layer_indx(2).indx)) sem(WN2(layer_indx(3).indx)) sem(WN2(layer_indx(4).indx))], 'ko')
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
bar([1 2], [nanmean(mi_sp(rs2)),  nanmean(mi(rs2))], 'BarWidth', .3)
bar([1.3, 2.3], [ nanmean(mi_sp(fs2)) nanmean(mi(fs2))], 'BarWidth', .3)
errorbar([1,1.3,2,2.3], [nanmean(mi_sp(rs2)), nanmean(mi_sp(fs2)), nanmean(mi(rs2)), nanmean(mi(fs2))], [sem(mi_sp(rs2)),sem(mi_sp(fs2)),sem(mi(rs2)) sem(mi(fs2))])
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

VIP_MI_rs_evoked = (nanmean(WN1L(rs2,indx1),2)- nanmean(WN1(rs2,indx1),2))./(nanmean(WN1L(rs2,indx1),2) + nanmean(WN1(rs2,indx1),2));
VIP_MI_fs_evoked = (nanmean(WN1L(fs2,indx1),2)- nanmean(WN1(fs2,indx1),2))./(nanmean(WN1L(fs2,indx1),2) + nanmean(WN1(fs2,indx1),2));
VIP_MI_rs_sp = (nanmean(SP1L(rs2,indx1),2)- nanmean(SP1(rs2,indx1),2))./(nanmean(SP1L(rs2,indx1),2) + nanmean(SP1(rs2,indx1),2));
VIP_MI_fs_sp = (nanmean(SP1L(fs2,indx1),2)- nanmean(SP1(fs2,indx1),2))./(nanmean(SP1L(fs2,indx1),2) + nanmean(SP1(fs2,indx1),2));

RUN_MI_rs_evoked = (WN1(rs2,1) - WN1(rs2,4))./(WN1(rs2,1) + WN1(rs2,4));
RUN_MI_fs_evoked = (WN1(fs2,1) - WN1(fs2,4))./(WN1(fs2,1) + WN1(fs2,4));
RUN_MI_rs_sp = (SP1(rs2,1) - SP1(rs2,4))./(SP1(rs2,1) + SP1(rs2,4));
RUN_MI_fs_sp = (SP1(fs2,1) - SP1(fs2,4))./(SP1(fs2,1) + SP1(fs2,4));
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
plot(b(rs2),a(rs2),'k.')
hold on; plot(b(fs2),a(fs2),'b.')
legend('rs', 'fs')
xlabel('running MI spont'); ylabel('running MI evoked')

% evoked and spont modulation are correlated for pupil
subplot(2,1,2); hold on
a = (WN1(:,2)-WN1(:,3))./(WN1(:,2)+WN1(:,3));
b = (SP1(:,2)-SP1(:,3))./(SP1(:,2)+SP1(:,3));
[rho,p]=corr(b, a, 'Type', 'Spearman', 'rows', 'pairwise')
plot(b(rs2),a(rs2),'k.')
hold on; plot(b(fs2),a(fs2),'b.')
legend('rs', 'fs')
xlabel('pupil MI spont'); ylabel('pupil MI evoked')


% running
WN_rs = (WN1(rs2,1) - WN1(rs2,4))./(WN1(rs2,1) + WN1(rs2,4));
WN_fs = (WN1(fs2,1) - WN1(fs2,4))./(WN1(fs2,1) + WN1(fs2,4));
WN_rsL = (WN1L(rs2,4) - WN1(rs2,4))./(WN1L(rs2,4) + WN1(rs2,4));
WN_fsL = (WN1L(fs2,4) - WN1(fs2,4))./(WN1L(fs2,4) + WN1(fs2,4));
WN_rs_plus = (WN1L(rs2,1) - WN1(rs2,4))./(WN1L(rs2,1) + WN1(rs2,4));
WN_fs_plus = (WN1L(fs2,1) - WN1(fs2,4))./(WN1L(fs2,1) + WN1(fs2,4));
WN_exp_rs = sum([WN_rs WN_rsL],2);
WN_exp_fs = sum([WN_fs WN_fsL],2);
[p,h,zstat] = ranksum(WN_rs_plus, WN_exp_rs);
[p,h,zstat] = ranksum(WN_fs_plus, WN_exp_fs);

SP_rs = (SP1(rs2,1) - SP1(rs2,4))./(SP1(rs2,1) + SP1(rs2,4));
SP_fs = (SP1(fs2,1) - SP1(fs2,4))./(SP1(fs2,1) + SP1(fs2,4));
SP_rsL = (nanmean(SP1L(rs2,2:3),2) - nanmean(SP1(rs2,2:3),2))./(nanmean(SP1L(rs2,2:3),2) + nanmean(SP1(rs2,2:3),2));
SP_fsL = (nanmean(SP1L(fs2,2:3),2) - nanmean(SP1(fs2,2:3),2))./(nanmean(SP1L(fs2,2:3),2) + nanmean(SP1(fs2,2:3),2));
SP_rs_plus = (SP1L(rs2,1) - SP1(rs2,4))./(SP1L(rs2,1) + SP1(rs2,4));
SP_fs_plus = (SP1L(fs2,1) - SP1(fs2,4))./(SP1L(fs2,1) + SP1(fs2,4));
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

WN_rs = (WN1(rs2,2) - WN1(rs2,3))./(WN1(rs2,2) + WN1(rs2,3));
WN_fs = (WN1(fs2,2) - WN1(fs2,3))./(WN1(fs2,2) + WN1(fs2,3));
WN_rsL = (WN1L(rs2,3) - WN1(rs2,3))./(WN1L(rs2,3) + WN1(rs2,3));
WN_fsL = (WN1L(fs2,3) - WN1(fs2,3))./(WN1L(fs2,3) + WN1(fs2,3));
WN_rs_plus = (WN1L(rs2,2) - WN1(rs2,3))./(WN1L(rs2,2) + WN1(rs2,3));
WN_fs_plus = (WN1L(fs2,2) - WN1(fs2,3))./(WN1L(fs2,2) + WN1(fs2,3));
% WN_exp_rs = [nanmean(WN_rs) nanmean(WN_rsL)];
% WN_exp_fs = [nanmean(WN_fs) nanmean(WN_fsL)];
WN_exp_rs = sum([WN_rs WN_rsL],2);
WN_exp_fs = sum([WN_fs WN_fsL],2);

[p,h,zstat] = ranksum(WN_rs_plus, WN_exp_rs);
[p,h,zstat] = ranksum(WN_fs_plus, WN_exp_fs);

SP_rs = (SP1(rs2,2) - SP1(rs2,3))./(SP1(rs2,2) + SP1(rs2,3));
SP_fs = (SP1(fs2,2) - SP1(fs2,3))./(SP1(fs2,2) + SP1(fs2,3));
SP_rsL = (SP1L(rs2,3) - SP1(rs2,3))./(SP1L(rs2,3) + SP1(rs2,3));
SP_fsL = (SP1L(fs2,3) - SP1(fs2,3))./(SP1L(fs2,3) + SP1(fs2,3));
SP_rs_plus = (SP1L(rs2,2) - SP1(rs2,3))./(SP1L(rs2,2) + SP1(rs2,3));
SP_fs_plus = (SP1L(fs2,2) - SP1(fs2,3))./(SP1L(fs2,2) + SP1(fs2,3));
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


