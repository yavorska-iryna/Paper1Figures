% The same analysis as
% WNModulationIndex_ThreeStates_by_CorticalLayer_epistatic but with FR
% compute FR WN responses
% three states: sit + small pupil, sit + large pupil, running  + large
% pupil
% by each cortical layer
%'WNdataLaserOFF_threeStateswithOFF_newState.mat' - .55 thresh for pupil,
%Moff = moff + poff
% WNdataLaserOFF_newState.mat -.55 thresh for pupil, moff = moff only

clear; close all; dbstop if error
variables_dir = 'C:\Users\lab\Resilio Sync\Paper1Figures\code\variables';
cd(variables_dir)
load('CellsQualityStats.mat')

load('WNdataLaserOFF_newState.mat'); % 'WNdataLaserOFF_newState1.mat has longer OFF responses
data = WNdataLaserOFF;
load('WNdataLaserON_newState.mat');
data1 = WNdataLaserON;
MeanpSSon = nan(length(data),1); MeanpSSoff = nan(length(data),1); MeanmSSon = nan(length(data),1); MeanmSSoff = nan(length(data),1);  maxFRall =[];
MeanpSSonL = nan(length(data),1); MeanpSSoffL =nan(length(data),1); MeanmSSonL = nan(length(data),1); MeanmSSoffL = nan(length(data),1);
color1 = [0.7 0.7 0.7]; color2 = [0.2 0.2 0.2];
load('CellsQualityStats.mat')
nreps_check_running = 7; %number of repetitions in each condition for comparison
nreps_check_pupil = 8;
meanWN = nan(length(data),4); meanWNL = nan(length(data),4);
CL = {[0 300], [301 400], [401 599], [600 2000]};
depths = nan(length(data),1); fs = zeros(length(data),1); Rs = zeros(length(data),1);
dirsM =[]; cellsM = 0;  dirsP =[]; cellsP = 0; recs = []; cells=[];
meanWN1M = nan(length(data),1);
meanWN1P = nan(length(data),1);
meanWN1ML = nan(length(data),1);
meanWN1PL = nan(length(data),1);

meanSPM = nan(length(data),1);
meanSPP = nan(length(data),1);
meanSPML = nan(length(data),1);
meanSPPL = nan(length(data),1);

meanWN1 = nan(length(data),1);
meanWN1L= nan(length(data),1);
meanSP= nan(length(data),1);
meanSPL = nan(length(data),1);

cdVIP;load('Silence_DistanceCorr_dirs.mat')
for cc =1:length(data)
    if data(cc).dir<30
        try
            meanSpikeCount = nanmean([data(cc).SpikeCountWN data(cc).SpikeCountSS data1(cc).SpikeCountWN data1(cc).SpikeCountSS]);
        catch
            meanSpikeCount = NaN;
        end
        if meanSpikeCount > 2  && CellsQualityStats.SNR(cc)>.5 && CellsQualityStats.uQ(cc)>10
            % exclude cells with very low spikecount, they usually have very large effects
            Spont = [data(cc).mSSon; data(cc).mSSoff; data(cc).pSSon]; %data(cc).mSSon; data(cc).mSSoff]; %Spont = nanmean(Spont);
            ON = [data(cc).mNON_on; data(cc).mNON_off; data(cc).pNON_on; data1(cc).mNON_on; data1(cc).mNON_off; data1(cc).pNON_on]; % data(cc).mNON_on; data(cc).mNON_off]; %ON = nanmean(ON);
            [rs, h, stats] = ranksum( ON(:), Spont(:));
            evoked(cc,1) = h;
            zstats(cc,1) = stats.zval;
            Sust = [data(cc).mNSustained_on; data(cc).mNSustained_off; data(cc).pNSustained_on; data1(cc).mNSustained_on; data1(cc).mNSustained_off; data1(cc).pNSustained_on];% data(cc).mNSustained_on; data(cc).mNSustained_off];
            [rs, h, stats] = ranksum( Sust(:), Spont(:));
            evoked(cc,2) = h;
            zstats(cc,2) = stats.zval;
            OFF = [data(cc).mNOFF_on; data(cc).mNOFF_off]; %data(cc).pNOFF_on data(cc).mNOFF_on; data(cc).mNOFF_off];
            [rs, h, stats] = ranksum( OFF(:), Spont(:));
            evoked(cc,3) = h;
            zstats(cc,3) = stats.zval;
            
            [rs, h, stats] = ranksum( data(cc).WNresponse, data(cc).SSresponse);
            evoked(cc,4) = h;
            zstats(cc,4) = stats.zval;
            
            if data(cc).nrepsWNMon > nreps_check_running && data(cc).nrepsWNMoff > nreps_check_running % && data1(cc).nrepsWNMon > nreps_check_running && data1(cc).nrepsWNMoff > nreps_check_running
                meanWN(cc,1)= nanmean([nanmean(nanmean(data(cc).mNOFF_on))]);% nanmean(nanmean(data(cc).mNSustained_on)) nanmean(nanmean(data(cc).mNOFF_on))]);
                meanWN(cc,4)= nanmean([nanmean(nanmean(data(cc).mNOFF_off))]);% nanmean(nanmean(data(cc).mNSustained_off)) nanmean(nanmean(data(cc).mNOFF_off))]);
                meanWN1M(cc) = nanmean(nanmean(data(cc).WNresponse));
                dirsM = [dirsM data(cc).dir];
                cellsM = cellsM+1;
            end
            
            if data(cc).nrepsWNPon > nreps_check_pupil && data(cc).nrepsWNPoff > nreps_check_pupil %&& data1(cc).nrepsWNPon > nreps_check_pupil && data1(cc).nrepsWNPoff > nreps_check_pupil
                meanWN(cc,2)= nanmean(nanmean(data(cc).pNOFF_on));%nanmean([nanmean(nanmean(data(cc).pNON_on)) nanmean(nanmean(data(cc).pNSustained_on)) nanmean(nanmean(data(cc).pNOFF_on))]);
                meanWN(cc,3)= nanmean(nanmean(data(cc).pNOFF_off));%nanmean([nanmean(nanmean(data(cc).pNON_off))] nanmean(nanmean(data(cc).pNSustained_off)) nanmean(nanmean(data(cc).pNOFF_off))]);
                meanWN1P(cc) = nanmean(nanmean(data(cc).WNresponse));
                dirsP = [dirsP data(cc).dir];
                cellsP = cellsP+1;
            end
            
            if data1(cc).nrepsWNMon > nreps_check_running && data1(cc).nrepsWNMoff > nreps_check_running %&& data(cc).nrepsWNMon > nreps_check_running && data(cc).nrepsWNMoff > nreps_check_running
                meanWNL(cc,1)= nanmean(nanmean(data1(cc).mNOFF_on));%nanmean([nanmean(nanmean(data1(cc).mNON_on)) nanmean(nanmean(data1(cc).mNSustained_on)) nanmean(nanmean(data1(cc).mNOFF_on))]);
                meanWNL(cc,4)= nanmean(nanmean(data1(cc).mNOFF_off));%nanmean([nanmean(nanmean(data1(cc).mNON_off)) nanmean(nanmean(data1(cc).mNSustained_off)) nanmean(nanmean(data1(cc).mNOFF_off))]);
                meanWN1ML(cc) = nanmean(nanmean(data1(cc).WNresponse));
            end
            
            if data1(cc).nrepsWNPon > nreps_check_pupil && data1(cc).nrepsWNPoff > nreps_check_pupil% && data(cc).nrepsWNPon > nreps_check_pupil && data(cc).nrepsWNPoff > nreps_check_pupil
                meanWNL(cc,2)= nanmean(nanmean(data1(cc).pNOFF_on));%nanmean([nanmean(nanmean(data1(cc).pNON_on)) nanmean(nanmean(data1(cc).pNSustained_on)) nanmean(nanmean(data1(cc).pNOFF_on))]);
                meanWNL(cc,3)= nanmean(nanmean(data1(cc).pNOFF_off));%nanmean([nanmean(nanmean(data1(cc).pNON_off)) nanmean(nanmean(data1(cc).pNSustained_off)) nanmean(nanmean(data1(cc).pNOFF_off))]);
                meanWN1PL(cc) = nanmean(nanmean(data1(cc).WNresponse));
            end
            
            
            % spont activity
             if data(cc).nrepsSSMon > nreps_check_running && data(cc).nrepsSSMoff > nreps_check_running && data1(cc).nrepsSSMon > nreps_check_running && data1(cc).nrepsSSMoff > nreps_check_running
                MeanmSSon(cc) = nanmean(nanmean(data(cc).mSSon)); MeanmSSoff(cc) = nanmean(nanmean(data(cc).mSSoff));
                meanSPM(cc) = nanmean(nanmean((data(cc).SSresponse)));
            end

            if data(cc).nrepsSSPon > nreps_check_pupil && data(cc).nrepsSSPoff > nreps_check_pupil && data1(cc).nrepsSSPon > nreps_check_pupil && data1(cc).nrepsSSPoff > nreps_check_pupil
                MeanpSSon(cc) =  nanmean(nanmean(data(cc).pSSon)); MeanpSSoff(cc) =  nanmean(nanmean(data(cc).pSSoff));
                meanSPP(cc) = nanmean(nanmean((data(cc).SSresponse)));
                
            end
            if data1(cc).nrepsSSMon > nreps_check_running && data1(cc).nrepsSSMoff > nreps_check_running && data(cc).nrepsSSMon > nreps_check_running && data(cc).nrepsSSMoff > nreps_check_running
                MeanmSSonL(cc) = nanmean(nanmean(data1(cc).mSSon)); MeanmSSoffL(cc) = nanmean(nanmean(data1(cc).mSSoff));
                meanSPML(cc) = nanmean(nanmean((data1(cc).SSresponse)));
            end
            
            if data1(cc).nrepsSSPon > nreps_check_pupil && data1(cc).nrepsSSPoff > nreps_check_pupil && data(cc).nrepsSSPon > nreps_check_pupil && data(cc).nrepsSSPoff > nreps_check_pupil
                MeanpSSonL(cc) = nanmean(nanmean(data1(cc).pSSon)); MeanpSSoffL(cc) = nanmean(nanmean(data1(cc).pSSoff));
                meanSPPL(cc) = nanmean(nanmean((data1(cc).SSresponse)));
            end
           
            % all cells with laser
            meanWN1(cc) = nanmean(nanmean(data(cc).WNresponse));
            meanWN1L(cc) = nanmean(nanmean(data1(cc).WNresponse));
            meanSP(cc) = nanmean(nanmean((data(cc).SSresponse)));
            meanSPL(cc) = nanmean(nanmean((data1(cc).SSresponse)));
            
        else
            evoked(cc,:,:,:) = [0 0 0 0];
            zstats(cc, :,:,:) = [NaN NaN NaN NaN];
            
        end
        
        if isempty(data(cc).maxFR) || data(cc).maxFR == 0
            data(cc).maxFR= NaN;
        end
        maxFRall = [maxFRall data(cc).maxFR];
        try
            depths(cc) = data(cc).depth;
            if data(cc).width < 0.79
                fs(cc) = 1;
            else
                Rs(cc) = 1;
            end
        catch
            fs(cc) = 0;
            Rs(cc) = 0;
        end
        
    end
    recs = [recs data(cc).dir];
    cells = [cells data(cc).cell];
end
allDirsM = unique(dirsM);
allDirsP = unique(dirsP);
evoked2 = evoked(:,3) & zstats(:,3)>0; 

SP = [MeanmSSon MeanpSSon MeanpSSoff MeanmSSoff];
SPL = [MeanmSSonL MeanpSSonL MeanpSSoffL MeanmSSoffL];
meanSP1 = meanSP(evoked2); %all evoked sp
meanSPL1 = meanSPL(evoked2); %laser
meanWN2 = meanWN1(evoked2); % all wn
meanWNL2 = meanWN1L(evoked2);

meanWN2M = meanWN1M(evoked2); %moves
meanWN2P = meanWN1P(evoked2);
meanWN2ML = meanWN1ML(evoked2);
meanWN2PL = meanWN1PL(evoked2);

meanSP2P = meanSPP(evoked2);
meanSP2M = meanSPM(evoked2);
meanSP2ML = meanSPML(evoked2);
meanSP2PL = meanSPPL(evoked2);


% run_mi =  (WN1(:,1) -WN1(:,4))./(WN1(:,1) + WN1(:,4));
% laser_mi = (meanWNL2 - meanWN2)./ (meanWNL2 +meanWN2);
% run_mi_sp =  (SP1(:,1) -SP1(:,4))./(SP1(:,1) + SP1(:,4));
% laser_mi_sp = (meanSPL1 - meanSP1)./ (meanSPL1 +meanSP1);



% state dependent data
WN1 = meanWN(evoked2,1:4);
WN1L = meanWNL(evoked2,1:4);
SP1 = SP(evoked2,:);
SP1L = SPL(evoked2,:);
depths1 = depths(evoked2);

fs1 = fs(evoked2); % only evoked fs and rs
rs1 = Rs(evoked2);
rs1 = logical(rs1); fs1 = logical(fs1);
recs1 = recs(evoked2); cells1 = cells(evoked2);

% modulation_indx1_rs = (WN1(rs1,1) -WN1(rs1,4))./(WN1(rs1,1) + WN1(rs1,4)); %running
% modulation_indx2_rs = (WN1(rs1,2) -WN1(rs1,3))./(WN1(rs1,2) + WN1(rs1,3)); %large pupil + sit
% modulation_indx1_fs = (WN1(fs1,1) -WN1(fs1,4))./(WN1(fs1,1) + WN1(fs1,4)); %running
% modulation_indx2_fs = (WN1(fs1,2) -WN1(fs1,3))./(WN1(fs1,2) + WN1(fs1,3)); %large pupil + sit
% 
% modulation_indx1_sp_rs = (SP1(rs1,1) -SP1(rs1,4))./(SP1(rs1,1) + SP1(rs1,4)); %running
% modulation_indx2_sp_rs = (SP1(rs1,2) -SP1(rs1,3))./(SP1(rs1,2) + SP1(rs1,3)); %large pupil + sit
% modulation_indx1_sp_fs = (SP1(fs1,1) -SP1(fs1,4))./(SP1(fs1,1) + SP1(fs1,4)); %running
% modulation_indx2_sp_fs = (SP1(fs1,2) -SP1(fs1,3))./(SP1(fs1,2) + SP1(fs1,3)); %large pupil + sit

firing_rate1_rs = WN1(rs1,1) - WN1(rs1,4); %running
firing_rate2_rs = WN1(rs1,2) - WN1(rs1,3); %pupil
firing_rate1_fs =  WN1(fs1,1) - WN1(fs1,4); %running
firing_rate2_fs = WN1(fs1,2) - WN1(fs1,3); %pupil

firing_rate1_sp_rs = SP1(rs1,1) - SP1(rs1,4); %running
firing_rate2_sp_rs = WN1(rs1,2) - SP1(rs1,3); %pupil
firing_rate1_sp_fs =  SP1(fs1,1) - SP1(fs1,4); %running
firing_rate2_sp_fs = SP1(fs1,2) -  SP1(fs1,3); %pupil

figure; 
subplot(1,2,1); hold on
plot(firing_rate1_rs, firing_rate1_sp_rs, 'k.')
plot(firing_rate1_fs, firing_rate1_sp_fs, 'b.')
xlabel('evoked FR change'); ylabel('spont FR change')
plot([0 0], [min(firing_rate1_sp_rs) max(firing_rate1_sp_rs)], 'r--');
plot([min(firing_rate1_rs) max(firing_rate1_rs)], [0 0], 'r--');
legend({'rs cells', 'fs cells'})
title('running')
subplot(1,2,2); hold on
plot(firing_rate2_rs, firing_rate2_sp_rs, 'k.')
plot(firing_rate2_fs, firing_rate2_sp_fs, 'b.')
xlabel('evoked FR change'); ylabel('spont FR change')
plot([0 0], [min(firing_rate2_sp_rs) max(firing_rate2_sp_rs)], 'r--');
plot([min(firing_rate2_rs) max(firing_rate2_rs)], [0 0], 'r--');
legend({'rs cells', 'fs cells'})
title('pupil')

% state means by RS FS WN SP
figure(50); hold on
bar([0.8 1.8], [nanmean(firing_rate1_sp_rs) nanmean(firing_rate2_sp_rs)], 'BarWidth', .1)
bar([0.9 1.9], [nanmean(firing_rate1_rs) nanmean(firing_rate2_rs)], 'BarWidth', .1)
bar([1  2 ], [ nanmean(firing_rate1_sp_fs) nanmean(firing_rate2_sp_fs ) ], 'BarWidth', .1)
bar([1.1 2.1], [nanmean(firing_rate1_fs) nanmean(firing_rate2_fs)], 'BarWidth', .1)
xticks([1:2])
xticklabels({'run', 'sit + large pupil'})
ylabel('Firing Rate Change')
errorbar( [0.8 0.9 1 1.1 1.8 1.9 2 2.1],[nanmean(firing_rate1_sp_rs) nanmean(firing_rate1_rs) nanmean(firing_rate1_sp_fs) nanmean(firing_rate1_fs)  ...
    nanmean(firing_rate2_sp_rs) nanmean(firing_rate2_rs) nanmean(firing_rate2_sp_fs) nanmean(firing_rate2_fs)], ...
    [sem(firing_rate1_sp_rs) sem(firing_rate1_rs) sem(firing_rate1_sp_fs) sem(firing_rate1_fs)  ...
    sem(firing_rate2_sp_rs) sem(firing_rate2_rs) sem(firing_rate2_sp_fs) sem(firing_rate2_fs)])
legend({'spont rs', 'evoked rs', 'spont fs', 'evoked fs'})
%%
recs1= (recs(evoked2));
recs_rs= (recs1(rs1));
recs_fs=(recs1(fs1));
figure; plot( recs_rs,firing_rate2_sp_rs, 'ko');
nanindx = find(isnan(firing_rate2_sp_rs)==1);
x=firing_rate2_sp_rs;
[p,tbl1,stats] = kruskalwallis(x, recs_rs);
c = multcompare(stats);
title('Spont RS')
x=firing_rate2_rs;
[p,tbl1,stats] = kruskalwallis(x, recs_rs);
c = multcompare(stats);
title('Evoked RS')
%[m,s]=grpstats(x,recs_rs,{'mean','sem'})
x=firing_rate2_sp_fs;
[p,tbl1,stats] = kruskalwallis(x, recs_fs);
c = multcompare(stats);
title('Spont FS')
%[m,s]=grpstats(x,recs_rs,{'mean','sem'})
r =unique(recs_rs);
%%
FIRING_RATE1 = []; FIRING_RATE2 =[];  FIRING_RATE1L = []; FIRING_RATE2L = [];
for cl = 1:length(CL)
    layer = CL{cl};
    indx = find(depths1 >layer(1) & depths1 <layer(2));
    
    firing_rate1 = WN1(indx,1) -WN1(indx,4); %running
    firing_rate2 = WN1(indx,2) -WN1(indx,3); %large pupil + sit
    
    firing_rate1L = WN1L(indx,1) -WN1L(indx,3); %running + laser
    firing_rate2L = WN1L(indx,2) -WN1L(indx,3); %large pupil + sit +laser
    
    firing_rate1Laser = nanmean(WN1L(indx,2:2:4),2) - nanmean(WN1(indx,2:2:4),2); %running + laser
    firing_rate2Laser = nanmean(WN1L(indx,2:3),2) - nanmean(WN1(indx,2:3),2); %pupil + laser
    
    fs2 = fs1(indx); fs2 = logical(fs2);
    rs2 = rs1(indx); rs2 = logical(rs2);
    n_layer_evoked_rsfs(cl,:) = [sum(rs2), sum(fs2)];
    
    recs2 = recs1(indx);
    cells2 = cells1(indx);
    %stats
    FIRING_RATE1 = [FIRING_RATE1; firing_rate1 ones(length(indx),1)*cl];
    FIRING_RATE2 = [FIRING_RATE2; firing_rate2 ones(length(indx),1)*cl];
    FIRING_RATE1L = [FIRING_RATE1L; firing_rate1L ones(length(indx),1)*cl];
    FIRING_RATE2L = [FIRING_RATE2L; firing_rate2L ones(length(indx),1)*cl];
    %
    meanMI1_Laser(cl) = nanmean(firing_rate1Laser);
    meanMI2_Laser(cl) = nanmean(firing_rate2Laser);
    semMI1_Laser(cl) = sem(firing_rate1Laser);
    semMI2_Laser(cl) = sem(firing_rate2Laser);
    
    meanMI1(cl) = nanmean(firing_rate1);
    n_layer_evoked(cl,1) = length(find(~isnan(firing_rate1)==1));
    n_layer_evokedL(cl) = length(find(~isnan(firing_rate1L)==1));
    meanMI2(cl) = nanmean(firing_rate2);
    n_layer_evoked(cl,2) = length(find(~isnan(firing_rate2)==1));
    semMI1(cl) = sem(firing_rate1);
    semMI2(cl) = sem(firing_rate2);
    meanMI1_fs(cl) = nanmean(firing_rate1(fs2));
    meanMI2_fs(cl) = nanmean(firing_rate2(fs2));
    semMI1_fs(cl) = sem(firing_rate1(fs2));
    semMI2_fs(cl) = sem(firing_rate2(fs2));
    meanMI1_rs(cl) = nanmean(firing_rate1(rs2));
    meanMI2_rs(cl) = nanmean(firing_rate2(rs2));
    semMI1_rs(cl) = sem(firing_rate1(rs2));
    semMI2_rs(cl) = sem(firing_rate2(rs2));
    
    % addative effect
    FR_evoked_plus = WN1L(indx,1) - WN1(indx,4);
    meanMI_evoked_plus(cl) = nanmean(FR_evoked_plus);
    semMI_evoked_plus(cl) = sem(FR_evoked_plus);
    n_layer_evoked_plus(cl) = length(find(~isnan(FR_evoked_plus)==1));
    FR_evoked_plus = WN1L(indx,2) - WN1(indx,3);
    meanMI2_evoked_plus(cl) = nanmean(FR_evoked_plus);
    semMI2_evoked_plus(cl) = sem(FR_evoked_plus);
    
    %laser on
    meanMI1L(cl) = nanmean(firing_rate1L);
    meanMI2L(cl) = nanmean(firing_rate2L);
    semMI1L(cl) = sem(firing_rate1L);
    semMI2L(cl) = sem(firing_rate2L);
    meanMI1_fsL(cl) = nanmean(firing_rate1L(fs2));
    meanMI2_fsL(cl) = nanmean(firing_rate2L(fs2));
    semMI1_fsL(cl) = sem(firing_rate1L(fs2));
    semMI2_fsL(cl) = sem(firing_rate2L(fs2));
    meanMI1_rsL(cl) = nanmean(firing_rate1L(rs2));
    meanMI2_rsL(cl) = nanmean(firing_rate2L(rs2));
    semMI1_rsL(cl) = sem(firing_rate1L(rs2));
    semMI2_rsL(cl) = sem(firing_rate2L(rs2));
    
    % spont act
    indx = find(depths1 >layer(1) & depths1 <layer(2));
    fs2 = fs1(indx); fs2 = logical(fs2);
    rs2 = rs1(indx); rs2 = logical(rs2);
    firing_rate1_sp = (SP1(indx,1) -SP1(indx,4)); %running + laser
    firing_rate2_sp = (SP1(indx,2) - SP1(indx,3)); %large pupil + sit +laser
    firing_rate1_spL = (SP1L(indx,1) -SP1L(indx,4)); %running + laser
    firing_rate2_spL = (SP1L(indx,2) - SP1L(indx,3)); %large pupil + sit +laser
    firing_rate1_spLaser = (nanmean(SP1L(indx,2:2:4),2) - nanmean(SP1(indx,2:2:4),2)); %running + laser
    firing_rate2_spLaser = (nanmean(SP1L(indx,2:3),2) - nanmean(SP1(indx,2:3),2)); %running + laser
    
    MI_sp_plus = (SP1L(indx,1) - SP1(indx,4))./(SP1L(indx,1) + SP1(indx,4));
    meanMI_sp_plus(cl) = nanmean(MI_sp_plus);
    semMI_sp_plus(cl) = sem(MI_sp_plus);
    n_layer_sp_plus(cl) = length(find(~isnan(MI_sp_plus)==1));
    MI2_sp_plus = (SP1L(indx,2) - SP1(indx,3))./(SP1L(indx,2) + SP1(indx,3));
    meanMI2_sp_plus(cl) = nanmean(MI_sp_plus);
    semMI2_sp_plus(cl) = sem(MI_sp_plus);
    
    meanMI1_sp(cl) = nanmean(firing_rate1_sp);
    n_layer_sp(cl,1) = length(find(~isnan(firing_rate1_sp)==1)); %number of cells in each layer
    n_layer_spL(cl) = length(find(~isnan(firing_rate1_spL)==1));
    meanMI1_spL(cl) = nanmean(firing_rate1_spL);
    meanMI1_spLaser(cl) = nanmean(firing_rate1_spLaser);
    semMI1_spLaser(cl) = sem(firing_rate1_spLaser);
    meanMI2_spLaser(cl) = nanmean(firing_rate2_spLaser);
    semMI2_spLaser(cl) = sem(firing_rate2_spLaser);
    
    meanMI2_sp(cl) = nanmean(firing_rate2_sp);
    n_layer_sp(cl,2) = length(find(~isnan(firing_rate2_sp)==1)); %number of cells in each layer
    meanMI2_spL(cl) = nanmean(firing_rate2_spL);
    semMI1_sp(cl) = sem(firing_rate1_sp);
    semMI2_sp(cl) = sem(firing_rate2_sp);
    semMI1_spL(cl) = sem(firing_rate1_spL);
    semMI2_spL(cl) = sem(firing_rate2_spL);
    
    meanMI1_sp_fs(cl) = nanmean(firing_rate1_sp(fs2));
    meanMI1_spL_fs(cl) = nanmean(firing_rate1_spL(fs2));
    meanMI2_sp_fs(cl) = nanmean(firing_rate2_sp(fs2));
    meanMI2_spL_fs(cl) = nanmean(firing_rate2_spL(fs2));
    semMI1_sp_fs(cl) = sem(firing_rate1_sp(fs2));
    semMI2_sp_fs(cl) = sem(firing_rate2_sp(fs2));
    semMI1_spL_fs(cl) = sem(firing_rate1_spL(fs2));
    semMI2_spL_fs(cl) = sem(firing_rate2_spL(fs2));
    
    meanMI1_sp_rs(cl) = nanmean(firing_rate1_sp(rs2));
    meanMI1_spL_rs(cl) = nanmean(firing_rate1_spL(rs2));
    meanMI2_sp_rs(cl) = nanmean(firing_rate2_sp(rs2));
    meanMI2_spL_rs(cl) = nanmean(firing_rate2_spL(rs2));
    semMI1_sp_rs(cl) = sem(firing_rate1_sp(rs2));
    semMI2_sp_rs(cl) = sem(firing_rate2_sp(rs2));
    semMI1_spL_rs(cl) = sem(firing_rate1_spL(rs2));
    semMI2_spL_rs(cl) = sem(firing_rate2_spL(rs2));
    
end

% state by layer without layer
figure(101); subplot(2,1,1); hold on
errorbar([1:4], meanMI1, semMI1, 'ko');
errorbar([1.2:4.2], meanMI2, semMI2, 'ro');
xlabel('Cortical Layer')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Firing Rate Change'); title('WN response')
% [p,tbl1,stats] = kruskalwallis(mi, layers);
% c = multcompare(stats);
subplot(2,1,2); hold on
errorbar([1:4], meanMI1_sp, semMI1_sp, 'ko');
errorbar([1.2:4.2], meanMI2_sp, semMI2_sp, 'ro');

xlabel('Cortical Layer')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Firing Rate Change'); title('Spont activity');

%State by layer with laser
figure(101); subplot(2,1,2); hold on
errorbar([1.1:4.1], meanMI1_spL, semMI1_spL, 'k>');
errorbar([1.3:4.3], meanMI2_spL, semMI2_spL, 'b>');
plot([0 5], [0 0], 'k--')
legend({ 'run laser off', 'sit + large pupil laser off', 'run laser on', 'sit + large pupil laser on'})
xlabel('Cortical Layer')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Firing Rate Change'); title('Spont activity - VIP laser on')
subplot(2,1,1); hold on
errorbar([1.1:4.1], meanMI1L, semMI1L, 'k>');
errorbar([1.3:4.3], meanMI2L, semMI2L, 'b>');
plot([0 5], [0 0], 'k--')
legend({ 'run laser off', 'sit + large pupil laser off', 'run laser on', 'sit + large pupil laser on'})
xlabel('Cortical Layer')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Firing Rate Change'); title('Evoked activity - VIP laser on')


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
ylabel('Firing Rate Change'); title('Spont Activity')
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
ylabel('Firing Rate Change'); title('WN responses')


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
ylabel('Firing Rate Change - Run'); title('Evoked responses')
figure(103); subplot(2,1,2); hold on
errorbar([1:4], meanMI1_sp, semMI1_sp, 'ko');
errorbar([1.1:4.1], meanMI1_spLaser, semMI1_spLaser, 'bo');
errorbar([1.2:4.2], meanMI_sp_plus, semMI_sp_plus, 'go');
plot([0 5], [0 0], 'k--')
legend({ 'run', 'laser' 'run + laser' })
xlabel('Cortical Layer')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Firing Rate Change - Run'); title('Spont responses')
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
ylabel('Firing Rate Change - Pupil'); title('Evoked responses')
subplot(2,1,2); hold on
errorbar([1:4], meanMI2_sp, semMI2_sp, 'ko');
errorbar([1.1:4.1], meanMI2_spLaser, semMI2_spLaser, 'bo');
errorbar([1.2:4.2], meanMI2_sp_plus, semMI2_sp_plus, 'go');
plot([0 5], [0 0], 'k--')
legend({ 'pupil', 'laser' 'pupil + laser' }); xlabel('Cortical Layer')
xticks([1:4]); xlim([0 5])
xticklabels({'2/3', '4', '5', '6'});
ylabel('Firing Rate Change - Pupil'); title('Spont responses')

fs1 = logical(fs1); rs1 = logical(rs1);
MI1_rs = (WN1(rs1,1) - WN1(rs1,4));
MI2_rs = (WN1(rs1,2) - WN1(rs1,3));
MI1_fs = (WN1(fs1,1) - WN1(fs1,4));
MI2_fs = (WN1(fs1,2) - WN1(fs1,3));

SP1_rs = (SP1(rs1,1) - SP1(rs1,4));
SP2_rs = (SP1(rs1,2) - SP1(rs1,3));
SP1_fs = (SP1(fs1,1) - SP1(fs1,4));
SP2_fs = (SP1(fs1,2) - SP1(fs1,3));

%laser VIP
MI1_rsL = (WN1L(rs1,1) - WN1L(rs1,4));
MI2_rsL = (WN1L(rs1,2) - WN1L(rs1,3));
MI1_fsL = (WN1L(fs1,1) - WN1L(fs1,4));
MI2_fsL = (WN1L(fs1,2) - WN1L(fs1,3));
SP1_rsL = (SP1L(rs1,1) - SP1L(rs1,4));
SP2_rsL = (SP1L(rs1,2) - SP1L(rs1,3));
SP1_fsL = (SP1L(fs1,1) - SP1L(fs1,4));
SP2_fsL = (SP1L(fs1,2) - SP1L(fs1,3));


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
ylabel('Firing Rate Change')

firing_rate1_rs = (WN1L(rs1,1) -WN1(rs1,1)); %running
firing_rate2_rs = (WN1L(rs1,2) -WN1(rs1,2)); %large pupil + sit
modulation_indx3_rs = (WN1L(rs1,4) -WN1(rs1,4)); %small pupil + sit
firing_rate1_fs = (WN1L(fs1,1) -WN1(fs1,1)); %running
firing_rate2_fs = (WN1L(fs1,2) -WN1(fs1,2)); %large pupil + sit
modulation_indx3_fs = (WN1L(fs1,4) -WN1(fs1,4)); %small pupil + sit

firing_rate1_sp_rs = (SP1L(rs1,1) -SP1(rs1,1)); %running
firing_rate2_sp_rs = (SP1L(rs1,2) -SP1(rs1,2)); %large pupil + sit
modulation_indx3_sp_rs = (SP1L(rs1,4) -SP1(rs1,4)); %small pupil + sit
firing_rate1_sp_fs = (SP1L(fs1,1) -SP1(fs1,1)); %running
firing_rate2_sp_fs = (SP1L(fs1,2) -SP1(fs1,2)); %large pupil + sit
modulation_indx3_sp_fs = (SP1L(fs1,4) -SP1(fs1,3)); %small pupil + sit
% VIP effect in ech state
figure(210); hold on
bar([.8  1.8  2.8 ], [nanmean(modulation_indx3_sp_rs) nanmean(firing_rate2_sp_rs) nanmean(firing_rate1_sp_rs)], 'BarWidth', .1)
bar([.9  1.9  2.9 ], [nanmean(modulation_indx3_sp_fs) nanmean(firing_rate2_sp_fs) nanmean(firing_rate1_sp_fs)], 'BarWidth', .1)
bar([1  2  3 ], [nanmean(modulation_indx3_rs) nanmean(firing_rate2_rs) nanmean(firing_rate1_rs)], 'BarWidth', .1)
bar([1.1  2.1 3.1 ], [nanmean(modulation_indx3_fs) nanmean(firing_rate2_fs) nanmean(firing_rate1_fs)], 'BarWidth', .1)

errorbar([.8 .9 1. 1.1 1.8 1.9 2 2.1 2.8 2.9 3 3.1], [nanmean(modulation_indx3_sp_rs) nanmean(modulation_indx3_sp_fs) nanmean(modulation_indx3_rs) nanmean(modulation_indx3_fs) ...
    nanmean(firing_rate2_sp_rs) nanmean(firing_rate2_sp_fs)  nanmean(firing_rate2_rs)  nanmean(firing_rate2_fs) ...
    nanmean(firing_rate1_sp_rs) nanmean(firing_rate1_sp_fs) nanmean(firing_rate1_rs) nanmean(firing_rate1_fs)],...
    [sem(modulation_indx3_sp_rs) sem(modulation_indx3_sp_fs) sem(modulation_indx3_rs) sem(modulation_indx3_fs) ...
    sem(firing_rate2_sp_rs) sem(firing_rate2_sp_fs)  sem(firing_rate2_rs)  sem(firing_rate2_fs) ...
    sem(firing_rate1_sp_rs) sem(firing_rate1_sp_fs) sem(firing_rate1_rs) sem(firing_rate1_fs)])
xticks([1:3])
xticklabels({'sit + small pupil', 'sit + large pupil', 'run'})
ylabel('Firing Rate Change'); title('VIP activation')
legend({'spont rs', 'spont fs', 'evoked rs', 'evoked fs'})

figure(213); all_MI=[];
mi_cl = []; mi_sp_cl = [];
for cl = 1:length(CL)
    layer = CL{cl};
    indx = find(depths1 >layer(1) & depths1 <layer(2));
    fs2 = fs1(indx); fs2 = logical(fs2);
    rs2 = rs1(indx); rs2 = logical(rs2);
    layer_indx(cl).indx = indx;
    
    mi_laser_m = (meanWN2ML(indx) - meanWN2M(indx));
    mi_laser_p = (meanWN2PL(indx) - meanWN2P(indx));
    mi_laser = (meanWNL2(indx) - meanWN2(indx));
    mi_run = (WN1(indx,1) - WN1(indx,4));
    
    meanMI(cl) = nanmean(mi_laser);
    semMI(cl) = sem(mi_laser);
    
    meanMI_laser_m(cl) = nanmean(mi_laser_m);
    semMI_laser_m(cl) = sem(mi_laser_m);
    
    meanMI_laser_p(cl) = nanmean(mi_laser_p);
    semMI_laser_p(cl) = sem(mi_laser_p);
    
    nonnan_indx = find(~isnan(mi_laser)==1);
    
    fs_layer(cl) = sum(fs2(nonnan_indx));
    rs_layer(cl) = sum(rs2(nonnan_indx));
    
    n_layer_VIP(cl,1) = sum(~isnan(mi_laser));
    mi_cl = [ mi_cl; mi_laser ones(length(mi_laser),1)*cl fs2 rs2];
    
    
    meanMI_rs(cl) = nanmean(mi_laser(rs2));
    semMI_rs(cl) = sem(mi_laser(rs2));
    meanMI_fs(cl) = nanmean(mi_laser(fs2));
    semMI_fs(cl) = sem(mi_laser(fs2));
    
    firing_rate1 = WN1L(indx,1) -WN1(indx,1); %running
    firing_rate2 = WN1L(indx,2) -WN1(indx,2); %large pupil + sit
    firing_rate3 = WN1L(indx,4) -WN1(indx,4);
    
    
    meanMI1(cl) = nanmean(firing_rate1);
    semMI1(cl) = sem(firing_rate1);
    meanMI2(cl) = nanmean(firing_rate2);
    semMI2(cl) = sem(firing_rate2);
    meanMI3(cl) = nanmean(firing_rate3);
    semMI3(cl) = sem(firing_rate3);
    
    
    %spont
    mi_sp_run = (SP1(indx,1) - SP1(indx,4));
    mi_sp_laser = (meanSPL1(indx) - meanSP1(indx));
    mi_sp_laser_m = (meanSP2ML(indx) - meanSP2M(indx));
    mi_sp_laser_p = (meanSP2PL(indx) - meanSP2P(indx));
    
    all_MI = [all_MI; mi_laser mi_sp_laser mi_laser_m mi_sp_laser_m mi_laser_p mi_sp_laser_p recs1(indx)' cells1(indx)' ones(length(indx),1)*cl fs2 rs2];
    
    n_layer_VIP(cl,2) = sum(~isnan(mi_sp_laser));
    mi_sp_cl = [ mi_sp_cl; mi_sp_laser ones(length(mi_sp_laser),1)*cl fs2 rs2];
    nonnan_indx = find(~isnan(mi_sp_laser)==1);
    fs_sp_layer(cl) = sum(fs2(nonnan_indx));
    rs_sp_layer(cl) = sum(rs2(nonnan_indx));
    
    
    meanMI_sp(cl) = nanmean(mi_sp_laser);
    semMI_sp(cl) = sem(mi_sp_laser);
    
    meanMI_sp_laser_m(cl) = nanmean(mi_sp_laser_m);
    semMI_sp_laser_m(cl) = sem(mi_sp_laser_m);
    meanMI_sp_laser_p(cl) = nanmean(mi_sp_laser_p);
    semMI_sp_laser_p(cl) = sem(mi_sp_laser_p);
    
    
    meanMI_sp_rs(cl) = nanmean(mi_sp_laser(rs2));
    semMI_sp_rs(cl) = sem(mi_sp_laser(rs2));
    meanMI_sp_fs(cl) = nanmean(mi_sp_laser(fs2));
    semMI_sp_fs(cl) = sem(mi_sp_laser(fs2));
    
    firing_rate1_sp = (SP1L(indx,1) -SP1(indx,1)); %running
    firing_rate2_sp = (SP1L(indx,2) -SP1(indx,2)); %large pupil + sit
    firing_rate3_sp = (SP1L(indx,4) -SP1(indx,4)); %small pupil pupil + sit
    
    meanMI1_sp(cl) = nanmean(firing_rate1_sp);
    semMI1_sp(cl) = sem(firing_rate1_sp);
    meanMI2_sp(cl) = nanmean(firing_rate2_sp);
    semMI2_sp(cl) = sem(firing_rate2_sp);
    meanMI3_sp(cl) = nanmean(firing_rate3_sp);
    semMI3_sp(cl) = sem(firing_rate3_sp);
end
figure; hold on
errorbar([1:4], [nanmean(meanWN2(layer_indx(1).indx)) nanmean(meanWN2(layer_indx(2).indx)) nanmean(meanWN2(layer_indx(3).indx)) nanmean(meanWN2(layer_indx(4).indx))], ...
    [sem(meanWN2(layer_indx(1).indx)) sem(meanWN2(layer_indx(2).indx)) sem(meanWN2(layer_indx(3).indx)) sem(meanWN2(layer_indx(4).indx))], 'ro')
xticks([1:4]); xticklabels({'2/3', '4' '5' '6'})
xlabel('Cortical Layers')
ylabel('mean FR (Hz), evoked and spont')

SP2 = meanSP1;
errorbar([1:4], [nanmean(SP2(layer_indx(1).indx)) nanmean(SP2(layer_indx(2).indx)) nanmean(SP2(layer_indx(3).indx)) nanmean(SP2(layer_indx(4).indx))], ...
    [sem(SP2(layer_indx(1).indx)) sem(SP2(layer_indx(2).indx)) sem(SP2(layer_indx(3).indx)) sem(SP2(layer_indx(4).indx))], 'ko')
xticks([1:4]); xticklabels({'2/3', '4' '5' '6'})
xlabel('Cortical Layers')

mi_laser = (meanWNL2 - meanWN2)./(meanWNL2 + meanWN2);
mi_sp_laser = (meanSPL1 - meanSP1)./(meanSPL1 + meanSP1);

mi_laser_m = (meanWN2ML - meanWN2M)./(meanWN2ML + meanWN2M);
mi_sp_laser_m = (meanSP2ML - meanSP2M)./(meanSP2ML + meanSP2M);

mi_laser_p = (meanWN2PL - meanWN2P)./(meanWN2PL + meanWN2P);
mi_sp_laser_p = (meanSP2PL - meanSP2P)./(meanSP2PL + meanSP2P);

%mean Laser Effect
figure; hold on
bar([1,2], [nanmean(mi_sp_laser), nanmean(mi_laser)])
errorbar([1,2], [nanmean(mi_sp_laser), nanmean(mi_laser)], [sem(mi_sp_laser), sem(mi_laser)])
xticks([1 2])
xticklabels({'spont','evoked'})
ylabel('VIP MI (mean/sem)')

figure; hold on
bar([1 1.1 1.2 2 2.1 2.2], [nanmean(mi_sp_laser), nanmean(mi_sp_laser_m) nanmean(mi_sp_laser_p) nanmean(mi_laser) nanmean(mi_laser_m) nanmean(mi_laser_p)])
errorbar([1 1.1 1.2 2 2.1 2.2], [nanmean(mi_sp_laser), nanmean(mi_sp_laser_m) nanmean(mi_sp_laser_p) nanmean(mi_laser) nanmean(mi_laser_m) nanmean(mi_laser_p)], ...
    [sem(mi_sp_laser), sem(mi_sp_laser_m) sem(mi_sp_laser_p) sem(mi_laser) sem(mi_laser_m) sem(mi_laser_p)])
xticks([1 2])
xticklabels({'spont','evoked'})
ylabel('VIP MI (mean/sem)')

%mean laser effect by cell type
figure; hold on
bar([1 2], [nanmean(mi_sp_laser(rs1)),  nanmean(mi_laser(rs1))], 'BarWidth', .3)
bar([1.3, 2.3], [ nanmean(mi_sp_laser(fs1)) nanmean(mi_laser(fs1))], 'BarWidth', .3)
errorbar([1,1.3,2,2.3], [nanmean(mi_sp_laser(rs1)), nanmean(mi_sp_laser(fs1)), nanmean(mi_laser(rs1)), nanmean(mi_laser(fs1))], [sem(mi_sp_laser(rs1)),sem(mi_sp_laser(fs1)),sem(mi_laser(rs1)) sem(mi_laser(fs1))])
xticks([1.15 2.15])
xticklabels({'spont','evoked'})
ylabel('VIP MI (mean/sem)')
legend({'rs', 'fs'})

figure; hold on
bar([1 1.1 1.2 2 2.1 2.2], [nanmean(mi_sp_laser(rs1)), nanmean(mi_sp_laser_m(rs1)) nanmean(mi_sp_laser_p(rs1)) nanmean(mi_laser(rs1)) nanmean(mi_laser_m(rs1)) nanmean(mi_laser_p(rs1))])
errorbar([1 1.1 1.2 2 2.1 2.2], [nanmean(mi_sp_laser(rs1)), nanmean(mi_sp_laser_m(rs1)) nanmean(mi_sp_laser_p(rs1)) nanmean(mi_laser(rs1)) nanmean(mi_laser_m(rs1)) nanmean(mi_laser_p(rs1))], ...
    [sem(mi_sp_laser(rs1)), sem(mi_sp_laser_m(rs1)) sem(mi_sp_laser_p(rs1)) sem(mi_laser(rs1)) sem(mi_laser_m(rs1)) sem(mi_laser_p(rs1))])
xticks([1 2])
xticklabels({'spont','evoked'})
ylabel('VIP MI (mean/sem)')
title('VIP modulation, rs, all, run, sit trials')
figure; hold on
bar([1 1.1 1.2 2 2.1 2.2], [nanmean(mi_sp_laser(fs1)), nanmean(mi_sp_laser_m(fs1)) nanmean(mi_sp_laser_p(fs1)) nanmean(mi_laser(fs1)) nanmean(mi_laser_m(fs1)) nanmean(mi_laser_p(fs1))])
errorbar([1 1.1 1.2 2 2.1 2.2], [nanmean(mi_sp_laser(fs1)), nanmean(mi_sp_laser_m(fs1)) nanmean(mi_sp_laser_p(fs1)) nanmean(mi_laser(fs1)) nanmean(mi_laser_m(rs1)) nanmean(mi_laser_p(fs1))], ...
    [sem(mi_sp_laser(rs1)), sem(mi_sp_laser_m(fs1)) sem(mi_sp_laser_p(fs1)) sem(mi_laser(fs1)) sem(mi_laser_m(fs1)) sem(mi_laser_p(fs1))])
xticks([1 2])
xticklabels({'spont','evoked'})
ylabel('VIP MI (mean/sem)')
title('VIP modulation, fs, all, run, sit trials')


%VIP effect by cell type and layer WN + SP1
figure(211); hold on
errorbar([0.9 1.9 2.9 3.9], [meanMI_sp_rs], [semMI_sp_rs], 'ko')
errorbar([1 2 3 4], [meanMI_sp_fs], [semMI_sp_fs], 'k*')
errorbar([1.1 2.1 3.1 4.1], [meanMI_rs], [semMI_rs], 'ro')
errorbar([1.2 2.2 3.2 4.2], [meanMI_fs], [semMI_fs], 'r*')
xticks([1:4]); xticklabels({'2/3', '4' '5' '6'})
xlabel('Cortical Layers')
ylabel('Firing Rate Change - VIP'); title('WN responses')
plot([0 5], [0 0], 'k--')
legend({'spont rs', 'spont fs', 'evoked rs', 'evoked fs'})

figure(211); hold on
errorbar([0.9 1.9 2.9 3.9], [meanMI_sp_rs], [semMI_sp_rs], 'ko')
errorbar([1 2 3 4], [meanMI_sp_fs], [semMI_sp_fs], 'k*')
errorbar([1.1 2.1 3.1 4.1], [meanMI_rs], [semMI_rs], 'ro')
errorbar([1.2 2.2 3.2 4.2], [meanMI_fs], [semMI_fs], 'r*')
xticks([1:4]); xticklabels({'2/3', '4' '5' '6'})
xlabel('Cortical Layers')
ylabel('Firing Rate Change - VIP'); title('WN responses')
plot([0 5], [0 0], 'k--')
legend({'spont rs', 'spont fs', 'evoked rs', 'evoked fs'})


%VIP effect by layer across cells
figure(212); hold on
errorbar([0.9 1.9 2.9 3.9], [meanMI_sp], [semMI_sp], 'ko')
errorbar([1 2 3 4], [meanMI], [semMI], 'bo')
xticks([1:4]); xticklabels({'2/3', '4' '5' '6'})
xlabel('Cortical Layers')
ylabel('Firing Rate Change - VIP');
plot([0 5], [0 0], 'k--')
legend({'spont', 'evoked'})

% indx1 = 2:2:4;
% indx2 = 3:4;
indx1 = 4;

VIP_MI_rs_evoked = (nanmean(WN1L(rs1,indx1),2)- nanmean(WN1(rs1,indx1),2));
VIP_MI_fs_evoked = (nanmean(WN1L(fs1,indx1),2)- nanmean(WN1(fs1,indx1),2));
VIP_MI_rs_sp = (nanmean(SP1L(rs1,indx1),2)- nanmean(SP1(rs1,indx1),2));
VIP_MI_fs_sp = (nanmean(SP1L(fs1,indx1),2)- nanmean(SP1(fs1,indx1),2));

RUN_MI_rs_evoked = (WN1(rs1,1) - WN1(rs1,4));
RUN_MI_fs_evoked = (WN1(fs1,1) - WN1(fs1,4));
RUN_MI_rs_sp = (SP1(rs1,1) - SP1(rs1,4));
RUN_MI_fs_sp = (SP1(fs1,1) - SP1(fs1,4));
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
xlabel('Firing Rate Change - VIP')
ylabel('Firing Rate Change - Run')
title('Evoked Activity')
legend({'rs', 'fs'})

subplot(2,1,2); hold on;
plot(VIP_MI_rs_sp, RUN_MI_rs_sp, 'k.');
plot(VIP_MI_fs_sp, RUN_MI_fs_sp, 'b.');
lsline;
xlabel('Firing Rate Change - VIP')
ylabel('Firing Rate Change - Run')
title('Spont Activity')
legend({'rs', 'fs'}) 

figure(217); hold on
subplot(2,1,1); hold on;
plot([VIP_MI_rs_evoked; VIP_MI_fs_evoked], [RUN_MI_rs_evoked; RUN_MI_fs_evoked], 'k.');
lsline;
xlabel('Firing Rate Change - VIP')
ylabel('Firing Rate Change - Run')
title('Evoked Activity')

subplot(2,1,2); hold on;
plot([VIP_MI_rs_sp; VIP_MI_fs_sp], [RUN_MI_rs_sp; RUN_MI_fs_sp], 'k.');
lsline;
xlabel('Firing Rate Change - VIP')
ylabel('Firing Rate Change - Run')
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
WN_rs = (WN1(rs1,1) - WN1(rs1,4));
WN_fs = (WN1(fs1,1) - WN1(fs1,4));
WN2 =  (WN1(:,1) - WN1(:,4));

WN_rsL = (WN1L(rs1,4) - WN1(rs1,4));
WN_fsL = (WN1L(fs1,4) - WN1(fs1,4));
WN2L = (WN1L(:,4) - WN1(:,4));

WN_rs_plus = (WN1L(rs1,1) - WN1(rs1,4));
WN_fs_plus = (WN1L(fs1,1) - WN1(fs1,4));
WN2_plus = (WN1L(:,1) - WN1(:,4));

WN_exp_rs = sum([WN_rs WN_rsL],2);
WN_exp_fs = sum([WN_fs WN_fsL],2);
WN2_exp = sum([WN2 WN2L],2);

[p,h,zstat] = ranksum(WN_rs_plus, WN_exp_rs);
[p,h,zstat] = ranksum(WN_fs_plus, WN_exp_fs);

SP_rs = (SP1(rs1,1) - SP1(rs1,4));
SP_fs = (SP1(fs1,1) - SP1(fs1,4));
SP2 = (SP1(:,1) - SP1(:,4));

SP_rsL = (nanmean(SP1L(rs1,2:3),2) - nanmean(SP1(rs1,2:3),2));
SP_fsL = (nanmean(SP1L(fs1,2:3),2) - nanmean(SP1(fs1,2:3),2));
SP2L = (nanmean(SP1L(:,2:3),2) - nanmean(SP1(:,2:3),2));

SP_rs_plus = (SP1L(rs1,1) - SP1(rs1,4));
SP_fs_plus = (SP1L(fs1,1) - SP1(fs1,4));
SP2_plus = (SP1L(:,1) - SP1(:,4));

SP_exp_rs = sum([SP_rs SP_rsL],2);
SP_exp_fs = sum([SP_fs SP_fsL],2);
SP2_exp = sum([SP2 SP2L],2);

[p,h,zstat] = ranksum(SP_rs_plus, SP_exp_rs);
[p,h,zstat] = ranksum(SP_fs_plus, SP_exp_fs);

figure; subplot(2,1,1);hold on
plot([0 0], [max(abs([WN2; WN2L]))*-1 max(abs([WN2; WN2L]))], 'k--')
plot([min([WN2; WN2L]) max([WN2; WN2L])], [0 0], 'k--')
plot(WN_rs, WN_rsL, 'ko')
plot(WN_fs, WN_fsL, 'go')
ylim([max(abs([WN2; WN2L]))*-1 max(abs([WN2; WN2L]))])
xlim([max(abs([WN2; WN2L]))*-1 max(abs([WN2; WN2L]))])
[rho,p]=corr(WN2, WN2L, 'Type', 'Spearman', 'rows', 'pairwise');
title(sprintf('Evoked, rho=%.4f,p=%.4f ', rho, p))
xlabel('run - VIP'); ylabel('VIP - run'); lsline

subplot(2,1,2);hold on
plot([0 0], [max(abs([SP2; SP2L]))*-1 max(abs([SP2; SP2L]))], 'k--')
plot([min([SP2; SP2L]) max([SP2; SP2L])], [0 0], 'k--')
plot(SP_rs, SP_rsL, 'ko')
plot(SP_fs, SP_fsL, 'go')
ylim([max(abs([WN2; WN2L]))*-1 max(abs([WN2; WN2L]))])
xlim([max(abs([WN2; WN2L]))*-1 max(abs([WN2; WN2L]))])
xlabel('run - VIP'); ylabel('VIP - run'); lsline
[rho,p]=corr(SP2, SP2L, 'Type', 'Spearman', 'rows', 'pairwise');
title(sprintf('Spont,  rho=%.4f,p=%.4f ', rho, p))

figure; subplot(2,1,1);hold on
plot(WN_rs_plus, WN_exp_rs, 'ko'); 
plot(WN_fs_plus, WN_exp_fs, 'go');
xlabel('VIP + run'); ylabel('Expected'); lsline
[rho,p]=corr(WN2_plus, WN2_exp, 'Type', 'Spearman', 'rows', 'pairwise');
title(sprintf('Evoked, rho=%.4f,p=%.4f ', rho, p))
subplot(2,1,2);hold on
plot(SP_rs_plus, SP_exp_rs, 'ko'); 
plot(SP_fs_plus, SP_exp_fs, 'go');
xlabel('VIP + run'); ylabel('Expected'); lsline
[rho, p]=corr(SP2_plus, SP2_exp, 'Type', 'Spearman', 'rows', 'pairwise');
title(sprintf('Spont, rho=%.4f,p=%.4f', rho, p))


numCells = [sum(~isnan(SP_rs)) sum(~isnan(SP_rsL)) sum(~isnan(SP_rs_plus)) sum(~isnan(SP_exp_rs));sum(~isnan(SP_fs)) sum(~isnan(SP_fsL)) sum(~isnan(SP_fs_plus)) sum(~isnan(SP_exp_fs))];

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
ylabel('Firing Rate Change')

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

WN_rs = (WN1(rs1,2) - WN1(rs1,3));
WN_fs = (WN1(fs1,2) - WN1(fs1,3));
WN_rsL = (WN1L(rs1,3) - WN1(rs1,3));
WN_fsL = (WN1L(fs1,3) - WN1(fs1,3));
WN_rs_plus = (WN1L(rs1,2) - WN1(rs1,3));
WN_fs_plus = (WN1L(fs1,2) - WN1(fs1,3));
% WN_exp_rs = [nanmean(WN_rs) nanmean(WN_rsL)];
% WN_exp_fs = [nanmean(WN_fs) nanmean(WN_fsL)];
WN_exp_rs = sum([WN_rs WN_rsL],2);
WN_exp_fs = sum([WN_fs WN_fsL],2);

[p,h,zstat] = ranksum(WN_rs_plus, WN_exp_rs);
[p,h,zstat] = ranksum(WN_fs_plus, WN_exp_fs);

SP_rs = (SP1(rs1,2) - SP1(rs1,3));
SP_fs = (SP1(fs1,2) - SP1(fs1,3));
SP_rsL = (SP1L(rs1,3) - SP1(rs1,3));
SP_fsL = (SP1L(fs1,3) - SP1(fs1,3));
SP_rs_plus = (SP1L(rs1,2) - SP1(rs1,3));
SP_fs_plus = (SP1L(fs1,2) - SP1(fs1,3));
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
ylabel('Firing Rate Change')

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
