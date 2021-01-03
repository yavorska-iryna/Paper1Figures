
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
meanON = nan(length(data),2, 21); % On response only (0 - 100 ms)
meanONL = nan(length(data),2, 21);
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
    if data(cc).dir < 30
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
            
            ON = [data(cc).mNON_on; data(cc).mNON_off]; 
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
                meanON(cc,1,:) = nanmean(data(cc).mNON_on);%-nanmean(nanmean(data(cc).mNson));
                meanON(cc,2,:) = nanmean(data(cc).mNON_off);%-nanmean(nanmean(data(cc).mNsoff));
                meanPreStim(cc,1) = nanmean(nanmean(data(cc).mNson));
                meanPreStim(cc,2) = nanmean(nanmean(data(cc).mNsoff));
                WNdirsM = [WNdirsM data(cc).dir]; % collect directory (rec#) for all recordings included in the analysis to make sure N is ok.
                WNcellsM = WNcellsM + 1; % count number of cells included in the analysis
            end
            
            % WN laser on trials
            if data1(cc).nrepsWNMon > nreps_check_running && data1(cc).nrepsWNMoff > nreps_check_running % && data(cc).nrepsWNMon > nreps_check_running && data(cc).nrepsWNMoff > nreps_check_running
                meanONL(cc,1,:) =nanmean(data1(cc).mNON_on);%-nanmean(nanmean(data1(cc).mNson));
                meanONL(cc,2,:) =nanmean(data1(cc).mNON_off);%-nanmean(nanmean(data1(cc).mNsoff));
                meanPreStimL(cc,1) = nanmean(nanmean(data1(cc).mNson));
                meanPreStimL(cc,2) = nanmean(nanmean(data1(cc).mNsoff));
            end
            
            % silent sound laser off
            if data(cc).nrepsSSMon > nreps_check_running && data(cc).nrepsSSMoff > nreps_check_running %&& data1(cc).nrepsSSMon > nreps_check_running && data1(cc).nrepsSSMoff > nreps_check_running
                SP(cc,1) = nanmean(nanmean(data(cc).mSSon));
                SP(cc,2) = nanmean(nanmean(data(cc).mSSoff));
                meanSilentSound(cc,2) = nanmean(nanmean(data(cc).mSSoff));
                meanSilentSound(cc,1) = nanmean(nanmean(data(cc).mSSon));
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
       try
            reps(cc,:,:,:) = [data(cc).nrepsWNMon data(cc).nrepsWNMoff size(data1(40).mNON_on,1) size(data1(40).mNON_off,1)];
        catch
            reps(cc,:,:,:) = [NaN NaN NaN NaN];
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

% get only cells with On responses (evoked 1 == On, zstat 1 >0 == activated)
evoked1 = logical(evoked(:,1) & zstats(:,1)>0); 
reps1 = reps(evoked1,:);
rs1 = logical(Rs(evoked1));
fs1= logical(fs(evoked1));
depths1 = depths(evoked1);

% sitting vs running
evoked_sitting = squeeze(meanON(evoked1,2,:));
evoked_running =squeeze(meanON(evoked1,1,:));

% subtract spont?
%evoked_sitting = evoked_sitting - (SP(evoked1,2));
%evoked_running = evoked_running - (SP(evoked1,1));

%evoked_sitting = evoked_sitting - (meanPreStim(evoked1,2));
%evoked_running = evoked_running - (meanPreStim(evoked1,1));


%sitting and running + laser
evoked_sitting_laser =squeeze(meanONL(evoked1,2,:));
evoked_running_laser =squeeze(meanONL(evoked1,1,:));

% smooth responses
sigma = 3;
num_cells = length(evoked_sitting);
sitting_st = [];
sitting_max_t = [];
for i =1:num_cells
    st = evoked_sitting(i,:);
   %st = st./reps1(i,2);
    %[x, st] = GaussSmooth(st, sigma, [0 100]);
    %sitting_st(i,:) =st./reps1(i,2);
    sitting_st(i,:) = smooth(st, sigma);
    try
        sitting_max_t(i) = find(st == max(st));
    catch
        sitting_max_t(i) = NaN;
    end
end

running_st = [];
running_max_t =[];
for i =1:num_cells
    st = evoked_running(i,:);
    %st = st./reps1(i,1);
    %[x, st] = GaussSmooth(st, sigma, [0 100]);
    %running_st(i,:) =st./reps1(i,1);
    running_st(i,:) = smooth(st, sigma);
    try
        running_max_t(i) = find(st == max(st));
    catch
        running_max_t(i) = NaN;
    end
end


% plot mean responses
figure; hold on
plot(nanmean(sitting_st)')
plot(nanmean(running_st)')
xticks([0,5,10,15,20])
xlim([0 20])
xticklabels({'0', '25', '50', '75', '100'})
legend('sitting', 'running' )
ylabel('Firing Rate (Hz)'); xlabel('time (ms)')
[P,H,STATS] = signrank(running_max_t, sitting_max_t)

figure; hold on
plot(nanmean(sitting_st(rs1,:))')
plot(nanmean(running_st(rs1,:))')
[P,H,STATS] = signrank(running_max_t(rs1), sitting_max_t(rs1))

plot(nanmean(sitting_st(fs1,:))')
plot(nanmean(running_st(fs1,:))')
xticks([0,5,10,15,20])
xlim([0 20])
xticklabels({'0', '25', '50', '75', '100'})
legend('sitting RS', 'running RS', 'sitting FS', 'running FS' )
ylabel('Firing Rate (Hz)'); xlabel('time (ms)')
[P,H,STATS] = signrank(running_max_t(fs1), sitting_max_t(fs1))

% break it down by layer

layer1 = nan(length(depths1),1);
for l = 1:length(CL)
    layer = CL{l};
    for d = 1:length(depths1)
        if depths1(d) > layer(1) & depths1(d) < layer(2)
            layer1(d) = l;
        end
    end
end

figure; hold on
layers = {'2/3', '4', '5', '6'}
for cl = 1:4
    subplot(1,4, cl); hold on
    indx = find(layer1 == cl); % indices of cells within this layer
    fs2 = fs1(indx); fs2 = logical(fs2); % fast spiking cells in this layer
    rs2 = rs1(indx); rs2 = logical(rs2); % regular spiking cells in this layer
    
    st = sitting_st(indx,:);
    plot(nanmean(st(rs2,:))', 'k-')
    plot(nanmean(st(fs2,:))', 'g-')
    
    st = running_st(indx,:);
    plot(nanmean(st(rs2,:))', 'k--')
    plot(nanmean(st(fs2,:))', 'g--')
    xticks([0,5,10,15,20])
    xlim([0 20])
    xticklabels({'0', '25', '50', '75', '100'})
    legend('sitting RS', 'sitting FS', 'running RS', 'running FS' )
    ylabel('Firing Rate (Hz)'); xlabel('time (ms)')
    title(sprintf('layer %s', layers{cl}))
end

%% the same for laser off and on 

sitting_st_laser = [];
sitting_max_t_laser = [];
for i =1:num_cells
    st = evoked_sitting_laser(i,:);
   %st = st./reps1(i,2);
    %[x, st] = GaussSmooth(st, sigma, [0 100]);
    %sitting_st(i,:) =st./reps1(i,2);
    sitting_st_laser(i,:) = smooth(st, sigma);
    try
        sitting_max_t_laser(i) = find(st == max(st));
    catch
        sitting_max_t_laser(i) = NaN;
    end
end

% plot mean responses
figure; hold on
plot(nanmean(sitting_st)', 'k')
plot(nanmean(sitting_st_laser)', 'c')
xticks([0,5,10,15,20])
xlim([0 20])
xticklabels({'0', '25', '50', '75', '100'})
legend('sitting', 'running' )
ylabel('Firing Rate (Hz)'); xlabel('time (ms)')
[P,H,STATS] = signrank(sitting_max_t_laser, sitting_max_t)


figure; hold on
plot(nanmean(sitting_st(rs1,:))', 'k')
plot(nanmean(sitting_st_laser(rs1,:))', 'c')
[P,H,STATS] = signrank(sitting_max_t_laser(rs1), sitting_max_t(rs1))

plot(nanmean(sitting_st(fs1,:))', 'b--')
plot(nanmean(sitting_st_laser(fs1,:))', 'c--')
xticks([0,5,10,15,20])
xlim([0 20])
xticklabels({'0', '25', '50', '75', '100'})
legend('sitting RS laser off', 'sitting RS laser on', 'sitting FS laser off', 'sitting FS laser on' )
ylabel('Firing Rate (Hz)'); xlabel('time (ms)')
[P,H,STATS] = signrank(sitting_max_t_laser(fs1), sitting_max_t(fs1))


figure; hold on
layers = {'2/3', '4', '5', '6'}
for cl = 1:4
    subplot(1,4, cl); hold on
    indx = find(layer1 == cl); % indices of cells within this layer
    fs2 = fs1(indx); fs2 = logical(fs2); % fast spiking cells in this layer
    rs2 = rs1(indx); rs2 = logical(rs2); % regular spiking cells in this layer
    
    st = sitting_st(indx,:);
    plot(nanmean(st(rs2,:))', 'k-')
    plot(nanmean(st(fs2,:))', 'k--')
    
    st = sitting_st_laser(indx,:);
    plot(nanmean(st(rs2,:))', 'c-')
    plot(nanmean(st(fs2,:))', 'c--')
    xticks([0,5,10,15,20])
    xlim([0 20])
    xticklabels({'0', '25', '50', '75', '100'})
    legend('RS laser off', 'FS laser off', 'RS laser on', 'FS laser on' )
    ylabel('Firing Rate (Hz)'); xlabel('time (ms)')
    title(sprintf('layer %s', layers{cl}))
end
