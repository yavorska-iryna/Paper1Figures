% Calculate Distance Correlation values for silent stimuli (VIP and PV)
clear; warning off; %close all

load('E:\djmaus-data\iraira\VIP_analysis\variables\SilencedirsVIP.mat')
DIRS=silence;
load('E:\djmaus-data\iraira\PV_analysis\variables\PVSilence.mat')
DIRS=[DIRS silence];
binwidth = [0.01 0.025 .05 0.1 .2 .4 .8 1.6 3.2 6.4 12.8]; % bins
cc = 0; %total cell counter 
POS =nan(length(DIRS),length(binwidth));
for d = 1:length(DIRS)
    try
        cd(DIRS{d})
        load('dirs.mat')
        load('notebook.mat')
        sp = loadKSdir(dirs{1}); 
        good_cells = find(sp.cgs == 2);
        load([dirs{1} '\RecLengths.mat'])
        load([dirs{1} '\maxChans.mat'])
        [maxChans1, maxChans1_indx] = sort(maxChans);
        load([dirs{1} '\WFs_stats.mat'])
        load('moves_trace1.mat')
        %calculate average speed
        m = abs(moves_trace1- median(moves_trace1));
        mM(d) = mean(m(m>0.5));
        
        laxis1= nan(length(moves_trace1),1);
        this_dir=pwd;
        silence_dir_indx = find(strcmp(dirs, this_dir));
        if silence_dir_indx == 1
            start = 0;
        else
            start = L(silence_dir_indx - 1);
        end
        recording_duration = L(silence_dir_indx) - start;
        good_cells1 = good_cells(maxChans1_indx); %reorganize cells by their probe location
        width = total_width(maxChans1_indx);
        shank1_labels = {}; shank2_labels = {};
        
        for b = 1: length(binwidth)
            ST1 = []; ST2 = []; ST_all =[]; data =[];
            inc1 = 0; inc2 = 0;
            switch1 = find(maxChans1 > 33);
            
            if isempty(switch1)
                switch1= 0;
            end
            
            for c = 1:length(good_cells1)
                cc = cc+1;
                if width(c)<.71
                    col = 'b';
                else
                    col = 'k';
                end
                
                if maxChans1(c) < 33
                    shank = 1;
                    shank1_labels = {shank1_labels num2str(maxChans1(c))};
                    inc1 = inc1+1;
                    offset = inc1;
                elseif maxChans1(c) > 32
                    shank = 2;
                    shank2_labels = {shank1_labels num2str(maxChans1(c)-32)};
                    inc2 = inc2 +1;
                    offset = inc2;
                end
                
                st1=sp.st(sp.clu == sp.cids(good_cells1(c)));
                st1 = st1(st1> start & st1 < L(silence_dir_indx));
                st1 = st1 - start;
                [N,x] = hist(st1, 0:binwidth(b):recording_duration);
                N=N/binwidth(b);
                N=smooth(N,10); N = N';
                N=N/max(N);
                data=[data; N];
                
                if shank == 1
                    ST1 = [ST1 N'];
                elseif shank == 2
                    ST2 = [ST2 N'];
                end
                
                if c == 1 || c == switch1(1) %switching to shank2
                    [M1,P1] = resampleTraces(N, moves_trace1, laxis1);
                end
                
            end % cells
            
            pos = 0; nn = 20;
            for n = 1:nn
                M1_shuffled = Shuffle(M1);
                try
                    DC_running_shuffled(d,b,1) = distcorr(ST1, M1_shuffled);
                    DC_running(d,b,1) = distcorr(ST1, M1);
                catch
                    DC_running_shuffled(d,b,1) = NaN;
                    DC_running(d,b,1) = NaN;
                end
                
                try
                    DC_running_shuffled(d,b,2) = distcorr(ST2, M1_shuffled);
                    DC_running(d,b,2) = distcorr(ST2, M1);
                    
                catch
                    DC_running_shuffled(d,b,2) = NaN;
                    DC_running(d,b,2) = NaN;
                end
                
                % repeat shuffle n times to see how many times it's
                % above non random dc value at 100 ms, 0.12
               ShuffledDC = nanmean(DC_running_shuffled(d,b,:));
               if ShuffledDC >= nanmean(DC_running(d,b,:) -DC_running_shuffled(d,b,:))
                    pos = pos +1;
               end
                
            end % number of shuffling
            POS(d,b) = pos;
            
        end %bin
        mouseID(d) =  str2num(nb.mouseID);
    end %try
end %dir

zero_indx = find(DC_running(:,1,1)==0);
DC_running(zero_indx,:,1) = NaN;
zero_indx = find(DC_running(:,1,2)==0);
DC_running(zero_indx,:,2) = NaN;

cd('C:\Users\lab\Documents\GitHub\Paper1Figures\code\variables') %cd('C:\Users\lab\Resi\Paper1Figures\code\variables');
save('Silence_DistanceCorr_final4.mat', 'DC_running', 'DC_running_shuffled','mouseID', 'POS')
DistCorr = DC_running - DC_running_shuffled;
figure; hold on
means  =  nanmean(squeeze(nanmean(DistCorr))');
SEMs = sem(nanmean(DistCorr,3));
errorbar([1:length(means)], means, SEMs, 'ko-')
ylabel('Distance Corr')
xlabel('bin time (sec)')
for i = 1:length(binwidth)
    time_labels{i} = num2str(binwidth(i));
end
xticklabels(time_labels)
xlim([0 length(means)+1])

figure; hold on
plot(DistCorr(:,:,1)', 'o-')
xticklabels(time_labels)
xlim([0 length(means)+1])
ylabel('Distance Corr')
xlabel('bin time (sec)')
figure; hold on
plot(DistCorr(:,:,2)', 'o-')
xticklabels(time_labels)
xlim([0 length(means)+1])
ylabel('Distance Corr')
xlabel('bin time (sec)')