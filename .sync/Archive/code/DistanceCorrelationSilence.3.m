% Calculate Distance Correlation values for silent stimuli (VIP and PV)
clear; warning off; %close all

load('E:\djmaus-data\iraira\VIP_analysis\variables\SilencedirsVIP.mat')
DIRS=silence;
load('E:\djmaus-data\iraira\PV_analysis\variables\PVSilence.mat')
DIRS=[DIRS silence];
binwidth = [.05 0.1 .2 .4 .8 1.6];
time_interval = 10/binwidth;

cc = 0; %total cell counter
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
        try
            load('pupil_long_axis_normalized.mat')
        catch
            %normalized_pupil
            cd(DIRS{d})
            load('pupil_long_axis_normalized.mat')
        end
        load('moves_trace1.mat')
        this_dir=pwd;
        silence_dir_indx = find(strcmp(dirs, this_dir));
        if silence_dir_indx == 1
            start = 0;
        else
            start = L(silence_dir_indx - 1);
        end
        recording_duration = L(silence_dir_indx) -start;
        good_cells1 = good_cells(maxChans1_indx); %reorganize cells by their probe location
        width = width(maxChans1_indx);
        shank1_labels = {}; shank2_labels = {};
        inc1 = 0; inc2 = 0;
        ST1 = []; ST2 = []; ST_all =[]; data =[];
        switch1 = find(maxChans1>33);
        if isempty(switch1)
            switch1= 0;
        end
        for b = 1: length(binwidth)
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
                l = length(st1);
                
                subplot(1,2,shank); hold on
                [N,x] = hist(st1, 0:binwidth(b):recording_duration);
                N=N/binwidth(b);
                N=smooth(N,10); N = N';
                N=N/max(N);
                data=[data; N];
                
                if c == 1 || c == switch1(1) %switching to shank2
                    [M1,P1] = resampleTraces(N, moves_trace1, laxis1);
                    if exist('Lasertrace')
                        [L1,S1] = resampleTraces(N, Lasertrace, Stimtrace);
                    end
                    M1=M1/max(M1);
                    moves_on = find(abs(M1)>.1);
                    moves_off = find(abs(M1)<.09);
                    pupil_on = find(P1>0.55);
                    pupil_off = find(P1<0.54);
                    movesRatio(d) = length(moves_on)/length(M1);
                    try
                        pupilRatio(d) = length(pupil_on)/length(P1);
                    catch
                        pupilRatio(d) = NaN;
                    end
                    
                    if shank == 1
                        ST1 = [ST1 N'];
                    elseif shank == 2
                        ST2 = [ST2 N'];
                    end
                    
                    %DC
                    %             try
                    NumCells(d,b,1)=size(ST1,2);
                    %                 DC_pupil(d,1) = distcorr(ST1, P1);
                    DC_running(d,b,1) = distcorr(ST1, M1);
                    %
                    %             catch
                    %                 DC_pupil(d,1) = NaN;
                    DC_running(d,b,1) = NaN;
                    %             end
                    %
                    %             try
                    %                 NumCells(d,2)=size(ST2,2);
                    %                 DC_pupil(d,2) = distcorr(ST2, P1);
                    DC_running(d,b,2) = distcorr(ST2, M1);
                    %             catch
                    %                 DC_pupil(d,2) = NaN;
                    %                 DC_running(d,2) = NaN;
                end
                mouseID(d) =  str2num(nb.mouseID);
                title(nb.mouseID)
            end %cell
        end %bin
    end %try
end %dir
%cdVIP; save('Silence_DistanceCorr_dirs_shanks.mat', 'DC_pupil', 'DC_running', 'mouseID', 'NumCells')


%check tunning properties
