%06.30.19 % collect data for response analysis
clear; close all
tic
cdPV
load('allWNdirs.mat');
DIRS = WNDIRS;

cdPV
load('allPINPdirs.mat') % buffer  cells (if WN was not recorded)
load('PinpedPVcellIDs.mat')
load('PinpedPVcellIDs_good_quality.mat')
save_dir = 'C:\Users\lab\Resilio Sync\Paper1Figures\code\variables';

ONresponse_region = [0:100];
Sustainedresponse_region = [100:600];
OFFresponse_region = [600:700];
spont_region = [-300:0]; % prior to the stimulus
SS_region = [-25:650];
laserStart =[]; laserWidth=[]; TR = [];
pupil_thresholds = [.55 .55 .50 .55 .55 .55 .55 .6 .6 .6 .55 .5 .45 .35 .55 ...
    .55 .45 .45 .55 .40 .45 .5 .55 .45 .5 .55 .35 .5 .4 .55 .4 .5 .47 .4 .55 .55 .55 .55 .55 .55 .55 .55 .55 .55 .55 .55];
stimulus = 'WNOFF';
cc=0; 
for d = 1:length(DIRS)
    if ~isempty(DIRS{d})
        cd(DIRS{d})
        load('dirs.mat')
        sp = loadKSdir(dirs{1});
        good_cells = find(sp.cgs == 2);
        %getStatesTC(pupil_thresholds(d));
        load('StateEvents.mat'); binwidth = 5 ;%m
        load('Events.mat'); load('notebook.mat')
        load([dirs{1} '\WFs_stats.mat']);
        load([dirs{1} '\probeDepth.mat'])
        for c = 1 : length(good_cells)
            cc=cc+1; pv = 0; pv_good = 0;
            outfilename = sprintf('outPSTH_ch-1c%d.mat', good_cells(c));
            t_filename = sprintf('ch-1_simpleclust_%d.t', good_cells(c));
            %             PlotTC_PSTH_single(pwd, t_filename)
            load(outfilename);
            
            for pv_indx = 1: length(PV_all) % mark pv cell 
                if strcmp(PV_all(pv_indx).dir, dirs{1}) && PV_all(pv_indx).id == out.KiloSort_ID
                    pv = 1;
                end
            end
            
            for pv_indx = 1: length(PV) % mark pv cell 
                if strcmp(PV(pv_indx).dir, dirs{1}) && PV(pv_indx).id == out.KiloSort_ID
                    pv_good = 1;
                end
            end
            
            nrepsWN =0;
            nrepsWNPon = 0; nrepsWNPoff = 0; nrepsSSPon = 0; nrepsSSPoff = 0;
            nrepsWNMon = 0; nrepsWNMoff = 0; nrepsSSMon = 0; nrepsSSMoff = 0;
            
            pNON_on = []; pNON_off = []; mNON_on = []; mNON_off = [];
            pNSustained_on = []; pNSustained_off = []; mNSustained_on = []; mNSustained_off = [];
            pNOFF_on = []; pNOFF_off = []; mNOFF_on = []; mNOFF_off = [];
            
            pWNon =[]; pWNoff =[]; mWNon =[]; mWNoff =[];
            pSSon =[]; pSSoff =[]; mSSon =[]; mSSoff =[];
            
            pNson = []; pNsoff =[]; mNson =[]; mNsoff =[];
            pNSS_on =[]; pNSS_off =[]; mNSS_on =[]; mNSS_off =[];
            
            nrepsWNPM = 0; pmNON = []; pmNSustained = []; pmNOFF = [];
            spikeCountPMON =[]; spikeCountPMOFF =[];  pmWN =[]; pmNs = [];
            nrepsSSPM = 0; pmNSS = []; pmSS =[]; mST_off=[]; mST_on=[];
            
            
            spikeCountPonON = []; spikeCountPoffON =[]; spikeCountPonOFF =[]; spikeCountPoffOFF =[];
            spikeCountMonON = []; spikeCountMoffON =[]; spikeCountMonOFF =[]; spikeCountMoffOFF =[];
            Ns=[]; Nr =[];
            
            if strcmp(stimulus, 'WNOFF')
                f = 1; a = 2; M1= squeeze(out.M1OFF(f,a,1,:)); mM1 = squeeze(out.mM1OFF(f,a)); laser = 0; response = 1;
                SS = out.SilentSoundOFF; mSS = out.mSilentSoundOFF; nrepsM1 = out.nrepsOFF(2); nrepsSS = out.nreps_ssOFF;
                if nrepsM1 > 1 && nrepsSS > 1 ; trials_recorded = 1; else; trials_recorded = 0; end
                
            elseif strcmp(stimulus, 'WNON')
                f = 1; a = 2; M1= squeeze(out.M1ON(f,a,1,:)); mM1 = squeeze(out.mM1ON(f,a)); laser = 1;response = 1;
                SS = out.SilentSoundON; mSS = out.mSilentSoundON; nrepsM1 = out.nrepsON(2); nrepsSS = out.nreps_ssON;
                if nrepsM1 > 1 && nrepsSS > 1 ; trials_recorded = 1; else; trials_recorded = 0; end
            end
            
            %stimulus evoked
            st = mM1.spiketimes;
            
            stON = st(st> ONresponse_region(1) & st<ONresponse_region(end));
            stSustained = st(st> Sustainedresponse_region(1) & st<Sustainedresponse_region(end));
            stOFF = st(st> OFFresponse_region(1) & st<OFFresponse_region(end));
            [NON, xON] = hist(stON, ONresponse_region(1):binwidth:ONresponse_region(end));
            [NSustained, xSustained] = hist(stSustained, Sustainedresponse_region(1):binwidth:Sustainedresponse_region(end));
            [NOFF, xOFF] = hist(stOFF, OFFresponse_region(1):binwidth:OFFresponse_region(end));
            
            st1 = st(st> 0 & st<600);
            N1 = hist(st1,[0:binwidth:600]);
            N1 = N1./nrepsM1;
            N1 = 1000*N1./binwidth;
            
            st2 =  mSS.spiketimes;
            st2 = st2(st2>0 & st2<600);
            N2 = hist(st2,[0:binwidth:600]);
            N2 = N2./nrepsM1;
            N2 = 1000*N2./binwidth;
            
            
            NON = NON./nrepsM1;
            NON = 1000*NON./binwidth;
            
            NSustained = NSustained./nrepsM1;
            NSustained = 1000*NSustained./binwidth;
            
            NOFF = NOFF./nrepsM1;
            NOFF = 1000*NOFF./binwidth;
            
            %prestimulus period
            st1 = st(st> spont_region(1) & st<spont_region(end));
            [Ns, xs] = hist(st1, spont_region(1):binwidth:spont_region(end));
            Ns = Ns./nrepsM1;
            Ns = 1000*Ns./binwidth;
            
            %silent stimulus
            st = mSS.spiketimes;
            
            stSS = st(st> SS_region(1) & st<SS_region(end));
            [NSS, xs] = hist(stSS, SS_region(1):binwidth:SS_region(end));
            NSS = NSS./nrepsSS;
            NSS = 1000*NSS./binwidth;
            
            maxFR = max([mean(NON), mean(NSustained), mean(NOFF), mean(Ns), mean(NSS)]);
            
            nrepsSS = 0;
            % single trials
            for i = 1:length(Events)
                if Events(i).laser==laser  %&& Events(i).LaserOnOff == laser % comment it out during off trials
                    if strcmp(Events(i).type, 'whitenoise')
                        
                        nrepsWN=nrepsWN+1;
                        st = M1(nrepsWN).spiketimes;
                        st2 = st;
                        SpikeCountWN(nrepsWN) = length(st);
                        stON = st(st> ONresponse_region(1) & st<ONresponse_region(end));
                        stSustained = st(st> Sustainedresponse_region(1) & st<Sustainedresponse_region(end));
                        stOFF = st(st> OFFresponse_region(1) & st<OFFresponse_region(end));
                        
                        [NON, xON] = hist(stON, ONresponse_region(1):binwidth:ONresponse_region(end));
                        [NSustained, xSustained] = hist(stSustained, Sustainedresponse_region(1):binwidth:Sustainedresponse_region(end));
                        [NOFF, xOFF] = hist(stOFF, OFFresponse_region(1):binwidth:OFFresponse_region(end));
                        
                        NON = 1000*NON./binwidth;
                        NSustained = 1000*NSustained./binwidth;
                        NOFF = 1000*NOFF./binwidth;
                        
                        sts = st(st> spont_region(1) & st<spont_region(end));
                        [Ns, xs] = hist(sts, spont_region(1):binwidth:spont_region(end));
                        Ns = 1000*Ns./binwidth;
                        
                        st = st(st>0 & st<600);
                        [Nr, xr] = hist(st, 0:binwidth:600);
                        Nr = 1000*Nr./binwidth;
                        
                        if pupil_indx(i) && motion_indx(i) == 0;
                            nrepsWNPon = nrepsWNPon+1;
                            pNON_on(nrepsWNPon,:) = NON;
                            pNSustained_on(nrepsWNPon,:) = NSustained;
                            pNOFF_on(nrepsWNPon,:) = NOFF;
                            spikeCountPonON(nrepsWNPon) = length(stON);
                            spikeCountPonOFF(nrepsWNPon) = length(stOFF);
                            pWNon(nrepsWNPon,:) = Nr;
                            pNson(nrepsWNPon,:) = Ns;
                        elseif pupil_indx(i) == 0 && motion_indx(i) == 0;
                            nrepsWNPoff = nrepsWNPoff+1;
                            pNON_off(nrepsWNPoff,:) = NON;
                            pNSustained_off(nrepsWNPoff,:) = NSustained;
                            pNOFF_off(nrepsWNPoff,:) = NOFF;
                            spikeCountPoffON(nrepsWNPoff) = length(stON);
                            spikeCountPoffOFF(nrepsWNPoff) = length(stOFF);
                            pWNoff(nrepsWNPoff,:) = Nr;
                            pNsoff(nrepsWNPoff,:) = Ns;
                        end
                        
                        if motion_indx(i);
                            nrepsWNMon = nrepsWNMon+1;
                            mNON_on(nrepsWNMon,:) = NON;
                            mNSustained_on(nrepsWNMon,:) = NSustained;
                            mNOFF_on(nrepsWNMon,:) = NOFF;
                            spikeCountMonON(nrepsWNMon) = length(stON);
                            spikeCountMonOFF(nrepsWNMon) = length(stOFF);
                            mWNon(nrepsWNMon,:) = Nr;
                            mNson(nrepsWNMon,:) = Ns;
                            mEvokedON_on(nrepsWNMon,:) = nanmean(NON)-nanmean(Ns);
                            mST_on = [mST_on st2];
                        elseif motion_indx(i)==0 %&& pupil_indx(i)==0;
                            nrepsWNMoff = nrepsWNMoff+1;
                            mNON_off(nrepsWNMoff,:) = NON;
                            mNSustained_off(nrepsWNMoff,:) = NSustained;
                            mNOFF_off(nrepsWNMoff,:) = NOFF;
                            spikeCountMoffON(nrepsWNMoff) = length(stON);
                            spikeCountMoffOFF(nrepsWNMoff) = length(stOFF);
                            mWNoff(nrepsWNMoff,:) = Nr;
                            mNsoff(nrepsWNMoff,:) = Ns;
                            mEvokedON_off(nrepsWNMoff,:) = nanmean(NON)-nanmean(Ns);
                            mST_off = [mST_off st2];
                        end
                        
                        %pupil without motion
                        if pupil_indx(i)==1 && motion_indx(i)==0;
                            nrepsWNPM = nrepsWNPM+1;
                            pmNON(nrepsWNPM,:) = NON;
                            pmNSustained(nrepsWNPM,:) = NSustained;
                            pmNOFF(nrepsWNPM,:) = NOFF;
                            spikeCountPMON(nrepsWNPM) = length(stON);
                            spikeCountPMOFF(nrepsWNPM) = length(stOFF);
                            pmWN(nrepsWNPM,:) = Nr;
                            pmNs(nrepsWNPM,:) = Ns;
                        end
                        
                        
                    elseif strcmp(Events(i).type, 'silentsound')
                        
                        nrepsSS=nrepsSS+1;
                        st = SS(nrepsSS).spiketimes;
                        SpikeCountSS(nrepsSS) = length(st);
                        stSS = st(st> SS_region(1) & st<SS_region(end));
                        [NSS, xSS] = hist(stSS,SS_region(1):binwidth:SS_region(end));
                        NSS = 1000*NSS./binwidth;
                        
                        stSS = st(st> out.xlimits(1) & st<out.xlimits(end));
                        NS = hist(stSS,out.xlimits(1):binwidth:out.xlimits(end));
                        NS = 1000*NS./binwidth;
                        
                        if pupil_indx(i) && motion_indx(i) == 0;
                            nrepsSSPon = nrepsSSPon+1;
                            pNSS_on(nrepsSSPon,:) = NSS;
                            pSSon(nrepsSSPon,:) = NS;
                        elseif pupil_indx(i) == 0 && motion_indx(i) == 0;
                            nrepsSSPoff = nrepsSSPoff+1;
                            pNSS_off(nrepsSSPoff,:) = NSS;
                            pSSoff(nrepsSSPoff,:) = NS;
                        end
                        
                        if motion_indx(i) == 1 ;
                            nrepsSSMon = nrepsSSMon+1;
                            mNSS_on(nrepsSSMon,:) = NSS;
                            mSSon(nrepsSSMon,:) = NS;
                        elseif motion_indx(i) == 0 %&& pupil_indx(i)== 0; 
                            nrepsSSMoff = nrepsSSMoff+1;
                            mNSS_off(nrepsSSMoff,:) = NSS;
                            mSSoff(nrepsSSMoff,:) = NS;
                        end
                        %pupil without motion
                        if pupil_indx(i)==1 && motion_indx(i)==0
                            nrepsSSPM = nrepsSSPM +1;
                            pmNSS(nrepsSSPM,:) = NSS; %within laser duration timescale
                            pmSS(nrepsSSPM,:) = NS; %full SS interval
                        end
                        
                        
                    end %WN /SS
                end % laser
                
                % plot check
                
            end %events
            
            % save it
            WNdata(cc).dir = d;
            WNdata(cc).cell = good_cells(c);
            
            WNdata(cc).nrepsWNPon =nrepsWNPon; % number of repetitions in each condition, 
            WNdata(cc).nrepsWNPoff =nrepsWNPoff;
            WNdata(cc).nrepsWNMon =nrepsWNMon;
            WNdata(cc).nrepsWNMoff =nrepsWNMoff;
            WNdata(cc).nrepsSSPon =nrepsSSPon;
            WNdata(cc).nrepsSSPoff =nrepsSSPoff;
            WNdata(cc).nrepsSSMon =nrepsSSMon;
            WNdata(cc).nrepsSSMoff =nrepsSSMoff;
            
            WNdata(cc).pNON_on =pNON_on; % pupil on, ON response (0-100)
            WNdata(cc).pNON_off =pNON_off;  % pupil off, ON response
            WNdata(cc).pNSustained_on =pNSustained_on; %pupil on sustained response (100-600)
            WNdata(cc).pNSustained_off =pNSustained_off;
            WNdata(cc).pNOFF_on =pNOFF_on; %pupil on OFF response 
            WNdata(cc).pNOFF_off =pNOFF_off;
            
            WNdata(cc).mNON_on =mNON_on;
            WNdata(cc).mNON_off =mNON_off;
            WNdata(cc).mNSustained_on =mNSustained_on;
            WNdata(cc).mNSustained_off =mNSustained_off;
            WNdata(cc).mNOFF_on =mNOFF_on;
            WNdata(cc).mNOFF_off =mNOFF_off;
            
            WNdata(cc).pNSS_on =pNSS_on;
            WNdata(cc).pNSS_off =pNSS_off;
            WNdata(cc).mNSS_on =mNSS_on;
            WNdata(cc).mNSS_off =mNSS_off;
            
            WNdata(cc).pWNon =pWNon;
            WNdata(cc).pWNoff =pWNoff;
            WNdata(cc).mWNon =mWNon;
            WNdata(cc).mWNoff =mWNoff;
            
            WNdata(cc).pSSon =pSSon;
            WNdata(cc).pSSoff =pSSoff;
            WNdata(cc).mSSon =mSSon;
            WNdata(cc).mSSoff =mSSoff;
            
            WNdata(cc).pNson =pNson;
            WNdata(cc).pNsoff =pNsoff;
            WNdata(cc).mNson =mNson;
            WNdata(cc).mNsoff =mNsoff;
            
            WNdata(cc).spikeCountPonON =spikeCountPonON;
            WNdata(cc).spikeCountPonOFF =spikeCountPonOFF;
            WNdata(cc).spikeCountPoffON =spikeCountPoffON;
            WNdata(cc).spikeCountPoffOFF =spikeCountPoffOFF;
            
            WNdata(cc).spikeCountMonON =spikeCountMonON;
            WNdata(cc).spikeCountMonOFF =spikeCountMonOFF;
            WNdata(cc).spikeCountMoffON =spikeCountMoffON;
            WNdata(cc).spikeCountMoffOFF =spikeCountMoffOFF;
            WNdata(cc).x =xr;
            WNdata(cc).binwidth = binwidth;
            WNdata(cc).maxFR = maxFR;
            
            WNdata(cc).nrepsWNPM = nrepsWNPM;
            WNdata(cc).pmNON = pmNON;
            WNdata(cc).pmNSustained = pmNSustained;
            WNdata(cc).pmNOFF = pmNOFF;
            WNdata(cc).spikeCountPMON = spikeCountPMON;
            WNdata(cc).spikeCountPMOFF = spikeCountPMOFF;
            WNdata(cc).pmWN = pmWN;
            WNdata(cc).pmNs = pmNs;
            
            WNdata(cc).laserStart = out.LaserStart;
            WNdata(cc).laserWidth  = out.LaserWidth;
            WNdata(cc).SpikeCountSS =SpikeCountSS;
            WNdata(cc).SpikeCountWN =SpikeCountWN;
            WNdata(cc).trials_recorded = trials_recorded; % where there laser off or laser on trials?
            WNdata(cc).durs = out.durs;
            WNdata(cc).width = total_width(c);
            WNdata(cc).endslope = endslope(c);
            WNdata(cc).pv = pv;
            WNdata(cc).pv_good = pv_good;
            WNdata(cc).depth = cds(c);
            WNdata(cc).WNresponse = N1;
            WNdata(cc).SSresponse = N2;
%             WNdata(cc). mEvokedON_off = mEvokedON_off;
%             WNdata(cc).mEvokedON_on = mEvokedON_on;
            WNdata(cc).mST_off = mST_off;
            WNdata(cc).mST_on = mST_on;
        end % cells
        close all
    else
        cd(PINPDIRS{d})
        load('dirs.mat')
        sp = loadKSdir(dirs{1});
        good_cells = find(sp.cgs == 2);
        
        for c = 1:length(good_cells)
            cc=cc+1;
            WNdata(cc).dir = d;
            WNdata(cc).cell = good_cells(c);
            WNdata(cc).x =NaN;
        end
        
    end
end
cdPV
WNdataLaserOFF = WNdata;
cd(save_dir)
save('WNdataLaserOFF.mat', 'WNdataLaserOFF')
toc

%01.20.20 ira three states : running, sitting + small pupil, sitting +
%large pupil. 


%go to WNresponsePlots to plot the results
% 'WNdataLaserOFF_threeSTates.mat' - sit+small pupil, sit+large pupil
% (pupilon), run+large pupil (runon) pupil off and run off should be the
% same
