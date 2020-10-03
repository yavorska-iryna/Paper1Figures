clear;
cdPV; load('PVPINPdirs.mat')
DIRS = PINPdirs;
cdVIP; load('PINPdirsVIP.mat')
DIRS = [DIRS PINPdirs];
rec = 0;
for D = 1:length(DIRS)
    cd(DIRS{D})
    load('dirs.mat')
    
    for d = 1:length(dirs)
        cd(dirs{d})
        try
            getStatesTC;
        end
        try
            load('pupil_long_axis_normalized.mat') %normalized to range
            load('pupil_long_axis_normalized2.mat') %normalized to max
        catch
            normalized_pupil;
            cd(dirs{d})
            load('pupil_long_axis_normalized2.mat')
        end
        
        load('moves_trace1.mat')
        
        if sum(laxis1)~=0 && sum(isnan(laxis1))~=length(laxis1)
            rec = rec+1;
            [M1, P1] = resampleTraces(laxis1, moves_trace1, laxis1);
            [~, P2] = resampleTraces(laxis1, moves_trace1, laxis2);
            M1 = M1/max(abs(M1));
            moves_indx = logical( abs(M1)>median(M1)*2); %covert to logical - moves
            
            % normalized to range
            [p1CountAll,x] = hist(P1,[0:.05:1]);
            [p1CountMoves,x] = hist(P1(moves_indx),[0:.05:1]);
            %normalized to max
            [p2CountAll,x] = hist(P2,[0:.05:1]);
            [p2CountMoves,x] = hist(P2(moves_indx),[0:.05:1]);
            
            prob1 = p1CountMoves./p1CountAll;
            prob2 = p2CountMoves./p2CountAll;
            
%             figure(D);hold on
%             plot(x, prob1);
%             
%             figure(100+D); hold on
%             plot(x, prob2);
%             
            PROB1(rec,:) = prob1;
            PROB2(rec,:) = prob2;
            p1SizeMoves(rec,:) = p1CountMoves;
            p2SizeMoves(rec,:) = p2CountMoves;
            p1SizeAll(rec,:) = p1CountAll;
            p2SizeAll(rec,:) = p2CountAll;
        else
            %fprintf('no pupil')
        end
        
        
    end
end

% plot
%probability of movement
figure;
shadedErrorBar(x,nanmean(PROB1),sem(PROB1),'lineprops','k-');
xlabel('pupil size, normalized by range')
ylabel('probability of movement')
figure;
shadedErrorBar(x,nanmean(PROB2),sem(PROB2),'lineprops','k-');
ylabel('probability of movement')
xlabel('pupil size, normalized by max')

%pupil size during movement
figure;
shadedErrorBar(x,nanmean(p1SizeMoves),sem(p1SizeMoves),'lineprops','k-');
xlabel('pupil size, normalized by range')
ylabel('# of movement events')
figure;
shadedErrorBar(x,nanmean(p2SizeMoves),sem(p2SizeMoves),'lineprops','k-');
ylabel('# of movement events')
xlabel('pupil size, normalized by max')



%pupil size during movement and all pupil sizes
figure; hold on
maxPupilMean = max(nanmean(p1SizeAll));
p1 = plot(x,nanmean(p1SizeAll)./maxPupilMean, 'r');
shadedErrorBar(x,nanmean(p1SizeAll)./maxPupilMean,sem(p1SizeAll)./maxPupilMean,'lineprops','r-');
p2 = plot(x,nanmean(p1SizeMoves)./maxPupilMean, 'k');
shadedErrorBar(x,nanmean(p1SizeMoves)./maxPupilMean,sem(p1SizeMoves)./maxPupilMean,'lineprops','k-');
xlabel('pupil size, normalized by range')
legend([p1 p2], {'pupil', 'pupil + movement'})

figure; hold on
maxPupilMean = max(nanmean(p2SizeAll));
p1 = plot(x,nanmean(p2SizeAll)./maxPupilMean, 'r');
p2 = plot(x,nanmean(p2SizeMoves)./maxPupilMean, 'k');
shadedErrorBar(x,nanmean(p2SizeAll)./maxPupilMean,sem(p2SizeAll)./maxPupilMean,'lineprops','r-');
shadedErrorBar(x,nanmean(p2SizeMoves)./maxPupilMean,sem(p2SizeMoves)./maxPupilMean,'lineprops','k-');
legend([p1 p2], {'pupil', 'pupil + movement'})
ylabel('# of pupil events normalized by max number')
xlabel('pupil size, normalized by max')

% still vs moves pupil sizes
p1SizeStill = p1SizeAll-p1SizeMoves;
figure; hold on
sumPupilMean = sum(nanmean(p1SizeAll));
p1 = plot(x,(nanmean(p1SizeStill)./sumPupilMean).*100, 'r');
shadedErrorBar(x,(nanmean(p1SizeStill)./sumPupilMean).*100,(sem(p1SizeStill)./sumPupilMean).*100,'lineprops','r-');
p2 = plot(x,(nanmean(p1SizeMoves)./sumPupilMean)*100, 'k-');
shadedErrorBar(x,(nanmean(p1SizeMoves)./sumPupilMean).*100,(sem(p1SizeMoves)./sumPupilMean)*100,'lineprops','k-');
xlabel('pupil size, normalized by range')
legend([p1 p2], {'still', 'moving'})
ylabel('% of recorded pupil sizes')

p2SizeStill = p2SizeAll-p2SizeMoves;
figure; hold on
sumPupilMean = sum(nanmean(p2SizeAll));
p1 = plot(x,(nanmean(p2SizeStill)./sumPupilMean).*100, 'r');
shadedErrorBar(x,(nanmean(p2SizeStill)./sumPupilMean).*100,(sem(p2SizeStill)./sumPupilMean).*100,'lineprops','r-');
p2 = plot(x,(nanmean(p2SizeMoves)./sumPupilMean)*100, 'k-');
shadedErrorBar(x,(nanmean(p2SizeMoves)./sumPupilMean).*100,(sem(p2SizeMoves)./sumPupilMean)*100,'lineprops','k-');
xlabel('pupil size, normalized by max')
legend([p1 p2], {'still', 'moving'})
ylabel('% of recorded pupil sizes')