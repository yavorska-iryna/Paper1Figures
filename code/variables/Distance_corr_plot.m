load('Silence_DistanceCorr_M1shuffled.mat')
load('Silence_DistanceCorr_paper.mat')

y = DC_running - DC_running_shuffled;

% replace 0 with NaN, 0 mean no running in a trace
for i = size(y,2)
    for j = size(y,3)
        indx = find(y(:,i,j)==0)
        y(indx,i,j) = NaN;
    end
end

binwidth = [.05 0.1 .2 .4 .8 1.6 5 10 20 40];

figure; hold on
means  =  nanmean(squeeze(nanmean(y))');
SEMs = nanmean(squeeze(sem(y))');
plot([1:length(means)], means, 'ko');
errorbar([1:length(means)], means, SEMs)
ylabel('Distance Corr')
xlabel('bin time (sec)')
for i = 1:length(binwidth)
    time_labels{i} = num2str(binwidth(i));
end
xticklabels(time_labels)
xlim([0 length(means)+1])

figure; hold on
stds = nanmean(squeeze(nanstd(y))');
plot([1:length(means)], means./stds, 'ko');
ylabel('Coefficient of variation')
xlabel('bin time (sec)')
for i = 1:length(binwidth)
    time_labels{i} = num2str(binwidth(i));
end
xticklabels(time_labels)
xlim([0 length(means)+1])
