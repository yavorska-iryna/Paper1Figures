cdVIP
load('WNdirsVIP.mat')
DIRS = WNdirs;
cdPV
load('PVWNdirs.mat')
DIRS = [DIRS WNdirs];
close all

behavior = [];
for d = 1:length(DIRS)
    if ~isempty(DIRS{d})
        cd(DIRS{d})
        [pupil_indx, motion_indx, pupil_trace, motion_trace, meanM2] = getStatesTC(0.5,d);
        behavior = [behavior; motion_indx meanM2];
    end
end

behavior;
motion_indx = logical(behavior(:,1));
nanmean(behavior(motion_indx,2))
sem(behavior(motion_indx,2))
nanmean(behavior(~motion_indx,2))

