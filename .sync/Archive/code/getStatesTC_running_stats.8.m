cdVIP
load('WNdirsVIP.mat')
DIRS = WNdirs;
cdPV
load('PVWNdirs.mat')
DIRS = [DIRS WNdirs];

behavior = [];
for d = 1:length(DIRS)
    if ~isempty(DIRS)
        cd(DIRS{d})
        [pupil_indx, motion_indx, pupil_trace, motion_trace, meanM2] = getStatesTC(0.5,d);
        behavior = [behavior; motion_indx meanM2];
    end
end
