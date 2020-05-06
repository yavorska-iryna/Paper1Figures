behavior = [];
for d = 1:length(DIRS)
    if ~empty(DIRS)
        cd(DIRS)
        [pupil_indx, motion_indx, pupil_trace, motion_trace, meanM2] = getStatesTC(thresh,d);
        behavior = [behavior motion_indx meanM2];
    end
end
