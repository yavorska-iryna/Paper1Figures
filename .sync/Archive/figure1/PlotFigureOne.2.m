clear; close all
cdPV
load('allWNdirs.mat')
DIRS = WNDIRS;
for d = 1:length(DIRS)
    if ~isempty(DIRS{d})
        cd(DIRS{d})
        load('dirs.mat')
        sp = loadKSdir(dirs{1});
        good_cells = find(sp.cgs == 2);
        spike_trace = ;
        pupil_trace = ;
        running_trace = ;
        stimulus_trace = ;
        laser_trace = ;
    end
end
