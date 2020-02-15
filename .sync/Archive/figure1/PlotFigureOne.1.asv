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
         try
            [scaledtrace1, timestamps, info] =load_open_ephys_data(sprintf('105_CH%d.continuous', 21));
        catch
            [scaledtrace1, timestamps, info] =load_open_ephys_data(sprintf('114_CH%d.continuous', 21));
         end
        load('moves_trace1.mat');
        load('pupil_long_axis_normalized.mat');
        [pupil_trace, running_trace] = resampleTraces(1:length(scaledtrace1)/30e3, laxis1, moves_trace1);
       [stimulus_trace, ~, ~] =load_open_ephys_data(sprintf('105_ADC%d.continuous', 3));
       [laser_trace, ~, ~] =load_open_ephys_data(sprintf('105_ADC%d.continuous', 3));

        stimulus_trace = ;
        laser_trace = ;
    end
end
