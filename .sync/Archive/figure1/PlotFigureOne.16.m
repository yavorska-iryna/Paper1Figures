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
        
        try
            [stimulus_trace, ~, ~] =load_open_ephys_data(sprintf('105_ADC%d.continuous', 1));
            [laser_trace, ~, ~] =load_open_ephys_data(sprintf('105_ADC%d.continuous', 3));
        catch
            [stimulus_trace, ~, ~] =load_open_ephys_data(sprintf('114_ADC%d.continuous', 1));
            [laser_trace, ~, ~] =load_open_ephys_data(sprintf('114_ADC%d.continuous', 3));
        end
        
        %normalize data for the plot
        stimulus_trace = stimulus_trace./(max(abs(stimulus_trace)));
        laser_trace = laser_trace./(max(abs(laser_trace)));
        if max(abs(moves_trace1))~=0
            running_trace = running_trace./max(abs(running_trace));
        end
        
        stimulus_indx = find(stimulus_trace > mean(stimulus_trace));
        
        
        %plot the figure
        figure; hold on;
        plot(running_trace, 'r')
        plot(pupil_trace+1, 'b')
        plot(stimulus_indx, ones(length(stimulus_indx),1)*2, 'm')
    end
end
