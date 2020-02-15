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
            [stimulus_trace, ~, ~] =load_open_ephys_data(sprintf('105_ADC%d.continuous', 1));
            [laser_trace, ~, ~] =load_open_ephys_data(sprintf('105_ADC%d.continuous', 3));
        catch
            [stimulus_trace, ~, ~] =load_open_ephys_data(sprintf('114_ADC%d.continuous', 1));
            [laser_trace, ~, ~] =load_open_ephys_data(sprintf('114_ADC%d.continuous', 3));
        end
        load('moves_trace1.mat');
        load('pupil_long_axis_normalized.mat');
         duration = length(stimulus_trace)/30e3;
        [pupil_trace, running_trace] = resampleTraces(1:duration/30e3, laxis1, moves_trace1);
        

        
        %normalize data for the plot
        stimulus_trace = stimulus_trace./(max(abs(stimulus_trace)));
        laser_trace = laser_trace./(max(abs(laser_trace)));
        if max(abs(moves_trace1))~=0
            running_trace = running_trace./max(abs(running_trace));
        end
        
        stimulus_indx = find(stimulus_trace > 0);
        
        
        %plot the figure
        figure; hold on;
        plot(running_trace, 'r')
        plot(pupil_trace+1, 'b')
        laser_start = stimlog(1).LaserStart/1000;
         laser_dur = stimlog(1).LaserStart/1000;
        for i = 1:length(Events)
            if strcmp(Events(i).type, 'whitenoise')
                start = Events(i). soundcard_trigger_timestamp_sec;
                dur = Evens(i).duration/1000;
                plot([start start+dur] [2 2], 'm-')
            end
            if Events(i).laser == 1
                
            end
        end
        
        plot(stimulus_indx, ones(length(stimulus_indx),1)*2, 'm')
    end
end
