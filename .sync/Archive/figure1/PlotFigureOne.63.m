clear; close all
cdPV
load('allWNdirs.mat')
DIRS = WNDIRS;
binwidth = 0.025;
for d = 10:10; %length(DIRS)
    if ~isempty(DIRS{d})
        cd(DIRS{d})
        load('dirs.mat'); load('notebook.mat'); load('Events.mat')
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
        figure; hold on; offset = 0; spike_traces = [];
         for c = 1:length(good_cells)
            st1 = sp.st(sp.clu==sp.cids(good_cells(c)));
            [st,recLength] = getSpikes(st1, DIRS{d});
            [ST, x] = hist(st, 0:binwidth:recLength);
            spike_trace = ST/binwidth; %convert to Hz
            spike_traces(c,:) = spike_trace./max(spike_trace);
            offset = offset+1;
            plot(x, smooth(spike_traces(c,:),5)+offset, 'k')
        end
        
        [pupil_trace, running_trace] = resampleTraces(1:recLength, laxis1, moves_trace1);
       
        %normalize data for the plot
        stimulus_trace = stimulus_trace./(max(abs(stimulus_trace)));
        laser_trace = laser_trace./(max(abs(laser_trace)));
        if max(abs(moves_trace1))~=0
            running_trace = running_trace./max(abs(running_trace));
        end
        
        %plot the figure
        figure; hold on;
        plot(smooth(running_trace,5), 'r')
        plot(pupil_trace*.7+.6, 'b')
        laser_start = stimlog(1).LaserStart/1000;
        laser_width = stimlog(1).LaserWidth/1000;
        for i = 1:length(Events)
             start = Events(i). soundcard_trigger_timestamp_sec;
             dur = Events(i).duration/1000;
            if strcmp(Events(i).type, 'whitenoise')
                plot([start start+dur], [1.5 1.5], 'm-', 'LineWidth', 4)
            else
                 plot([start start+dur], [1.5 1.5], '-', 'LineWidth', 4, 'Color', [.8 .8 .8])
            end
            if Events(i).laser == 1
                plot([start+laser_start start+laser_width], [1.45 1.45], 'c-', 'LineWidth', 4)
            end
        end %Events
        plot(x, smooth(spike_traces(46,:),5)+2.2, 'k')
        xlim([0 recLength])
        xlabel('time (sec)')
    end
end %dirs
