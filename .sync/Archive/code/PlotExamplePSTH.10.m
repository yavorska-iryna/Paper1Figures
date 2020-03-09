clear; close all; dbstop if error
variables_dir = 'C:\Users\lab\Resilio Sync\Paper1Figures\code\variables';

cd(variables_dir);
load('WNdataLaserOFF.mat'); % 'WNdataLaserOFF.mat (OFF1 -  longer responses), dynamic pupil threshold, but varified.
data = WNdataLaserOFF;
load('WNdataLaserON.mat');
data1 = WNdataLaserON;
load('CellsQualityStats.mat')
load('evoked_indx_epistatic.mat')
load('example_cells_epistatic.mat')

% good examples: 387
sigma = 10;
for cc = 1:length(cell_number1)
    c =  cell_number1(cc);
    if ~isnan(c)
        st = data(c).mST_off;
        [x, sitting_st] = GaussSmooth(st, sigma, [-300 900]);
        sitting_st = sitting_st./data(c). nrepsWNMoff;
        st = data(c).mST_on;
        [x, running_st] = GaussSmooth(st, sigma, [-300 900]);
        running_st = running_st./data(c). nrepsWNMon;
        st = data1(c).mST_off;
        [x, sitting_st_laser_on] = GaussSmooth(st, sigma, [-300 900]);
        sitting_st_laser_on = sitting_st_laser_on./data1(c). nrepsWNMoff;
        st = data1(c).mST_on;
        [x, running_st_laser_on] = GaussSmooth(st, sigma, [-300 900]);
        running_st_laser_on = running_st_laser_on./data1(c). nrepsWNMon;
        predicted_trace = nanmean([running_st' sitting_st_laser_on'], 2);
        figure; hold on; plot(x,sitting_st,'k', 'LineWidth', 2)
        plot(x,running_st, 'k--',  'LineWidth', 2)
        plot(x,sitting_st_laser_on, 'c',  'LineWidth', 2)
        plot(x,running_st_laser_on, 'c--',  'LineWidth', 2)
        plot(x,predicted_trace, 'r',  'LineWidth', 2)
        plot([0 0], [-5 150], 'm--')
        plot([-50 -50], [-10 150], 'c--')
        plot([-50 750], [-10 -10], ' c' , 'LineWidth', 7)
        plot([0 600], [-5 -5], ' m',  'LineWidth', 7)
        legend(' sitting laser off', ' running laser off', ' sitting laser on', ' running laser on', ' predicted')
        ylabel('Firing Rate (Hz)'); xlabel('time (ms)')
        xlim([-60 150]); ylim([-11 100])
    end
end