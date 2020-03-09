clear; close all; dbstop if error
variables_dir = 'C:\Users\lab\Resilio Sync\Paper1Figures\code\variables';

cd(variables_dir);
load('WNdataLaserOFF.mat'); % 'WNdataLaserOFF.mat (OFF1 -  longer responses), dynamic pupil threshold, but varified.
data = WNdataLaserOFF;
load('WNdataLaserON.mat');
data1 = WNdataLaserON;
load('CellsQualityStats.mat')
load ('C:\Users\lab\Resilio Sync\Paper1Figures\code\colorblind_colormap.mat'); % use color blind friendly colors
load('evoked_indx_epistatic.mat')


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
        figure; hold on; plot(x,sitting_st,'k')
        plot(x,running_st, 'k--')
        plot(x,sitting_st_laser_on, 'c')
        plot(x,running_st_laser_on, 'c--')
    end
end