% collect statistics about each cell (SNR, spike width, etc)
%ira 01.21.2020

cdPV
load('allPINPdirs.mat')
DIRS = PINPDIRS;
indx = 1;
for d = 1: length(DIRS)
    if ~isempty(DIRS{d})
        cd(DIRS{d})
        load('dirs.mat')
        cd(dirs{1}) %cd to master dir
        load('WFs_stats.mat')
        load('SNR.mat')
        load('clusterQualityMetrics.mat')
        good_cells=find(cgs ==2);
        uQ = uQ(good_cells); uQ = uQ(:); SNR = SNR(:); cids = cids(good_cells); cids = cids(:); total_width = total_width(:); endslope = endslope(:);
        CellsQualityStats_dir = table(uQ, SNR, cids, total_width, endslope);
        if d ==1
            CellsQualityStats =CellsQualityStats_dir;
        else
            CellsQualityStats = [CellsQualityStats; CellsQualityStats_dir];
        end
    end
end
save_dir = 'C:\Users\lab\Resilio Sync\Paper1Figures\code\variables';
cd(save_dir)
save('CellsQualityStats.mat', 'CellsQualityStats')
