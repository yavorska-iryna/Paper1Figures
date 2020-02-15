function getCellsStats
cdPV
load('allPINPdirs.mat')
DIRS = PINPDIRS;

for d = 1: length(DIRS)
    if ~isempty(DIRS{d})
        cd(DIRS{d})
        load('dirs.mat')
        cd(dirs{1}) %cd to master dir
        load('WFs_stats.mat')
        load('SNR.mat')
        load('clusterQualityMetrics.mat')
        good_cells=find(cgs ==2);
        
        CellsQualityStats = table(uQ(good_cells), SNR, cids(good_cells), total_width, endslope);
        
    end
end