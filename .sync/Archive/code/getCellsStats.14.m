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
        uQ = uQ(:); SNR = SNR(:); cids = cids(:); total_width = total_width(:); endslope = endslope(:);
        CellsQualityStats_dir = table(uQ(good_cells), SNR, cids(good_cells), total_width, endslope);
        CellsQualityStats = [CellsQualityStats CellsQualityStats_dir];
        
    end
end