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
        CellsQualityStats(indx : (indx+length(good_cells)-1)) = table(uQ(good_cells), SNR, cids(good_cells), total_width, endslope);
        indx = indx+length(good_cells);
        
    end
end