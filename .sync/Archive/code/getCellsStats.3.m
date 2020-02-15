function getCellsStats
cdPV
load('allPINPdirs.mat')
DIRS = PINPDIRS;

for d = 1: length(DIRS)
    if ~isempty(DIRS{d})
        cd(DIRS{d})
        load('dirs.mat')
        cd(dirs{1}) %cd to master dir
        
        
        
    end
end