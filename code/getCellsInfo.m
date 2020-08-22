function  [WF, SNR, uQ] = getCellsInfo(directory)

% collects mean wafeform, SNR and uQ from directories

cd(directory)
load('dirs.mat')
cd(dirs{1})
load('WFs.mat')
load('SNR.mat')
SNR = SNR';
load('clusterQualityMetrics.mat', 'uQ', 'cgs')
uQ =  uQ(cgs ==2);


    
    

