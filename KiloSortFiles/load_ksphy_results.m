clear

addpath(genpath('/users/kacie/documents/matlab/npy-matlab-master/')) % used to get kilo/phy results into matlab
addpath(genpath('/users/kacie/documents/matlab/spikes-master'))      % from cortex lab on github 


dates = {'190326_B_p01'}; 

for d = 1:length(dates)
    
clearvars -except dates d 

myKsDir     = ['/Volumes/Drobo2/data/neurophys/KiloSort-ed/LGNgrouped/' dates{d} ];

filename    = [myKsDir '/spike_clusters.npy']; 

cluster_dat = readNPY(filename); 

filename    = [myKsDir '/spike_times.npy'];

times_dat   = readNPY(filename); 

% sp        = loadKSdir(myKsDir);

[~,~,~,ch]  = ksDriftmap(myKsDir); % ch is ordered index not label on probe (i think). 

% get channel and probe info and filenames:  
clusters            = double(unique(cluster_dat));
good_               = zeros(length(clusters),1); 
clusterMap          = nan(length(clusters),3);

filename            = [myKsDir '/cluster_groups.csv'];

if exist(filename,'file')
id_dat              = readtable(filename);
    for c  = 1:length(clusters)
    good_(c)        = strcmp(id_dat.group(id_dat.cluster_id == clusters(c)),'good');
    end
else
   good_            = ones(length(clusters),1); 
end

for c    = 1:length(clusters)
    if clusters(c) == 0, continue, end
    
    chan            = double(mode(ch(cluster_dat == clusters(c))));
    clusterMap(c,:) = [clusters(c) chan good_(c)]; 
    
end

load([myKsDir '/chanMap.mat'],'probetype','chanMapLabels','nprobes'); 
load([myKsDir '/concatInfo.mat']); 

ss.clusterMap    = clusterMap;
ss.probetype     = probetype; 
ss.nprobes       = nprobes; 
ss.chanLabels    = chanMapLabels; 
ss.ftps          = ftps; 
ss.headers       = allheads; 
ss.fnames        = fnames; 
ss.spikeClusters = double(cluster_dat); 
ss.spikeTimes    = double(times_dat); 

% split into separate files: 

for f = 1:length(ss.fnames)
   
    ss.(['x' fnames{f}]).clusterMap = clusterMap; 
    I = ss.spikeTimes >= ftps(f,1)...
        & ss.spikeTimes <= ftps(f,2); 
    
    ss.(['x' fnames{f}]).('spikeTimes')      = ss.('spikeTimes')(I); 
    if f > 1
         ss.(['x' fnames{f}]).('spikeTimes') = ss.(['x' fnames{f}]).('spikeTimes')- ftps(f-1,2);    
    end
    ss.(['x' fnames{f}]).('spikeClusters')   = ss.('spikeClusters')(I); 
    
    
end

% get waveforms 

gwfparams.dataDir       = myKsDir;                       % KiloSort/Phy output folder
gwfparams.fileName      = [myKsDir(end-11:end) '.dat'];  % .dat file containing the raw
gwfparams.dataType      = 'int16';                       % Data type of .dat file (this should be BP filtered)
gwfparams.nCh           = size(ss.chanLabels,1);         % Number of channels that were streamed to disk in .dat file
gwfparams.wfWin         = [-40 40];                      % Number of samples before and after spiketime to include in waveform
gwfparams.nWf           = 3000;                          % Number of waveforms per unit to pull out
gwfparams.spikeTimes    = times_dat;                     % Vector of cluster spike times (in samples) same length as .spikeClusters
gwfparams.spikeClusters = cluster_dat;                   % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes
ss.wf                   = getWaveForms(gwfparams);       % .dat is not HF filtered ! 
ss.wf.wfWin             = gwfparams.wfWin;  


save(['/volumes/bigdrobo2/drobo2/data/doughek/' myKsDir(end-11:end) '_ss.mat'],'ss','-v7.3'); 
end

