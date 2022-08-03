% testrun
clear; close all;

% define file
BRdatafile = '160128_I_dotmapping001';
datafile = sprintf('/Volumes/Drobo2/USERS/Michele/KiloSort-ed Data/%s/%s.ns6.dat',BRdatafile,BRdatafile);
fpath = fileparts(datafile);

% set confguration
load(fullfile(fpath, 'chanMap.mat'),'nprobes')
run('/Users/coxm/Documents/ephys-analysis/KiloSort Utils/config.m');

% preprocess 
[rez, DATA, uproj] = preprocessData(ops); % preprocess data and extract spikes for initialization

%error('stop here to before fitting templates')

% sort!
rez                = fitTemplates(rez, DATA, uproj);  % fit templates iteratively
rez                = fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)

% save python results file for Phy
rezToPhy(rez, fpath);

%error('stop here to check in PHY before auto merge')

% AUTO MERGES
rez = merge_posthoc2(rez);

% save python results file for Phy
rezToPhy(rez, fpath);

% save and clean up
save(fullfile(fpath,  'rez.mat'), 'rez', '-v7.3');

% % remove temporary file
delete(ops.fproc);
%%

out = KiloSort2SpikeStruct(rez);
%%
uClusters = unique(out.spikeClusters);
for c = 1:length(uClusters)
    clust = uClusters(c);
    I = out.spikeClusters == clust;
    
    wave = out.spikeWaves(:,:,I);
    wave = mean(wave,3);
    wave = bsxfun(@minus,wave,wave(:,1));
    offset = 0.5 * max(range(wave,1));
    wave = bsxfun(@plus,wave,[1:size(wave,1)]' * offset);
    
    
    figure;
    subplot(2,2,[1 3])
    plot(wave')
    axis tight
    title(clust)
    legend(fliplr(out.chanIDs),'Location','BestOutside')
    box off
    
    subplot(2,2,2)
    scatter(out.peakSpikeCh(I),out.peakWaveCh(I));
    axis tight; hold on
    plot(xlim,xlim,'k');
    title('peakSpikeCh x peakWaveCh');
    box off
    
    subplot(2,2,4)
    scatter(out.peakSpikeCh(I),out.peakTempCh(I));
    axis tight; hold on
    plot(xlim,xlim,'k');
    title('peakSpikeCh x peakTempCh');
    box off


    
   

end





