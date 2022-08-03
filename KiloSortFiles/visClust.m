function h = visClust(ss,uClusters)
% MAC, Dec 2016

% ss is the output of KiloSort2SpikeStruct(rez);
% uClusters is an optional input of cluster id numbers
% without uClusters speficied, code plots all clusters. 

if nargin < 2
    uClusters = unique(ss.spikeClusters);
end

h = [];
for c = 1:length(uClusters)
    h(c) = figure;
    
    clust = uClusters(c);
    J      = ss.clusterMap(:,1) == clust;
    good   = ss.clusterMap(J,3);
    eidx   = ss.clusterMap(J,2);
    if eidx > 0
        elabel = ss.chanIDs{eidx};
    else
        elabel = 'n/a';
    end
     
    I = ss.spikeClusters == clust;
    spks = ss.spikeTimes(I);
    
    wave = ss.spikeWaves(:,:,I);
    wave = mean(wave,3);
    wave = bsxfun(@minus,wave,wave(:,1));
    offset = 0.5 * max(range(wave,1)) * [1:size(wave,1)]';
    wave = bsxfun(@plus,wave,offset);
    offset = sort(offset);
    
    subplot(3,2,[1 3 5])
    plot(ss.spikeWavesTM,wave')
    axis tight
    set(gca,'YTick',offset(1:5:end),'YTickLabel',ss.chanIDs(1:5:end),...
        'Box','off','TickDir','out');
    title(sprintf('Clust %u Waveform\nel = %s   pass = %u',clust,elabel,good),'interpreter','none')
    xlabel('Time (ms)')
    if mod(length(ss.chanIDs),22) || mod(length(ss.chanIDs),24)
        set(gca,'YDir','reverse')
    end
    
    % corelation of time v. number spks/sec
    subplot(3,2,2); cla
    tm = spks(1):ss.Fs:spks(end); % 1 second bins
    ct = histc(spks,tm);
    tm(end) = []; ct(end) = [];    
    [r, p] = corrcoef(tm,ct);
    r = r(2); p = p(2);
    bar(tm ./ ss.Fs,ct,'histc');
    box off; axis tight; 
    title(sprintf('Number of Spike over Time\n r = %0.2f    p = %0.3f',r,p));
    ylabel('spk/s');
    
    subplot(3,2,4); cla
    scatter(ss.peakSpikeCh(I),ss.peakWaveCh(I)); hold on
    lim = [min([ss.peakSpikeCh(I);ss.peakWaveCh(I)]) max([ss.peakSpikeCh(I);ss.peakWaveCh(I)])];
    axis equal;
    set(gca,...
        'YTick',[lim(1):lim(2)],'YTickLabel',ss.chanIDs(lim(1):lim(2)),...
        'XTick',[lim(1):lim(2)],'XTickLabel',ss.chanIDs(lim(1):lim(2)),...
        'Box','off','TickDir','out');
    plot(lim,lim,'k');
    title('peakSpikeCh x peakWaveCh');

    
    subplot(3,2,6)
    scatter(ss.peakSpikeCh(I),ss.peakTempCh(I)); hold on
    lim = [min([ss.peakSpikeCh(I);ss.peakWaveCh(I)]) max([ss.peakSpikeCh(I);ss.peakWaveCh(I)])];
    axis equal;
    set(gca,...
        'YTick',[lim(1):lim(2)],'YTickLabel',ss.chanIDs(lim(1):lim(2)),...
        'XTick',[lim(1):lim(2)],'XTickLabel',ss.chanIDs(lim(1):lim(2)),...
        'Box','off','TickDir','out');
    plot(lim,lim,'k');
    title('peakSpikeCh x peakTempCh');
    
    
    
end

% uClusters = ss.clusterMap(:,1);
% tmW = ss.spikeWavesTM;
% WAVE = ss.spikeWaves;
% for c = 1:length(uClusters);
% 
% clust     = uClusters(c);
% I         = ss.spikeClusters == clust;
% ss.clusterMap(c,:)
% 
% clf
% subplot(1,2,1)
% clear wave
% wave = mean(WAVE(:,:,I),3);
% wave = bsxfun(@minus,wave,wave(:,1));
% offset = 0.5 * max(range(wave,1)) * [1:size(wave,1)]';
% wave = bsxfun(@plus,wave,offset);
% offset = sort(offset);
% plot(tmW,wave'); hold on
% axis tight
% set(gca,'YTick',offset(1:5:end),'YTickLabel',1:5:length(offset),...
%     'Box','off','TickDir','out');
% title(sprintf('clust = %u\nMean WAVE',clust))
% 
% subplot(1,2,2)
% clear wave stmap width span_*
% wave = mean(WAVE(:,:,I),3);
% wave = bsxfun(@minus,wave,wave(:,1));
% wave(wave>0) = 0;
% wave  = abs(wave);
% stmap = wave > 0.3*max(max(wave));
% imagesc(stmap); set(gca,'Ydir','normal'); title(sprintf('clust = %u\nThresholded Map',clust))
% pause
% end