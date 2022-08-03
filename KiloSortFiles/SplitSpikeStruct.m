
function ss_split = SplitSpikeStruct(ss,fnames,ftps,fun_version,flag_save2disk,flag_savespikeWaves)

% KiloSort Utils
% splits up concatenated data
% Jan 2017
% MAC

rcrit = 0.4;

if nargin < 5
    flag_save2disk = true;
    flag_savespikeWaves = false; 
elseif nargin == 5
     flag_savespikeWaves = false; 
end

if isempty(fun_version)
    fun_version = 1; 
end

if ischar(ss)
    temp = fileparts(ss);
    if isempty(temp)
        fdir = sprintf('/Volumes/Drobo2/DATA/NEUROPHYS/KiloSort-ed/concatBR_tuning/%s/',rez);
    else
        fdir = ss;
    end
    clear temp ss
    load([fdir filesep 'ss.mat'],'-mat');
else
    fdir = fileparts(ss.fbinary);
end

if nargin == 1 || isempty(fnames) || isempty(ftps)
    load([fdir filesep 'concatInfo.mat'],'fnames','ftps')
end

% check each clusters's stabilty over time
uClusters  = unique(ss.spikeClusters);
stableClust  = NaN(size(uClusters))';
CorrResult = NaN(length(size(uClusters)),2);
for c = 1:length(uClusters)
    
    clust = uClusters(c);
    I = ss.spikeClusters == clust;
    spks = ss.spikeTimes(I);
    
    % get time v. number spks/sec
    tm = spks(1):ss.Fs:spks(end); % 1 second bins
    ct = histc(spks,tm);
    tm(end) = []; ct(end) = [];
    
    
    % corelation of time v. number spks/sec
    [r, p] = corrcoef(tm,ct);
    r = r(2); p = p(2);
    
    %         figure; plot(tm,ct,'o')
    %     title(sprintf('clust = %u\n r = %0.2f    p = %0.3f',clust,r,p))
    
    % exclude it corelation is greater than 0.4
    if abs(r) <= rcrit
        stableClust(c) = 1;
    else
        stableClust(c) = 0;
    end
    CorrResult(c,:) = [r p];
    
end

fields = {'clusterMap','crit','chanIDs','timestamp','fbinary','Fs'};
clear ss_split

ss_split.fnames = fnames;
for j = 1:length(fields)
    
    switch fields{j}
        case 'clusterMap'
            ss_split.(fields{j}) = [ss.(fields{j}) stableClust];
        case 'crit'
            ss_split.(fields{j}) = [ss.(fields{j}) rcrit];
        case 'timestamp'
            ss_split.(fields{j}) = [ss.(fields{j}) now];
        otherwise
            ss_split.(fields{j}) = ss.(fields{j});
    end
end

for i = 1:length(fnames)
    
    I =   ss.spikeTimes >=  ftps(i,1)...
        & ss.spikeTimes <=  ftps(i,2);
    
    ss_split.(['x' fnames{i}]).('spikeTimes')    =  ss.('spikeTimes')(I); 
    
    if i > 1
        ss_split.(['x' fnames{i}]).('spikeTimes')    =  ss_split.(['x' fnames{i}]).('spikeTimes') - ftps(i-1,2); %correct times so that they are in reference to the original, single file
    end
    
     ss_split.(['x' fnames{i}]).('spikeClusters') = ss.('spikeClusters')(I);
    if flag_savespikeWaves & fun_version == 2
        ss_split.(['x' fnames{i}]).('spikeWaves')    = ss.('spikeWaves');
        ss_split.(['x' fnames{i}]).('spikeWavesSTD') = ss.('spikeWavesSTD');
        %ss_split.(['x' fnames{i}]).('spikeWavesCI')  = ss.('spikeWavesCI');
        ss_split.(['x' fnames{i}]).('spikeWavesN')   = ss.('spikeWavesN');
        ss_split.(['x' fnames{i}]).('spikeWavesTM')  = ss.('spikeWavesTM');
    else
        ss_split.(['x' fnames{i}]).('spikeWaves') = ss.('spikeWaves')(:,:,I);
    end
end

if flag_save2disk
    % save to disk for future use
    fdir  = fileparts(ss.fbinary);
    fname = [fdir filesep 'ss_split.mat'];
    save(fname,'ss_split');
end