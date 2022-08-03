function ss = KiloSort2SpikeStruct(rez,flag_save2disk,fun_version)
% MAC, Dec 2016 / Jan 2017

crit = [0.3 0.9 3]; % [%mag, ms, chan])

if nargin < 2
    flag_save2disk = false;
end

if nargin < 3
    fun_version = 1;
end

% rez can be the actual var or the directory contining the rez var or the BRdatafile name
if ischar(rez)
    temp = fileparts(rez);
    if isempty(temp)
        fdir = sprintf('/Volumes/Drobo2/DATA/NEUROPHYS/KiloSort-ed/%s/',rez);
    else
        fdir = rez;
    end
    clear temp1 temp2 rez
    load([fdir filesep 'rez.mat'],'-mat');
else
    fdir = rez.ops.root;
end

%%
switch fun_version
    
    case 1
        
        % extract info from rez
        spikeTimes     = rez.st3(:,1)';
        spikeClusters  = 1+rez.st3(:,5)';
        spikeTemplates = rez.st3(:,2)';
        
        % get timing and channel info
        win = [-30:30];
        tmW  = win ./ rez.ops.fs * 1000;
        NchanTOT = rez.ops.NchanTOT;
        
        
        % get raw data around spiketimes
        % also get peak template channel
        fid = fopen(rez.ops.fbinary, 'r');
        
        peakSpikeCh = zeros(size(spikeClusters));
        WAVE = zeros(NchanTOT,numel(win),numel(spikeTimes));
        TEMP = zeros(sum(rez.connected),size(rez.dWU,1),numel(spikeTimes));
        for i = 1:length(spikeTimes)
            
            % binary start bit and duration
            st  = 2* (spikeTimes(i) -1 + win(1)) * NchanTOT;
            dr  = numel(win);
            if st < 0
                % spike to close to begining of file
                continue
            end
            
            % wave data from binary (only read in snips)
            clear dat
            fseek(fid,st,'bof');
            dat = fread(fid,[NchanTOT dr],'int16');
            
            if size(dat,2) < dr && (i > length(spikeTimes)-5)
                % spike is too close to end of file
                continue
            end
            
            WAVE(:,:,i) = dat;
            
            % template
            spktemplate = rez.dWU(:,:,spikeTemplates(i));
            TEMP(:,:,i) = spktemplate';
            
            % channel with the largest template for this spike
            [~,ch]    = max(range(spktemplate,1));
            peakSpikeCh(i) = ch;
            
        end
        fclose(fid);
        
        % organize data with chanMap, remove unconnected channels
        WAVE = WAVE(rez.ops.chanMap(rez.connected),:,:);
        
        % find channel index with maximum amplitude template/WAVE for each cluster
        % similar to peakSpikeCh, but for each cluster
        peakTempCh = zeros(size(spikeClusters));
        peakWaveCh = zeros(size(spikeClusters));
        spikeCh    = zeros(size(spikeClusters));
        spikeRank  = zeros(size(spikeClusters));
        
        uClusters  = unique(spikeClusters);
        clusterMap = zeros(length(uClusters),3); % cluster ID, Best Channel, Good Unit
        
        for c = 1:length(uClusters)
            clust     = uClusters(c);
            I         = spikeClusters == clust;
            
            % get template data, find max
            templates =  unique(spikeTemplates(I));
            t = squeeze(range((rez.dWU(:,:,templates)),1));
            m = max(max(t));
            if any(size(t) ==1)
                chidxT = find(t == m);
            else
                [chidxT, ~] = find(t == m);
            end
            peakTempCh(I) = chidxT;
            
            % determin cluster electrode from template data
            kcoords = rez.ops.kcoords;
            if numel(unique(kcoords)) > 1
                k = kcoords(chidxT);
            else
                k = unique(kcoords);
            end
            
            % get wave data
            clear chidx
            w = range(mean(WAVE(:,:,I),3),2);
            w(kcoords ~= k)  = 0;
            [~,chidxW]       = max(w);
            peakWaveCh(I)    = chidxW;
            
            % map cluster to a channel
            final_chidx = 0;
            if abs(diff([chidxW chidxT])) <= 1
                final_chidx = chidxW;
            else
                possible_chidx = intersect(...
                    unique(peakSpikeCh(I)),...
                    find(kcoords == k));
                N = histc(peakSpikeCh(I),possible_chidx);
                [~,idx]=max(N);
                chidx = possible_chidx(idx);
                if abs(diff([chidxW chidx])) <= 1 ...
                        || abs(diff([chidxT chidx])) <= 1
                    final_chidx = chidx;
                end
            end
            spikeCh(I) = final_chidx;
            
            % determin spatiotemporal extent of wave
            goodunit = 0;
            if final_chidx > 0
                clear wave stmap width span*
                wave = mean(WAVE(:,:,I),3);
                wave(kcoords ~= k,:) = 0;
                wave = bsxfun(@minus,wave,wave(:,1));
                wave(wave>0) = 0;
                wave  = abs(wave);
                stmap = wave >= crit(1)*max(max(wave));
                pix0  = sub2ind(size(wave),final_chidx, find(tmW == 0));
                % use imagetoobox
                if any(any(stmap))
                    cc    = bwconncomp(stmap);
                    idx   = cellfun(@(x) any(x == pix0),cc.PixelIdxList);
                    if any(idx)
                        stats = regionprops(cc,'BoundingBox','Centroid');
                        width = (stats(idx).BoundingBox(3:4)) ;
                        span  = width .* [(1 ./ rez.ops.fs * 1000) 1];
                        goodunit = span(1) <= crit(2) && span(2) <= crit(3);
                    end
                end
            end
            spikeRank(I) = goodunit;
            
            % save clusterMap
            clusterMap(c,1) = clust;
            clusterMap(c,2) = final_chidx;
            clusterMap(c,3) = goodunit;
        end
        
        ss.clusterMap    = clusterMap;
        ss.spikeTimes    = spikeTimes;
        ss.spikeClusters = spikeClusters;
        ss.spikeChannels = spikeCh; % best localized channel
        ss.spikeAutoRank = spikeRank;
        ss.spikeWaves    = WAVE;
        ss.spikeWavesTM  = tmW;
        ss.spikeTemps    = TEMP;
        
        % other info re: localization
        ss.peakSpikeCh   = peakSpikeCh; % defined per spike, based on template
        ss.peakTempCh    = peakTempCh;  % defined per cluster, based on template
        ss.peakWaveCh    = peakWaveCh;  % defined per cluster, based on wave data
        
        % stuff for recordkeeping
        ss.crit          = crit;
        ss.chanIDs       = rez.ops.chanMapLabels(rez.connected);
        ss.fbinary       = rez.ops.fbinary;
        [~,BRdatafile,~] = fileparts(rez.ops.fbinary);
        ss.BRdatafile    = BRdatafile(1:end-4);
        ss.Fs            = rez.ops.fs;
        ss.timestamp     = now;
        
        
    case 2
         
        addpath('/volumes/drobo/users/kacie/code/nbanalysis/FilterM_20Jul2011')
        % get timing and channel info
        win            = [-30:30];
        tmW            = win ./ rez.ops.fs * 1000;
        NchanTOT       = rez.ops.NchanTOT;
        spikeTimes     = rez.st3(:,1)';
        spikeTemplates = rez.st3(:,2)';
        spikeClusters  = 1+rez.st3(:,5)';
        uClusters      = unique(spikeClusters);
        clusterMap     = zeros(length(uClusters),3); % cluster ID, Best Channel, Good Unit
        spikeCh        = zeros(size(spikeClusters));
        spikeRank      = zeros(size(spikeClusters));
        
        hpc       = 1000;
        nyq       = rez.ops.fs./2;
        hWn       = hpc/nyq;
        [bwb,bwa] = butter(2,hWn,'high');
  
        for c = 1:length(uClusters)
            
            fprintf('%u clust of %u clusters\n\n',c,length(uClusters))
            clear peakWaveCh peakSpikeCh clustIDs clustTimes clust WAVE TEMP peakTempCh  peakTempCh spktemplate w  
            
            
            clust       = uClusters(c);
            clustIDs    = find(spikeClusters == clust);
            clustTimes  = spikeTimes(spikeClusters == clust);
            peakWaveCh  = zeros(size(clustIDs));
            peakSpikeCh = zeros(size(clustIDs));
            fid         = fopen(rez.ops.fbinary, 'r');
            
            WAVE        = zeros(NchanTOT,numel(win),numel(clustTimes));
            TEMP        = zeros(sum(rez.connected),size(rez.dWU,1),numel(clustTimes));
 
            
            for i = 1:length(clustTimes)
                
                % binary start bit and duration
                st  = 2* (clustTimes(i) -1 + win(1)) * NchanTOT;
                dr  = numel(win);
                if st < 0
                    % spike to close to begining of file
                    continue
                end
                
                % wave data from binary (only read in snips)
                clear dat
                fseek(fid,st,'bof');
                dat = fread(fid,[NchanTOT dr],'int16');
                
                if size(dat,2) < dr && (i > length(clustTimes)-20)
                    % spike is too close to end of file
                    continue
                end
               
                WAVE(:,:,i) = dat;
                
                % template
                spktemplate = rez.dWU(:,:,spikeTemplates(clustIDs(i)));
                TEMP(:,:,i) = spktemplate';
     
            end
            fclose(fid);

           % organize data with chanMap, remove unconnected channels//high pass filter data:
            [~,peakSpikeCh] = max(range(TEMP,2));
            peakSpikeCh     = squeeze(peakSpikeCh); 
            WAVE            = WAVE(rez.ops.chanMap(rez.connected),:,:);
            h_WAVE          = permute(WAVE,[2 1 3]); 

            WAVE            = FiltFiltM(bwb,bwa,WAVE,2); 
            spkN(c)         = size(WAVE,3);


            % find channel index with maximum amplitude template/WAVE for this cluster
            % similar to peakSpikeCh, but for each cluster
            peakTempCh = zeros(size(clustIDs));
            
            % get template data, find max
            templates =  unique(spikeTemplates(clustIDs));
            t = squeeze(range((rez.dWU(:,:,templates)),1));
            m = max(max(t));
            if any(size(t) ==1)
                chidxT = find(t == m);
            else
                [chidxT, ~] = find(t == m);
            end
            peakTempCh(:) = chidxT;
  
            % determine cluster electrode from template data
            kcoords = rez.ops.kcoords;
            if numel(unique(kcoords)) > 1
                k = kcoords(chidxT);
            else
                k = unique(kcoords);
            end
     
            % get wave data
            clear chidx
            w                    = range(nanmean(WAVE,3),2);
            w(kcoords ~= k)      = 0;
            [~,chidxW]           = max(w);
            peakWaveCh(:)        = chidxW;
      
            % map cluster to a channel
            final_chidx = 0;
            if abs(diff([chidxW chidxT])) <= 1
                final_chidx = chidxW;
            else
                possible_chidx = intersect(...
                    unique(peakSpikeCh),...
                    find(kcoords == k));
                N = histc(peakSpikeCh,possible_chidx);
                [~,idx]=max(N);
                chidx = possible_chidx(idx);
                if abs(diff([chidxW chidx])) <= 1 ...
                        || abs(diff([chidxT chidx])) <= 1
                    final_chidx = chidx;
                end
            end
            spikeCh(clustIDs) = final_chidx;
 
            % determin spatiotemporal extent of wave
            goodunit = 0;
            if final_chidx > 0
                clear wave stmap width span*
                wave = mean(WAVE,3);
                wave(kcoords ~= k,:) = 0;
                wave = bsxfun(@minus,wave,wave(:,1));
                wave(wave>0) = 0;
                wave  = abs(wave);
                stmap = wave >= crit(1)*max(max(wave));
                pix0  = sub2ind(size(wave),final_chidx, find(tmW == 0));
                % use imagetoobox
                if any(any(stmap))
                    cc    = bwconncomp(stmap);
                    idx   = cellfun(@(x) any(x == pix0),cc.PixelIdxList);
                    if any(idx)
                        stats = regionprops(cc,'BoundingBox','Centroid');
                        width = (stats(idx).BoundingBox(3:4)) ;
                        span  = width .* [(1 ./ rez.ops.fs * 1000) 1];
                        goodunit = span(1) <= crit(2) && span(2) <= crit(3);
                    end
                end
            end
            spikeRank(clustIDs) = goodunit;
            
            % save clusterMap
            clusterMap(c,1) = clust;
            clusterMap(c,2) = final_chidx;
            clusterMap(c,3) = goodunit;
            
            if any(final_chidx)
                save_spikes{c} = squeeze(WAVE(final_chidx,:,:));
            else
                save_spikes{c} = [];
            end
            
        end
  
        
        ss.clusterMap    = clusterMap;
        ss.spikeTimes    = spikeTimes;
        ss.spikeClusters = spikeClusters;
        ss.spikeChannels = spikeCh;        % best localized channel
        ss.spikeAutoRank = spikeRank;
        ss.spikeWaves    = save_spikes; 
        
        % other info re: localization
        ss.peakSpikeCh   = peakSpikeCh; % defined per spike, based on template
        ss.peakWaveCh    = peakWaveCh;  % defined per cluster, based on wave data
        %ss.peakTempCh   = peakTempCh;  % defined per cluster, based on template
                
        % stuff for recordkeeping
        ss.crit          = crit;
        ss.chanIDs       = rez.ops.chanMapLabels(rez.connected);
        ss.fbinary       = rez.ops.fbinary;
        [~,BRdatafile,~] = fileparts(rez.ops.fbinary);
        ss.BRdatafile    = BRdatafile(1:end-4);
        ss.Fs            = rez.ops.fs;
        ss.timestamp     = now;
        
        ss.spikeALL = nan(61,length(ss.spikeClusters));
        fake_wave   = nan(61,1);
        for c = 1:length(uClusters) 
            clear cidx; 
            cidx = find(ss.spikeClusters == uClusters(c)); 
            if clusterMap(c,2) == 0
                ss.spikeALL(:,cidx) =  repmat(fake_wave,[1 length(cidx)]); 
            else
                ss.spikeALL(:,cidx) =  ss.spikeWaves{c}; 
            end
        end
        
    case 3
  
        addpath('/volumes/drobo/users/kacie/code/nbanalysis/FilterM_20Jul2011')
        % get timing and channel info
        win            = -150:150;
        tmW            = win ./ rez.ops.fs * 1000;
        NchanTOT       = rez.ops.NchanTOT;
        spikeTimes     = rez.st3(:,1)';
        spikeTemplates = rez.st3(:,2)';
        spikeClusters  = 1+rez.st3(:,5)';
        uClusters      = unique(spikeClusters);
        clusterMap     = zeros(length(uClusters),3); % cluster ID, Best Channel, Good Unit
        spikeCh        = zeros(size(spikeClusters));
        spikeRank      = zeros(size(spikeClusters));
        
        hpc       = 1000;
        nyq       = rez.ops.fs./2;
        hWn       = hpc/nyq;
        [bwb,bwa] = butter(2,hWn,'high');
  
        % subsample waveforms and localize cluster to a channel first 
        for c = 1:length(uClusters)
            fprintf('%u clust of %u clusters\n\n',c,length(uClusters))
            clear peakWaveCh peakSpikeCh clustIDs clustTimes clust WAVE TEMP peakTempCh  peakTempCh spktemplate w  
            
            
            clust       = uClusters(c);
            clustIDs    = find(spikeClusters == clust);
            clustTimes  = spikeTimes(spikeClusters == clust);
            peakWaveCh  = zeros(size(clustIDs));
            peakSpikeCh = zeros(size(clustIDs));
            fid         = fopen(rez.ops.fbinary, 'r');
            
          
            TEMP        = zeros(sum(rez.connected),size(rez.dWU,1),numel(clustTimes));
            mxsamp      = min(numel(clustTimes),10000); 
            subsampidx  = floor(linspace(1,numel(clustTimes),mxsamp));  
            subsamp     = length(subsampidx); 
  
            somewaves   = zeros(NchanTOT,numel(win),subsamp);
            ct = 0; 
            for i = 1:length(clustTimes)
                
                % binary start bit and duration
                st  = 2* (clustTimes(i) -1 + win(1)) * NchanTOT;
                dr  = numel(win);
                if st < 0
                    % spike to close to begining of file
                    continue
                end
                
                % wave data from binary (only read in snips)
                clear dat
                fseek(fid,st,'bof');
                dat = fread(fid,[NchanTOT dr],'int16');
                
                if size(dat,2) < dr && (i > length(clustTimes)-20)
                    % spike is too close to end of file
                    continue
                end
                % template
                spktemplate = rez.dWU(:,:,spikeTemplates(clustIDs(i)));
                TEMP(:,:,i) = spktemplate';
                
                if any(ismember(subsampidx,i))
                    ct = ct + 1; 
                somewaves(:,:,ct) = dat;
                end
     
            end
            fclose(fid);
            
            somewaves       = somewaves(rez.ops.chanMap(rez.connected),:,:);
            somewaves       = FiltFiltM(bwb,bwa,somewaves,2); 
            
            % organize data with chanMap, remove unconnected channels//high pass filter data:
            [~,peakSpikeCh] = max(range(TEMP,2));
            peakSpikeCh     = squeeze(peakSpikeCh); 
            % find channel index with maximum amplitude template/WAVE for this cluster
            % similar to peakSpikeCh, but for each cluster
            peakTempCh = zeros(size(clustIDs));
            
            % get template data, find max
            templates =  unique(spikeTemplates(clustIDs));
            t = squeeze(range((rez.dWU(:,:,templates)),1));
            m = max(max(t));
            if any(size(t) ==1)
                chidxT = find(t == m);
            else
                [chidxT, ~] = find(t == m);
            end
            peakTempCh(:) = chidxT;
  
            % determine cluster electrode from template data
            kcoords = rez.ops.kcoords;
            if numel(unique(kcoords)) > 1
                k = kcoords(chidxT);
            else
                k = unique(kcoords);
            end
            
             % get wave data
            clear chidx
            w                    = range(nanmean(somewaves,3),2);
            w(kcoords ~= k)      = 0;
            [~,chidxW]           = max(w);
            peakWaveCh(:)        = chidxW;
      
            % map cluster to a channel
            final_chidx = 0;
            if abs(diff([chidxW chidxT])) <= 1
                final_chidx = chidxW;
            else
                possible_chidx = intersect(...
                    unique(peakSpikeCh),...
                    find(kcoords == k));
                N = histc(peakSpikeCh,possible_chidx);
                [~,idx]=max(N);
                chidx = possible_chidx(idx);
                if abs(diff([chidxW chidx])) <= 1 ...
                        || abs(diff([chidxT chidx])) <= 1
                    final_chidx = chidx;
                end
            end
            
            clusterMap(c,1) = clust;
            clusterMap(c,2) = final_chidx;
            
        end
        

        for c = 1:length(uClusters)
            
            fprintf('%u clust of %u clusters\n\n',c,length(uClusters))
            clear peakWaveCh peakSpikeCh clustIDs clustTimes clust WAVE TEMP peakTempCh  peakTempCh spktemplate w  
            
            fid         = fopen(rez.ops.fbinary, 'r');
            
            clust       = uClusters(c);
            clustIDs    = find(spikeClusters == clust);
            clustTimes  = spikeTimes(spikeClusters == clust);
            fid         = fopen(rez.ops.fbinary, 'r');
            
            WAVE        = zeros(numel(clustTimes),numel(win));
                  
            for i = 1:length(clustTimes)
                
                % binary start bit and duration
                st  = 2* (clustTimes(i) -1 + win(1)) * NchanTOT;
                dr  = numel(win);
                if st < 0
                    % spike to close to begining of file
                    continue
                end
                
                % wave data from binary (only read in snips)
                clear dat
                fseek(fid,st,'bof');
                dat = fread(fid,[NchanTOT dr],'int16');
                
                if size(dat,2) < dr && (i > length(clustTimes)-20)
                    % spike is too close to end of file
                    continue
                end
                if clusterMap(c,2) == 0 
                 
                else
                clear hwave 
                hwave      = dat(rez.ops.chanMap(rez.connected) == clusterMap(c,2),:); 
                WAVE(i,:)  =  hwave;   
                end
          
            end
            fclose(fid);


            WAVE            = FiltFiltM(bwb,bwa,WAVE,2); 
            spkN(c)         = size(WAVE,1);
%
 
            % determine spatiotemporal extent of wave
            final_chidx = clusterMap(c,2);
            goodunit = 0;
            %             if final_chidx > 0
            %                 clear wave stmap width span*
            %                 wave = mean(WAVE,2);
            %                 wave(kcoords ~= k,:) = 0;
            %                 wave = bsxfun(@minus,wave,wave(:,1));
            %                 wave(wave>0) = 0;
            %                 wave  = abs(wave);
            %                 stmap = wave >= crit(1)*max(max(wave));
            %                 pix0  = sub2ind(size(wave),final_chidx, find(tmW == 0));
            %                 % use imagetoobox
            %                 if any(any(stmap))
            %                     cc    = bwconncomp(stmap);
            %                     idx   = cellfun(@(x) any(x == pix0),cc.PixelIdxList);
            %                     if any(idx)
            %                         stats = regionprops(cc,'BoundingBox','Centroid');
            %                         width = (stats(idx).BoundingBox(3:4)) ;
            %                         span  = width .* [(1 ./ rez.ops.fs * 1000) 1];
            %                         goodunit = span(1) <= crit(2) && span(2) <= crit(3);
            %                     end
            %                 end
            %             end
            clear wave stmap widths stats cc 
            crit  = 0.25;
            win   = -150:150;
            wave  = nanmean(WAVE,1);
            wave  = bsxfun(@minus,wave,mean(wave));
            wave(wave > 0) = 0;
            wave  = abs(wave);
            stmap = wave >= crit(1) * max(max(wave));
            %use imagetoolbox
            if any(any(stmap))
                cc     = bwconncomp(stmap);
                stats  = regionprops(cc,'BoundingBox','Centroid');
                widths = cell2mat({stats.BoundingBox}');
                if (length(stats) < 3) && ...
                        (all(widths(:,3) > 2 & widths(:,3) < 15))
                    goodunit = 1;
                end
            end
            spikeRank(clustIDs) = goodunit;
       
            clusterMap(c,3) = goodunit;
            
            if any(final_chidx)
                save_spikes{c} = WAVE;
            else
                save_spikes{c} = [];
            end

            % save clusterMap

            clusterMap(c,3) = goodunit;
            
        end
        
        
    
        ss.clusterMap    = clusterMap;
        ss.spikeTimes    = spikeTimes;
        ss.spikeClusters = spikeClusters;
        ss.spikeChannels = spikeCh;        % best localized channel
        ss.spikeAutoRank = spikeRank;
        ss.spikeWaves    = save_spikes; 
       
        % stuff for recordkeeping
        ss.crit          = crit;
        ss.chanIDs       = rez.ops.chanMapLabels(rez.connected);
        ss.fbinary       = rez.ops.fbinary;
        [~,BRdatafile,~] = fileparts(rez.ops.fbinary);
        ss.BRdatafile    = BRdatafile(1:end-4);
        ss.Fs            = rez.ops.fs;
        ss.timestamp     = now;
        
        ss.spikeALL = nan(length(win),length(ss.spikeClusters));
        fake_wave   = nan(length(win),1);
        for c = 1:length(uClusters) 
            clear cidx; 
            cidx = find(ss.spikeClusters == uClusters(c)); 
            if clusterMap(c,2) == 0
                ss.spikeALL(:,cidx) =  repmat(fake_wave,[1 length(cidx)]); 
            else
                ss.spikeALL(:,cidx) =  ss.spikeWaves{c}'; 
            end
        end
        
end


if flag_save2disk
    % save to disk for future use
    save([fdir filesep 'ss.mat'],'ss','-v7.3');
end



