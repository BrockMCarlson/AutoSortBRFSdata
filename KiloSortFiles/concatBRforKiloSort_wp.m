function [datafile, sucess] = concatBRforKiloSort(filelist,targetpath,flag_checkforexisting,depth)

if  nargin < 3
    flag_checkforexisting = false;
end
if nargin < 2 || isempty(targetpath)
    targetpath = '/Volumes/Drobo2/DATA/NEUROPHYS/KiloSort-ed/concatBR_default'; 
end


addpath(genpath('/volumes/Drobo/Users/Kacie/Code/fNPMK/NPMK'));

sucess = 1; 
%% setup

fnames = cell(size(filelist)); 
fpath = cell(size(filelist)); 
brdrname = [];
for i = 1:length(filelist)
    
    [pathstr,BRdatafile,~] = fileparts(filelist{i});
    
    if isempty(pathstr) 
        if isempty(brdrname)
            locations = {...
                'Drobo2','022';...
                'Drobo','022';...
                'Drobo','021'};
            for j = 1:size(locations, 1)
                temp = sprintf('/Volumes/%s/DATA/NEUROPHYS/rig%s/%s/',locations{j,1},locations{j,2},BRdatafile(1:8));
                if  exist([temp BRdatafile '.nev'],'file');
                    brdrname = temp;
                    break
                end
            end
        end
        pathstr = brdrname; 
    end
    
    fpath{i} = pathstr;
    fnames{i} = BRdatafile; 
end
    
header = unique(cellfun(@(x) x(1:8), fnames,'UniformOutput',0));
if length(header) > 1
    sucess = 0; 
    fprintf('\nfiles have diffrent headers\n') 
    return
end
datpath = sprintf('%s/%s',targetpath,header{1});
if ~exist(datpath,'dir') 
    mkdir(datpath)
end

% add depth label to file to differentiate it from files recorded at diff depths 
if nargin == 4
datafile = [datpath filesep header{1} depth '.dat'];
else
    datafile = [datpath filesep header{1} '.dat'];
end
if flag_checkforexisting && exist(datafile,'file')
    sucess = NaN;
    fprintf('\n found on disk %s',datafile)
    return
end

%% sort files in time
fields   = {'DateTime'};
clear allheads 
for e = 1:length(fnames)
   
    clear ns6file ns6head
    ns6file  = [fpath{e} filesep  fnames{e} '.ns6'];
    ns6head  = openNSx(ns6file,'noread');
    if isempty(ns6head)
        sucess = 0; 
        fprintf('\ncould not find ns6file: %s\n',ns6file)
        return
    end
    
    allheads.(fields{1}){1,e} = ns6head.MetaTags.(fields{1});
    
end
TIMES = datenum(allheads.DateTime,'dd-mmm-yyyy HH:MM:SS');
[~,I]=sort(TIMES,'ascend');
fnames = fnames(I);
fpath = fpath(I); 

%% load headers and check

% check if this is single electrode or multielectrode data:
for e = 1:length(fnames)
    
    clear ns6file ns6head labels 
    ns6file  = [fpath{e} filesep  fnames{e} '.ns6'];
    ns6head  = openNSx(ns6file,'noread');
    if isempty(ns6head)
        sucess = 0;
        fprintf('\ncould not find ns6file: %s\n',ns6file)
        return
    end
    
    labels = {ns6head.ElectrodesInfo.Label};
    if length(cell2mat(strfind(labels,'eA'))) > 1 | length(cell2mat(strfind(labels,'eB'))) > 1 | ...
            length(cell2mat(strfind(labels,'eC'))) > 1 | length(cell2mat(strfind(labels,'eD'))) > 1
        recordingtype = 'multi';
    else
        recordingtype = 'single';
    end
    
end

fields   = {'SamplingLabel','ChannelCount','SamplingFreq','TimeRes','DateTime','FileSpec','DataPoints'};

if strcmp(recordingtype,'single')
    
     clear allheads
    for e = 1:length(fnames)
        
        clear ns6file ns6head
        ns6file  = [fpath{e} filesep  fnames{e} '.ns6'];
        ns6head  = openNSx(ns6file,'noread');
        if isempty(ns6head)
            sucess = 0;
            fprintf('\ncould not find ns6file: %s\n',ns6file)
            return
        end
        
        clear allheads 
        padchanN = 3; 
        % cycle through pertenent fields
        for f = 1:length(fields)
            if strcmp(fields{f},'ChannelCount')
                allheads.(fields{f})(1,e) = ns6head.MetaTags.(fields{f}) + padchanN; 
            elseif iscell(ns6head.MetaTags.(fields{f}))
                allheads.(fields{f}){1,e} = ns6head.MetaTags.(fields{f}){1};
            elseif ischar(ns6head.MetaTags.(fields{f}))
                allheads.(fields{f}){1,e} = ns6head.MetaTags.(fields{f});
            else
                allheads.(fields{f})(1,e) = ns6head.MetaTags.(fields{f});
            end
        end
        
        % remake labels with white noise channels: 
        banks  = {ns6head.ElectrodesInfo.ConnectorBank}; 
        bank   = strcat('e',banks(cellfun('isempty',strfind(banks,'E')))); 
        h_labels = {ns6head.ElectrodesInfo.Label}'; 
        eID    = find(~cellfun('isempty',strfind(h_labels,bank))); 
        oID    = find(cellfun('isempty',strfind(h_labels,bank))); 
        labels{1} = h_labels{eID}(1:4); 
        for i = 2:padchanN+1
            labels{i} = sprintf('%s%02u',bank{1},i); 
        end
        ct = 1; 
        for i = (padchanN+1):(padchanN)+length(oID)
            ct = ct + 1; 
            labels{i} = h_labels{ct}; 
        end
        
        h_chanids  = double(ns6head.MetaTags.ChannelID);
        chanids    = nan(length(labels),1); 
        chanids(1) = h_chanids(eID); 
        chanids(padchanN+1:end) = h_chanids(oID); 
        
        

        allheads.('PausedFile')(1,e) = ns6head.RawData.PausedFile;
        allheads.('ElectrodeLabel')(:,e) = labels;
        allheads.('ChannelID')(:,e) = double(chanids);
        allheads.('NS6File'){e,1} = ns6file;
        
    end
    
else
    
    clear allheads
    for e = 1:length(fnames)
        
        clear ns6file ns6head
        ns6file  = [fpath{e} filesep  fnames{e} '.ns6'];
        ns6head  = openNSx(ns6file,'noread');
        if isempty(ns6head)
            sucess = 0;
            fprintf('\ncould not find ns6file: %s\n',ns6file)
            return
        end
        
        % cycle through pertenent fields
        for f = 1:length(fields)
            if iscell(ns6head.MetaTags.(fields{f}))
                allheads.(fields{f}){1,e} = ns6head.MetaTags.(fields{f}){1};
            elseif ischar(ns6head.MetaTags.(fields{f}))
                allheads.(fields{f}){1,e} = ns6head.MetaTags.(fields{f});
            else
                allheads.(fields{f})(1,e) = ns6head.MetaTags.(fields{f});
            end
        end
        
        allheads.('PausedFile')(1,e) = ns6head.RawData.PausedFile;
        allheads.('ChannelID')(:,e) = double(ns6head.MetaTags.ChannelID);
        allheads.('ElectrodeLabel')(:,e) = {ns6head.ElectrodesInfo.Label}';
        allheads.('NS6File'){e,1} = ns6file;
        
    end 

end

pass = true; 
for f = 1:length(fields)
    switch fields{f}
        case {'ChannelCount' 'SamplingFreq' 'TimeRes','SamplingLabel','FileSpec'}
           if length(unique(allheads.(fields{f}))) > 1
               pass = false; 
           end
           
        case 'ChannelID'
            if any(any(diff(allheads.ChannelID,[],2)))
                pass = false;
            end
            
        case 'ElectrodeLabel'
           ref = allheads.ElectrodeLabel(:,1); 
           for e = 2:length(fnames)
               if ~all(strcmp(ref,allheads.ElectrodeLabel(:,e)))
                   pass = false;
               end
           end
    end
end

if ~pass
    sucess = 0; 
    fprintf('\ncannot concatenate files\n')
    return
end


%% extract and save chan map 
ns6file  = [fpath{1} filesep  fnames{1} '.ns6'];
disp('Extracting header info and checking file size...')
out = NSHead2MapFile(ns6file,'single',allheads);
save([datpath '/chanMap.mat'],'-struct','out')
clear ns6file out 

%% main concatenation 

nch  = length({ns6head.ElectrodesInfo.Label}); 
data = [];
ftps = nan(length(fnames),2); 

for e = 1:length(fnames)
    clear ln; 
    fprintf('\nConcatenating file %u of %u...\n',e,length(fnames))

    ns6file  = [fpath{e} filesep  fnames{e} '.ns6'];
    x  = openNSxHL(ns6file);
    ln = length(x) / nch;
    
    % NEED TO INSERT DATA CHECKS HERE
    if strcmp(recordingtype,'single')
        pad  = [randn([(padchanN-1)*ln,1])]; 
        data = [data; cat(1,x(1:ln),pad,x(ln+1:end))]; 
        ln   = size(data,1); 
    else
        data = [data; x];
    end
    

    % also get ftps and fnames for later 
    if e == 1
        ftps(e,1) = 1;
        ftps(e,2) = ln;
    else
        ftps(e,1) = ftps(e-1,2) + 1;
        ftps(e,2) = ftps(e,1) + ln - 1;
    end

end

if strcmp(recordingtype,'multi')
if ~isequal(diff(ftps,[],2)+1,allheads.DataPoints')
    sucess = 0; 
    fprintf('\nissue with data size\n')
    return
end
end

%% Opening the output file for saving
FIDw = fopen(datafile, 'w+', 'ieee-le');

% Writing data into file
fprintf('\nSaving concatenated file...')
fwrite(FIDw, data, 'int16');
fclose(FIDw);
fprintf('done!\n')

% save concat info
save([datpath filesep 'concatInfo.mat'],'fnames','ftps','allheads','-mat')
