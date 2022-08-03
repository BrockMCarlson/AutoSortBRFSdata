function [datafile, sucess,fnames,ftps,allheads] = concatBRforKiloSort(filelist,targetpath,flag_checkforexisting,depth)

if  nargin < 3
    flag_checkforexisting = false;
end
if nargin < 2 || isempty(targetpath)
    targetpath = '/Volumes/Drobo2/DATA/NEUROPHYS/KiloSort-ed/concatBR_default'; 
end


addpath(genpath('/volumes/Drobo/Users/Kacie/Code/fNPMK/NPMK'));

sucess = 1; 
fnames = []; 
ftps = []; 
allheads = []; 
%% setup

fnames = cell(size(filelist)); 
fpath  = cell(size(filelist)); 
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
datpath = sprintf('%s/%s%s',targetpath,header{1},depth);
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
    load([targetpath '/' header{1}(1:8) '/concatInfo.mat'])
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
TIMES  = datenum(allheads.DateTime,'dd-mmm-yyyy HH:MM:SS');
[~,I]  =sort(TIMES,'ascend');
fnames = fnames(I);
fpath  = fpath(I); 

%% load headers and check

fields   = {'SamplingLabel','ChannelCount','SamplingFreq','TimeRes','DateTime','FileSpec','DataPoints'};
clear allheads  
allheads.PauseNoDataPoints = cell(1,3);  
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
       if strcmp(fields{f},'DataPoints') && ns6head.RawData.PausedFile
           allheads.(fields{f})(1,e) = sum(ns6head.MetaTags.(fields{f}));
       elseif iscell(ns6head.MetaTags.(fields{f})) 
           allheads.(fields{f}){1,e} = ns6head.MetaTags.(fields{f}){1};
       elseif ischar(ns6head.MetaTags.(fields{f}))
           allheads.(fields{f}){1,e} = ns6head.MetaTags.(fields{f});
       else
           allheads.(fields{f})(1,e) = ns6head.MetaTags.(fields{f}); 
       end
   end
 
   allheads.('PausedFile')(1,e) = ns6head.RawData.PausedFile;
   allheads.('ChannelID')(:,e)  = double(ns6head.MetaTags.ChannelID);
   allheads.('ElectrodeLabel')(:,e) = {ns6head.ElectrodesInfo.Label}';
   allheads.('NS6File'){e,1} = ns6file; 
  
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
out = NSHead2MapFile(ns6file);
save([datpath '/chanMap.mat'],'-struct','out')
clear ns6file out 

%% main concatenation 

nch  = allheads.ChannelCount(1);
data = [];
ftps = nan(length(fnames),2); 

for e = 1:length(fnames)
    clear ln x ns6file segfiles 
    fprintf('\nConcatenating file %u of %u...\n',e,length(fnames))
     fprintf('\n%s\n',fnames{e})
    if ~ allheads.PausedFile(e)
        ns6file  = [fpath{e} filesep  fnames{e} '.ns6'];
        x        = openNSxHL(ns6file);
        % NEED TO INSERT DATA CHECKS HERE
        data = [data; x];
      
        ln = length(x) / nch;
    else
        ns6file  = [fpath{e} filesep  fnames{e} '.ns6'];
        % if the data haven't been segmented, segment
        if ~exist([ns6file(1:end-4) '-s001.ns6'],'file')
            splitNSxPauses(ns6file)
        end
        % find segmented files
        sln       = 0; 
        segfiles = dir([fpath{e} filesep  fnames{e} '-s*' '.ns6']);
        for s = 1:length(segfiles)
            clear ns6file x 
            ns6file  = [fpath{e} filesep  segfiles(s).name];
            x        = openNSxHL(ns6file);
            data     = [data; x];
            
            % also get ftps and fnames for later
            sln = sln + length(x) / nch; 
        end
        ln = sln; 
    end
   
        if e == 1
            ftps(e,1) = 1;
            ftps(e,2) = ln;
        else
            ftps(e,1) = ftps(e-1,2) + 1;
            ftps(e,2) = ftps(e,1) + ln - 1;
        end
    
end



if ~isequal(diff(ftps,[],2)+1,allheads.DataPoints')
    sucess = 0; 
    fprintf('\nissue with data size\n')
    return
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

