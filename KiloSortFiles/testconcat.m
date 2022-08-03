clear; close all
header = '160905_E';
files = {'dotmapping001','dotmapping002','rfori001'};

%% setup kilo sort dir

datpath = sprintf('/Volumes/Drobo2/USERS/Michele/testconcat/%s',header);
cd(datpath)

%% load headers and check
fields   = {'SamplingLabel','ChannelCount','SamplingFreq','TimeRes','DateTime','FileSpec','DataPoints'};
clear allheads  
for e = 1:length(files)
   
    clear ns6file ns6head
    ns6file  = [fpath header '_' files{e} '.ns6'];
    ns6head  = openNSx(ns6file,'noread');
   
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
           for e = 2:length(files)
               if ~all(strcmp(ref,allheads.ElectrodeLabel(:,e)))
                   pass = false;
               end
           end
    end
end

if ~pass
    error('cannot concatenate files')
end

TIMES = datenum(allheads.DateTime,'dd-mmm-yyyy HH:MM:SS');
[~,I]=sort(TIMES,'ascend');
files = files(I);

nch = allheads.ChannelCount(1);

%% chan map [WILL NEED TO CHECK THAT ALL FILES HAVE SAME CONFIG]
fpath = sprintf('/Volumes/%s/DATA/NEUROPHYS/rig%s/%s/%s','Drobo2','022',header);
ns6file = [fpath header '_' files{1} '.ns6'];
disp('Extracting header info and checking file size...')
out = NSHead2MapFile(ns6file);
save([datpath '/chanMap.mat'],'-struct','out')

%%
data = []; ftps = nan(length(files),2); 
for e = 1:length(files)
    ns6file = [fpath header '_' files{e} '.ns6'];
    x  = openNSxHL(ns6file);
    % NEED TO INSERT DATA CHECKS HERE
    data = [data; x];
    
    % also get ftps and fnames for later 
    ln = length(x) / nch;
    if e == 1
        ftps(e,1) = 1;
        ftps(e,2) = ln;
    else
        ftps(e,1) = ftps(e-1,2) + 1;
        ftps(e,2) = ftps(e,1) + ln - 1;
    end
    [~,BRdatafile,~] = fileparts(ns6file);
    fnames{e} = BRdatafile; 
end

if ~isequal(diff(ftps,[],2)+1,allheads.DataPoints')
    error('issue with data size')
end

% Opening the output file for saving
datafile = [datpath filesep header '.dat'];
FIDw = fopen(datafile, 'w+', 'ieee-le');

% Writing data into file
disp('Saving concatenated file...')
fwrite(FIDw, data, 'int16');
fclose(FIDw);
disp('done!')

% save concat info
save([datpath filesep header '.cat'],'fnames','ftps','allheads','-mat')


%% basicKiloSort
rez = basicKiloSort(datafile);

%% rez 2 ss
ss = KiloSort2SpikeStruct(rez,0);
sss = SplitSpikeStruct(ss,fnames,ftps,0); 