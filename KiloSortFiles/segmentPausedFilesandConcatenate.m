function segmentPausedFilesandConcatenate(filelist)

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


for e = 1:length(fnames)
    clear ns6file ns6head
    ns6file  = [fpath{e} filesep  fnames{e} '.ns6'];
    ns6head  = openNSx(ns6file,'noread');
    if isempty(ns6head)
        fprintf('\ncould not find ns6file: %s\n',ns6file)
        return
    end
    
    if ns6head.RawData.PausedFile
        if ~exist([ns6file(1:end-4) '-s001.ns6'],'file')
            splitNSxPauses(ns6file)
        end
        [splitfiles] = dir([fpath{i} fnames{e} '-s' '*']);
        
    end
end

