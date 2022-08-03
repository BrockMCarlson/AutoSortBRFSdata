clear 

addpath('/volumes/drobo/Users/Kacie/Code/fNPMK/NPMK');
addpath('/volumes/drobo/Users/Kacie/Code/fNPMK/NPMK/NSx Utilities');

tic
expIDs = {'mcosinteroc'};
% locations = {'Drobo','021';'Drobo2','021';...
%             'Drobo','022';'Drobo2','022'};
locations = {'LaCie'}; 
        
targetdir = '/Volumes/LaCie/KiloSort-ed_KD/';
flag_checkforexisting = true;

ct = 0; FILES = {};
for j = 1:size(locations, 1)
    

    %fpath = sprintf('/Volumes/%s/DATA/NEUROPHYS/rig%s/%s',locations{j,1},locations{j,2});

    fpath = '/volumes/LaCie/DATA_KD/';
    
    listing = dir(fpath);
    listing(find(~cellfun('isempty',strfind({listing.name},'_T')))) = []; 
   
    for l = 1:length(listing)
        
        if listing(l).isdir
            
            for e = 1:length(expIDs)
                
                datafiles = dir(sprintf('%s/%s/**%s**.ns6',fpath,listing(l).name,expIDs{e}));
                
                if isempty(datafiles)
                    continue
                end
                
                % exclude known problem sessions
                n = datenum(listing(l).name(1:6),'yymmdd');
                if strcmp('I',listing(l).name(end)) 
                    % some cherry picking of I34 sessions to avoid LGN
                    if n < datenum('160128','yymmdd')
                        continue
                    elseif any(strcmp(...
                            {'160226_I','160307_I','160308_I'},...
                            listing(l).name))
                        % DO NOT SKIP THESE DATES
                        % but these sessions had injections, so early files
                        % are the ones to look at, see Kacie's notes / DataLog
                    elseif n > datenum('160215','yymmdd')
                        continue
                    end
                elseif strcmp('151203_E',listing(l).name)...
                        ||strcmp('151223_E',listing(l).name)...
                        % quick exclusion of E sessions with bad chanel maps
                    continue
                end
                
                for d = 1:length(datafiles)
                    BRdatafile = datafiles(d).name(1:end-4);
                    ns6file = sprintf('%s/%s/%s.ns6',fpath,listing(l).name,BRdatafile);
                    datfile = sprintf('%s/%s/%s.ns6.dat',targetdir, BRdatafile,BRdatafile);
                    [datpath,~,~] = fileparts(datfile);
                    
                    % extract header info into channel map format,save in target directoy
                    if ~(flag_checkforexisting && exist([datpath '/chanMap.mat'],'file'))
                        
                        clear out
                        disp('Extracting header info and checking file size...')
                        out = NSHead2MapFile(ns6file);
                        
                        % check file length, exclude files less than 1.5 min long
                        dur = out.NS_Header.MetaTags.DataDurationSec / 60;
                        if dur < 1.5 
                            disp('Too short, skipping')
                            continue
                        end
                        
                        mkdir(datpath);
                        save([datpath '/chanMap.mat'],'-struct','out')
                                                
                    end

                    if ~(flag_checkforexisting && exist(datfile,'file'))
                        
                        % create binary data file
                        NSxToHL(ns6file);
                        %
                        % copy and save data to 'targetdir' destination
                        disp('copying...')
                        copyfile(...
                            [ns6file '.dat'],...
                            datfile);
                        delete([ns6file '.dat']); % BE VERY CAREFUL WITH THIS COMMAND! DO NOT ACCIDENTLY DELETE .NS6, MUST HAVE .dat 
                        toc
                        disp('done!')
                        
                    end
                    
                end
            end
        end
    end
end



