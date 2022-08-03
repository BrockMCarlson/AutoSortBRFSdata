% setup LGN kilsort

clear
addpath('/volumes/drobo/Users/Kacie/Code/fNPMK/NPMK');
addpath('/volumes/drobo/Users/Kacie/Code/fNPMK/NPMK/NSx Utilities');
addpath(genpath('/users/kacie/documents/Github'));
rng('default');
flag_checkforexisting  = 0;
DataList               = importDataList('LGN-V1 I34');
targetpath             = '/Volumes/Drobo2/DATA/NEUROPHYS/KiloSort-ed/LGNgrouped';
bigdrobopath           = '/volumes/BigDrobo2/Drobo2/data/';
dates                  = unique(DataList.Datestr(:));
missing                = []; 
specific_date          = find(ismember(dates,{'190326'}));

for d =  specific_date' 
    
    
    fprintf('d : %u\n',d); 
    
    clear these_idfilelist
    these_id               =  find(ismember(DataList.Datestr,dates{d}));
    for p =[these_id]'
        
        clearvars -except DataList these_id p d targetpath flag_checkforexisting dates monocdates  bigdrobopath missing extHDpath
        
        %try
         depth = sprintf('_p%02u',sum(ismember(DataList.Datestr(1:p-1),DataList.Datestr(p))) + 1);
        if DataList.Drobo(p) == 1
            drname = sprintf('/volumes/drobo/data/neurophys/rig%03u/%s_%s/',DataList.Rig(p),DataList.Datestr{p},DataList.Monkey{p});
        elseif DataList.Drobo(p) == 2 & (datenum(dates{d},'yymmdd') < datenum('05/05/2018','mm/dd/yyyy'))
            drname = sprintf('/volumes/drobo2/data/neurophys/rig%03u/%s_%s/',DataList.Rig(p),DataList.Datestr{p},DataList.Monkey{p});
            
        else
            drname = sprintf('%s/rig%03u/%s_%s/',bigdrobopath,DataList.Rig(p),DataList.Datestr{p},DataList.Monkey{p});
        end
        
            
            cd(drname);
            file_head = {'posdisparitydrft','brfs','conedrft','color','cinterocdrft','rforidrft','darkrest','evpL','evpR','rfori','evp','disparitydrft','coneinterocdrft','colorflicker','bwflicker'};
            ct = 0;
        
            for each = 1:length(file_head)
                
                clear BRdatafile ns6file datfile datpath out dur rez ss filenum
                if ~isfield(DataList,file_head{each})
                    continue
                end
                if isempty(cell2mat(DataList.(file_head{each})(p)))
                    continue
                end
                
                filenum = cell2mat(DataList.(file_head{each})(p));
                for ff = 1:length(filenum)
                    ct = ct + 1;
                    BRdatafile = strcat(DataList.Datestr{p},'_',DataList.Monkey{p},sprintf('_%s%03u',file_head{each},filenum(ff)));
                    ns6file    = sprintf('%s/%s.ns6',drname,BRdatafile);
                    filelist{ct} = ns6file;
                end
                
            end

            % check for file pauses. if there are file pauses, segment
           
            [datafile, success,fnames,ftps,allheads] = concatBRforKiloSort(filelist,targetpath,flag_checkforexisting,depth);

            
            if success == 1 || isnan(success)
                [~,pfname] = fileparts(datafile);
                % basicKiloSort
                rez = basicKiloSort(datafile);
                if ~exist(strcat(targetpath,'/',pfname(1:8),depth,'/'))
                    mkdir(strcat(targetpath,'/',pfname(1:8),depth,'/'))
                end
                save(strcat(targetpath,'/',pfname(1:8),depth,'/rez.mat'),'rez','-v7.3');
            end

            %rez 2 ss
            ss  = KiloSort2SpikeStruct(rez,0,3);
            saveoutpath  = '/volumes/drobo2/data/neurophys/doughek/kiloout/LGNgrouped_waves/';
            
            if ~exist(saveoutpath), mkdir(saveoutpath), end;
            save(strcat(saveoutpath,pfname,'_kilo','_presplit_ss','.mat'),'ss','-v7.3');
            save(strcat(saveoutpath,pfname,'_kilo','_concatInfo_ss','.mat'),'fnames','ftps','allheads','-v7.3');
    end
    
end
     
% figure,
% for c = 1:length(save_spikes)
%     if ~isempty(save_spikes{c})
%          plot(win,nanmean(save_spikes{c},1)); hold on; 
%     end
% end
%%
% for d = 1 %:length(dates)
%     fprintf('d : %u\n',d); 
%     
%     clear these_idfilelist ss pfname depth 
%     these_id               =  find(ismember(DataList.Datestr,dates{d}));
%     for p = [these_id]'
%         
%         clearvars -except DataList theseid p d targetpath flag_checkforexisting dates monocdates  bigdrobopath missing extHDpath
%         
%         %try
%         if DataList.Drobo(p) == 1
%             drname = sprintf('/volumes/drobo/data/neurophys/rig%03u/%s_%s/',DataList.Rig(p),DataList.Datestr{p},DataList.Monkey{p});
%         elseif DataList.Drobo(p) == 2 & (datenum(dates{d},'yymmdd') < datenum('05/05/2018','mm/dd/yyyy'))
%             drname = sprintf('/volumes/drobo2/data/neurophys/rig%03u/%s_%s/',DataList.Rig(p),DataList.Datestr{p},DataList.Monkey{p});
%             
%         else
%             drname = sprintf('%s/rig%03u/%s_%s/',bigdrobopath,DataList.Rig(p),DataList.Datestr{p},DataList.Monkey{p});
%         end
%     end
%     
%      load([targetpath '/' DataList.Datestr{p} '_' DataList.Monkey{p} '/rez.mat']); 
%      
% 
%      
%      
%      ss  = KiloSort2SpikeStruct(rez,0,3);
%      ss.clusterMap(:,3)
%      saveoutpath  = '/volumes/drobo2/data/neurophys/doughek/kiloout/LGNgrouped_waves/';
%      depth = sprintf('_p%02u',sum(ismember(DataList.Datestr(1:p-1),DataList.Datestr(p))) + 1);
%      pfname       = [DataList.Datestr{p} '_' DataList.Monkey{p} depth]; 
%      save(strcat(saveoutpath,pfname,'_kilo','_presplit_ss','.mat'),'ss','-v7.3');
% end

