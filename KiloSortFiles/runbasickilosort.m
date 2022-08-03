targetdir = '/Volumes/Drobo2/DATA/NEUROPHYS/KiloSort-ed/';
list = dir(targetdir);
flag_checkforexisting = true;

erct = 0;
for i = 1:length(list)
    if list(i).isdir
        clear rezfile fullbinaryfile
        
        rezfile = sprintf('%s/%s/rez.mat',targetdir,list(i).name);
        if flag_checkforexisting && exist(rezfile,'file')
            continue
        end
         
        
        fullbinaryfile = sprintf('%s/%s/%s.ns6.dat',targetdir,list(i).name,list(i).name);
        if exist(fullbinaryfile,'file')
            fprintf('\n SORTING : %s\n',fullbinaryfile);
            try
            [~]=basicKiloSort(fullbinaryfile);
            catch err
                fprintf('\nERROR: %s\n',err.message)
                erct = erct +1
                ERR{erct,1} = fullbinaryfile;
                ERR{erct,2} = err;
            end
        else
            fprintf('\n NO FILE : %s\n',fullbinaryfile);
        end
    end
end
save([targetdir 'matlab.mat'],'ERR')
erct
ERR