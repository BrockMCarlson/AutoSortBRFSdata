% ConvertAutoSortOutputsToNevStructure
% BMC 8/8/2022

clear


global CODEDIR
CODEDIR = 'C:\Users\Brock\Documents\MATLAB\GitHub\AutoSortBRFSdata';
cd(CODEDIR)

autoSortOutput = 'E:\ppnev_OutputsFromAutoSort';
outputDirForFormatedPpnevFiles = 'E:\ppnev_FinalFileOutputs';
cd(autoSortOutput)
list = dir(autoSortOutput);
fileList = {list.name}.';
fileList = fileList(3:end,1);



for i = 1:length(fileList)

    cd(autoSortOutput)
    filename = fileList{i};
    disp(filename)
    S = load(filename);
    N = fieldnames(S);

    idx_spk = strncmp(N,'spk',3);
    TimeStamp = []; 
    varNumber_spk = find(idx_spk);
    for k = 1:length(varNumber_spk)
        spkFieldName = N{varNumber_spk(k)};
        spk = S.(spkFieldName);
        TimeStamp = [TimeStamp spk]; 
    end

    idx_elec = strncmp(N,'elec',4);
    Electrode = []; 
    varNumber_elec = find(idx_elec);
    for l = 1:length(varNumber_elec)
        elecFieldName = N{varNumber_elec(l)};
        elec = S.(elecFieldName);
        Electrode = [Electrode elec]; 
    end

    Unit = ones(size(TimeStamp));
    [TimeStamp,idx] = sort(TimeStamp);
    Electrode = Electrode(idx);
    
    % load NEV
    clear NEV ppNEV
    tebaDir = 'T:\Brock - backups\Backup - WD harddrive - 220311\all BRFS';
    nevFileName = strcat(tebaDir,filesep,filename(1:8),filesep,filename(1:end-4),'.nev');
    NEV = openNEV(nevFileName,'noread','nomat','nosave');
    
    % organize ppNEV structure for outpu
    ppNEV.MetaTags             = NEV.MetaTags; 
    ppNEV.ElectrodesInfo       = NEV.ElectrodesInfo; 
    ppNEV.Data.SerialDigitalIO = NEV.Data.SerialDigitalIO; 
    clear NEV
    
    ppNEV.Data.Spikes.TimeStamp = TimeStamp; 
    ppNEV.Data.Spikes.Electrode = Electrode; 
    ppNEV.Data.Spikes.Unit = Unit; 
    
    ppNEV.filename = filename; 
    ppNEV.multiplier = 1.5; %BMC 8/9/2022 run multiplier


    ppNEV.timestamp = now; 

    cd(outputDirForFormatedPpnevFiles)
    fileSaveName = strcat(filename(1:end-4),'.ppnev');
    save(fileSaveName,'ppNEV','-v7.3')


end



