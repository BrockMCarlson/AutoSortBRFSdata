function Directories = AutoSortSetup(user)

% helper function 
% MAC, Feb 2020
global Directories
 

switch user
    
    case {'BrockWork'}
        Directories.NS6DIR   = 'T:\Brock - backups\Backup - WD harddrive - 220311\1 brfs ns6 files\';
        Directories.RIGDIR   = 'T:\Brock - backups\Backup - WD harddrive - 220311\all BRFS\';
        Directories.CODEDIR  = 'C:\Users\Brock\Documents\MATLAB\GitHub\GitHub\AutoSortBRFSdata';
        Directories.IDXDIR   = 'T:\Brock - backups\Backup - WD harddrive - 220311\5 diIDX dir\';
        Directories.STIMDIR  = 'T:\diSTIM - adaptdcos&CRF\STIM\';
        Directories.OUTDIR   = 'C:\Users\Brock\Box\JOV brfs submission\formattedDataOutputs\methodsFigOutput\';
        cd(Directories.CODEDIR)

    case {'BrockHome'}
        Directories.NS6DIR   = 'T:\Brock - backups\Backup - WD harddrive - 220311\1 brfs ns6 files\';
        Directories.RIGDIR   = 'D:\all BRFS\';
        Directories.CODEDIR  = 'C:\Users\Brock Carlson\Documents\GitHub\AutoSortBRFSdata';
        Directories.IDXDIR   = 'D:\5 diIDX dir';
        Directories.STIMDIR  = 'T:\diSTIM - adaptdcos&CRF\STIM\';
        Directories.OUTDIR   = 'C:\Users\Brock Carlson\Box\JOV brfs submission\formattedDataOutputs\methodsFigOutput\';
        cd(Directories.CODEDIR)
        
end



           
        
end
