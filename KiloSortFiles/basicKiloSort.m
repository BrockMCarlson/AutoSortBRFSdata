function rez = basicKiloSort(datafile)
% MAC, Dec 2016


% extract root sorting path
fpath = fileparts(datafile);

% direct to KiloSort Toolbox
addpath('/Users/kacie/documents/GitHub/npy-matlab');
addpath(genpath('/Users/kacie/documents/GitHub/KiloSort'));
addpath(genpath('/Users/kacie/documents/GitHub/KiloSortUtils'));

% set confguration
load(fullfile(fpath, 'chanMap.mat'),'nprobes')
run('/Users/kacie/documents/GitHub/KiloSortUtils/config.m');

% preprocess 
[rez, DATA, uproj] = preprocessData(ops); % preprocess data and extract spikes for initialization

if size(DATA,2) == 1
    top_pad = randn([size(DATA,1) 15 size(DATA,3)]);
    bot_pad = randn([size(DATA,1) 15 size(DATA,3)]);
    DATA = cat(2,top_pad,DATA,bot_pad); 
end

% sort!
rez                = fitTemplates(rez, DATA, uproj);  % fit templates iteratively
rez                = fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)

% Auto Merge
rez = merge_posthoc2(rez);

% save python results file for Phy
rezToPhy(rez, fpath);

% save rez var with data
save(fullfile(fpath,  'rez.mat'), 'rez', '-v7.3');

% remove temporary file
delete(ops.fproc);

% Close all FIDS
fclose('all');