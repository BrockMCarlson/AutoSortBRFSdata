% setupAutoSortForBRFS
% BMC 8/3/22


%% Find Unique dates (to make sure you have all BRFS files ya dummy)
cd('T:\Brock - backups\Backup - WD harddrive - 220311\1 brfs ns6 files')
holder = dir;
clear fullName dateName
for i = 3:length(holder)
    fullName{i,1} = holder(i).name;
    dateName(i,:) = holder(i).name(1:8);
end
uniqueDates = unique(dateName,'rows');

%%
inputDataDirectory = 'E:\brfsNs6FilesToAutoSort';
outputPpnevDirectory = 'E:\ppnev_OutputsFromAutoSort';
cd(inputDataDirectory)
list = dir(inputDataDirectory);
fileList = {list.name}.';
fileList = fileList(3:end,1);

for i = 1:length(fileList)
    cd(inputDataDirectory)
    filename = fileList{i};
    multiplier = 1.5; % default is 2.2
    [ppNEV, WAVES] = offlineBRAutoSort(filename,multiplier);

    cd(outputPpnevDirectory)
    saveName = strcat(filename(1:end-4),'.ppnev');
    save(saveName,'ppNEV','WAVES')
end

cd(outputPpnevDirectory)
save('testPpnevOutput.mat','ppNEV','WAVES')