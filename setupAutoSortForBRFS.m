% setupAutoSortForBRFS
% BMC 8/3/22
clear

%% Setup
global CODEDIR
CODEDIR = 'C:\Users\Brock Carlson\Documents\GitHub\AutoSortBRFSdata';
cd(CODEDIR)

inputDataDirectory = 'D:\all BRFS';
outputPpnevDirectory = 'E:\AutoSortOutput - HOME';

%% Find Unique dates 
clear uniqueDates

% 


% allBRFS
cd('D:\all BRFS')
holder = dir;
clear fullName dateName
for i = 3:length(holder)
    fullName{i,1} = holder(i).name;
    dateName(i,:) = holder(i).name(1:8);
end
uniqueDates.allBRFS = unique(dateName,'rows');

% diSTIM - adaptdcos&CRF\AutoSort-ed
cd('T:\diSTIM - adaptdcos&CRF\AutoSort-ed')
holder = dir;
clear fullName dateName
for i = 3:length(holder)
    fullName{i,1} = holder(i).name;
    dateName(i,:) = holder(i).name(1:8);
end
uniqueDates.currentFiles = unique(dateName,'rows');

test = ismember(uniqueDates.allBRFS,uniqueDates.currentFiles,'rows');
missingDates = uniqueDates.allBRFS()

%% Make list of all rfori and brfs files
% What info do I want? -- 
% I would like a list of all ns6 files to autosort. I do NOT need the same
% list of nev files, because I can pull that out wiht: 
% % [pathstr,name,ext] = fileparts(filename); 
% % filename = [pathstr filesep name]; 
% % NS6_header = openNSx([filename '.ns6'],'noread');
% % NEV = openNEV(strcat(filename,'.nev'),'noread','nomat','nosave');



cd(inputDataDirectory)
list1 = dir(inputDataDirectory);
count = 0;

for i = 3:length(list1)
    folderName = strcat(list1(i).folder,filesep,list1(i).name);
    cd(folderName)
    %BRFS
    list2_brfs  = dir('*brfs*.ns6');
    for j = 1:length(list2_brfs)
        count = count + 1;
        fullFileList{count,1} = strcat(list2_brfs(j).folder, filesep, list2_brfs(j).name);
    end

    %rfori
    list2_rfori = dir('*rfori*.ns6');
    for k = 1:length(list2_rfori)
        count = count + 1;
        fullFileList{count,1} = strcat(list2_rfori(k).folder, filesep, list2_rfori(k).name);
    end
end



%% run list
cd(inputDataDirectory)
for i = 1:length(fullFileList)
    filename = fullFileList{i};
    disp(strcat('fileNum_',string(i),'_of_',string(length(fullFileList)),'--',filename))
    multiplier = 2.2; % default is 2.2
    [ppNEV] = offlineBRAutoSort(filename,multiplier,outputPpnevDirectory); 
    %save
    cd(outputPpnevDirectory)
    [pathstr,name,ext] = fileparts(filename); 
    saveName = strcat(name,'.ppnev');
    save(saveName,'ppNEV')
end




