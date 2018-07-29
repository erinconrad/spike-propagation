function [ignoreChLabels] = ...
    findChsToIgnore(pt,whichPt,chNames)


[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;


%% Load electrode file
electrodeFile = pt(whichPt).electrode_labels;

if isempty(electrodeFile) ==1
    error('Warning, cannot find electrode file for pt %s\n',pt(whichPt).name);
end

fileID = fopen([electrodeFolder,electrodeFile]);
out=textscan(fileID, '%s', 'whitespace',',');
out = out{1};
nChans = length(out)/8;
fclose(fileID);

%% Get the names of the channels in the file
fileNames = {};
for i = 1:nChans
    tempname = out{(i-1)*8+5};
    fileNames = [fileNames,tempname];
    
end

ignoreChLabels = {};
for i = 1:length(chNames)
    [Lia,~] = ismember(chNames{i},fileNames);
    if Lia == 0
       ignoreChLabels = [ignoreChLabels,chNames{i}]; 
    end
    
end


end





