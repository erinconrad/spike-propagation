function electrodeData = chanLocUseGdf(unignoredChLabels,electrodeFile)

%{
This file makes .mat structures with electrode data. Importantly, it uses
the gdf file to separate the channels that I ignored in the spike detection
script from those that I did not ignore, and then makes channels that are
aligned with the gdf data.

This REQUIRES that you have already calculated spike times and made a gdf
file. Use getSpikeTimes.m to make gdf files.

%}




%% which patient
if nargin  == 0
    mainFolder = '/Users/erinconrad/Desktop/residency stuff/R25/actual work/';
    fname = 'HUP078_T1_19971218_electrode_labels.csv';
    electrodeFile = [mainFolder,'data/','electrodeData/',fname];
    % get patient name
    C = strsplit(fname,'_');
    electrodeData.ptName = C{1};
    ptname = electrodeData.ptName;

    spikePaths
end



%% Other files


if nargin == 0
%% Load gdf file
    load([gdfFolder,electrodeData.ptName,'_gdf.mat'])
end

%% Load electrode file

fileID = fopen(electrodeFile);

out=textscan(fileID, '%s', 'whitespace',',');
out = out{1};
nChans = length(out)/8;

%% make electrode data struct

% this will only have the unignored channels
electrodeData.electrodes(length(unignoredChLabels),1) = struct();

% this will have labels of the channel I am ignoring
electrodeData.ignoreChs = {};

for i = 1:nChans % this loops through the electrode file channels
    
    %% Decide whether to include the channel
    
    % Get label of the current channel
    tempname = out{(i-1)*8+5};
    ignoreCh = 0;
    
    % loop through all the channels I do NOT want to ignore
    for j = 1:length(unignoredChLabels)
        
       % if the current channel matches one of the designated unignore channels
       if strcmp(tempname,unignoredChLabels{j}) == 1
      
           % fill up the electrode data in the struct such that the index
           % is the same as the index of the channel in the gdf file
           electrodeData.electrodes(j).x = str2double(out{(i-1)*8+2});
           electrodeData.electrodes(j).y = str2double(out{(i-1)*8+3});
           electrodeData.electrodes(j).z = str2double(out{(i-1)*8+4});
           electrodeData.electrodes(j).xyz = [electrodeData.electrodes(j).x,...
           electrodeData.electrodes(j).y,electrodeData.electrodes(j).z];
           electrodeData.electrodes(j).name = out{(i-1)*8+5};
           electrodeData.electrodes(j).type = out{(i-1)*8+6};
           ignoreCh = 0;
           break
       else
           ignoreCh = 1;
       end
        
    end
    
    if ignoreCh == 1
        % add it to the list of channels to ignore
        electrodeData.ignoreChs = [electrodeData.ignoreChs,out{(i-1)*8+5}];
    end
    
 
    
end

electrodeData.unignoredChs = cell(length(electrodeData.electrodes),1);

for i = 1:length(electrodeData.electrodes)
    electrodeData.unignoredChs{i} = electrodeData.electrodes(i).name;
    
end

% make a dump of the locations of the correctly indexed channels
electrodeData.locs = zeros(length(electrodeData.electrodes),4);
for i = 1:length(electrodeData.electrodes)
    electrodeData.locs(i,:) = [i,electrodeData.electrodes(i).xyz];
    
end

%% Save struct
if nargin == 0
    save([electrodeFolder,electrodeData.ptName,'_chanLoc.mat'], 'electrodeData')
end

end
