clear

%% To do
% Add the ieeg.org file name

%% File names
outputFile = 'desiredTimes.mat';
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
%electrodeFile = [electrodeFolder,csvFile];
p1 = genpath(scriptFolder);
addpath(p1);
ptInfo = loadjson(jsonfile);
ptnames = fieldnames(ptInfo.PATIENTS);
gdfFolder = 'gdf/';

%% Desired times
totalTime = 3600*24; % Look 12 hours before the sz and 12 hours after
chunkTime = 3600; % Save it in one hour chunks
nchunks = ceil(totalTime/chunkTime);
window = 3600;
overlap = 600;

%% Loop through the patients in the json file
for i = 1:length(ptnames)
    info = ptInfo.PATIENTS.(ptnames{i});
    pt(i).name = ptnames{i};
    pt(i).ignore_electrodes = info.IGNORE_ELECTRODES;
    pt(i).ieeg_name = ieegNames(pt(i).name);
    pt(i).chLocationFile = [resultsFolder,'chLocations/',pt(i).name,'_chLocations.mat'];
    
    % Get electrode label file names
    origname = info.ELECTRODE_LABELS;
    temp = strsplit(origname,'/');
    fname = temp{end};
    pt(i).electrode_labels = [electrodeFolder,fname];
    
    % Get seizures
    szs = fieldnames(info.Events.Ictal);
    
    for j = 1:length(szs)
       sz = info.Events.Ictal.(szs{j});
       
       if isfield(sz,'SeizureEEC') == 0
           continue
       end
       
       pt(i).sz(j).onset = sz.SeizureEEC;
       pt(i).sz(j).offset = sz.SeizureEnd;
       
       initialTime = max(pt(i).sz(j).onset - totalTime/2,1);
       
       % This will break if it's too close to the end of the file
       finalTime = pt(i).sz(j).onset + totalTime/2;
       
       % Initialize run times and file names
       pt(i).sz(j).runTimes = zeros(nchunks,2);
       pt(i).sz(j).chunkFiles = cell(nchunks,1);
       
       % Create the times
       for k = 1:nchunks
           startTime = initialTime + (k-1)*chunkTime;
           endTime = min(finalTime, startTime + chunkTime);
           pt(i).sz(j).runTimes(k,:) = [startTime, endTime];
           pt(i).sz(j).chunkFiles{k} = ...
               [resultsFolder,gdfFolder,pt(i).name,'_sz_',sprintf('%d',j),'_times_',sprintf('%d',startTime),...
               '-',sprintf('%d',endTime),'.mat'];
       end
       
       % Also create the times for doing the moving window to calculate RL,
       % which will be done much later
       timesRL = zeros(totalTime/overlap + 1 - window/overlap,2);
       timesRL(1,:) = [0 window];
       pt(i).sz(j).blockRL(1).times = [0 window];
       for t = 2:size(timesRL,1)
           timesRL(t,:) = [timesRL(t-1,1)+overlap timesRL(t-1,2) + overlap];
           pt(i).sz(j).blockRL(t).times = timesRL(t,:);
       end
       pt(i).sz(j).timesRL = timesRL;

    end
  
end


save([resultsFolder,outputFile],'pt');



