%% getDesiredRunTimes

% This file, meant to be ported to Borel, loops through all the patients in
% the json file containing basic patient info, and fills up a structure
% with desired times to calculate spikes, as well as output file locations
% for the future spike files

clear

%% File names
outputFile = 'desiredTimes.mat';
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
%electrodeFile = [electrodeFolder,csvFile];
p1 = genpath(scriptFolder);
addpath(p1);
ptInfo = loadjson(jsonfile);
ptnames = fieldnames(ptInfo.PATIENTS);

%% Desired times
totalTime = 3600*24; % Look 12 hours before the sz and 12 hours after
chunkTime = 2000; % Save it in 2000 s chunks (ieeg crashes if request more than this)
nchunks = ceil(totalTime/chunkTime);
window = 3600; % For spatial organization calculation, calculate SO over a one hour window
overlap = 600; % Allow 10 minutes of overlap between the SO windows


%% Loop through the patients in the json file
for i = 1:length(ptnames)
    info = ptInfo.PATIENTS.(ptnames{i});
    
    % Get basic info
    pt(i).name = ptnames{i};
    %pt(i).ignore_electrodes = info.IGNORE_ELECTRODES;
    [pt(i).ieeg_name,electrodeFile,pt(i).thresh,pt(i).dmin,pt(i).elecNotes] = ieegAndElectodeNames(pt(i).name);
    pt(i).chLocationFile = [pt(i).name,'_chLocations.mat'];
    pt(i).sz_onset = info.SeizureOnset;

    % Get electrode label file names
    pt(i).electrode_labels = electrodeFile;
    
    % Get seizures
    szs = fieldnames(info.Events.Ictal);
    
    whichSz = 0;
    
    for j = 1:length(szs)
       sz = info.Events.Ictal.(szs{j});
       
       if isfield(sz,'SeizureEEC') == 0
           continue
       end
       
       skipSz = 0;
       
       % If the seizure is the same time as a prior, skip this one
       if j > 1
           for k = 1:whichSz-1
               if abs(pt(i).sz(k).onset - sz.SeizureEEC) < 1
                  skipSz = 1;
                   
               end
               
           end
           
       end
       
       if skipSz == 1
           continue
       end
       
       whichSz = whichSz + 1;
       
       % Get seizure onset and offset
       pt(i).sz(whichSz).onset = sz.SeizureEEC;
       pt(i).sz(whichSz).offset = sz.SeizureEnd;
       pt(i).sz(whichSz).electrodes = sz.SEIZURE_ONSET_ELECTRODES;
       
       % If the seizure onset is before the prior seizure onset, switch
       % positions
       if j > 1
           if pt(i).sz(whichSz).onset < pt(i).sz(whichSz-1).onset
               firstSz = [pt(i).sz(whichSz).onset pt(i).sz(whichSz).offset];
               secondSz = [pt(i).sz(whichSz-1).onset pt(i).sz(whichSz-1).offset];
               pt(i).sz(whichSz-1).onset = firstSz(1);
               pt(i).sz(whichSz-1).offset = firstSz(2);
               
               pt(i).sz(whichSz).onset = secondSz(1);
               pt(i).sz(whichSz).offset = secondSz(2);
               
           end
           
       end
       
    end
    
end

for i = 1:length(ptnames)
    
    % Loop through seizures
    for j = 1:length(pt(i).sz)
       
       sztimes =  [pt(i).sz(j).onset,pt(i).sz(j).offset];
       
       % Start detecting spikes 12 hours before the seizure onset
       initialTime = max(pt(i).sz(j).onset - totalTime/2,1);
       
       % If it is not the first seizure, can start detecting just after the
       % run of the prior seizure (removes redundant data collection)
       if j > 1
           initialTime = max(pt(i).sz(j).onset - totalTime/2,...
               pt(i).sz(j-1).runTimes(end,2)); 
       end
       
       % This will break if it's too close to the end of the file
       finalTime = pt(i).sz(j).onset + totalTime/2;
       
       % If it's not the first seizure, nchunks might be smaller
       nchunks = ceil((finalTime - initialTime)/chunkTime);
       
       % Initialize run times and file names
       pt(i).sz(j).runTimes = zeros(nchunks,2);
       pt(i).sz(j).chunkFiles = cell(nchunks,1);
       pt(i).sz(j).EKGchunkFiles = cell(nchunks,1);

       
       % Create the times (in 2000 second chunks) over which we will detect
       % spikes
       for k = 1:nchunks
           startTime = initialTime + (k-1)*chunkTime;
           endTime = min(finalTime, startTime + chunkTime);
           pt(i).sz(j).runTimes(k,:) = [startTime, endTime];
           
           % define the output file
           pt(i).sz(j).chunkFiles{k} = ...
               [pt(i).name,'_sz_',sprintf('%d',j),'_times_',sprintf('%d',startTime),...
               '-',sprintf('%d',endTime),'.mat'];
           
           pt(i).sz(j).EKGchunkFiles{k} = ...
               [pt(i).name,'_sz_',sprintf('%d',j),'_times_',sprintf('%d',startTime),...
               '-',sprintf('%d',endTime),'_ekg_','.mat'];
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
       
       if size(pt(i).sz(j).runTimes,1) == 0
           continue;
       end
           
       % Add additional things to say if it's during a seizure
       for k = 1:size(timesRL,1)
           blocktimes = pt(i).sz(j).blockRL(k).times + pt(i).sz(j).runTimes(1,1);
           pt(i).sz(j).blockRL(k).adjustedTimes = blocktimes;
           
           % if the start of the block is before the end of the seizure,
           % and the end of the block is after the start of the seizure,
           % then the block includes some seizure time
           if blocktimes(1) < sztimes(2) && blocktimes(2) > sztimes(1)
               pt(i).sz(j).blockRL(k).ictal = 1;
           else
               pt(i).sz(j).blockRL(k).ictal = 0;
           end
           
           
       end

    end
    
    % Add a date stamp
    pt(i).timestamp = datetime('now');
  
end


save([resultsFolder,'ptStructs/',outputFile],'pt');



