%% getDesiredRunTimes

% This file, meant to be ported to Borel, loops through all the patients in
% the json file containing basic patient info, and fills up a structure
% with desired times to calculate spikes, as well as output file locations
% for the future spike files

clear

%% File names
outputFile = 'long_times.mat';
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
%electrodeFile = [electrodeFolder,csvFile];
p1 = genpath(scriptFolder);
addpath(p1);
ptInfo = loadjson(jsonfile);
ptnames = fieldnames(ptInfo.PATIENTS);

% Do trickery
ptnames(strcmp(ptnames,'HUP111B')==1) = [];
ptnames(strcmp(ptnames,'Study004')==1) = [];

%load([resultsFolder,'ptStructs/',outputFile])

%% Desired times
% Look 12 hours before the sz and 12 hours after
chunkTime = 2000; % Save it in 2000 s chunks (ieeg crashes if request more than this)
totalTime = 3600*24;

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
       
       % If the seizure number is 1000, skip seizure
       if sz.SeizureEEC == -1
           skipSz = 1;
       end
       
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


%% Get run times for each patient
for i = 1:length(ptnames)
    
    if isempty(pt(i).sz) ==1, continue; end
    
    % If the time between seizures is less than 100 hours, run the whole time
    if pt(i).sz(end).offset -  pt(i).sz(1).onset < 3600*200
        pt(i).allTimes = [max(0,pt(i).sz(1).onset - 3600*12),...
            pt(i).sz(end).offset + 3600*12];
    
    % If it's longer than this, then break into chunks
    else
        allTimes = [];
    
        % Loop through seizures
        for j = 1:length(pt(i).sz)

           sztimes =  [pt(i).sz(j).onset,pt(i).sz(j).offset];

           % Start detecting spikes 12 hours before the seizure onset
           initialTime = max(pt(i).sz(j).onset - totalTime/2,1);

           % If it is not the first seizure, can start detecting just after the
           % run of the prior seizure (removes redundant data collection)
           if j > 1
               initialTime = max(pt(i).sz(j).onset - totalTime/2,...
                   allTimes(j-1,2)); 
           end

           % This will break if it's too close to the end of the file
           finalTime = pt(i).sz(j).onset + totalTime/2;


           allTimes = [allTimes;initialTime finalTime];
         end
         pt(i).allTimes = allTimes;
         
    end
end

for i = 1:length(pt)
    
   % Get sz times 
   pt(i).newSzTimes = [];
   for j = 1:length(pt(i).sz)
       pt(i).newSzTimes(j,:) = [pt(i).sz(j).onset,pt(i).sz(j).offset];
   end 
    
   if isempty(pt(i).allTimes) == 1, continue; end
   
   pt(i).runTimes = [];
   pt(i).chunkFiles = {};
   
   for j = 1:size(pt(i).allTimes,1)
       initialTime = pt(i).allTimes(j,1);
       finalTime = pt(i).allTimes(j,2);
       totalTime = finalTime - initialTime;
       
       nchunks = ceil(totalTime/chunkTime);
       
       for k = 1:nchunks
           startTime = initialTime + (k-1)*chunkTime;
           endTime = min(finalTime, startTime + chunkTime);
           pt(i).runTimes = [pt(i).runTimes;startTime, endTime];
           
           % define the output file
           pt(i).chunkFiles = [pt(i).chunkFiles;...
               [pt(i).name,'_subset_',sprintf('%d',j),'_times_',sprintf('%d',startTime),...
               '-',sprintf('%d',endTime),'.mat']];
           
           
       end
       
   end
    
    
    
end
       
       
 

save([resultsFolder,'ptStructs/',outputFile],'pt');



