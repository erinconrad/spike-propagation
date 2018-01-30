%{

This is the main code to run the pipeline to 1) detect spikes, 2) detect
spike sequences, and 3) calculate the spatial organization of these spike
sequences at various time points prior to seizures.

It first finds the times of seizures and breaks the time intervals before
the seizures into chunks. Then for each chunk, it goes through and detects
spikes (calling getSpikeTimes). Then for each chunk, it looks through those
spikes and detects spike sequences and calculates the spatial organization
(calling mainSequences).

At some point, this could be expanded to loop through multiple patients.

This requires that I have a file with electrode data for the patient
already made

%}


function Patient = main
tic
%% Parameters

% data name (for ieeg.org)
dataName = 'HUP78_phaseII-Annotations';   

% The patient name with format as used in the json file
ptname = 'HUP078';

% The number of the patient
pt = 78;

% The sampling rate of the EEG data
fs = 512;

% How many seconds you want per block. Max allowable appears to be 2000, or
% possibly less
sPerBlock = 500;

% How many blocks you want to compare before the seizure
nblocks = 10;


%% Get paths and load seizure info and channel info
[electrodeFile,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
p1 = genpath(scriptFolder);
addpath(p1);
ptInfo = loadjson(jsonfile);


%load([ptname,'_chanLoc.mat'])
Patient(pt).seizures = ptInfo.PATIENTS.(ptname).Events.Ictal;
szNames = fieldnames(Patient(pt).seizures);

%% Run the getSpikes script once as a dummy run just to produce a file of electrode locations
[~,electrodeData] = getSpikeTimes(0,dataName,electrodeFile,ptInfo,pwfile,1);

%% Define seizure onset and offset times for each seizure
for i = 1:length(fieldnames(Patient(pt).seizures))
    Patient(pt).sz(i).onset = Patient(pt).seizures.(szNames{i}).SeizureEEC;
    Patient(pt).sz(i).offset = Patient(pt).seizures.(szNames{i}).SeizureEnd;
end

%% Define the start and stop times of each block prior to the seizure

% Loop through all the seizures
for i = 1:length(Patient(pt).sz)
    
    
    
    % Skip the seizure if it's too close to the start of the data
    if Patient(pt).sz(i).onset < nblocks*sPerBlock+2
        continue
    end

    % Skip the seizure if it's too close to the prior seizure
    if i~=1
        if Patient(pt).sz(i).onset - Patient(pt).sz(i-1).offset < nblocks*sPerBlock+1
           continue 
        end
    end
    
    % The initial time of the first block for the seizure is the seizure
    % onset time minus the number of blocks x time per block
    initialTime = Patient(pt).sz(i).onset-nblocks*sPerBlock-1;
    
    % Loop through the blocks
    for j = 1:nblocks
        
        % Advance the time to the next block start
        newStartTime = initialTime + (j-1)*sPerBlock+1;
        
        % Define the run times
        Patient(pt).sz(i).runTimes(j,1:2) = [newStartTime, newStartTime+sPerBlock-1];
    end
    
end

%% Do the full analysis on each block

% Loop through all seizures
for i = 1:length(Patient(pt).sz)
    
    fprintf('Doing seizure %d of %d\n',i,length(Patient(pt).sz));
    
    % Skip if we're not running anything for the seizure
   if isempty(Patient(pt).sz(i).runTimes) == 0
       
       % Loop through all blocks
       for j = 1:length(Patient(pt).sz(i).runTimes)
           fprintf('Doing block %d of %d in seizure %d of %d\n',...
               j,length(Patient(pt).sz(i).runTimes),i,length(Patient(pt).sz));
           
          mm = java.lang.Runtime.getRuntime.maxMemory;
          tm = java.lang.Runtime.getRuntime.totalMemory;
          fm = java.lang.Runtime.getRuntime.freeMemory;
          
          fprintf('Max memory =  %1.2e, total memory = %1.2e, free memory = %1.2e\n',...
              mm, tm, fm);
           
           % Establish start and stop times
           desiredTimes = Patient(pt).sz(i).runTimes(j,1:2);
           
           % calculate gdf (spike times and locations) for the block
           fprintf('Detecting spikes\n');
           [gdf,~] = getSpikeTimes(desiredTimes,dataName,electrodeFile,ptInfo,pwfile,0);
           
           % Get spike sequences and spatial organization for the block
           fprintf('Detecting sequences and calculating spatial organization\n');
           Patient(pt).sz(i).block(j).data = mainSequences(gdf,electrodeData);
       end
   end
    
end

%% Put the spatial organizations together in a single array for each seizure
for i = 1:length(Patient(pt).sz)
   if isempty(Patient(pt).sz(i).runTimes) == 0
       for j = 1:length(Patient(pt).sz(i).runTimes)
           Patient(pt).sz(i).spatialOrg(j) = ...
               Patient(pt).sz(i).block(j).data.spatialOrg;
       end
   end
    
end

save([resultsFolder,'Patient.mat'],'Patient');
toc
end

