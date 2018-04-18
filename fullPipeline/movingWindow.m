function P = movingWindow(gdf)

% make it skip seizure??
% This function takes a broad block of time before and after a seizure,
% detects spikes and then spike sequences, and then it goes through and
% uses a moving window to calculate the spatial organization of the spike
% sequences in each window, so that I can then plot the changing spatial
% organization over time

%% Parameters to change every time

% remove depths?
remove_electrodes = 0;
remove_type = {'D'};

% do cleaning step?
doClean = 0;
ss_thresh     = 15;                   
tt_thresh     = 0.015;  

% some distance between channels for channel weights. 
dmin = 15; 

% threshold for Marsh lab algorithm
tmul = 13;

% Output file name to save
outputName = 'HUP80_movingWindow_sz1_24hrs_noremove_noclean.mat';
%outputName = 'HUP78_oneMinBlocks.mat';
%outputName = 'HUP107_movingWindow_test.mat';

outputGDF = 'HUP080_sz1_gdf.mat';
%outputGDF = 'HUP107_sz1_gdf.mat';

% data name (for ieeg.org)
%dataName = 'HUP107_phaseII';
dataName = 'HUP80_phaseII';
%dataName = 'HUP78_phaseII-Annotations';  

% CSV file with electrode locations
%csvFile = 'HUP107_T1_19991219_electrode_labels.csv';
csvFile = 'HUP080_T1_19991213_electrode_labels.csv';
%csvFile = 'HUP078_T1_19971218_electrode_labels.csv';

% The patient name with format as used in the json file
%ptname = 'HUP107';
ptname = 'HUP080';
%ptname = 'HUP078';

% The number of the patient
%pt = 107;
pt = 80;
%pt = 78;

% How many seconds you want per block (the total chunk of time)
sPerBlock = 3600*24; %12 hours before the seizure and 12 hours after the seizure

% How much overlap for the moving window of recruitment latencies
window = 600; % 600 s = 10 minutes

% How much time to calculate recruitment latency in each window
timeRL = 3600; % 3600 s = 1 hour (less than this and the stats probably would not be very robust)

% I break the block into chunks to detect spikes, because calling iEEG.org
% with more than this causes Java permgen memory errors
chunkSize = 1000;

% parameters that the spike detector needs
vanleer = 0;
vtime = 0;
keepEKG = 0;
ignore = 1;
funnyname = 0;

%% Get paths and load seizure info and channel info
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
electrodeFile = [electrodeFolder,csvFile];
p1 = genpath(scriptFolder);
addpath(p1);
ptInfo = loadjson(jsonfile);

Patient(pt).seizures = ptInfo.PATIENTS.(ptname).Events.Ictal;
szNames = fieldnames(Patient(pt).seizures);

%% Load EEG data info
% calling this with 0 and 0 means I will just get basic info like sampling
% rate and channel labels
data = getiEEGData(dataName,0,0,pwfile);  
fs = data.fs;

%% Run the getSpikes script once as a dummy run just to produce a file of electrode locations
dummyRun = 1;
outputData = 0;
[~,electrodeData,~] = getSpikeTimes(0,ptname,dataName,electrodeFile,ptInfo,pwfile,...
    dummyRun,vanleer,vtime,outputData,keepEKG,ignore,funnyname,tmul);
Patient(pt).electrodeData = electrodeData;

%% Define seizure onset and offset times for each seizure
for i = 1:length(fieldnames(Patient(pt).seizures))
    Patient(pt).sz(i).onset = Patient(pt).seizures.(szNames{i}).SeizureEEC;
    Patient(pt).sz(i).offset = Patient(pt).seizures.(szNames{i}).SeizureEnd;
end

%% Define the start and stop times of each block prior to the seizure

% Loop through all the seizures
for i = 1:1%length(Patient(pt).sz)
    
    
    % Skip the seizure if it's too close to the start of the data
    if Patient(pt).sz(i).onset < sPerBlock/2+2
     %   continue
    end
    
    % The initial time of the first block for the seizure is the seizure
    % onset time minus the number of blocks x time per block/2    
    initialTime = max(Patient(pt).sz(i).onset-sPerBlock/2,1);
    
    
    % Define the run times
    Patient(pt).sz(i).runTimes(1:2) = [initialTime, initialTime+sPerBlock-1];
    Patient(pt).sz(i).szTimes = [Patient(pt).sz(i).onset Patient(pt).sz(i).offset];
    
    
end

%% Detect spikes

% Skip spike detection if I inputted a gdf
if nargin ==0
    
% Loop through all seizures
for i = 1:1%length(Patient(pt).sz)
    tic
    fprintf('Doing seizure %d of %d\n',i,length(Patient(pt).sz));
    
    % Skip if we're not running anything for the seizure
   if isempty(Patient(pt).sz(i).runTimes) == 0
       
   
       % Establish start and stop times
       desiredTimes = Patient(pt).sz(i).runTimes(1:2);

       % Break it up into chunks to avoid permgen memory errors
       nchunks = ceil(sPerBlock/chunkSize);
       gdf = [];
       for k = 1:nchunks 
           tic
           fprintf('Doing chunk %d of %d\n',k,nchunks);

           currTimes(1) = desiredTimes(1) + (k-1)*chunkSize;
           currTimes(2) = min(desiredTimes(1) + sPerBlock - 1, currTimes(1) + chunkSize-1);

           %% calculate gdf (spike times and locations) for the chunk in the block
           fprintf('Detecting spikes\n');
           dummyRun = 0;
           [gdft,~,~] = getSpikeTimes(currTimes,ptname,dataName,electrodeFile,ptInfo,pwfile,...
               dummyRun,vanleer,vtime,outputData,keepEKG,ignore,funnyname,tmul);

           if isempty(gdft) == 0
               % Adjust the times based on what chunk it is
               gdft(:,2) = gdft(:,2)+currTimes(1)-desiredTimes(1);

               % Add the spikes found in this chunk to all the spikes in
               % the block
               gdf = [gdf;gdft];
           end
           fprintf('Chunk %d took %1.1f\n',k,toc);
       end

        
       fprintf('It took %1.1f seconds to do seizure %d\n',toc,i);
   end
    
end

end



%% Save the gdf so that we can skip the above in the future
save([resultsFolder,outputGDF],'gdf');

if nargin == 0
    P = 0;
    
else
    
%% Remove depth electrodes
if remove_electrodes == 1
    
    gdf = removeChs(gdf,electrodeData,remove_type);
    
end

%% Get spike sequences for the block
for i = 1:length(Patient(pt).sz)
    
    if isempty(Patient(pt).sz(i).runTimes) == 1
        continue
    end
    
    Patient(pt).sz(i).stats.nspikes = size(gdf,1);
    Patient(pt).sz(i).stats.spikefreq = size(gdf,1)/sPerBlock;
    
    fprintf('Detecting sequences\n');
    Patient(pt).sz(i).data = mainSequences(gdf,electrodeData, fs);

    Patient(pt).sz(i).stats.nseqs = size(Patient(pt).sz(i).data.sequences,2)/2;
    Patient(pt).sz(i).stats.seqfreq = Patient(pt).sz(i).stats.nseqs/sPerBlock;
    
    % Do cleaning step
    if doClean == 1
        allSeq = Patient(pt).sz(i).data.sequences;
        iclean = spt_seqclust(Patient(pt).sz(i).data.xyChan,allSeq,ss_thresh,tt_thresh);
        
        % Get an array of the iclean indices and the subsequent ones going up by 2
        % per sequence to help index allseq
        y = zeros(length(iclean)*2, 1); 
        y(1:2:end-1)=(iclean-1)*2+1; 
        y(2:2:end)=(iclean-1)*2+2;
        
        cleanSeq = allSeq(:,y);
        Patient(pt).sz(i).data.oldSequences = allSeq;
        Patient(pt).sz(i).data.sequences = cleanSeq;
        Patient(pt).sz(i).data.iclean = iclean;
        
    end
    
end

%% Create an array of times to calculate recruitment latency and spatial org
for i = 1:length(Patient(pt).sz)
    times = zeros(sPerBlock/window + 1 - timeRL/window,2);
    times(1,:) = [0 timeRL];
    Patient(pt).sz(i).tblock(1).times = [0 timeRL];
    for j = 2:size(times,1)
        times(j,:) = [times(j-1,1)+window times(j-1,2) + window];
        Patient(pt).sz(i).tblock(j).times = times(j,:);
    end
   Patient(pt).sz(i).timesRL = times; 
end


%% Loop through the times and find the corresponding sequence indices
for i = 1:length(Patient(pt).sz)
    
    if isempty(Patient(pt).sz(i).data) == 1
        continue
    end
    
    % Get the times for the seizure
    times = Patient(pt).sz(i).timesRL;
    
    nseq = size(Patient(pt).sz(i).data.sequences,2)/2;
    
    % Loop through the times
    for j = 1:size(times,1)
        
        % initialize an empty array of sequence indices for that time
        Patient(pt).sz(i).tblock(j).sIdx = [];
        
        % Loop through the sequences
        for s = 1:nseq
            
            % get the correct column
            col = s*2;
            
            % Get the time of the first spike in the sequence
            seqtime = Patient(pt).sz(i).data.sequences(1,col);
            
            % if the spike falls between the times
            if seqtime >= times(j,1) && seqtime <= times(j,2)
                
                % add that sequence index to the array for that time block
                Patient(pt).sz(i).tblock(j).sIdx =...
                    [Patient(pt).sz(i).tblock(j).sIdx,s];
                
            % if the spike comes before the start of the time
            elseif seqtime < times(j,1)
                % move to the next sequence
                continue
                
            % if the spike comes after the end of the time
            elseif seqtime >times(j,2)
                
                % break out of the while loop, thus moving onto the next
                % time block
                break
            end
        end
    end
    
end

% Now i need to go through each time block, and for each block grab the
% appropriate sequences and calculate recruitment latencies and spatial
% orgs

%% Calculate recruitment latencies and spatial orgs for each block
for i = 1:length(Patient(pt).sz)
   if isempty(Patient(pt).sz(i).data) == 1
        continue
    end 
    
   sequences = Patient(pt).sz(i).data.sequences;
   for j = 1:length(Patient(pt).sz(i).tblock)
        sIdx = Patient(pt).sz(i).tblock(j).sIdx;
        
        newseq = [];
        
        % Get the columns from the sIdx
        for k = 1:length(sIdx)
            newseq = [newseq,sequences(:,k*2-1:k*2)];
        end
         
        [recruitmentLatencySingle,spikeCount] = ...
            getRecruitmentLatency(newseq,Patient(pt).sz(i).data.xyChan);
        
        [Patient(pt).sz(i).tblock(j).avgRecruitmentLat,...
           Patient(pt).sz(i).tblock(j).spatialOrg] = ...
           getSpatialOrg(recruitmentLatencySingle,Patient(pt).sz(i).data.xyChan,dmin);
       
       
   end 
end

P = Patient;
save([resultsFolder,outputName],'P');

end

end