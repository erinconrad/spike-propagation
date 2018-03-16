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

This requires that you have the following data:
- the EEG data (on ieeg.org); you need the name of the file you want
- the csv file with electrode locations for the patient
- the json file with 1) the seizure times for the patient, and 2) identity
of which electrodes to ignore for the patient

%}


function Patient = main
tic
%% Parameters to change every time
% Cleaning parameters

% spatial and temporaly thresholds used for cleaning algorithm. It finds all
% points that fall within ss_thresh of the reference point and tt_thresh of
% the reference point (similar time to lead spike). Then, the matched point
% with the shortest Euclidean distance from the reference point was used to
% compute the similarity score 1-d_match/ss_thresh. In the paper ss_thresh
% was 1.5 cm (although they used 2 D measurements) and tt_thresh was 15 ms.
% Sam's code uses 1.5 for ss_thresh and 3 for tt_thresh. However, I assume
% the 3 for tt_thresh is because they do everything in indices, not
% seconds. I am not sure if the 1:1 conversion here between cm and points
% is correct. It looks like my locations of points are in mm because
% distance between 2 electrodes in my coordinates is about 10 units, and
% the interelectrode distance should be 10 mm. So 1 unit = 1 mm.
ss_thresh     = 15;    % PROBABLY WRONG                
tt_thresh     = 0.015;  

% some distance between channels for channel weights. In the paper this was
% 1.5 cm so I need to figure out the conversion. I think for me 1 unit is 1 mm. MUST CHECK.
dmin = 15; 

% Output file name to save
outputName = 'HUP80_allsz_5blocks_icandpreic_1000s_cleanall2.mat';
%outputName = 'HUP78_oneMinBlocks.mat';

% data name (for ieeg.org)

dataName = 'HUP80_phaseII';
%dataName = 'HUP78_phaseII-Annotations';  

% CSV file with electrode locations
csvFile = 'HUP080_T1_19991213_electrode_labels.csv';
%csvFile = 'HUP078_T1_19971218_electrode_labels.csv';

% The patient name with format as used in the json file
ptname = 'HUP080';
%ptname = 'HUP078';

% The number of the patient
pt = 80;
%pt = 78;

% How many seconds you want per block. Max allowable appears to be 2000, or
% possibly less
sPerBlock = 500;

% Include ictal period?
includeIc = 0;

% How many blocks you want to compare before the seizure
nblocks = 1;

% Remove EKG artifact?
rmEKGArtifact = 0;
prox = 0.01; %10 ms
indexToms = 0;

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
[~,electrodeData,~] = getSpikeTimes(0,dataName,electrodeFile,ptInfo,pwfile,1,0,0,0,0,1,0);

%% Define seizure onset and offset times for each seizure
for i = 1:length(fieldnames(Patient(pt).seizures))
    Patient(pt).sz(i).onset = Patient(pt).seizures.(szNames{i}).SeizureEEC;
    Patient(pt).sz(i).offset = Patient(pt).seizures.(szNames{i}).SeizureEnd;
end

%% Define the start and stop times of each block prior to the seizure

% Loop through all the seizures
for i = 1:1%length(Patient(pt).sz)
    
    
    
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
    % onset time minus the number of blocks x time per block (nblocks -1
    % makes it include the ictal period)
    if includeIc == 1
        initialTime = Patient(pt).sz(i).onset-(nblocks-1)*sPerBlock-1;
    elseif includeIc == 0
        initialTime = Patient(pt).sz(i).onset-(nblocks)*sPerBlock-1;
    end
    
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
       for j = 1:size(Patient(pt).sz(i).runTimes,1)
           tic
           fprintf('Doing block %d of %d in seizure %d of %d\n',...
               j,length(Patient(pt).sz(i).runTimes),i,length(Patient(pt).sz));
           
          
           % Establish start and stop times
           desiredTimes = Patient(pt).sz(i).runTimes(j,1:2);
           
           %% calculate gdf (spike times and locations) for the block
           fprintf('Detecting spikes\n');
           [gdf,~,~] = getSpikeTimes(desiredTimes,dataName,electrodeFile,ptInfo,pwfile,0,0,0,0,0,1,0);
           size(gdf)
           
           %% EKG artifact removal
           if rmEKGArtifact == 1
               % calculate gdf and values of EKG channels
               [gdfEKG,~,~] = getSpikeTimes(desiredTimes,dataName,electrodeFile,...
                   ptInfo,pwfile,0,0,0,0,1);

               % remove spikes that occur too close to EKG channel spikes
                gdf = removeEKGArtifact(gdf,gdfEKG,prox);
           end
           
           %% Get spike sequences and spatial organization for the block
           fprintf('Detecting sequences and calculating spatial organization\n');
           Patient(pt).sz(i).block(j).data = mainSequences(gdf,electrodeData, fs);
           toc
       end
   end
    
end


%% concatenate all sequences for the purpose of cleaning

% The array with all sequences
allseq = [];

% Nx3 array where each row is a sequence, the first column is the seizure
% number, the 2nd is the block number, and the 3rd is the number of the
% sequence in that seizure and block
block_seiz_ids = [];

for i = 1:length(Patient(pt).sz)
   for j = 1:length(Patient(pt).sz(i).block)
      
       % initialize cleanseq as empty
      Patient(pt).sz(i).block(j).data.cleanseq = [];
      
      % get the sequences for the current seizure and block
      currseq =  Patient(pt).sz(i).block(j).data.sequences;
      
      % pad the shorter sequence set with zeros
      diff = size(currseq,1) - size(allseq,1);
      if diff < 0
          currseq = vertcat(currseq,zeros(abs(diff),size(currseq,2)));
      elseif diff > 0
          allseq = vertcat(allseq, zeros(diff,size(allseq,2)));
      end         
          
     
      % concatenate
      allseq = [allseq, currseq];
      
      % fill up matrix of identities of block and seizure number
      
      % the number of sequences is the number of columns in the current
      % sequence divided by 2
      nseq = size(currseq,2)/2;
      
      % The number of rows is the number of sequences in this
      % seizure/block, and then the first column is just the seizure
      % number, the 2nd is the block number, and then the third goes from 1
      % to the number of sequences
      toadd = [ones(nseq,2).*[i j],(1:nseq)'];
      block_seiz_ids = [block_seiz_ids; toadd];
       
   end
    
end

%% Clean

% Get the indices of the clean sequences
iclean =  spt_seqclust(Patient(pt).sz(1).block(1).data.xyChan,allseq,ss_thresh,tt_thresh);

% Get an array of the iclean indices and the subsequent ones going up by 2
% per sequence to help index allseq
y = zeros(length(iclean)*2, 1); 
y(1:2:end-1)=(iclean-1)*2+1; 
y(2:2:end)=(iclean-1)*2+2;

% Keep the clean ones
allseq = allseq(:,y);
keptseqs = block_seiz_ids(iclean,:);

%% Rebuild structure

% loop through the indices of the kept sequences
for i = 1:size(keptseqs,1)
    
    % getting the corresponding index in allseq
    allSeqIdx = (i-1)*2+1;
    
    % get the sequence
    seqData = allseq(:,allSeqIdx:allSeqIdx+1);
    
    % add it to the clean sequence data for the correct seizure and block
    Patient(pt).sz(keptseqs(i,1)).block(keptseqs(i,2)).data.cleanseq = ...
        [Patient(pt).sz(keptseqs(i,1)).block(keptseqs(i,2)).data.cleanseq,...
        seqData];   
end

%% Get recruitment latency and spatial organization
for i = 1:length(Patient(pt).sz)
   for j = 1:length(Patient(pt).sz(i).block)
       
       % Remove rows of all zeros
       Patient(pt).sz(i).block(j).data.cleanseq =...
           Patient(pt).sz(i).block(j).data.cleanseq(any(...
           Patient(pt).sz(i).block(j).data.cleanseq,2),:);
       
       % Get recruitment latency
        if size(Patient(pt).sz(i).block(j).data.cleanseq,2) >2
            [recruitmentLatencySingle,spikeCount] = ...
                getRecruitmentLatency(Patient(pt).sz(i).block(j).data.cleanseq,...
                Patient(pt).sz(i).block(j).data.xyChan);
            
            
            [Patient(pt).sz(i).block(j).data.avgRecruitmentLat,...
            Patient(pt).sz(i).block(j).data.spatialOrg] =...
            getSpatialOrg(recruitmentLatencySingle,Patient(pt).sz(i).block(j).data.xyChan,indexToms,dmin);
        else
            Patient(pt).sz(i).block(j).data.avgRecruitmentLat = nan;...
                Patient(pt).sz(i).block(j).data.spatialOrg = nan;
        end
       
   end
end

%% Put the spatial organizations together in a single array for each seizure
for i = 1:length(Patient(pt).sz)
   if isempty(Patient(pt).sz(i).runTimes) == 0
       for j = 1:size(Patient(pt).sz(i).runTimes,1)
           Patient(pt).sz(i).spatialOrg(j) = ...
               Patient(pt).sz(i).block(j).data.spatialOrg;
       end
       
       
   end
    
end

%save([resultsFolder,outputName],'Patient');
toc
end

