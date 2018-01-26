
%{
TO DO:
- [] try source localization instead? I may be falsely assuming that
wherever I see the spike is where the actual source of the spike is.
Instead it could be that there is a spike below the surface that I am
seeing at multiple locations on the cortex at the same time
- [] think about the fact that so many spikes seem to occur at the same
time. Think about the time steps and confirm that I am not losing sig figs.
It looks like the closest time steps are 0.005 seconds, so sampling at 200
Hz, which is what Sam's paper says.
- [] compare my code to Sam's
- [] go through the spatial propagation code, especially the spatial
organization code, and check to make sure no glaring errors
- [] correct some of the hard coded numbers to make sure physiologic
- [] decide if I want to include this sequence similarity matrix de-noising
step that they used in the Sam Tomlinson paper
- [] decide if I want to do bootstrapping
- [] can I optimize spike detection algorithm for the higher sampling rate?

- [x] come up with nice way to visualize sequences

---------------------------------------------------------------------------
    Erin comments and code

- This is altered from code from Sam Tomlinson, written in 2014
Tomlinson, Samuel B., et al. "Spatiotemporal Mapping of Interictal Spike Propagation: A Novel Methodology Applied to Pediatric Intracranial EEG Recordings." Frontiers in neurology 7 (2016).

---------------------------------------------------------------------------

    Variable descriptions

Patient:
- this is the main output variable
- Patient.xyChan contains the xyz coordinates of all the channels
- Patient.gdf contains the channel location and time of each spike
- Patient.sequences is a 1xn cell, where n is the number of segments the
data is divided into. Within each cell, there are nChannel cells, each of
which contains an mx(2*s) array. s is the number of spike sequences
starting in that channel (and it is 2*s wide because the first row contains
the channel number and the second row contains the time) and m is the
number of steps in the spike sequence.
%}

%% parameters
ptname = 'HUP078';

% max time of candidate spike from first spike in sequence (50 ms in paper)
t1 = .05; 

% max time from preceding spike (15 ms in paper)
t2 = .015; 

% minimum sequence length (how many spikes per sequence, 5 in the paper)
minSeqLength = 5; 

% maximum percentage of a sequence that is ties (when the spike hits
% multiple channels at the same time. I don't think this is made explicit
% in the paper. NOTE!!! I am not using this currently as the latest dropbox
% code took out this requirement. 
maxPercTies = 0.7; 

indexToms = 1; % Conversion between indices and ms (I am not using this right now, which concerns me)

dmin = 10.5; % some distance between channels for channel weights

% This is used when two channels are not within an allowable distance of
% each others, but if the 2 channels spike together this frequently, then
% I allow the jump regardless. This was 0.05 in the paper.
minConcurrentFreq = 0.05;

% The allowable distance between 2 channels of successive spikes. This is
% new in my methods because the Tomlinson paper instead used a
% "partitioning" approach and required that successive spikes be in
% adjacent partitions. Need to think about this number.
maxDist = 4;


interval = 10000;  
% Hardcode number of spikes/seg


%% Load data files
spikePaths
load([dataLoc,'gdf/',ptname,'_gdf.mat']);
load([dataLoc,'electrodeData/',ptname,'_chanLoc.mat']);


%% Initialize patient variables
Patient                = struct;

% nx4 array, where n is the number of unignored channels, 1st column is
% index, and 2nd-4th columns are x,y,z positions
Patient.xyChan         = electrodeData.locs;
Patient.gdf            = gdf; 

% Prepare to loop through gdf 
nspikes           = size(gdf,1);
segments          = 1:interval:nspikes;

%{
WOAH what is this for - I think they were picking random segments or
something?
p                 = randperm(length(segments)-1);
p                 = p(1:10);                % Hardcode # of segments 
segments          = segments(p);  
          %}

chans             = Patient.xyChan(:,1); 

%% Get sequences

for s = 1:length(segments)
    
    start_row = segments(s);
    end_row = min(start_row + interval,length(gdf));
    disp((s/length(segments))*100)
    
    % Extract spike sequences from the gdf  
    [Patient.sequences{1,s},Patient.discarded{1,s}] = ...
        getSequencesErin(gdf(start_row:end_row,:), Patient.xyChan,...
        t1, t2, minSeqLength, maxPercTies, minConcurrentFreq, maxDist);
    
    % Store IEEG indices
    Patient.indices{1,s} = [start_row, end_row];
        
end


%% Erin Code - get recruitment latency and spatial organization
Patient = getRecruitmentLatency(Patient);
Patient = getSpatialOrg(Patient,indexToms,dmin);

if 1 == 0 % NEED TO FIGURE THIS OUT
% Get statistics from spike sequences
[Patient.statistics,Patient.centroids,Patient.mpv,Patient.repeats,...
    Patient.network]= get_stats(Patient, chans, interval, nRandomNetworks); % Erin now passes interval and nRandomNetworks
end

% get and print total number of sequences remaining and percent discarded
totalSeq = 0;
totalDisc = 0;
for i = 1:size(Patient.discarded)
   totalSeq = totalSeq + Patient.discarded{i}.remaining;
   totalDisc = totalDisc +  Patient.discarded{i}.total;
end

percentDisc = totalDisc/(totalSeq+totalDisc)*100;

fprintf('There are %d total sequences remaining. %1.1f percent were discarded for violating rules\n',...
    totalSeq,percentDisc);

%Output patient data 
save([resultLoc,ptname,'_sequences.mat'],'Patient')  
    