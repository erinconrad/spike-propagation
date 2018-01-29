function Patient = mainSequences(gdf,electrodeData)

%{

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

% If I don't pass a gdf, use the one written here
if nargin == 0
    filename = 'gdf267000-268000.mat';
    ptname = 'HUP078';

end

% max time of candidate spike from first spike in sequence (50 ms in paper)
t1 = .05; 

% max time from preceding spike (15 ms in paper)
t2 = .015; 

% minimum sequence length (how many spikes per sequence, 5 in the paper)
minSeqLength = 5; 

% This is my criteria for potentially throwing out a sequence if there are
% too many simultaneous spikes. If this is set to 0, I never throw any
% sequences out. If this is set to 1, I use the minimum number of
% non-simultaneous spikes criterion. If this is set to 2, I use the maximum
% percentage that is ties criterion.
uniquenessCheck = 0;

% minimum number of non-simultaneous spikes in a sequence
minUniqueSeqLength = 3;

% maximum percentage of a sequence that is ties (when the spike hits
% multiple channels at the same time). This was not a criterion in the
% Tomlinson paper. 
maxPercTies = 0.7; 

% some distance between channels for channel weights. In the paper this was
% 1.5 cm so I need to figure out the conversion.
dmin = 10.5; 

% This is used when two channels are not within an allowable distance of
% each others, but if the 2 channels spike together this frequently, then
% I allow the jump regardless. This was 0.05 in the paper.
minConcurrentFreqRel = 0.05;

% This is used when two channels are not within an allowable distance of
% each others, but if the 2 channels spike together this frequently, then
% I allow the jump regardless. This was 5 in the paper.
minConcurrentFreqAbs = 5;

% The allowable distance between 2 channels of successive spikes. This is
% new in my methods because the Tomlinson paper instead used a
% "partitioning" approach and required that successive spikes be in
% adjacent partitions. Need to think about this number.
maxDist = 20;

% How long the spike sequence needs to be after the spatial constraint
minSpikesCloseEnough = 3;

% I don't know what this does.
indexToms = 1;

% Cleaning parameters
ss_thresh     = 1.5;                      
tt_thresh     = 3;  


%% Load data files if I don't pass a gdf
if nargin == 0
    spikePaths
    load([dataLoc,'gdf/',filename]);
    load([dataLoc,'electrodeData/',ptname,'_chanLoc.mat']);
    addpath('/Users/erinconrad/Desktop/residency stuff/R25/actual work/scripts/my scripts/sequenceCleaning/');


end

%% Add path to the cleaning scripts

%% Initialize patient variables
Patient                = struct;

% nx4 array, where n is the number of unignored channels, 1st column is
% index, and 2nd-4th columns are x,y,z positions
Patient.xyChan         = electrodeData.locs;
Patient.gdf            = gdf; 

% Prepare to loop through gdf 
nspikes           = size(gdf,1);
chans             = Patient.xyChan(:,1); 


%% Get sequences
[Patient.sequences,Patient.discarded] = ...
    getSequencesErin(gdf, Patient.xyChan,...
    t1, t2, minSeqLength, uniquenessCheck, minUniqueSeqLength, ...
    maxPercTies, minConcurrentFreqAbs, minConcurrentFreqRel,...
    maxDist, minSpikesCloseEnough);

%% Do cleaning step
% At some point I should consider whether to do the cleaning step in each
% block (as I am doing now) or with all the data. If I do it with all the
% data then I will have it on n_blocks*n_seizures as much data, but it will
% be much more complicated to do

Patient.dirtysequences = Patient.sequences;
Patient =  spt_seqclust(Patient,ss_thresh,tt_thresh);
iclean = Patient.seq_track;
y = zeros(length(iclean)*2, 1); y(1:2:end-1)=(iclean-1)*2+1; y(2:2:end)=(iclean-1)*2+2;
Patient.sequences = Patient.sequences(:,y);

Patient.discarded.cleaning = -size(Patient.sequences,2)/2 + size(Patient.dirtysequences,2)/2;
Patient.discarded.total = Patient.discarded.total + Patient.discarded.cleaning;
Patient.discarded.remaining = Patient.discarded.origNum - Patient.discarded.total;


%% Get recruitment latency and spatial organization
Patient = getRecruitmentLatency(Patient);
Patient = getSpatialOrg(Patient,indexToms,dmin);



% get and print total number of sequences remaining and percent discarded
totalSeq = 0;
totalDisc = 0;

totalSeq = totalSeq + Patient.discarded.remaining;
totalDisc = totalDisc +  Patient.discarded.total;


percentDisc = totalDisc/(totalSeq+totalDisc)*100;
Patient.discarded.percentdiscarded = percentDisc;

%fprintf('There are %d total sequences remaining. %1.1f percent were discarded for violating rules\n',...
%    totalSeq,percentDisc);

%Output patient data 
if nargin==0
    save([resultLoc,ptname,'_sequences_noTieRestrictions.mat'],'Patient')  
end

end
    