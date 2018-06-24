function Patient = mainSequences(gdf,electrodeData,fs)

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
- Patient.sequences is an mx(2*s) array. s is the number of spike sequences
 (and it is 2*s wide because the first row contains
the channel number and the second row contains the time) and m is the
number of steps in the spike sequence.
%}

%% parameters
% The choice of these is at this point arbitrary and is very important for
% the final result

% If I don't pass a gdf, use the one written here
if nargin == 0
    filename = 'gdfTemp.mat';
    ptname = 'HUP078';
    fs =  512;

end

% max time of candidate spike from first spike in sequence (50 ms in paper)
t1 = .05; 

% max time from preceding spike (15 ms in paper)
t2 = .015; 

% minimum sequence length for a sequence to even initially be considered
% (how many spikes per sequence, 5 in the paper)
minSeqLength = 5;  %5

% This is my criteria for potentially throwing out a sequence if there are
% too many simultaneous spikes. If this is set to 0, I never throw any
% sequences out. If this is set to 1, I use the minimum number of
% non-simultaneous spikes criterion. If this is set to 2, I use the maximum
% percentage that is ties criterion. The paper did not have any check (so
% this was set to 0)
uniquenessCheck = 0;

% minimum number of non-simultaneous spikes in a sequence. (Only used if
% uniquenessCheck = 1.)
minUniqueSeqLength = 3;

% maximum percentage of a sequence that is ties (when the spike hits
% multiple channels at the same time). (Only used if uniquenessCheck = 2).
maxPercTies = 0.7; 

% some distance between channels for channel weights. In the paper this was
% 1.5 cm. For me 1 unit is 1 mm. 
dmin = 15; 

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
% adjacent partitions. I downloaded a paper Hirsh1991 (Synaptic physiology
% of horizontal connections in the cats visual cortex), which says that CV
% is about 0.3-1 m/s. Let's assume that the fastest it can go is 10 m/s.
% This is 10000 mm/s. I believe my spatial units are in mm and my time units
% are in s.
maxSpeed = 10000; %10,000 mm/s

% How long the spike sequence needs to be after the spatial constraint.
% This was 5 in Sam's code but not explicitly mentioned in the paper.
minSpikesCloseEnough = 5; %5

% I don't know what this does.
indexToms = 1;

cleaning = 0;

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
ss_thresh     = 15;              
tt_thresh     = 0.015;  


%% Load data files if I don't pass a gdf
if nargin == 0
    spikePaths
    load([dataLoc,'gdf/',filename]);
    load([dataLoc,'electrodeData/',ptname,'_chanLoc.mat']);
    addpath('/Users/erinconrad/Desktop/residency stuff/R25/actual work/scripts/my scripts/sequenceCleaning/');


end


%% Initialize patient variables
Patient                = struct;

% nx4 array, where n is the number of unignored channels, 1st column is
% index, and 2nd-4th columns are x,y,z positions
Patient.xyChan         = electrodeData.locs;


%% Get sequences
[sequences,Patient.discarded] = ...
    getSequencesErin(gdf, Patient.xyChan,...
    t1, t2, minSeqLength, uniquenessCheck, minUniqueSeqLength, ...
    maxPercTies, minConcurrentFreqAbs, minConcurrentFreqRel,...
    maxSpeed, fs, minSpikesCloseEnough);

%% Do cleaning step
% Now I do this at the end once all blocks are obtained

%{
dirtysequences = sequences;
if cleaning == 1
    
    if size(sequences,2) >2
    iclean =  spt_seqclust(Patient.xyChan,sequences,ss_thresh,tt_thresh);
    y = zeros(length(iclean)*2, 1); y(1:2:end-1)=(iclean-1)*2+1; y(2:2:end)=(iclean-1)*2+2;
    sequences = dirtysequences(:,y);
    end

    %% Get recruitment latency and spatial organization
    if size(sequences,2) >2
    [recruitmentLatencySingle,spikeCount] = getRecruitmentLatency(sequences,Patient.xyChan);
    [Patient.avgRecruitmentLat,Patient.spatialOrg] = getSpatialOrg(recruitmentLatencySingle,Patient.xyChan,indexToms,dmin);
    else
        Patient.avgRecruitmentLat = nan;Patient.spatialOrg = nan;
    end

    Patient.discarded.cleaning = -size(sequences,2)/2 + size(dirtysequences,2)/2;
    Patient.discarded.total = Patient.discarded.total + Patient.discarded.cleaning;

end
%}


Patient.discarded.remaining = Patient.discarded.origNum - Patient.discarded.total;
Patient.sequences = sequences;



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
    