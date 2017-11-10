
%{
TO DO:
- go through the spatial propagation code, especially the spatial
organization code, and check to make sure no glaring errors
- correct some of the hard coded numbers to make sure physiologic
- make some nice way to visualize sequences
- write in something to make it clear how many sequences are being removed
at each step
- decide if I want to include this sequence similarity matrix de-noising
step that they used in the Sam Tomlinson paper
- decide if I want to do bootstrapping

---------------------------------------------------------------------------
    Erin comments and code

- This is altered from code from Sam Tomlinson, written in 2014
Tomlinson, Samuel B., et al. "Spatiotemporal Mapping of Interictal Spike Propagation: A Novel Methodology Applied to Pediatric Intracranial EEG Recordings." Frontiers in neurology 7 (2016).

---------------------------------------------------------------------------
%}

%% parameters
ptname = 'HUP078';

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
interval          = 10000;                   % Hardcode number of spikes/seg
segments          = 1:interval:nspikes;

%{
WOAH what is this for - I think they were picking random segments or
something?
p                 = randperm(length(segments)-1);
p                 = p(1:10);                % Hardcode # of segments 
segments          = segments(p);  
          %}

chans             = Patient.xyChan(:,1); 

for s = 1:length(segments)
    
    start_row = segments(s);
    end_row = min(start_row + interval,length(gdf));
    disp((s/length(segments))*100)
    
    % Extract spike sequences from the gdf  
    [Patient.sequences{1,s},Patient.discarded{1,s}] = getSequencesErin(gdf(start_row:end_row,:), Patient.xyChan);
    
    % Store IEEG indices
    Patient.indices{1,s} = [start_row, end_row];
        
end


% Erin Code
Patient = getRecruitmentLatency(Patient);
Patient = getSpatialOrg(Patient);

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
    