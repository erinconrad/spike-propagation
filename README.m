%{

Pipeline for the Spatial Organization project

- written by Erin Conrad, 2017, University of Pennsylvania
- uses modified scripts from Sam Tomlinson to generate spike sequences
Tomlinson, Samuel B., et al. "Spatiotemporal Mapping of Interictal Spike Propagation: A Novel Methodology Applied to Pediatric Intracranial EEG Recordings." Frontiers in neurology 7 (2016).


To run the pipeline, navigate to forPorting and run doItAll. This assumes
you have:
- files in the expected locations with electrode locations
- json file with basic clinical data and electrodes to ignore
- a pw file for ieeg
- output folders for result

To visualize things, go to dataVisualization:
- to visualize spike detections, run portVisualizeSpikes
- to visualize sequence detections, run portVisualizeSequences
- to visualize gifs of individual sequences, run portVisualizeChLocations
- to visualize avg sequence paths, run portAvgPath
- to visuazlie spike frequency or spatial organization, run portChangingSO


Index:

MAIN FILES TO RUN THE PIPELINE (in forPorting):

- doItAll: calls the other main files in sequence to run the full pipeline
- getDesiredRunTimes: gets the files and times over which you want to run
the pipeline
- getElectrodeData: gets the electrode data for running the pipeline
- make_gdf: detects spikes (does very little on its own)
- portGetSpikes: called by make_gdf to get the data and then calls a
specific spike detector
- portGetSequences: detects sequences once you have spikes
- portGetSpatialOrg: calculates spatial organization once you have
sequences

-------------------------------------------

SPIKE DETECTION FILES (in detectSpikes):

- getSpikeTimes: similar to portGetSpikes but it can be called by various
other data visualization files
- getiEEGData: returns a chunk of iEEG data
- chParser: takes a string representing the name of an EEG channel and
returns a new string that agrees with iEEG naming conventions
- confirmFs: makes sure I have the sampling frequency included

Specific spike detectors (in specificDetectors in detectSpikes):

- spike_detector_hilbert_v16_nodownsample: Hoameng's spike detector
- spike_detector_Erin: Hoameng's spike detector with minor edits by Erin
- fspk2: Sam's spike detector
- fspk3: Sam's spike detector, edited by Erin so that the baseline
amplitude (used to determine if the spike crosses an amplitude threshold)
is calculated from a one-minute moving window rather than the whole
dataset)
- FindPeaks: a file called by fspk2 and fspk3 that finds peaks and troughs
in the EEG data
- eegfilt: a file called by fspk2 and fspk3 that filters the data using a
butterworth filter

--------------------------------------------

SEQUENCE DETECTION FILES (in getSequences):

- mainSequences: the main file called to detect sequences
- getSequencesErin: the main subFile called by mainSequences to get
sequences
- reorder_tiesErin: called by getSequencesErin to take the sequences
detected and reorder ties
- spatialConstraint: called by getSequencesErin to remove spikes that are
too far away
- decideIfIctal: decide if the sequences are ictal or interictal

Sequence cleaning files (in sequenceCleaning in getSequences):
- a bunch...

---------------------------------------

FILES TO GET SPATIAL INFORMATION FROM SEQUENCES (in spatialInfo):

- getRecruitmentLatency: get the recruitment latency for each sequence
- getSpatialOrg: takes a bunch of recruitment latencies and outputs the
average recruitment latency and the moran index
- getwij: calculates the weighting matrix
- moranStats: called by getSpatialOrg to get the moran index, and stats on
it. I got the stats using formulas from wikipedia...


%}