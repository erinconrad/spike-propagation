
%% Pipeline for the Spike location project
% written by Erin Conrad, 2019, University of Pennsylvania

%% OVERVIEW
%{
To run the full spike location analysis, you will need access to the
electrode locations, the clinical info, and the ieeg portal. Then, do the
following steps:

1) Detect spikes and spike sequences
2) Cluster spikes according to location
3) Analyze cluster changes
4) Alpha/delta analysis
5) Area of influence analysis

Note that the electrode locations and clinical info for the patients in the
Spike location project can be found in the example patient structure
example_pt.mat.

%}


%% RUN COMPLETE PIPELINE ON EXAMPLE DATA
%{
To illustrate the steps of this project, we included a single script
run_example_data, which completes the entire pipeline from spike detection 
through clustering by spike location for a 500 second sample of EEG data 
for a single patient (Study 029). To run this pipeline, load the .mat file
example_pt.mat. Then call the script run_example_data.m, located in the
folder prepare_example/. Because this contains such a short duration of
data, time series analyses cannot meaningfully be performed on it. However,
a longer duration of data, structured according to this data structure,
could be run with the same pipeline. Alternatively, analysis can be run on
ieeg data using the steps described above.

%}



%% 1: DETECT SPIKES AND SPIKE SEQUENCES (in forPorting)
%{
To do this, you will need to run the following scripts in the following
order (all in forPorting). It will output a pt structure with information
on the spike sequences.

- longRunTimes: gets the files and times over which you want to run
the pipeline

- long_electrode_data: gets the electrode data for running the pipeline

- long_make_gdf: detects spikes (does very little on its own, but calls
scripts in detectSpikes)

- long_sequences: detects sequences once you have spikes (calls scripts in
 getSequences)

%}

%% 2: CLUSTER SPIKES ACCORDING TO LOCATION (in analysis/clustering)
%{
To do this, you will need to run the following script (in
analysis/clustering)

- getClusters: this takes a pt struct with spike sequence info produced
from Step 1 and outputs a cluster structure and also plots a bunch of
sample sequences for each cluster.

After running this, you will need to manually visually examine each of the
sample sequences and then update the cluster struct with information about
which clusters are bad (those with >50% false sequences).

%}

%% 3: ANALYZE CLUSTER CHANGES (in analysis/clustering)
%{
This contains 2 analyses which can be done in any particular order.

- CNewStats: this takes a pt struct with spike
sequence info, a cluster structure with cluster info, and does statistical
tests to check for a change in cluster distribution over time, as well as
tests to compare the pre-ictal, inter-ictal, and post-ictal cluster
distributions. It also does some clinical correlations.

- CVar: this takes a pt struct with spike sequence
info and a cluster struct with cluster info and calculates the number of
hours needed to capture 80% of the variability in spike location.

%}

%% 4: Alpha/delta analysis
%{
This analysis first consists of calculating the alpha/delta ratio and then
correlating this with spike location info. And so the following steps must
be done in the order written:

- alphaDelta (in detectSpikes/alphaDelta): this will load a pt struct with
spike sequence info (although it will really only use the run times from
this struct) and it will grab ieeg data and calculate the alpha delta power
ratio in 2000 second bins). It calls innerAlphaDelta to do the actual work
(which is in the same folder). This will output a power struct with the
alpha delta power ratio.

- AD_AR (in analysis/alphaDelta): this take a pt struct with spike sequence
info, a cluster struct with cluster info, and a power struct with
alpha/delta power ratio and correlate spike location with alpha/delta
ratio.

- r_stats.r and r_lin.r (in analysis/r_stats): these r scripts take
data returned from AD_AR and perform a generalized linear model with
autoregressive errors for the analyses correlating alpha/delta ratio with
the proportion of spikes in the predominant cluster as well as the spike
rate.

- r_stats.m and r_rate_stats.m (in analysis/r_stats): these Matlab scripts
take the data returned from the R scripts and do aggregate level analyses
to determine if there is overall a significant correlation between sleep
and spike spatial distribution/spike rate

- glm_model.m (in analysis/alphaDelta): this Matlab script determines 
whether there is a relationship between the proportion of spikes in the 
predominant cluster and the alpha delta ratio, assuming a GLM model for the
relationship, rather than a GLARMA model (for which you would need to use R).

- ad_remove_spikes.m (in detectSpikes/confirm_ad): this script is the first
part of the analysis to test if the alpha/delta ratio is affected by
spikes. It recalculates the alpha/delta ratio removing spike times from the
calculation over each electrode.

- ad_check_analysis.m (in detectSpikes/confirm_ad): this is the second part
of the analysis to test the effect of spikes on the alpha/delta ratio. It
takes the new and old alpha-delta ratios (keeping and removing spikes) and
performs a Pearson correlation coefficient between the two.

%}

%% 5: Area of influence analysis
%{
This is a single script in analysis/influence:

- CInfluence: this take a pt struct, a cluster struct, and calculates
distances of various electrodes of interest, including the area with the
largest area of influence, from the nearest SOZ.

%}

%% Dependency 1: SPIKE DETECTION FILES (in detectSpikes)

%{
- getSpikesSimple: this is called by long_make_gdf (in Step 1) and grabs
data for spike detection, calls the specific spike detector, and does some
artifact removal
- getiEEGData: returns a chunk of iEEG data
- chParser: takes a string representing the name of an EEG channel and
returns a new string that agrees with iEEG naming conventions
- confirmFs: makes sure I have the sampling frequency included

Specific spike detectors (in detectSpikes/specificDetectors):
- fspk6: this is the spike detector I used. Confusingly, it is called when
whichDetector is 7.

%}

%% Dependency 2: SEQUENCE DETECTION FILES (in getSequences)
%{

- mainSequences: the main file called to detect sequences
- getSequencesErin: the main subFile called by mainSequences to get
sequences
- reorder_tiesErin: called by getSequencesErin to take the sequences
detected and reorder ties
- spatialConstraint: called by getSequencesErin to remove spikes that are
too far away
- decideIfIctal: decide if the sequences are ictal or interictal

%}



