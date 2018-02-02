%{

Pipeline for the Spatial Organization project

- written by Erin Conrad, 2017, University of Pennsylvania
- uses scripts from Sam Tomlinson (to generate spike sequences) and Radek
Janca (to detect spikes)

Tomlinson, Samuel B., et al. "Spatiotemporal Mapping of Interictal Spike Propagation: A Novel Methodology Applied to Pediatric Intracranial EEG Recordings." Frontiers in neurology 7 (2016).

Janca, Radek, et al. "Detection of interictal epileptiform discharges using signal envelope distribution modelling: application to epileptic and non-epileptic intracranial recordings." Brain topography 28.1 (2015): 172-183.


To run this pipeline, go to the folder "fullPipeline" and run main. In
main, you should change the patient name, the output file name, the
ieeg.org data name, and the electrode file name, and the patient number, as
well as the sampling rate of the EEG data. The pipeline works as follows:

1) it gets the seizure onset and offset times for each seizure that the
patient has
2) from this, it defines the start and stop times for each block we want to
detect sequences over prior to each seizure
3) it detects spikes in each block
4) it detects sequences of spikes in each block
5) it calculates the spatial organization of the spike sequences for each
block

Non-script files needed to run this pipeline, in addition to knowing what
file you want to look up on ieeg.org, include the json file with clinical
info (seizure times and electrodes to ignore), as well as the csv file with
electrode locations. These csv files can be found on borel in
/gdrive/public/USERS/lkini/aim3/results/PATIENT NAME/aim1

This code tends to fail when requesting larger chunks of data, but only on
Borel. It raises PermGen memory errors in the IEEG database code,
suggesting that somehow on Borel this code is creating a memory leak

%}