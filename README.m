%{

Pipeline for the spatial organization project

- written by Erin Conrad, 2017, University of Pennsylvania
- uses scripts from Sam Tomlinson (to generate spike sequences) and Radek
Janca

Tomlinson, Samuel B., et al. "Spatiotemporal Mapping of Interictal Spike Propagation: A Novel Methodology Applied to Pediatric Intracranial EEG Recordings." Frontiers in neurology 7 (2016).

Janca, Radek, et al. "Detection of interictal epileptiform discharges using signal envelope distribution modelling: application to epileptic and non-epileptic intracranial recordings." Brain topography 28.1 (2015): 172-183.




1: Detect spikes and output them to a gdf file
- to do this, go to detectSpikes and run getSpikeTimes. 
- need to input the file name you want to load from the iEEG portal
- this will create a gdf file with spike times
- this ignores channels that the json file tells you to ignore 

2: Make a structure containing the channel locations for the unignored
channels
- to do this, go to make ChannelStructs and run chanLocUseGdf
- this requires that the gdf file with spike times already exists because
it uses this to ensure that the channel indices in this struct will be the
same as the channel indices in the gdf file
- this also requires that you have a csv file with electrode locations

3: Combine the spike times from 1) and the channel locations from 2) to
calculate spike propagation information
- to do this, go to newSpatialOrgSubscripts and run main
- this takes the gdf file and the channel location file you created and
finds spike sequences and saves the sequences and spatial organization


%}