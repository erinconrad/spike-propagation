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

%}