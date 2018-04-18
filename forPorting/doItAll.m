%% doItAll

% This script just calls all of the subscripts, in order

clear

% loop through all patients and get the correct run times and file output
% locations
getDesiredRunTimes

% Detect all spikes, output spike detections in 2000 second chunks in a gdf
% file
make_gdf

% Calculate sequences from spikes
portGetSequences

% Calculate the spatial organization
portGetSpatialOrg