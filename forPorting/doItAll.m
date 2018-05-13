%% doItAll

% This script just calls all of the subscripts, in order

clear

% loop through all patients and get the correct run times and file output
% locations
fprintf('Getting run times\n');
getDesiredRunTimes

% get electrode data for the patients
fprintf('Getting electrode data\n');
getElectrodeData

% Detect all spikes, output spike detections in 2000 second chunks in a gdf
% file
fprintf('Detecting spikes\n');
make_gdf

% Calculate sequences from spikes
fprintf('Detecting sequences\n');
portGetSequences

% Calculate the spatial organization
fprintf('Calculating stats on sequences\n');
portGetSpatialOrg