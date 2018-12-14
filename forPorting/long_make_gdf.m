%% make_gdf

% This script takes a patient structure with information about what times
% to look for spikes over, and then this actually detects those spikes and
% outputs them to a gdf file


clear

%% Parameters

% Should I re-run the spike detection and overwrite gdf file if it already
% exists?
overwrite = 0; 

% Should we try to merge the patient structure with an existing, incomplete
% patient structure?
merge = 1;

%% File names
% starting file (unless we are merging with an existing gdf file)
timeFile = 'long_electrode_data.mat'; 

% output file
newptfile = 'long_gdf.mat';

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
p1 = genpath(scriptFolder);
addpath(p1);
gdfFolder = [resultsFolder,'long_gdf/'];
chLocationsFolder = 'chLocations/';


%% Load file with filenames and run times
if merge == 1 && exist([resultsFolder,'ptStructs/',newptfile],'file') ~= 0
    load([resultsFolder,'ptStructs/',newptfile]);
else
    load([resultsFolder,'ptStructs/',timeFile]);
end

%% Loop through patients, szs, run times
% which patients to run over
for i = [1,2,5,6,7,10,11,13,14,15,16,21,23]%[3 4 8 9 12 17 18 19 20 22 24 25 27 30 31]
    
    % load the info about which detector, the tmul, and the absthresh
    thresh =  pt(i).thresh;

    if thresh.whichDetector ~= 7
        error('look!');
    end
    
    % make patient directory
    mkdir([gdfFolder,pt(i).name,'/']);
    
    % Get ieeg.org data name
    dataName =  pt(i).ieeg_name;
    if isempty(dataName) == 1
        continue
    end
    
    % get electrode file
    electrodeFile = pt(i).electrode_labels;
    if isempty(electrodeFile) ==1
        continue
    end
    
    % Skip it if I haven't entered a desired tmul
    if isempty(pt(i).thresh.tmul) == 1
        fprintf('Warning, no threshold info entered for %s\n',pt(i).name);
        continue
    end
    
    % Skip if no run times
    if isfield(pt(i),'runTimes') == 0
        fprintf('Warning, no run times entered for %s\n',pt(i).name);
        continue
    end
    
    % Loop through run times and run spike detector
    for k = 1:size(pt(i).runTimes,1)

        fprintf('Doing chunk %d of %d for patient %d \n',...
            k,size(pt(i).runTimes,1),i);
    
        desiredTimes = [pt(i).runTimes(k,:)];

       % Skip the run if it already exists and we're not overwriting
        if exist([gdfFolder,pt(i).name,'/',pt(i).chunkFiles{k}],'file') ~= 0
            if overwrite == 0
                fprintf('File %s already found, skipping\n',[pt(i).chunkFiles{k}]);
                continue
            end
        end

        % Run the spike detector
        tic
        [gdf,extraOutput] = getSpikesSimple(pt,i,desiredTimes,thresh.whichDetector,thresh,0);
        toc

        % Save gdf file
        removed = extraOutput.removed;
        save([gdfFolder,pt(i).name,'/',pt(i).chunkFiles{k}],...
            'gdf','thresh','removed');

        % Resave pt file now that I have fs
        save([resultsFolder,'ptStructs/',newptfile],'pt');

    end
   
end

