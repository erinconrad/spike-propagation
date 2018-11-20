%% make_gdf

% This script takes a patient structure with information about what times
% to look for spikes over, and then this actually detects those spikes and
% outputs them to a gdf file


clear

%% Parameters

% 4 is the fspk3 detector, which is my edited version of the Marsh lab
% detector, edited to use a moving window of one minute over which I
% calculate the baseline amplitude for the relative amplitude threshold

% 5 is fspk4, which is like fspk3 but it uses different parameters if
% looking at depth electrodes
whichDetector = 7;

% Should I re-run the spike detection and overwrite gdf file if it already
% exists?
overwrite = 0; 

% Should we try to merge the patient structure with an existing, incomplete
% patient structure?
merge = 1;

%% File names
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
%electrodeFile = [electrodeFolder,csvFile];
p1 = genpath(scriptFolder);
addpath(p1);
timeFile = 'long_electrode_data.mat'; 
gdfFolder = [resultsFolder,'long_gdf/'];
mkdir(gdfFolder)
chLocationsFolder = 'chLocations/';
newptfile = 'long_gdf.mat';

%% Load file with filenames and run times
if merge == 1 && exist([resultsFolder,'ptStructs/',newptfile],'file') ~= 0
    load([resultsFolder,'ptStructs/',newptfile]);
else
    load([resultsFolder,'ptStructs/',timeFile]);
end

%% Loop through patients, szs, run times
for i = [22,24,25]%:length(pt) % STILL NEED TO DO 17
    pt(i).thresh.whichDetector = whichDetector;
    thresh =  pt(i).thresh;
    
    
    mkdir([resultsFolder,'long_gdf/',pt(i).name]);
    
    dataName =  pt(i).ieeg_name;
    if isempty(dataName) == 1
        continue
    end
    
    electrodeFile = pt(i).electrode_labels;
    if isempty(electrodeFile) ==1
        continue
    end
    
    % Skip it if I haven't entered a desired tmul
    if isempty(pt(i).thresh.tmul) == 1
        continue
    end
    
    if isfield(pt(i),'runTimes') == 0
        continue
    end
    
    
    
    for k = 1:size(pt(i).runTimes,1)

        fprintf('Doing chunk %d of %d for patient %d \n',...
            k,size(pt(i).runTimes,1),i);

    
        desiredTimes = [pt(i).runTimes(k,:)];


       
        if exist([gdfFolder,pt(i).name,'/',pt(i).chunkFiles{k}],'file') ~= 0
            if overwrite == 0
                fprintf('File %s already found, skipping\n',[pt(i).chunkFiles{k}]);
                continue
            end
        end



        % Run the spike detector

        tic
        [gdf,extraOutput] = getSpikesSimple(pt,i,desiredTimes,whichDetector,thresh,0);
        toc



        vanleer = extraOutput.vanleer;
        removed = extraOutput.removed;
        % Save gdf file
        save([gdfFolder,pt(i).name,'/',pt(i).chunkFiles{k}],'gdf','removed','thresh');


        % Resave pt file now that I have fs
        save([resultsFolder,'ptStructs/',newptfile],'pt');

    end
   
end

