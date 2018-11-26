%% ictal_gdf

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
whichDetector = 8;

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
gdfFolder = [resultsFolder,'ictal_gdf/'];
chLocationsFolder = 'chLocations/';
newptfile = 'ptIctalGDF.mat';

%% Load file with filenames and run times
if merge == 1 && exist([resultsFolder,'ptStructs/',newptfile],'file') ~= 0
    fprintf('Adding current run to existing file %s\n',newptfile);
    load([resultsFolder,'ptStructs/',newptfile]);
else
    load([resultsFolder,'ptStructs/',timeFile]);
end

%% Loop through patients, szs, run times
for i = 1:length(pt)
    
    [~,tmul,absthresh] = icChsToIgnore(pt(i).name);
    
    pt(i).ictal_thresh.whichDetector = whichDetector;
    
    if isempty(tmul) == 1 || isempty(absthresh) == 1
        continue
    end
    
    if isempty(tmul) == 0
    pt(i).ictal_thresh.tmul = tmul;
    end
    
    if isempty(absthresh) == 0
    pt(i).ictal_thresh.absthresh = absthresh;
    end
    
    thresh =  pt(i).ictal_thresh;
    
    
    mkdir([gdfFolder,pt(i).name]);
    
    dataName =  pt(i).ieeg_name;
    if isempty(dataName) == 1
        continue
    end
    
    electrodeFile = pt(i).electrode_labels;
    if isempty(electrodeFile) ==1
        continue
    end
    
    if isempty(pt(i).ictal_thresh) == 1
        fprintf('Warning, ictal thresh empty for %s, skipping\n',pt(i).name);
    end
    
    % Skip it if I haven't entered a desired tmul
    if isempty(pt(i).ictal_thresh.tmul) == 1
        fprintf('Warning, no tmul for %s, skipping\n',pt(i).name);
        continue
    end
    
    
    for j = 1:size(pt(i).newSzTimes,1)
       

            fprintf('Doing %s seizure %d\n',...
                pt(i).name,j);
            
          
            % Define the desired times to be the seizure times
            % need to think about restricting further
            desiredTimes = [pt(i).newSzTimes(j,1) pt(i).newSzTimes(j,2)];
       
            filename = sprintf('%s_sz_%d_ictal.mat',...
                pt(i).name,j);
            
            pt(i).sz(j).ictal_fn = filename;
            
            if exist([gdfFolder,pt(i).name,'/',filename],'file') ~= 0
                if overwrite == 0
                    fprintf('File %s already found, skipping\n',filename);
                    continue
                end
            end
            
            % Run the spike detector
            
            
            [gdf,extraOutput] = getSpikesSimple(pt,i,desiredTimes,whichDetector,thresh,0);
            
            
          
                
        
            
            % Save gdf file
            save([gdfFolder,pt(i).name,'/',filename],'gdf','thresh','desiredTimes');
            
           
            % Resave pt file now that I have fs
            save([resultsFolder,'ptStructs/',newptfile],'pt');
            
       
    end
end



