clear


%% File names
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
%electrodeFile = [electrodeFolder,csvFile];
p1 = genpath(scriptFolder);
addpath(p1);
ptWithFs = 'ptWithfs.mat';
gdfFolder = 'gdf/';
chLocationsFolder = 'chLocations/';
ptWithSeq = 'ptWithSeq.mat';

%% Load file with filenames and run times
load([resultsFolder,ptWithFs]);

%% Loop through patients and seizures
for i = 1:length(pt)
    
    if exist(pt(i).chLocationFile,'file') == 0
        continue
    end
    
    for j = 1:length(pt(i).sz)
        if isfield(pt(i).sz,'runTimes') == 0
            continue
        end
        
        if isfield(pt,'fs') == 0
            continue
        end
       
        
        
        gdf_all = [];
        
        for k = 1:length(pt(i).sz(j).chunkFiles)
            
            if exist([pt(i).sz(j).chunkFiles{k}],'file') == 0
                continue
            end
            
            % Load gdf file
            load([pt(i).sz(j).chunkFiles{k}]);
            
            % Load chLocations file
            load([pt(i).chLocationFile]);
            
            if isempty(gdf) == 1
                continue
            end
            
            gdf(:,2) = gdf(:,2) + pt(i).sz(j).runTimes(k,1) - pt(i).sz(j).runTimes(1,1);
            gdf_all = [gdf_all;gdf];
            
        end
        
        % Now that you have all the spikes for the desired patient and
        % seizure, calculate sequences
        pt(i).sz(j).stats.nspikes = size(gdf_all,1);
        pt(i).sz(j).data = mainSequences(gdf_all,electrodeData, pt(i).fs);
        pt(i).sz(j).stats.nseqs = size(pt(i).sz(j).data.sequences,2)/2;
            
    end
    
    
end


save([resultsFolder,ptWithSeq],'pt');

