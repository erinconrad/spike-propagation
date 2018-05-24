clear

%% Remove EKG artifact
rmEKG = 1;
prox = 0.02; % NEED TO CHANGE

%% Remove depth electrodes
rmDepth = 1;
rmType = 'D';

%% File names
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
%electrodeFile = [electrodeFolder,csvFile];
p1 = genpath(scriptFolder);
addpath(p1);
ptWithFs = 'ptPostGDF.mat';
gdfFolder = [resultsFolder,'gdf/'];
chLocationsFolder = 'chLocations/';
ptWithSeq = 'ptWithSeq.mat';

%% Load file with filenames and run times
load([resultsFolder,'ptStructs/',ptWithFs]);

%% Loop through patients and seizures
for i = 1:length(pt)
    
    % Get electrode data
    electrodeData =  pt(i).electrodeData;
    
    for j = 1:length(pt(i).sz)
        if isfield(pt(i).sz,'runTimes') == 0
            continue
        end
        
        if isempty(pt(i).electrodeData) == 1
            continue
        end
        
        if isfield(pt,'fs') == 0
            continue
        end
        
        if isempty(pt(i).fs) == 1
            pt(i).fs = 512;
            fprintf('Warning, no fs for patient %d, assuming it is 512 Hz\n',i);
        end
       
        
        
        gdf_all = [];
        gdf_ekg_all = [];
        empty_all = [];
        noise_all = [];
        
        for k = 1:length(pt(i).sz(j).chunkFiles)
            
            if exist([gdfFolder,pt(i).name,'/',pt(i).sz(j).chunkFiles{k}],'file') == 0
                continue
            end
            
            % Load gdf file
            load([gdfFolder,pt(i).name,'/',pt(i).sz(j).chunkFiles{k}]);
            
            if isempty(gdf) == 1
                continue
            end
            
            % Load gdf ekg file
            load([gdfFolder,pt(i).name,'/',pt(i).sz(j).EKGchunkFiles{k}]);
            
            % remove EKG artifact
            if rmEKG == 1
                gdf = removeEKGArtifact(gdf,gdf_ekg,prox);
            end
            
            gdf(:,2) = gdf(:,2) + pt(i).sz(j).runTimes(k,1) - pt(i).sz(j).runTimes(1,1);
            bad.noise(:,1) = bad.noise(:,1) + pt(i).sz(j).runTimes(k,1) - pt(i).sz(j).runTimes(1,1);
            bad.empty = bad.empty + pt(i).sz(j).runTimes(k,1) - pt(i).sz(j).runTimes(1,1);
            
            if isempty(gdf_ekg) == 0
                gdf_ekg(:,2) = gdf_ekg(:,2) + pt(i).sz(j).runTimes(k,1) - pt(i).sz(j).runTimes(1,1);
            end
            
            gdf_all = [gdf_all;gdf];
            gdf_ekg_all = [gdf_ekg_all;gdf_ekg];
            empty_all = [empty_all;bad.empty];
            noise_all = [noise_all;bad.noise];
            
        end
        
        if rmDepth == 1
            gdf_all = removeChs(gdf_all,electrodeData,rmType);
        end
        
        % Now that you have all the spikes for the desired patient and
        % seizure, calculate sequences
        pt(i).sz(j).stats.nspikes = size(gdf_all,1);
        pt(i).sz(j).data = mainSequences(gdf_all,electrodeData, pt(i).fs);
        pt(i).sz(j).stats.nseqs = size(pt(i).sz(j).data.sequences,2)/2;
        pt(i).sz(j).stats.empty = empty_all;
        pt(i).sz(j).stats.noise = noise_all;
        pt(i).sz(j).stats.noise_ictal = zeros(size(pt(i).sz(j).stats.noise,1),1);
        
        % for noise function, get if it's ictal
        for t = 1:size(pt(i).sz(j).stats.noise,1)
            time = pt(i).sz(j).stats.noise(t,1) + pt(i).sz(j).runTimes(1,1);
            if time >= pt(i).sz(j).onset && ...
                    time <= pt(i).sz(j).offset
                pt(i).sz(j).stats.noise_ictal(t) = 1;
            end
            
        end
        
        % Add EKG times
        pt(i).sz(j).ekg = gdf_ekg_all;
            
    end
    
    
end


save([resultsFolder,'ptStructs/',ptWithSeq],'pt');

