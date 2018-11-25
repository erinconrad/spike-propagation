clear

%% Remove EKG artifact
rmEKG = 0;


% how close the EKG spike can be from the EEG spike to throw out the EEG
% spike
prox = 0.02; %20 ms

% Remove depth electrodes
rmDepth = 0;
rmType = 'D';

% remove noisy electrodes? I don't think I need to because I remove noisy
% times for each electrode
rmNoisy = 0;


%% File names
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
%electrodeFile = [electrodeFolder,csvFile];
p1 = genpath(scriptFolder);
addpath(p1);
ptWithFs = 'ptIctalGDF.mat';
gdfFolder = [resultsFolder,'ictal_gdf/'];
chLocationsFolder = 'chLocations/';
ptWithSeq = 'ptIctalSeq.mat';

%% Load file with filenames and run times
load([resultsFolder,'ptStructs/',ptWithFs]);


%% Loop through patients and seizures
for i = 1:length(pt)
    
    % Get electrode data
    electrodeData =  pt(i).electrodeData;
    
    if isempty(pt(i).ictal_thresh) == 1
        continue;
    end
    
    if isempty(pt(i).ictal_thresh.tmul) == 1
        continue
    end
    
    for j = 1:length(pt(i).sz)
       
        if isempty(pt(i).electrodeData) == 1
            continue
        end
        
        if isfield(pt,'fs') == 0
            continue
        end
        
        
        gdf_all = [];
        %{
        empty_all = [];
        
        %}
        
        if isfield(pt(i).sz(j),'ictal_fn') == 0
            continue
        end
        
        filename = pt(i).sz(j).ictal_fn;
            
        if exist([gdfFolder,pt(i).name,'/',filename],'file') == 0
            continue
        end

        % Load gdf file
        load([gdfFolder,pt(i).name,'/',filename]);

        if isempty(gdf) == 1
            continue
        end

        if thresh.tmul ~= pt(i).ictal_thresh.tmul || ...
                thresh.absthresh ~= pt(i).ictal_thresh.absthresh
            error('Error, the thresholds in the gdf file are different from those in the pt struct\n');

        end

        % Sort by times
        times = gdf(:,2);
        chs = gdf(:,1);
        [times,I] = sort(times);
        chs = chs(I);
        gdf_sorted = [chs,times];

        if isequal(gdf,gdf_sorted) == 0
            fprintf('warning, gdf is not appropriately sorted\n');
        end
        gdf = gdf_sorted;


       
        gdf_all = [gdf_all;gdf];


        
        
        if rmDepth == 1
            gdf_all = removeChs(gdf_all,electrodeData,rmType);
        end
        
        % remove spikes that came from noisy times
        %{
        if rmNoisy == 1
            [gdf_all,rmNoisy,rmEmpty] = removeNoisy(gdf_all,noise_all,empty_all);
        end
        %}
        
         
        
        % Now that you have all the spikes for the desired patient and
        % seizure, calculate sequences
        pt(i).sz(j).stats.nspikes = size(gdf_all,1);
        
       if size(gdf_all,1) > 0
        pt(i).sz(j).data = mainSequences(gdf_all,electrodeData, pt(i).fs);
        pt(i).sz(j).stats.nseqs = size(pt(i).sz(j).data.sequences,2)/2;
        % Make fancy new matrix for sequences
        pt(i).sz(j).seq_matrix = ...
            makeSeqMatrix(pt(i).sz(j).data.sequences,length(pt(i).channels));
       else
           pt(i).sz(j).data.sequences = [];
           pt(i).sz(j).seq_matrix = [];
       end
        
       
        
        
      
        % vanleer
        %pt(i).sz(j).vanleer = vanleer_all;
    end
    
    
    
end

%scatter(linspace(pt(5).sz(1).runTimes(1,1),pt(5).sz(1).runTimes(end,2),length(noise_all)),noise_all(:,1));


save([resultsFolder,'ptStructs/',ptWithSeq],'pt');



