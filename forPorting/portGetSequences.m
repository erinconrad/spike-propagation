clear

%% Remove EKG artifact
rmEKG = 0;

%% Do vanleer
doVanleer = 0;

% how close the EKG spike can be from the EEG spike to throw out the EEG
% spike
prox = 0.02; %20 ms

% Remove depth electrodes
rmDepth = 1;
rmType = 'D';

% remove noisy electrodes? I don't think I need to because I remove noisy
% times for each electrode
rmNoisy = 0;

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
    
    if isempty(pt(i).thresh.tmul) == 1
        continue
    end
    
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
        
        %{
        if isempty(pt(i).fs) == 1
            pt(i).fs = 512;
            fprintf('Warning, no fs for patient %d, assuming it is 512 Hz\n',i);
        end
       %}
        
        
        gdf_all = [];
        gdf_ekg_all = [];
        noise_all = [];
        vanleer_all.spike_times = [];
        vanleer_all.rms = [];
        vanleer_all.delay = [];
        %{
        empty_all = [];
        
        %}
        
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
            
            %gdf(:,2) = gdf(:,2) + pt(i).sz(j).runTimes(k,1) - pt(i).sz(j).runTimes(1,1);
            
            %{
            bad.noise(:,1) = bad.noise(:,1) + pt(i).sz(j).runTimes(k,1) - pt(i).sz(j).runTimes(1,1);
            bad.empty = bad.empty + pt(i).sz(j).runTimes(k,1) - pt(i).sz(j).runTimes(1,1);
            %}
            
            %{
            if isempty(gdf_ekg) == 0
                gdf_ekg(:,2) = gdf_ekg(:,2) + pt(i).sz(j).runTimes(k,1) - pt(i).sz(j).runTimes(1,1);
            end
            %}
            
            gdf_all = [gdf_all;gdf];
            gdf_ekg_all = [gdf_ekg_all;gdf_ekg];
            
            if doVanleer == 1
            vanleer_all.spike_times = [vanleer_all.spike_times;vanleer.spikeTimes];
            vanleer_all.delay = [vanleer_all.delay;vanleer.delay];
            vanleer_all.rms = [vanleer_all.rms;vanleer.rms];
            end
            
            
           % noise_all = [noise_all;noise];
            %{
            empty_all = [empty_all;bad.empty];
            
            %}
            
        end
        
        if rmDepth == 1
            gdf_all = removeChs(gdf_all,electrodeData,rmType);
        end
        
        % remove spikes that came from noisy times
        %{
        if rmNoisy == 1
            [gdf_all,rmNoisy,rmEmpty] = removeNoisy(gdf_all,noise_all,empty_all);
        end
        %}
        
        
        % Plot noise
        %{
        chans_noise = [10 15 20 25 30];
        y_offset = linspace(1,length(chans_noise)*4,length(chans_noise));
        x_loc = linspace(pt(i).sz(j).runTimes(1,1),...
            pt(i).sz(j).runTimes(1,1)+60*size(noise_all,1),size(noise_all,1));
        figure
        for chan = 1:length(chans_noise)
            col = zeros(size(noise_all,1),3);
            for nn = 1:size(noise_all,1)
                if noise_all(nn,chans_noise(chan)) == -1, col(nn,:) = [0 0 1]; else, col(nn,:) = [1 0 0]; end
            end

            scatter(x_loc,(noise_all(:,chans_noise(chan))+y_offset(chan)),100,col,'filled')
            hold on
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
        
        %% Cleaning
        %{
        dirty_seq = pt(i).sz(j).data.sequences;
        ss_thresh     = 15;              
        tt_thresh     = 0.015;  

        if size(dirty_seq,2) > 2
            
            iclean =  spt_seqclust(pt(i).electrodeData.locs,...
                dirty_seq,ss_thresh,tt_thresh);
            y = zeros(length(iclean)*2, 1); y(1:2:end-1)=(iclean-1)*2+1; y(2:2:end)=(iclean-1)*2+2;
            pt(i).sz(j).data.sequences = dirtysequences(:,y);
            
            
        end
        %}
        
        
        
        %{
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
        
        %}
        
        % Add EKG times
        pt(i).sz(j).ekg = gdf_ekg_all;
        
        
        
        
      
        % vanleer
        %pt(i).sz(j).vanleer = vanleer_all;
    end
    
    
    
end

%scatter(linspace(pt(5).sz(1).runTimes(1,1),pt(5).sz(1).runTimes(end,2),length(noise_all)),noise_all(:,1));

save([resultsFolder,'ptStructs/',ptWithSeq],'pt');

