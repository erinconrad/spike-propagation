clear

%% Remove EKG artifact
rmEKG = 0;

%% Do vanleer
doVanleer = 0;

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
ptWithFs = 'ptPostVanleerGDF.mat';
gdfFolder = [resultsFolder,'gdf_vanleer/'];
chLocationsFolder = 'chLocations/';
ptWithSeq = 'ptWithVanleerSeq.mat';

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
        vanleer_all.times = [];
        vanleer_all.rms = [];
        %{
        empty_all = [];
        
        %}
        
        for k = 1:length(pt(i).sz(j).chunkFiles)
            
            if exist([gdfFolder,pt(i).name,'/',pt(i).sz(j).chunkFiles{k}],'file') == 0
                continue
            end
            
            % Load gdf file
            load([gdfFolder,pt(i).name,'/',pt(i).sz(j).chunkFiles{k}]);
          
            
            if thresh.tmul ~= pt(i).thresh.tmul || ...
                    thresh.absthresh ~= pt(i).thresh.absthresh
                error('Error, the thresholds in the gdf file are different from those in the pt struct\n');
               
            end
            
            if isempty(vanleer) == 1
                continue
            end
            
          
            vanleer_all.times = [vanleer_all.times;vanleer.times];
            vanleer_all.rms = [vanleer_all.rms;vanleer.rms];
         
            
            
           % noise_all = [noise_all;noise];
            %{
            empty_all = [empty_all;bad.empty];
            
            %}
            
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
        
        pt(i).sz(j).seq_matrix = vanleerGetSeqs(vanleer_all);
      
        
       
        % vanleer
        %pt(i).sz(j).vanleer = vanleer_all;
    end
    
    
    
end

%scatter(linspace(pt(5).sz(1).runTimes(1,1),pt(5).sz(1).runTimes(end,2),length(noise_all)),noise_all(:,1));

save([resultsFolder,'ptStructs/',ptWithSeq],'pt');

