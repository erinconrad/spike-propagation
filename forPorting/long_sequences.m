%{

This script detects spike sequences once we have a bunch of spike
detections

%}


clear

%% Do vanleer
doVanleer = 0;

% Remove depth electrodes
rmDepth = 0;
rmType = 'D';

% Should we try to merge the patient structure with an existing, incomplete
% patient structure?
merge = 1;


%% File names
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
%electrodeFile = [electrodeFolder,csvFile];
p1 = genpath(scriptFolder);
addpath(p1);
ptWithFs = 'long_gdf.mat';
gdfFolder = [resultsFolder,'long_gdf/'];
chLocationsFolder = 'chLocations/';
ptWithSeq = 'long_seq.mat';


%% Load gdf or sequence file
if merge == 1 && exist([resultsFolder,'ptStructs/',ptWithSeq ],'file') ~= 0
    load([resultsFolder,'ptStructs/',ptWithSeq ]);
else
    load([resultsFolder,'ptStructs/',ptWithFs]);
end


%% Loop through patients and seizures
for i = 1:31
    
    fprintf('Doing %s\n',pt(i).name);
    
    % Get electrode data
    electrodeData =  pt(i).electrodeData;
    
    if isempty(pt(i).thresh.tmul) == 1
        continue
    end
    
    
    if isfield(pt(i),'runTimes') == 0
        continue
    end

    if isempty(pt(i).electrodeData) == 1
        continue
    end

    if isfield(pt,'fs') == 0
        continue
    end



    gdf_all = [];
    removed_all = [];


    for k = 1:length(pt(i).chunkFiles)

        if exist([gdfFolder,pt(i).name,'/',pt(i).chunkFiles{k}],'file') == 0
            continue
        end

        % Load gdf file
        load([gdfFolder,pt(i).name,'/',pt(i).chunkFiles{k}]);

        if isempty(gdf) == 1
            continue
        end

        if thresh.tmul ~= pt(i).thresh.tmul || ...
                thresh.absthresh ~= pt(i).thresh.absthresh
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
        removed_all = [removed_all;removed];

      
    end

    if rmDepth == 1
        gdf_all = removeChs(gdf_all,electrodeData,rmType);
    end

    % Get stats on spike rates per channel
    if isempty(gdf_all) == 0
        nchs = length(pt(i).channels);
        npch = zeros(nchs,1);
        for ich = 1:nchs
            npch(ich) = sum(gdf_all(:,1) == ich);
        end
        rate(i).npch = npch;
    end
    
    % Now that you have all the spikes for the desired patient and
    % seizure, calculate sequences
    pt(i).stats.nspikes = size(gdf_all,1);
    
    if size(gdf_all,2) >2, doMorph =1; else, doMorph = 0; end

   if size(gdf_all,1) > 0
    pt(i).data = mainSequences(gdf_all,electrodeData, pt(i).fs);
    pt(i).stats.nseqs = size(pt(i).data.sequences,2)/2;
    % Make fancy new matrix for sequences
    pt(i).seq_matrix = ...
        makeSeqMatrix(pt(i).data.sequences,length(pt(i).channels),doMorph);
   else
       pt(i).data.sequences = [];
       pt(i).seq_matrix = [];
   end
   
   pt(i).removed = removed_all;

    %% Cleaning (I don't do this for spike paper)
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



    % vanleer
    %pt(i).sz(j).vanleer = vanleer_all;
    
    
    
    
end

%scatter(linspace(pt(5).sz(1).runTimes(1,1),pt(5).sz(1).runTimes(end,2),length(noise_all)),noise_all(:,1));


%save([resultsFolder,'ptStructs/',ptWithSeq],'pt');
save([resultsFolder,'ptStructs/','rate.mat'],'rate');



