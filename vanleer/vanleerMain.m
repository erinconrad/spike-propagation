function [delay,rms,electrodeData] = vanleerMain

%{

This function calculates the spatiotemporal propagation patterns of spikes
using methods by Ann Vanleer published in:

Vanleer et al. Millimeter-scale epileptiform spike propagation
patterns and their relationship to seizures. J Neural Eng. 2016; 13(2).

Briefly, detects spikes and for each spike looks at a bin of time around
that spike. Within each bin it calculates the rms power of each channel and
the delay at which each channel has its power maximum in that bin. Using
these delay and rms power maps as features, it then does a PCA to reduce
the dimensionality of the features. Using this reduced set of features, it
then clusters spikes into one of 10 different clusters. It then compares
the proportion of spikes in each cluster that are ictal versus interictal
to see if the proportion differs among the clusters.


%}

%% Parameters to change every time

% Output file name to save
outputName = 'HUP80_test.mat';
%outputName = 'HUP78_oneMinBlocks.mat';

% data name (for ieeg.org)
dataName = 'HUP80_phaseII';
%dataName = 'HUP78_phaseII-Annotations';  

% CSV file with electrode locations
csvFile = 'HUP080_T1_19991213_electrode_labels.csv';
%csvFile = 'HUP078_T1_19971218_electrode_labels.csv';

% The patient name with format as used in the json file
ptname = 'HUP080';
%ptname = 'HUP078';

% The number of the patient
pt = 80;
%pt = 78;

% segment time (50 ms segments going from 2 ms before the spike detection
% to 48 ms after the spike detection)
vtime = [-0.002,0.048];

% number of permutations in permutation test
nboot = 1e4;

% Remove EKG artifact? This is not set up to be able to do this
rmEKGArtifact = 0;

%% Get paths and load seizure info and channel info
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
electrodeFile = [electrodeFolder,csvFile];
p1 = genpath(scriptFolder);
addpath(p1);
ptInfo = loadjson(jsonfile);


Patient(pt).seizures = ptInfo.PATIENTS.(ptname).Events.Ictal;
szNames = fieldnames(Patient(pt).seizures);

%% Load EEG data info
% calling this with 0 and 0 means I will just get basic info like sampling
% rate and channel labels
data = getiEEGData(dataName,0,0,pwfile);  
fs = data.fs;

%% Run the getSpikes script once as a dummy run just to produce a file of electrode locations
[~,electrodeData,~] = getSpikeTimes(0,dataName,electrodeFile,ptInfo,pwfile,1,0,vtime,0,0);

%% Define seizure onset and offset times for each seizure
for i = 1:length(fieldnames(Patient(pt).seizures))
    Patient(pt).sz(i).onset = Patient(pt).seizures.(szNames{i}).SeizureEEC;
    Patient(pt).sz(i).offset = Patient(pt).seizures.(szNames{i}).SeizureEnd;
    Patient(pt).sz(i).duration = Patient(pt).sz(i).offset - Patient(pt).sz(i).onset;
end

%% Define the start and stop times of each block
% For each seizure, I will have 2 blocks: an ictal block containing all the
% time in the seizure, and a pre-ictal block, containing an equivalent
% amount of time immediately prior to the seizure. I should consider
% potentially using more removed pre-ictal data, or possibly both.

% Loop through all the seizures
for i = 1:length(Patient(pt).sz)
    
    
    
    % Skip the seizure if it's too close to the start of the data
    if Patient(pt).sz(i).onset < Patient(pt).sz(i).duration+2
        continue
    end

    % Skip the seizure if it's too close to the prior seizure
    if i~=1
        if Patient(pt).sz(i).onset - Patient(pt).sz(i-1).offset < Patient(pt).sz(i).duration+1
           continue 
        end
    end
    
    % Go back in time the duration of the seizure to get the pre-ictal
    % initial time
    initialTime = Patient(pt).sz(i).onset-Patient(pt).sz(i).duration-1;
    
    
    Patient(pt).sz(i).preSzStart =  initialTime;
    Patient(pt).sz(i).preSzEnd = initialTime + Patient(pt).sz(i).duration;
      
    
end

%% Run the full analysis on each block

% initialize the rms array, the delay array, which seizure it's in array,
% and seizure or no seizure array
allrms = [];
alldelay = [];
whichsz = [];
szOrNot = [];

% Loop through all seizures
for i = 1:2%length(Patient(pt).sz)
    
    fprintf('Doing seizure %d of %d\n',i,length(Patient(pt).sz));
    
   
       % Loop through ictal and pre-ictal
       for j = 1:2
           
           if j == 1
               %seizure
               issz = 1;
               desiredTimes = [Patient(pt).sz(i).onset,Patient(pt).sz(i).offset];
           elseif j == 2
               %preseizure
               issz = 0;
               desiredTimes = [Patient(pt).sz(i).preSzStart,Patient(pt).sz(i).preSzEnd];

           end
          
           
           
           % Detect spikes and get rms power and delay maps for each block
           fprintf('Detecting spikes\n');
           [gdf,~,~] = getSpikeTimes(desiredTimes,dataName,electrodeFile,ptInfo,pwfile,0,1,vtime,0,0);
           
           
           
           % root mean square power
           allrms = [allrms;gdf.rms];
           
           % delays
           alldelay = [alldelay;gdf.delay];
           whichsz = [whichsz;i*ones(size(gdf.rms,1),1)];
           
           % an array saying if the spike is ictal or not
           szOrNot = [szOrNot;issz*ones(size(gdf.rms,1),1)];
           
      
        end
    
end


%% Use PCA to reduce the dimensionality

% concatenate the delay and rms power features
all_features = [alldelay,allrms];

% run PCA
[coeff,score,latent,tsquared,explained,mu] = pca(all_features);

% figure out how many dimensions I need to keep to explain 99% of the
% variance
sum_e = 0;
for i = 1:length(explained)
   sum_e = sum_e + explained(i);
   if sum_e > 99
       ndim = i;
       break
   end
end

new_coeff = coeff(:,1:ndim);
new_scores = score(:,1:ndim);

% Run a k-medoids algorithm to cluster the data, do it 30 times for fun
for i = 1:30
    [cluster_assignment{i},C{i},sumd{i},~,midx{i},~] = kmedoids(new_scores,10,'Distance','cityblock');
    
    metric(i) = sum(sumd{i});
end


% Take the results of the clustering algorithm that worked the best
[~,minidx] = min(metric);
cluster_assignment = cluster_assignment{minidx};
C = C{minidx};
sumd = sumd{minidx};
midx = midx{minidx};




%% Do bootstrap
stats =  vPermTest(cluster_assignment,szOrNot,10,nboot);


%% Things for plots
if 1 == 1

    % Reconstruct all_features with reduced dimensionality
    temp = new_scores*new_coeff';
    n_delay = temp(:,1:size(temp,2)/2); n_rms = temp(:,size(temp,2)/2+1:end);


    gold_delay = alldelay(midx,:); gold_rms = allrms(midx,:);

    % Plots

    for i = 1:size(gold_delay,1)
       plotDelays(gold_delay(i,:),gold_rms(i,:),electrodeData.locs); 

    end

    idx = cluster_assignment;
    figure;
    X=new_scores;
    plot(X(idx==1,1),X(idx==1,2),'r.','MarkerSize',7)
    hold on
    plot(X(idx==2,1),X(idx==2,2),'b.','MarkerSize',7)
    hold on
    plot(X(idx==3,1),X(idx==3,2),'g.','MarkerSize',7)
    plot(C(1:3,1),C(1:3,2),'co',...
         'MarkerSize',7,'LineWidth',1.5)
    legend('Cluster 4','Cluster 5','Cluster 6','Medoids',...
           'Location','NW');
    title('Cluster Assignments and Medoids');
    hold off
    end


end