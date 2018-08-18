clear

%% Remove depth electrodes
rmDepth = 1;
rmType = 'D';

%% File names
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
%electrodeFile = [electrodeFolder,csvFile];
p1 = genpath(scriptFolder);
addpath(p1);
ptWithFs = 'ptWithSeq.mat';
gdfFolder = [resultsFolder,'gdf/'];
chLocationsFolder = 'chLocations/';
ptVanleer = 'ptVanleer.mat';

% number of permutations in permutation test
nboot = 1e3;

%% Load file with filenames and run times
load([resultsFolder,'ptStructs/',ptWithFs]);

% initialize the rms array, the delay array, which seizure it's in array,
% and seizure or no seizure array
allrms = [];
alldelay = [];
whichsz = [];
szOrNot = [];

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
        
        for k = 1:length(pt(i).sz(j).chunkFiles)
            
            if exist([gdfFolder,pt(i).name,'/',pt(i).sz(j).chunkFiles{k}],'file') == 0
                continue
            end
            
            % Load gdf file
            load([gdfFolder,pt(i).name,'/',pt(i).sz(j).chunkFiles{k}]);
            
            if isempty(vanleer) == 1
                continue
            end
            
            % root mean square power
            allrms = [allrms;vanleer.rms];
           
            % delays
            alldelay = [alldelay;vanleer.delay];
            whichsz = [whichsz;j];

            % an array saying if the spike is ictal or not
            szOrNot = [szOrNot;vanleer.ictal];
            
        end
    end
    
    %% Get delays and power for seizures and not seizures
    pt(i).vanleer.alldelay_ic = alldelay(logical(szOrNot),:);
    pt(i).vanleer.alldelay_inter = alldelay(~logical(szOrNot),:);
    pt(i).vanleer.allrms_ic = allrms(logical(szOrNot),:);
    pt(i).vanleer.allrms_inter = allrms(~logical(szOrNot),:);

    pt(i).vanleer.avgdelay_ic = nanmean(pt(i).vanleer.alldelay_ic);
    pt(i).vanleer.avgdelay_inter = nanmean(pt(i).vanleer.alldelay_inter);
    pt(i).vanleer.avgrms_ic = nanmean(pt(i).vanleer.allrms_ic);
    pt(i).vanleer.avgrms_inter = nanmean(pt(i).vanleer.allrms_inter);

            
    %% Use PCA to reduce the dimensionality
    fprintf('Doing PCA\n');
    % concatenate the delay and rms power features
    all_features = [alldelay,allrms];

    % run PCA
    [coeff,score,latent,tsquared,explained,mu] = pca(all_features);

    % figure out how many dimensions I need to keep to explain 99% of the
    % variance
    sum_e = 0;
    for x = 1:length(explained)
       sum_e = sum_e + explained(x);
       if sum_e > 99
           ndim = x;
           break
       end
    end

    new_coeff = coeff(:,1:ndim);
    new_scores = score(:,1:ndim);

    % Run a k-medoids algorithm to cluster the data, do it 30 times for fun
    fprintf('Doing clustering\n');
    for x = 1:30
        [cluster_assignment{x},C{x},sumd{x},~,midx{x},~] = kmedoids(new_scores,10,'Distance','cityblock');

        metric(x) = sum(sumd{x});
    end


    % Take the results of the clustering algorithm that worked the best
    [~,minidx] = min(metric);
    cluster_assignment = cluster_assignment{minidx};
    C = C{minidx};
    sumd = sumd{minidx};
    midx = midx{minidx};




    %% Do bootstrap
    %stats =  vPermTest(cluster_assignment,szOrNot,10,nboot);
    
    

    %% Things for plots
    

    % Reconstruct all_features with reduced dimensionality
    temp = new_scores*new_coeff';
    n_delay = temp(:,1:size(temp,2)/2); n_rms = temp(:,size(temp,2)/2+1:end);


    gold_delay = alldelay(midx,:); gold_rms = allrms(midx,:);
    
    pt(i).vanleer.stats = stats;
    pt(i).vanleer.gold_delay = gold_delay;
    pt(i).vanleer.gold_rms = gold_rms;
    pt(i).vanleer.cluster_assignment = cluster_assignment;

    % Plots

    for t = 1:3
       plotDelays(gold_delay(t,:),gold_rms(t,:),pt(i).electrodeData.locs); 

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

save([resultsFolder,'ptStructs/',ptVanleer],'pt')

fprintf('end\n');

