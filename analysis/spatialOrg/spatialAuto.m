function spatialAuto(pt,whichPts)


% Parameters
removeTies = 1;
doPlots = 0;

[~,~,scriptFolder,resultsFolder,~] = fileLocations;
p1 = genpath(scriptFolder);
addpath(p1);

if isempty(whichPts) == 1
    for i = 1:length(pt)
        if isempty(pt(i).seq_matrix) == 0
            whichPts = [whichPts,i];
        end
    end
end

destFolder = [resultsFolder,'spatialAuto/'];
mkdir(destFolder)

for whichPt = whichPts
    %% Patient parameters
    fprintf('Doing %s\n',pt(whichPt).name);
    locs = pt(whichPt).electrodeData.locs(:,2:4);
    nchs = size(locs,1);
    szTimes = pt(whichPt).newSzTimes;
    soz = pt(whichPt).newSOZChs; 
    [~,~,~,dmin,~] = ieegAndElectodeNames(pt(whichPt).name);
    wij = getwij(locs,dmin);
    outcome(whichPt) = getOutcome(pt(whichPt).name);
    
    
    %% Get all sequences
    seq_matrix = pt(whichPt).seq_matrix;
    
    %% Remove ties
    if removeTies == 1
        keep = ones(size(seq_matrix,2),1);
        for s = 1:size(seq_matrix,2)
           curr_seq = seq_matrix(:,s);
           nonans = curr_seq(~isnan(curr_seq));
           norepeats = unique(nonans);
           if length(norepeats) < 0.5*length(nonans)
               keep(s) = 0;
           end
        end
        seq_matrix(:,keep==0) = [];
        fprintf(['%s had %d sequences (%1.2f of all sequences) deleted'...
        'for having >50 percent ties\n%d sequences remain\n'],...
        pt(whichPt).name,sum(keep == 0),sum(keep == 0)/length(keep),sum(keep==1));

    end
    
    %% Remove ictal sequences
    first_time = min(seq_matrix,[],1);
    t = (any(first_time >= (szTimes(:,1)-repmat(60,size(szTimes,1),1)) ...
        & first_time <= szTimes(:,2),2));
    seq_matrix(:,t) = [];
    fprintf('Removed %d ictal spikes \n',sum(t));
    
    %% Get average recruitment Latency
    RL = nanmean(seq_matrix-min(seq_matrix,[],1),2);
    
    %% Info on observed versus expected per channel
    n_spikes = sum(~isnan(seq_matrix(:)));
    exp_per_ch = n_spikes/nchs;
    actual_per_ch = sum(~isnan(seq_matrix),2);
    
    %% Calculate moran index (spatial autocorrelation)
    MI(whichPt) = moranStats(RL',wij,nchs);
    
    %% Plot RL
    if doPlots == 1
        figure
        scatter3(locs(:,1),locs(:,2),locs(:,3),actual_per_ch+1,'k');
        hold on
        scatter3(locs(~isnan(RL),1),locs(~isnan(RL),2),locs(~isnan(RL),3),...
            actual_per_ch(~isnan(RL)),RL(~isnan(RL)),'filled');
        colorbar
        title(sprintf('%s RL, MI is %1.1f',pt(whichPt).name,MI.I));
    end
    
    
end



% get all moran indices and outcomes
allMI = [];
allOutcome = [];
allHospital = [];
for whichPt = whichPts
   allMI = [allMI;MI(whichPt).I]; 
   allOutcome = [allOutcome;outcome(whichPt)];
   if strcmp(pt(whichPt).name(1),'H') == 1
       allHospital = [allHospital;1];
   elseif strcmp(pt(whichPt).name(1),'S') == 1
       allHospital = [allHospital;2];
   else
       error('what');
   end
end

%% Do stats correlating MI with outcome

% Spearman correlation coefficient (non parametric rank)
[rho,pval] = corr(allMI(~isnan(allMI)),allOutcome(~isnan(allMI)),...
    'Type','Spearman');

%% Plot
figure
scatter(allMI,allOutcome,100,'filled');
xlabel('Spatial autocorrelation (Moran Index)');
ylabel('Outcome (modified Engel)');
title(sprintf('Autocorrelation vs outcome, p = %1.2f',pval))






end