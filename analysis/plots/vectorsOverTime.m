function vectorsOverTime(pt,whichPts)

%% Parameters
sm_span = 100;
window = 3600;

[~,~,~,resultsFolder,~] = fileLocations;


for whichPt = whichPts

saveFolder = [resultsFolder,'plots/',pt(whichPt).name,'/','vectors/'];

% reinitialize the sequences
seq_all = {};
    

% Patient specific stuff
dmin = pt(whichPt).dmin;
nchs = length(pt(whichPt).channels);



% Get wij
xyChan = pt(whichPt).electrodeData.locs;
wij = getwij(xyChan,dmin);


%% First need to divide it up into multiple seizure chunks to avoid overlapping

% just look at first seizure
sequences = pt(whichPt).sz(1).seq_matrix;
sequences(sequences==0) = nan; % WHY ARE THERE ANY ZEROS?????
firstSpikes = min(sequences,[],1);  

% remove sequences occuring during the seizure
sz = [pt(whichPt).sz(1).onset pt(whichPt).sz(1).offset];
sequences(:,firstSpikes>=sz(1) & firstSpikes<=sz(2)) = [];
seq_all{1} = sequences;


% Loop through the other seizures
for j = 2:length(pt(whichPt).sz)
    szTimeLast = pt(whichPt).sz(j-1).onset;
    szTime = pt(whichPt).sz(j).onset;
    totalTime = pt(whichPt).sz(j).runTimes(end,2) - pt(whichPt).sz(j).runTimes(1,1);
    
    sequences = pt(whichPt).sz(j).seq_matrix;
    sequences(sequences==0) = nan; % WHY ARE THERE ANY ZEROS?????

    
    if szTime - szTimeLast > totalTime
        % If the seizure time is more than 24 hours after the last seizure,
        % put this in a new chunk for plotting purposes
        seq_all{end+1} = sequences;
    else
        % if it is less than 24 hours after the last seizure, then need to
        % add any spikes that occured after the last seizure run time to
        % this new chunk (and ignore spikes from this seizure that would
        % have already been captured in the first seizure run)
        keepAfter = pt(whichPt).sz(j-1).runTimes(end,2);
        firstSpikes = min(sequences,[],1);  
        seq_to_keep = sequences(:,firstSpikes > keepAfter);

        % remove sequences occuring during the seizure
        sz = [pt(whichPt).sz(j).onset pt(whichPt).sz(j).offset];
        firstSpikes = min(seq_to_keep,[],1);  
        seq_to_keep(:,firstSpikes>=sz(1) & firstSpikes<=sz(2)) = [];
        
        seq_all{end} = [seq_all{end},seq_to_keep];
    end

end


%% Recombine all sequences for the patient
all_seq_cat = [];
trackingNo = [];
for i = 1:length(seq_all)
    seq = seq_all{i};
    all_seq_cat = [all_seq_cat,seq];
    trackingNo = [trackingNo, i*ones(1,size(seq,2))];
end

all_times = min(all_seq_cat,[],1);


%% Get all of the vectors
[all_vecs,early,late] = (getVectors2(all_seq_cat,pt(whichPt).electrodeData));

%% Plot an example vector
%{
figure
x = xyChan(:,2);
y = xyChan(:,3);
z = xyChan(:,4);
scatter3(x,y,z,60,'markeredgecolor','k');
hold on
e=scatter3(early(500,1),early(500,2),early(500,3),...
        60,'g','filled');
l=scatter3(late(500,1),late(500,2),late(500,3),...
        60,'r','filled');
plot3([early(500,1) late(500,1)],...
        [early(500,2) late(500,2)],...
        [early(500,3) late(500,3)],...
        'k','LineWidth',2);
grid off
axis off
title('Example spike sequence');
legend([e l],{'Early ch mean','Late ch mean'});
set(gca,'FontSize',15);
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
set(gca,'ZTickLabel',[]);
set(gcf,'Position',[50 100 500 500])
saveas(gcf,[saveFolder,pt(whichPt).name,'example_seq.png']);
close(gcf)
%}

%% Test if the early locations are different from the late locations
% This basically tells me if there is an overall direction, that the
% vectors aren't completely randomly oriented

% MANOVA and Hotelling's give the same result

% MANOVA
whichOne = [ones(size(early,1),1);zeros(size(late,1),1)];
[d,p_diff_MANOVA,stats] = manova1([early;late],whichOne,0.05);
diff_real = norm(mean(late)-mean(early));

% Hotelling's
early_HT = [ones(length(early),1), early];
late_HT = [2*ones(length(late),1), late];
p_diff_HT = HotellingT2([early_HT;late_HT],0.05);

%% Plot histograms to see if vectors are normally distributed
% Need to be more or less normal to justify the above tests

% Early points
figure
set(gca,'FontSize',15)
whichAx(1) = 'x';
whichAx(2) = 'y';
whichAx(3) = 'z';

for i = 1:3
    
set(gca,'FontSize',15)
subplot(3,2,(i-1)*2+2)
histogram(normrnd(mean(early(:,i)),std(early(:,i)),size(early(:,i))),100)
title(sprintf('Random normal data for %s component of early mean',whichAx(i)))
set(gca,'FontSize',15)
xl = get(gca,'xlim');
    
    
set(gca,'FontSize',15)
subplot(3,2,(i-1)*2+1)
histogram(early(:,i),100)
title(sprintf('Real data for %s component of early mean',whichAx(i)))
set(gca,'xlim',xl);
set(gca,'FontSize',15)  


end

makePlotPretty
set(gcf,'Position',[50 100 900 800])
saveas(gcf,[saveFolder,pt(whichPt).name,'early_hist.png']);
close(gcf)

% Late
figure
set(gca,'FontSize',15)
whichAx(1) = 'x';
whichAx(2) = 'y';
whichAx(3) = 'z';
for i = 1:3
    
set(gca,'FontSize',15)
subplot(3,2,(i-1)*2+2)
histogram(normrnd(mean(late(:,i)),std(late(:,i)),size(late(:,i))),100)
title(sprintf('Random normal data for %s component of late mean',whichAx(i)))
set(gca,'FontSize',15)    
xl = get(gca,'xlim');

set(gca,'FontSize',15)   
set(gca,'FontSize',15)
subplot(3,2,(i-1)*2+1)
histogram(late(:,i),100)
title(sprintf('Real data for %s component of late mean',whichAx(i)))
set(gca,'xlim',xl);
set(gca,'FontSize',15)  

end

makePlotPretty;
set(gcf,'Position',[50 100 900 800])
saveas(gcf,[saveFolder,pt(whichPt).name,'late_hist.png']);
close(gcf)


%% Plot the vectors over time
% I am using a smoothing function
figure
for i = 1:length(seq_all)
    
    temp_all_times = all_times(:,trackingNo==i);
    temp_all_vecs = all_vecs(trackingNo==i,:);
    
    x=plot(temp_all_times/3600,smooth(temp_all_vecs(:,1),sm_span),'b','LineWidth',2);
    hold on
    y=plot(temp_all_times/3600,smooth(temp_all_vecs(:,2),sm_span),'r','LineWidth',2);
    z=plot(temp_all_times/3600,smooth(temp_all_vecs(:,3),sm_span),'g','LineWidth',2);
end

for j = 1:length(pt(whichPt).sz)
   yl = get(gca,'ylim'); 
   szOnset = pt(whichPt).sz(j).onset;
   sz = plot([szOnset szOnset]/3600,yl,'k--','LineWidth',2);
end

xlabel('Hour');
ylabel('Magnitude and sign of vector component (mm)');
legend([x,y,z,sz],{'x-component','y-component','z-component','seizure times'});

title(sprintf('Spike propagation vector for %s',pt(whichPt).name));
set(gca,'FontSize',15)

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

set(gcf,'Position',[50 100 1200 400])

mkdir(saveFolder)
saveas(gcf,[saveFolder,pt(whichPt).name,'vec_time.png']);
close(gcf)


%% Plot histograms to see if vectors are normally distributed
figure

set(gca,'FontSize',15)
whichAx(1) = 'x';
whichAx(2) = 'y';
whichAx(3) = 'z';

for i = 1:3

set(gca,'FontSize',15)
subplot(3,2,(i-1)*2+2)
histogram(normrnd(mean(all_vecs(:,i)),std(all_vecs(:,i)),size(all_vecs(:,i))),100)
title(sprintf('Random normal data for %s component of vector',whichAx(i)))
set(gca,'FontSize',15)    
xl = get(gca,'xlim');
    
set(gca,'FontSize',15)
subplot(3,2,(i-1)*2+1)
histogram(all_vecs(:,i),100)
title(sprintf('Real data for %s component of vector',whichAx(i)))
set(gca,'FontSize',15)
set(gca,'xlim',xl);


end
makePlotPretty;
set(gcf,'Position',[50 100 900 800])
saveas(gcf,[saveFolder,pt(whichPt).name,'hist.png']);
close(gcf)

%% Test if vectors are different for each different hour long chunk
nchunks = ceil((all_times(end)-all_times(1))/window);
times = zeros(nchunks,2);

% define chunks
for tt = 1:nchunks
     times(tt,:) = [(tt-1)*window + all_times(1),...
         tt*window + all_times(1)];
end

whichChunk = zeros(size(all_vecs,1),1);
for i = 1:size(all_times,2)
    for tt = 1:nchunks
        if all_times(i) >= times(tt,1) && all_times(i) <= times(tt,2)
            whichChunk(i) = tt;
        end
    end
end


[d,p_chunk_comparison,stats] = manova1(all_vecs,whichChunk,0.05);

fprintf(['For %s:\nEarly vs late channel position difference p = %1.2e\n'...
    'Change in vector over time p_x = %1.2e, p_y = %1.2e, p_z=%1.2e\n'],...
    pt(whichPt).name,p_diff_HT,p_chunk_comparison(1),...
    p_chunk_comparison(2),p_chunk_comparison(3));

end


