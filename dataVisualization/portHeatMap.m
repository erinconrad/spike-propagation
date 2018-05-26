function portHeatMap(P,pt,sz,whichToPlot)

%% Parameters
ptname = P(pt).name;
delay = 0.4; % time delay between steps

% output file name
[~,~,~,resultsFolder,~] = fileLocations;
filename = [sprintf('%s',ptname),'_sz_',sprintf('%d',sz),...
    '_heatmap.gif'];

% base size of electrodes
baseSizeElec = 20;

% base size of lines
baseSizeLine = 0.5;

% multiplication factor for electrodes
multElec = 3;

% multiplication factor for lines
multLines = 1;


%% Initialize channel arrays
% Get the channel locations for the patient
chLocs = P(pt).sz(sz).data.xyChan(:,2:4);
allChs = P(pt).sz(sz).data.xyChan(:,1);


allRec = [];
for block = 1:length(P(pt).sz(sz).blockRL)
    allRec = [allRec,P(pt).sz(sz).blockRL(block).avgRecruitmentLat];
end

crange = [nanmedian(allRec) - 2*nanstd(allRec),nanmedian(allRec) + 2*nanstd(allRec)];

%% Loop through blocks
for block = 1:length(P(pt).sz(sz).blockRL)
    
    chHits = zeros(size(allChs));
    pathHits = zeros(length(allChs),length(allChs));
    
    %% Get the sequence data
    % Get the columns from the sIdx
    sIdx = P(pt).sz(sz).blockRL(block).sIdx;
    sequences = P(pt).sz(sz).data.sequences;
    newseq = [];
    for k = 1:length(sIdx)
        newseq = [newseq,sequences(:,k*2-1:k*2)];
    end
    temp = newseq;

    % Remove the times, now this has as many columns as there are sequences
    temp = temp(:,1:2:end-1);


    %% Loop through the sequences to fill up the arrays

    for i = 1:size(temp,2)

        % Get a single sequence
       seq = temp(:,i); 
       seq = seq(any(seq~=0,2),:); % remove rows of zeros
       seq = round(seq);

       % Loop through the spikes in the sequence
       for j = 1:length(seq)

          % add a hit each time the channel is hit
          chHits(seq(j)) = chHits(seq(j)) + 1; 

          % add a path hit each time the path is hit
          if j ~=1

              % for all non first spikes, save the pair of the prior channel
              % and the current channel as a hit on the path
              pathHits(seq(j-1),seq(j)) = pathHits(seq(j-1),seq(j)) + 1;
          end

       end

    end

    %% plotting

    fig = figure;

    if whichToPlot == 1
        toPlot = chHits;
    elseif whichToPlot == 2
        toPlot = P(pt).sz(sz).blockRL(block).avgRecruitmentLat;
    end
    % make all channels a base color
    scatter3(chLocs(:,1),chLocs(:,2),chLocs(:,3),...
        baseSizeElec,toPlot,'filled')
    hold on

   
    grid off
    axis off
    title(['HUP ',sprintf('%d',pt),' average recruitment latency for block ',sprintf('%d',block)]);
    set(gca,'FontSize',15);
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    set(gca,'ZTickLabel',[]);
    colormap jet
    if whichToPlot == 2
        caxis(crange)
    end
    colorbar
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.4, 0.4, 0.3, 0.45]);

    % capture the figure as a frame in the gif
    F(block) = getframe(fig);
    im = frame2im(F(block));
    [imind,cm] = rgb2ind(im,256);
    
    if block == 1
        imwrite(imind,cm,[resultsFolder,'plots/',ptname,'/',filename],'gif', 'Loopcount',inf,'DelayTime',delay);
    else
        imwrite(imind,cm,[resultsFolder,'plots/',ptname,'/',filename],'gif','WriteMode','append','DelayTime',delay);
    end

    close(fig)

end



end