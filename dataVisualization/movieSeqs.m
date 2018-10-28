function movieSeqs(seqs,locs,info)

%% Parameters
circleSize = 60;
delay = 0.5; % time delay between steps
x = locs(:,1); y = locs(:,2); z = locs(:,3);
nseqs = size(seqs,2);
prows = 2;
columns =  nseqs/prows;

%% Get seqs in easy plotting order
curr_max = 0;
for i = 1:nseqs
   curr_seq = seqs(:,i); 
   ifull = zeros(size(curr_seq));
   non_nans = ~isnan(curr_seq);
   [~,~,idx] = unique(curr_seq(non_nans));
   ifull(non_nans) = idx;
   plot_order{i} = ifull;
   temp_max = max(ifull);
   if temp_max > curr_max, curr_max = temp_max; end
end

% Number of frames
n_frames = curr_max;


%% Plot base
fig = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.4, 0.95, 0.8]);

suptitle(sprintf('%d most representative sequences for cluster %d, %s',...
    nseqs,info.cluster,info.name));
set(gca,'FontSize',15);

for s = 1:nseqs
    whichcol = ceil(s/2);
    whichrow = mod(s+1,2)+1;
    subplot('Position',[(whichcol-1)*1/columns (prows-whichrow)*1/prows 1/columns 1/prows])
    scatter3(x,y,z,circleSize,'markeredgecolor','k');
    hold on
    grid off
    axis off
    
end




makeGif(fig,1,delay,info.outputFile)


%% Plot frames
for i = 1:n_frames

for s = 1:nseqs
    whichcol = ceil(s/2);
    whichrow = mod(s+1,2)+1;
    subplot('Position',[(whichcol-1)*1/columns (prows-whichrow)*1/prows 1/columns 1/prows]);
    hold on
    
    seq_plot_order = plot_order{s};
    if i < max(seq_plot_order)
        chs = find(seq_plot_order==i);
        scatter3(locs(chs,1),locs(chs,2),locs(chs,3),circleSize,'r','filled');
       
        
    end
    
    grid off
    axis off
    
    
    
end
makeGif(fig,i+1,delay,info.outputFile)
end

close(fig)
end