%{ 
visualizeChLocations

This code takes the Patient structure, which you can get by running "main"
in "newSpatialOrgSubscripts" and allows you to visualize a desired spike
sequence as a gif. It assumes equal time steps between each successive
spike (which is not actually true).

%}

function portVisualizeChLocations(P,pt,sz,s)

ptname = P(pt).name;
nseqs = size(P(pt).sz(sz).data.sequences,2)/2;

if isempty(s) == 1
    s = randsample(1:nseqs,1);
end

% output file name
filename = [sprintf('%s',ptname),'_sz_',sprintf('%d',sz),...
    '_sequence_',sprintf('%d',s),'.gif'];

[~,~,~,resultsFolder,~] = fileLocations;

% Get the channel locations for the patient
chLocs = P(pt).sz(sz).data.xyChan(:,2:4);

% Parameters
circleSize = 20;
delay = 0.4; % time delay between steps

temp = P(pt).sz(sz).data.sequences; % which segment and starting channel

seq = temp(:,(s-1)*2+1:(s-1)*2+2); % which sequence
seq = seq(any(seq~=0,2),:); % remove rows of zeros



% loop through all the spikes in the sequence
for iTime = 1:size(seq,1)
    
    
    fig = figure;
    
    % make all channels a base color
    scatter3(chLocs(:,1),chLocs(:,2),chLocs(:,3),...
        circleSize,'b','filled')
    hold on
    
    
    % loop through all the spikes up to the current spike (want the current
    % spike and all previous spikes to change color to highlight the path)
    for i = 1:iTime
       scatter3(chLocs(seq(i,1),1),chLocs(seq(i,1),2),chLocs(seq(i,1),3),...
           circleSize,'r','filled')
       % seq(i,1) is the location of the ith spike in the sequence. I am
       % making all the channels that have been activated so far turn red
       
        if i ~= 1
           plot3([chLocs(seq(i,1),1),chLocs(seq(i-1,1),1)],[chLocs(seq(i,1),2),chLocs(seq(i-1,1),2)],...
               [chLocs(seq(i,1),3),chLocs(seq(i-1,1),3)],'r');
           % I am also connecting the dots to more clearly highlight the
           % path
        end
    end
    
    grid off
    axis off
    
    % capture the figure as a frame in the gif
    F(iTime) = getframe(fig);
    im = frame2im(F(iTime));
    [imind,cm] = rgb2ind(im,256);
    
    if iTime == 1
        imwrite(imind,cm,[resultsFolder,'plots/',ptname,'/',filename],'gif', 'Loopcount',inf,'DelayTime',delay);
    else
        imwrite(imind,cm,[resultsFolder,'plots/',ptname,'/',filename],'gif','WriteMode','append','DelayTime',delay);
    end

    close(fig)
end

end