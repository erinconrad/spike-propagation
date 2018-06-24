function plotSpikeTimes(data,spikes,indices,chs)

% This function plots EEG data and detected spike times for desired indices
% and desired channels

% Establish colors
colors = {'b','r','g','k','m'};

% Loop through desired channels
for x = 1:length(chs)
    iCh = chs(x);
    
    % which color to plot
    whichColor = mod(x,length(colors));
    if whichColor == 0
        whichColor = 5; 
    end
    whichColor = colors{whichColor};
    
    % Make a y axis offset
    yOffset = (x-1)*300+300;
    
    % Get the time in seconds of the first desired index, relative to the
    % start of the data
    startTime = data.times(indices(1)) - data.times(1);
    
    % Get the time in seconds of the last desired index, relative to the
    % start of the data
    endTime = data.times(indices(end)) - data.times(1);
    
    % plot the data, where x axis is time in seconds relative to the start
    % of the data
    plot(data.times(indices)-ones(1,length(indices))*data.times(1),...
        data.values(indices,iCh)+yOffset*ones(1,length(indices)),whichColor);
    hold on
    
    % find the spikes that occurred in the desired channel and occurred before
    % the desired end time
    correctCHSpikesEarly = intersect(find(spikes.chan==iCh),find(spikes.pos < endTime));
    
    % Also find those that occurred after the desired start time
    whichSpikes = intersect(correctCHSpikesEarly,find(spikes.pos > startTime));
    
    
    % plot the spike times on top of the original plot
    scatter(spikes.pos(whichSpikes),(yOffset+100)*(ones(length(whichSpikes),1)),...
        whichColor,'filled')
end



end