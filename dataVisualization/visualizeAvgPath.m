%{ 
visualizeAvgPath

This makes a graph where the size of the electrode corresponds to how many
spikes went through there, and the width of the line connecting two
electrodes corresponds to how often that path was traversed.


%}

function visualizeAvgPath(P,sz,block)

%% Parameters

Patient = P(80).sz(sz).block(block).data;

% output file name
[~,~,~,resultsFolder,~] = fileLocations;
filename = ['HUP080_avgsequences','.png'];

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
chLocs = Patient.xyChan(:,2:4);
allChs = Patient.xyChan(:,1);
chHits = zeros(size(allChs));
pathHits = zeros(length(allChs),length(allChs));


%% Get the sequence data
% Change these lines each time. Which segment and which starting channel
% and which spike sequence to visualize
temp = Patient.cleanseq;%sequences; % which segment and starting channel

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

% make all channels a base color
scatter3(chLocs(:,1),chLocs(:,2),chLocs(:,3),...
    baseSizeElec,'b','filled')
hold on

% Loop through the channels
for i = 1:length(allChs)
    
    % if the channel has a hit
   if chHits(i) ~= 0
       
       % weight the dot by how many hits
      %circleSize = baseSizeElec + (chHits(i)-1)*multElec;
      
      % plot the dot
      %scatter3(chLocs(i,1) ,chLocs(i,2),chLocs(i,3),circleSize,'r','filled');
      
     
   end
   
   % loop through the other channels, looking for path hits
   for j = 1:length(allChs)
       
      if pathHits(i,j) ~= 0
         
          % weight the line by the number of path hits
          linesize = (pathHits(i,j)-1)*multLines+baseSizeLine;
          
          % plot the line
          plot3([chLocs(i,1),chLocs(j,1)],...
              [chLocs(i,2),chLocs(j,2)],...
               [chLocs(i,3),chLocs(j,3)],'r','lineWidth',linesize);
      end
   end
end


saveas(fig,[resultsFolder,filename]);

end
