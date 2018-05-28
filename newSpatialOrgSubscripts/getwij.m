function wij = getwij(xyChan,dmin)

nChannels = size(xyChan,1);

% Calculate dij, distances between channels and weights wij
dij = nan(nChannels,nChannels);
wij = nan(nChannels,nChannels);
for iChannel = 1:nChannels
   for jChannel = 1:nChannels
       if iChannel == jChannel
           dij(iChannel,jChannel) = nan;
           wij(iChannel,jChannel) = nan;
       else
           % calculate Euclidean distance between the two channels
           dij(iChannel,jChannel) =  sqrt(sum((xyChan(iChannel,2:end) - xyChan(jChannel,2:end)).^2));
           
           % if this distance is less than a minimum distance
           if dij(iChannel,jChannel) <= dmin
               
               % make the weight be 1/d
               wij(iChannel,jChannel) = 1/dij(iChannel,jChannel);
               
               if wij(iChannel,jChannel) == inf
                   wij(iChannel,jChannel) = 0;
                   fprintf('Warning, channels %d and %d have the same location, making spatial weight 0\n',iChannel,jChannel);
               end
           else
               
               % if they're further, set the weight to 0
               wij(iChannel,jChannel) = 0;
           end
       end
      
   end
    
end




end