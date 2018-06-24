function gdf = removeEKGArtifact(gdf,gdfEKG,prox)

tooClose = zeros(size(gdf,1),1);

% Loop through spikes
for i = 1:size(gdf,1)
    
    % Get spike time
    time = gdf(i,2);
    
    % Loop through EKG spikes
    for j = 1:size(gdfEKG,1)
        
        % if the spike occurs to close to the EKG spike
       if abs(gdf(i,2)-gdfEKG(j,2)) < prox 
           
           % Mark it to be discarded
           tooClose(i) = 1;
           
       end
        
    end


end

% Remove spikes too close to EKG spikes
gdf(tooClose==1,:) = [];