function newgdf = removeEKGArtifact(gdf,gdfEKG,prox)


newgdf = [];

for i = 1:size(gdf,1)
    time = gdf(i,2);
    tooClose = 0;
    for j = 1:size(gdfEKG,1)
       if abs(gdf(i,2)-gdfEKG(j,2)) < prox 
        tooClose = 1;
           
       end
        
    end
    if tooClose == 0
        newgdf = [newgdf;gdf(i,:)];
    end


end