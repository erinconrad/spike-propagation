function sa = areaConnecting(locs,doplot,all_locs)

shp = alphaShape(locs);
sa = surfaceArea(shp);

if doplot == 1
    figure
    scatter3(all_locs(:,1),all_locs(:,2),all_locs(:,3),100);
    hold on   
    scatter3(locs(:,1),locs(:,2),locs(:,3),100,'filled');
    
    cv = plot(shp);
    alpha(cv,0.05) 
    title(sprintf('Surface area: %1.2f',sa))
    pause
    close(gcf)
    
end

end