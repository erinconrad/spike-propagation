function sa = areaConnecting(locs,doplot)

shp = alphaShape(locs);
sa = surfaceArea(shp);

if doplot == 1
    figure
    scatter3(locs(:,1),locs(:,2),locs(:,3),100);
    hold on   
    cv = plot(shp);
    alpha(cv,0.05) 
    pause
    close(gcf)
    
end

end