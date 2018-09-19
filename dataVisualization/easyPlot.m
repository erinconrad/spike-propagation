function easyPlot(thing,xyChan)

dotsize = 100;
chLocs = xyChan(:,2:4);


figure
scatter3(chLocs(:,1),chLocs(:,2),chLocs(:,3),dotsize,thing,'filled','MarkerEdgeColor','k')

end