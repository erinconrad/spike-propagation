function plotDelays(delay,power,chLocs)
chLocs = chLocs(:,2:4);
dotsize = 70;


x = chLocs(:,1);
y = chLocs(:,2);
z = chLocs(:,3);


figure
subplot(2,1,1)
scatter3(x,y,z,dotsize,delay,'filled')

colormap jet
colorbar

subplot(2,1,2)
scatter3(x,y,z,dotsize,power,'filled')


%caxis([0 8])
colormap jet
colorbar

end