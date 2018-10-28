function brainImageOfAnything(thing,chLocs,info)

colormin = min(min(thing));
colormax = max(max(thing));


x = chLocs(:,1);
y = chLocs(:,2);
z = chLocs(:,3);
minxy = min(min(chLocs(:,1)),min(chLocs(:,2)));
maxxy = max(max(chLocs(:,1)),max(chLocs(:,2)));
rangexy = maxxy-minxy;
[xq,yq] = meshgrid(minxy:0.05*rangexy:maxxy);
z1 = griddata(x,y,z,xq,yq,'natural');
zcol = griddata(x,y,thing,xq,yq,'natural');
fig = figure;

%surf(xq,yq,z1,zcol)
scatter3(x,y,z,60,'k');
hold on
scatter3(x,y,z,60,thing,'filled');

end