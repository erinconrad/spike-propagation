function movieChangeVectors(early,vecs,chLocs,info)

% Parameters
circleSize = 60;
delay = 0.5; % time delay between steps


nsteps = size(early,1);

for iTime = 1:nsteps
   fig = figure;
   
   
    x = chLocs(:,1);
    y = chLocs(:,2);
    z = chLocs(:,3);
    minxy = min(min(chLocs(:,1)),min(chLocs(:,2)));
    maxxy = max(max(chLocs(:,1)),max(chLocs(:,2)));
    rangexy = maxxy-minxy;
    [xq,yq] = meshgrid(minxy:0.05*rangexy:maxxy);
    
   
    
    scatter3(chLocs(:,1),chLocs(:,2),chLocs(:,3),...
        circleSize,'markeredgecolor','k');
    
    hold on
     
   
    scatter3(early(iTime,1),early(iTime,2),early(iTime,3),...
        circleSize,'g','filled');
    scatter3(early(iTime,1)+vecs(iTime,1),early(iTime,2)+vecs(iTime,2),...
        early(iTime,3)+vecs(iTime,3),...
        circleSize,'r','filled');
    plot3([early(iTime,1) early(iTime,1)+vecs(iTime,1)],...
        [early(iTime,2) early(iTime,2)+vecs(iTime,2)],...
        [early(iTime,3) early(iTime,3)+vecs(iTime,3)],...
        'k','LineWidth',2);

    grid off
    axis off
    title(info.title{iTime});
    set(gca,'FontSize',15);
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    set(gca,'ZTickLabel',[]);
    %caxis([colormin,colormax])
    %colormap(flipud(gray))
    
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.4, 0.4, 0.3, 0.45]);
    
    % capture the figure as a frame in the gif
    F(iTime) = getframe(fig);
    im = frame2im(F(iTime));
    [imind,cm] = rgb2ind(im,256);
    
    if iTime == 1
        imwrite(imind,cm,info.save,'gif', 'Loopcount',inf,'DelayTime',delay);
    else
        imwrite(imind,cm,info.save,'gif','WriteMode','append','DelayTime',delay);
    end

    close(fig)
   
    
end

end