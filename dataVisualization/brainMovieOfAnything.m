function brainMovieOfAnything(thing,chLocs,info)

% Parameters
circleSize = 50;
delay = 0.4; % time delay between steps

colormin = min(min(thing));
colormax = max(max(thing));

nsteps = size(thing,1);
nchs = size(thing,2);

for iTime = 1:nsteps
   fig = figure;
   
    scatter3(chLocs(:,1),chLocs(:,2),chLocs(:,3),...
        circleSize,'markeredgecolor','k');
    
    hold on
     
   
    scatter3(chLocs(:,1),chLocs(:,2),chLocs(:,3),...
        circleSize,thing(iTime,:),'filled')

    grid off
    axis off
    title(info.title{iTime});
    set(gca,'FontSize',15);
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    set(gca,'ZTickLabel',[]);
    %caxis([colormin,colormax])
    colormap(flipud(gray))
    
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