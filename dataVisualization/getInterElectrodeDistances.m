function getInterElectrodeDistances(pt,i)



chLocs = pt(i).electrodeData.locs(:,2:4);
fig = figure;

while 1
    % plot channel locations
    scatter3(chLocs(:,1),chLocs(:,2),chLocs(:,3),...
        30,'b','filled')
    hold on
    grid off

    title(sprintf('%s',pt(i).name));
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    set(gca,'ZTickLabel',[]);

    % Ask the user to select a point for the start of the reference vector
    dcm_obj = datacursormode(fig);
    set(dcm_obj,'DisplayStyle','datatip',...
        'SnapToDataVertex','off','Enable','on')
    disp('Select the first channel, then press enter');
    pause
    c_info = getCursorInfo(dcm_obj);
    start = c_info.Position;
    scatter3(start(1),start(2),start(3),30,'g','filled');

    % Ask them to select the 2nd point
    set(dcm_obj,'DisplayStyle','datatip',...
        'SnapToDataVertex','off','Enable','on')
    disp('Select the second channel, then press enter');
    pause
    c_info = getCursorInfo(dcm_obj);
    finish = c_info.Position;
    scatter3(finish(1),finish(2),finish(3),30,'r','filled');

    % Plot the distance
    fprintf('The distance between the electrodes is %1.1f mm.\n',sqrt(sum((start-finish).^2)));

    % Save this vector in the electrode data
    pt(i).electrodeData.ref_vector = [start;finish];

    disp('Press enter to select additional channel pairs');
    pause
    hold off

end


end

