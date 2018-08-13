%% Load electrode data file
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
elecStructFolder = [resultsFolder,'ptStructs/allPtsElectrodeData/'];
structWithElecFile = [elecStructFolder,'ptWithElectrodeData.mat'];
load(structWithElecFile)


%% Loop through all patients
for i = 1:length(pt)

    if isfield(pt(i).electrodeData,'ref_vector') == 1
        fprintf('Already did this for patient %s, skipping...\n',...
            pt(i).name);
        continue
    end
    chLocs = pt(i).electrodeData.locs(:,2:4);
    fig = figure;

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
    disp('Select the first point in the reference vector, then press enter');
    pause
    c_info = getCursorInfo(dcm_obj);
    start = c_info.Position;
    scatter3(start(1),start(2),start(3),30,'g','filled');

    % Ask them to select the 2nd point
    set(dcm_obj,'DisplayStyle','datatip',...
        'SnapToDataVertex','off','Enable','on')
    disp('Select the second point in the reference vector, then press enter');
    pause
    c_info = getCursorInfo(dcm_obj);
    finish = c_info.Position;
    scatter3(finish(1),finish(2),finish(3),30,'r','filled');

    % Show the line between them
    plot3([start(1),finish(1)],[start(2),finish(2)],[start(3),finish(3)],'k','LineWidth',2);

    % Save this vector in the electrode data
    pt(i).electrodeData.ref_vector = [start;finish];

    disp('Press enter to move to the next patient');
    pause
    close(fig)
end

%newFile = 'ptWithElectodeDataRefVector.mat';
save(structWithElecFile);
