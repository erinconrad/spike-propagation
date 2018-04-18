function showElectrodeLocations(P)

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
ptInfo = loadjson(jsonfile);



% base size of electrodes
baseSizeElec = 20;

for i = 1:length(P)
    
    name = P(i).name;
    
    filename = [name,'_locations.png'];
    
    if exist([electrodeFolder,filename],'file') ~= 0
        fprintf('File %s already found, skipping\n',filename);
        continue
    end

    
   [ieeg_name,electrode_name] =  ieegAndElectodeNames(name);
   electrodeFile = [electrodeFolder,electrode_name];
   
   if isempty(ieeg_name) == 1 || isempty(electrode_name) == 1
       continue
   end
   
   data = getiEEGData(ieeg_name,0,0,pwfile);  
   dummyRun = 1;
   outputData = 0;
   [~,electrodeData,~] = getSpikeTimes(0,name,ieeg_name,electrodeFile,ptInfo,pwfile,...
    dummyRun,0,0,outputData,0,1,0);
   
   chLocs = electrodeData.locs(:,2:4);
   fig = figure;

    % make all channels a base color
    scatter3(chLocs(:,1),chLocs(:,2),chLocs(:,3),...
        baseSizeElec,'b','filled')
    grid off
    
    title([name,' electrode locations']);
    set(gca,'FontSize',15);
    
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    set(gca,'ZTickLabel',[]);
    
    saveas(fig,[electrodeFolder,filename]);
    close(fig)
    
end