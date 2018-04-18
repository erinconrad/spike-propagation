function spikeVerification(P)

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;


thresholds_to_try = [13 10 8];
times_to_try = [1 2 3];
soz = [0 1];

for i = 1:length(P)
    outputFolder = [resultsFolder,'spike verification/',P(i).name,'/'];

    if isempty(dir([outputFolder,'*.png'])) == 0
        fprintf('Patient %s presumably done, skipping\n',P(i).name);
        continue
    end
    
    if isempty(P(i).ieeg_name) == 1 || isempty(P(i).electrode_labels) == 1
        fprintf('Missing ieeg_name or electrode labels for patient %s, skipping\n',P(i).name);
        continue
    end
    
    for s = soz
    
        nchs = size(P(i).electrodeData.locs,1);
        chIds = randsample(nchs,10);
    for j = times_to_try
       
        for k = thresholds_to_try
            
            
            portVisualizeSpikes(P,i,1,j,s,chIds,k);
            
            
            
        end
        
    end
    
    end
    
    
end



end