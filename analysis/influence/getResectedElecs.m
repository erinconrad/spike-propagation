function pt = getResectedElecs(pt,cluster,whichPts)

%{
Add the identities of the electrodes overlying resected cortex to the
patient structure
%}

[electrodeFolder,~,scriptFolder,resultsFolder,~] = fileLocations;

if isempty(whichPts) == 1
    for i = 1:length(pt)
        if isempty(pt(i).seq_matrix) == 0
            if size(cluster(i).bad_cluster) < cluster(i).k
                whichPts = [whichPts,i];
            end
        end
    end
end

if isequal(whichPts,[1,4,6,7,8,9,12,14,15,16,17,18,19,20,22,24,25,27,30,31]) == 0
   % error('Warning, not doing correct patients!\n');
end

resecFolder = [electrodeFolder,'resectedElecs/for_erin/'];


%{
for whichPt = whichPts
    name = pt(whichPt).name;
    
    % Find appropriate file in resected folder
    listing = dir([resecFolder,name,'*']);
    
    if length(listing) ~= 1
        error('what\n');
    end
    
    % Open file
    T = readtable([resecFolder,listing.name],'Delimiter',',',...
        'ReadVariableNames',false);
    
    new_elecs = [];
    new_labels = {};
    
    for i = 1:size(T,1)
        
        % Get the label and add it to my cell array
        elec_label = T.Var2(i);
        new_labels = [new_labels;elec_label];
        
        found_it = 0;
        
        % Find the corresponding electrode number in my electrode field
        for j = 1:length(pt(whichPt).electrodeData.electrodes)
            
            
            if strcmp(elec_label,pt(whichPt).electrodeData.electrodes(j).name) == 1
                % Add the electrode number to my array
                new_elecs = [new_elecs;j];
                found_it = 1;
                
                % break out of inner loop
                break
                
            end

        end
        
        if found_it == 0
            error('Did not find electrode!\n');
        end
        
    end
    
    
end


%}

for whichPt = whichPts
    name = pt(whichPt).name;
    
    % Find appropriate file in resected folder
    listing = dir([resecFolder,name,'*']);
    
    if length(listing) ~= 1
        error('what\n');
    end
    
    % Open file
    T = readtable([resecFolder,listing.name],'Delimiter',',',...
        'ReadVariableNames',false);
    
    % Check that nums agree with labels
    bad_names = 0;
    for i = 1:size(T,1)
        elec_num = T.Var1(i);
        elec_label = T.Var2(i);
        
        % This disagreement happens sometimes because I define my electrode
        % numbers based on the order in ieeg not the order in the electrode
        % files
        
        if elec_num == 0
            bad_names = 1;
            break
        end
        
        if strcmp(elec_label,...
                pt(whichPt).electrodeData.electrodes(elec_num).name) == 0
            fprintf('Warning, disagreement between electrode labels for %s\n',...
                name);
            bad_names = 1;
            break
            
        end

    end
    
    if bad_names == 0
    
        % Add resected elecs to structure
        pt(whichPt).resecElecs = T.Var1;
        pt(whichPt).resecLabels = T.Var2;
    else
        % Need to find the right ones
        new_elecs = [];
        new_labels = {};
        for i = 1:size(T,1)
            elec_num = T.Var1(i);
            elec_label = T.Var2(i);
            
            for j = 1:length(pt(whichPt).electrodeData.electrodes)
                if strcmp(pt(whichPt).electrodeData.electrodes(j).name,...
                        elec_label) == 1
                    new_elecs = [new_elecs;j];
                    new_labels = [new_labels;elec_label];
                    
                    break
                end
                
                % If I never found the electrode, throw an error
                if j == length(pt(whichPt).electrodeData.electrodes)
                    error('What\n');
                end
                
            end
            
        end
        
        if length(new_elecs) ~= size(T,1)
            error('what\n');
        end
        
        pt(whichPt).resecElecs = new_elecs;
        pt(whichPt).resecLabels = new_labels;
        
    end
        
    
end

%}


end