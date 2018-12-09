function outcome = getOutcome(name)


%{

http://seizure.mgh.harvard.edu/engel-surgical-outcome-scale/
%}

[~,jsonfile,scriptFolder,~,~] = fileLocations;
p1 = genpath(scriptFolder);
addpath(p1);
ptInfo = loadjson(jsonfile);
ptnames = fieldnames(ptInfo.PATIENTS);

% Do trickery
ptnames(strcmp(ptnames,'HUP111B')==1) = [];
ptnames(strcmp(ptnames,'Study004')==1) = [];

% Get the outcome from the json file
for i = 1:length(ptnames)
    if strcmp(name,ptnames{i}) == 1
        info = ptInfo.PATIENTS.(ptnames{i});
        outcome = (info.Outcome);
    end
end

%% Rename the outcome

if strcmp(outcome(1),'I') == 1
    % They are using ILAE definitions :(
    if strcmp(outcome,'ILAE1') == 1
        % completeley seizure free (I am calling this IA)
        new_outcome = 1;
    elseif strcmp(outcome,'ILAE2') == 1
        % only auras (I am calling this IB)
        new_outcome = 1.25;
    elseif strcmp(outcome,'ILAE3') == 1
        % 1-3 sz days per year (I am calling this IIB)
        new_outcome = 2.25;
    elseif strcmp(outcome,'ILAE4') == 1
        % 4 sz days per year to 50% reduction baseline (I am calling this
        % IIIA)
        new_outcome = 3;
    elseif strcmp(outcome,'ILAE5') == 1
        % Less than 50% reduction of sz days (I am calling this IVA)
        new_outcome = 4;
    elseif strcmp(outcome,'ILAE6') == 1
        % More than 100% increase of sz days (I am calling this IVC)
        new_outcome = 4.5;
    end
    
    old_outcome = outcome;
    outcome = new_outcome;
        
    
else
    % They are using Engel definitions
    
    outcome = str2num(outcome);

    % Get mantissa
    mantissa = uint16((outcome - floor(outcome))*10);

    % Re-scale mantissa 
    if mantissa == 1
        new_mantissa = 0;
    elseif mantissa == 2
        new_mantissa = 2.5;
    elseif mantissa == 3
        new_mantissa = 5;
    elseif mantissa == 4
        new_mantissa = 7.5;
    end

    old_outcome = outcome;
    outcome = floor(outcome) + new_mantissa/10;

end



end