function pt = updateClinical(pt)

clinicalFile = 'Clinical.xls';
T = readtable(clinicalFile);

whichPts = 1:46;
%[1,4,6,7,8,9,12,14,15,16,17,18,19,20,22,24,25,27,30,31];

for whichPt = whichPts
    name = pt(whichPt).name;
    found = 0;
    for i = 1:size(T,1)
        if strcmp(name,T.name(i)) == 1
            found = 1;
            
            % Update clinical info
            pt(whichPt).clinical.outcome = T.outcome{i};
            pt(whichPt).clinical.pathology = T.path{i};
            pt(whichPt).clinical.ageOnset = T.ageOnset{i};
            pt(whichPt).clinical.ageSurgery = T.ageSurg{i};
            pt(whichPt).clinical.sex = T.sex{i};
            pt(whichPt).clinical.seizureOnset = T.soz{i};
                
            
        end
        
    end
    if found == 0
        error('Could not find %s\n',name);
    end
    
end




end