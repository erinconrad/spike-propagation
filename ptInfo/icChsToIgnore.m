function [elecs,tmul,absthresh] = icChsToIgnore(name)

if strcmp(name,'HUP068') ==1
    elecs = {};
    tmul = 5;
    absthresh = 300;
elseif strcmp(name,'HUP070') == 1
    elecs = {'LST2','LST1','LFP6','LFP1'};
elseif strcmp(name,'HUP078') == 1
    elecs = {'AMY1','AMY2'};
    tmul = 4;
    absthresh = 300;
elseif strcmp(name,'HUP086') == 1
    elecs = {'LMPI6','LMPI5','LMPI4','LMPI3',...
        'LMPI2','LMPI1','LAST1','LAST2','LAST3',...
        'LAST4'};
elseif strcmp(name,'HUP107') == 1
    elecs = {};
elseif strcmp(name,'HUP111A') == 1
    elecs = {};
elseif strcmp(name,'HUP116') == 1
    elecs = {};
elseif strcmp(name,'Study022') == 1
    elecs = {};
elseif strcmp(name,'Study028') == 1
    elecs = {};
end

end