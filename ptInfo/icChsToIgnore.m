function [elecs,tmul,absthresh] = icChsToIgnore(name)

tmul = [];
absthresh = [];

if strcmp(name,'HUP064') == 1
    elecs = {};
    tmul = 9; % try lowering tmul
    absthresh = 300;
elseif strcmp(name,'HUP065') ==1
    elecs = {'RPO3','RPO4','RPO5','RPO6','RST1',...
        'RST2','RST3','RST4','RST5','RST6'};
    tmul = 9;
    absthresh = 300;
elseif strcmp(name,'HUP068') ==1
    elecs = {};
    tmul = 5;
    absthresh = 300;
elseif strcmp(name,'HUP070') == 1
    elecs = {'LST2','LST1','LFP6','LFP1'};
    % SKIP HUP070, too short
elseif strcmp(name,'HUP073') == 1
    elecs = {};
    tmul = [];
    absthresh = [];
elseif strcmp(name,'HUP074') == 1
    elecs = {'HD2','FOP7','FOP8'};
    tmul = 6;
    absthresh = 300;
elseif strcmp(name,'HUP075') == 1
    elecs = {};
    tmul = [];
    absthresh = [];
elseif strcmp(name,'HUP078') == 1
    elecs = {'AMY1','AMY2'};
    tmul = 4;
    absthresh = 300;
elseif strcmp(name,'HUP080') == 1
    elecs = {'PST2','PST1','MST2','MST1'};
    tmul = 9;
    absthresh = 300;
elseif strcmp(name,'HUP082') == 1
    elecs = {};
    tmul = [];
    absthresh = [];
elseif strcmp(name,'HUP086') == 1
    elecs = {'LMPI6','LMPI5','LMPI4','LMPI3',...
        'LMPI2','LMPI1','LAST1','LAST2','LAST3',...
        'LAST4','LG41','LG48','RSPI6'};
    tmul = 9;
    absthresh = 300;
    %'LG33','LG34',};
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