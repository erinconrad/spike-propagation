function AEDs = getAEDs(name)

AEDs = {};

if strcmp(name,'HUP068') ==1
    AEDs = {'CBZ','ZNS','VPA'};
elseif strcmp(name,'HUP070') ==1
    AEDs = {'LTG','OXC'};
elseif strcmp(name,'HUP078') ==1
    AEDs = {'PHT','LEV','LTG'};
elseif strcmp(name,'HUP080') ==1
    AEDs = {'OXC'};
elseif strcmp(name,'HUP086') ==1
    AEDs = {'OXC','TPM'};
elseif strcmp(name,'HUP106') ==1
    AEDs = {'LTG','ZNS'};
elseif strcmp(name,'HUP107') ==1
    AEDs = {'LTG'};
elseif strcmp(name,'HUP111A') ==1
    AEDs = {'LCS','ZNS'};
elseif strcmp(name,'HUP116') ==1
    AEDs = {'OXC','TPM','ONF'};
elseif strcmp(name,'Study016') ==1
    AEDs = {};
elseif strcmp(name,'Study019') ==1
    AEDs = {};
elseif strcmp(name,'Study020') ==1
    AEDs = {};
elseif strcmp(name,'Study022') ==1
    AEDs = {};
elseif strcmp(name,'Study028') ==1
    AEDs = {};
elseif strcmp(name,'Study029') ==1
    AEDs = {};
end

end

