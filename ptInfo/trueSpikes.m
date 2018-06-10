%{

1) Tiny positive wave at start
2) Huge negative deflection
3) asymmetric rising
3) drops below the baseline
5) duration < 200 ms (spike < 70 ms)
6) slow wave following the discharge
(7) there should be a field)


%}


function validated = trueSpikes(name)

if strcmp(name,'HUP064') == 1
    validated.times = [10801.11; 10861.95; 10880.03; 10881.60;42138.89;...
         42074.98; 42076.13; 42120.25;42208.89; 42347.86;...
          42626.85; 42647.54;42847.97; 43025.92;86426.04;...
         86841.83;];
    validated.chs = {'LG2','LG3','LG4','LG10','LG17','LG25','LG34','LG52','LG54','LG59','LG3','LIH5','RFR3','LST1','LST2','LST3'};
elseif strcmp(name, 'HUP065') == 1
    validated = [];
elseif strcmp(name,'HUP068') ==1
    validated = [];
elseif strcmp(name,'HUP070') ==1
    validated = [];
elseif strcmp(name,'HUP073') == 1
    validated = [];
elseif strcmp(name,'HUP074') ==1
    validated = [];
elseif strcmp(name,'HUP075')==1
    validated = [];
elseif strcmp(name,'HUP078') == 1
    validated = [];
elseif strcmp(name,'HUP080') == 1
    validated = [];
elseif strcmp(name,'HUP082') == 1
    validated = [];
elseif strcmp(name,'HUP083') == 1
    validated = [];
elseif strcmp(name,'HUP086') == 1
    validated = [];
elseif strcmp(name,'HUP087') == 1
    validated = [];
elseif strcmp(name,'HUP088') == 1
    validated = [];
elseif strcmp(name,'HUP094') == 1
    validated = [];
elseif strcmp(name,'HUP105') == 1
    validated = [];
elseif strcmp(name,'HUP106') == 1
    validated = [];
elseif strcmp(name,'HUP107') == 1
    validated = [];
elseif strcmp(name,'HUP111') == 1
    validated = [];
elseif strcmp(name,'Study004') == 1
    validated = [];
elseif strcmp(name,'Study012') == 1
    validated = [];
elseif strcmp(name,'Study016') == 1
    validated = [];
elseif strcmp(name,'Study017') == 1
    validated = [];
elseif strcmp(name,'Study019') == 1
    validated = [];
elseif strcmp(name,'Study020') == 1
    validated = [];
elseif strcmp(name,'Study021') == 1
    validated = [];
elseif strcmp(name,'Study022') == 1
    validated = [];
elseif strcmp(name,'Study023') == 1
    validated = [];
elseif strcmp(name,'Study026') == 1
    validated = [];
elseif strcmp(name,'Study028') == 1
    validated = [];
elseif strcmp(name,'Study029') == 1
    validated = [];
elseif strcmp(name,'Study033') == 1
    validated = [];
elseif strcmp(name,'Study038') == 1
    validated = [];
elseif strcmp(name,'CHOP08') == 1
    validated = [];
elseif strcmp(name,'CHOP14') == 1
    validated = [];
elseif strcmp(name,'CHOP20') == 1
    validated = [];
elseif strcmp(name,'CHOP27') == 1
    validated = [];
elseif strcmp(name,'CHOP38') == 1
    validated = [];
elseif strcmp(name,'CHOP40') == 1
    validated = [];
elseif strcmp(name,'CHOP42') == 1
    validated = [];
elseif strcmp(name,'CHOP43') == 1
    validated = [];
elseif strcmp(name,'CHOP45') == 1
    validated = [];
elseif strcmp(name,'CHOP46') == 1
    validated = [];
elseif strcmp(name,'CHOP47') == 1
    validated = [];
end


end