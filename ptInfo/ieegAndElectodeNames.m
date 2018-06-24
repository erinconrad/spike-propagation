%% ieegAndElectrodeNames

% This (hand written) file takes the json-format patient name and outputs
% the corresponding ieeg.org data name and the electrode file location


function [ieeg_name,electrode_name,tmul,absthresh,ignore] =  ieegAndElectodeNames(name)

tmul = 13;
absthresh = 300;

if strcmp(name,'HUP064') == 1
    ieeg_name = 'HUP64_phaseII-Annotations';
    electrode_name = 'HUP064_T1_19991230_electrode_labels.csv';
    ignore = {"EKG1","EKG2","RR","Rate","X1","X2","X3","X4"};
elseif strcmp(name, 'HUP065') == 1
    ieeg_name = 'HUP65_phaseII-Annotations';
    electrode_name = 'HUP065_T1_19991230_electrode_labels.csv';
    ignore = {"X1","X2","EKG1","EKG2"};
elseif strcmp(name,'HUP068') ==1
    ieeg_name = [];
    electrode_name = 'HUP068_T1_19991230_electrode_labels.csv';
    ignore = {"EKG1","EKG2","X1","X2"};
elseif strcmp(name,'HUP070') ==1
    ieeg_name = [];
    electrode_name = 'HUP070_T1_19980321_electrode_labels.csv';
    ignore = {"EKG1","EKG2","Rate","RR"};
elseif strcmp(name,'HUP073') == 1
    ieeg_name = 'HUP73_phaseII-Annotations';
    electrode_name = 'HUP073_T1_19990708_electrode_labels.csv';
    ignore = {"EKG1","EKG2"};
elseif strcmp(name,'HUP074') ==1
    ieeg_name = 'HUP74_phaseII';
    electrode_name = 'HUP074_T1_19990327_electrode_labels.csv';
    ignore = {"EKG1","EKG2"};
elseif strcmp(name,'HUP075')==1
    ieeg_name = 'HUP75_phaseII';
    electrode_name = 'HUP075_T1_19990604_electrode_labels.csv';
    ignore = {"EKG1","EKG2","RR","Rate"};
elseif strcmp(name,'HUP078') == 1
    ieeg_name = 'HUP78_phaseII-Annotations';
    electrode_name = 'HUP078_T1_19971218_electrode_labels.csv';
    ignore = {"EKG1","EKG2"};
elseif strcmp(name,'HUP080') == 1
    ieeg_name = 'HUP80_phaseII';
    electrode_name = 'HUP080_T1_19991213_electrode_labels.csv';
    ignore = {"EKG1","EKG2"};
elseif strcmp(name,'HUP082') == 1
    ieeg_name = 'HUP82_phaseII';
    electrode_name = 'HUP082_T1_19991212_electrode_labels.csv';
    ignore = {"EKG1","EKG2"};
elseif strcmp(name,'HUP083') == 1
    ieeg_name = 'HUP83_phaseII';
    electrode_name = 'HUP083_T1_19991217_1_electrode_labels.csv';
    ignore = {"EKG1","EKG2"};
elseif strcmp(name,'HUP086') == 1
    ieeg_name = [];
    electrode_name = 'HUP086_T1_19990728_electrode_labels.csv';
    ignore = {"EKG1","EKG2"};
elseif strcmp(name,'HUP087') == 1
    ieeg_name = [];
    electrode_name = 'HUP087_T1_19991226_electrode_labels.csv';
    ignore = {"EKG1","EKG2"};
elseif strcmp(name,'HUP088') == 1
    ieeg_name = 'HUP88_phaseII';
    electrode_name = 'HUP088_T1_19991222_electrode_labels.csv';
    ignore = {"EKG1"};
elseif strcmp(name,'HUP094') == 1
    ieeg_name = [];
    electrode_name = 'HUP094_T1_19991225_electrode_labels.csv';
    ignore = {"EKG1","EKG2"};
elseif strcmp(name,'HUP105') == 1
    ieeg_name = 'HUP105_phaseII';
    electrode_name = 'HUP105_T1_19981114_electrode_labels.csv';
    ignore = {"EKG1"};
elseif strcmp(name,'HUP106') == 1
    ieeg_name = 'HUP106_phaseII';
    electrode_name = 'HUP106_T1_19990923_electrode_labels.csv';
    ignore = {"ECG1","ECG2"};
elseif strcmp(name,'HUP107') == 1
    ieeg_name = 'HUP107_phaseII';
    electrode_name = 'HUP107_T1_19991219_electrode_labels.csv';
    ignore = {"ECG1","ECG2"};
elseif strcmp(name,'HUP111') == 1
    ieeg_name = 'HUP111_phaseII_D01';
    electrode_name = 'HUP111_T1_19990822_electrode_labels.csv';
    ignore = {"EKG1","EKG2"};
elseif strcmp(name,'Study004') == 1
    ieeg_name = 'Study 004-2';
    electrode_name = 'Study004_electrode_labels.csv';
    ignore = {};
elseif strcmp(name,'Study012') == 1
    ieeg_name = 'Study 012-2';
    electrode_name = 'Study012_T1_19991231_electrode_labels.csv';
    ignore = {};
elseif strcmp(name,'Study016') == 1
    ieeg_name = 'Study 016';
    electrode_name = 'Study016_T1_19991230_electrode_labels.csv';
    ignore = {};
elseif strcmp(name,'Study017') == 1
    ieeg_name = 'Study 017';
    electrode_name = 'Study017_T1_19990919_electrode_labels.csv';
    ignore = {};
elseif strcmp(name,'Study019') == 1
    ieeg_name = 'Study 019';
    electrode_name = 'Study019_T1_19990811_electrode_labels.csv';
    ignore = {};
elseif strcmp(name,'Study020') == 1
    ieeg_name = 'Study 020';
    electrode_name = 'Study020_T1_19991230_electrode_labels.csv';
    ignore = {};
elseif strcmp(name,'Study021') == 1
    ieeg_name = 'Study 021';
    electrode_name = 'Study021_electrode_labels.csv';
    ignore = {};
elseif strcmp(name,'Study022') == 1
    ieeg_name = 'Study 022';
    electrode_name = 'Study022_T1_19990129_electrode_labels.csv';
    ignore = {};
elseif strcmp(name,'Study023') == 1
    ieeg_name = 'Study 023';
    electrode_name = 'Study023_electrode_labels.csv';
    ignore = {};
elseif strcmp(name,'Study026') == 1
    ieeg_name = 'Study 026';
    electrode_name = 'Study026_electrode_labels.csv';
    ignore = {};
elseif strcmp(name,'Study028') == 1
    ieeg_name = 'Study 028';
    electrode_name = 'Study028_T1_19991115_electrode_labels.csv';
    ignore = {};
elseif strcmp(name,'Study029') == 1
    ieeg_name = 'Study 029';
    electrode_name = 'Study029_T1_19980624_electrode_labels.csv';
    ignore = {};
elseif strcmp(name,'Study033') == 1
    ieeg_name = 'Study 033';
    electrode_name = 'Study033_electrode_labels.csv';
    ignore = {};
elseif strcmp(name,'Study038') == 1
    ieeg_name = 'Study 038';
    electrode_name = 'Study038_electrode_labels.csv';
    ignore = {};
elseif strcmp(name,'CHOP08') == 1
    ieeg_name = [];
    electrode_name = 'CHOP08_electrode_labels.csv';
    ignore = {};
elseif strcmp(name,'CHOP14') == 1
    ieeg_name = [];
    electrode_name = 'CHOP14_electrode_labels.csv';
    ignore = {};
elseif strcmp(name,'CHOP20') == 1
    ieeg_name = [];
    electrode_name = 'CHOP20_electrode_labels.csv';
    ignore = {};
elseif strcmp(name,'CHOP27') == 1
    ieeg_name = [];
    electrode_name = 'CHOP27_electrode_labels.csv';
    ignore = {};
elseif strcmp(name,'CHOP38') == 1
    ieeg_name = [];
    electrode_name = 'CHOP38_electrode_labels.csv';
    ignore = {};
elseif strcmp(name,'CHOP40') == 1
    ieeg_name = [];
    electrode_name = 'CHOP40_electrode_labels.csv';
    ignore = {};
elseif strcmp(name,'CHOP42') == 1
    ieeg_name = [];
    electrode_name = 'CHOP42_electrode_labels.csv';
    ignore = {};
elseif strcmp(name,'CHOP43') == 1
    ieeg_name = [];
    electrode_name = 'CHOP43_electrode_labels.csv';
    ignore = {};
elseif strcmp(name,'CHOP45') == 1
    ieeg_name = [];
    electrode_name = 'CHOP45_electrode_labels.csv';
    ignore = {};
elseif strcmp(name,'CHOP46') == 1
    ieeg_name = [];
    electrode_name = 'CHOP46_electrode_labels.csv';
    ignore = {};
elseif strcmp(name,'CHOP47') == 1
    ieeg_name = [];
    electrode_name = 'CHOP47_electrode_labels.csv';
    ignore = {};
end
    
    
    
end