function newStr = chParser(origStr)

% split it with spaces
C = strsplit(origStr,' ');

% I would expect all of the names to start with EEG
if strcmp(C{1},'EEG') == 0
    fprintf('Warning, there is something weird in the channel labels for channel %d in patient %s\n',i,dataName);
    C = C{1};

else
    % remove 'EEG '
    C = strrep(origStr,[C{1},' '],'');

    % Remove -Ref
    D = strsplit(C,'-');

    C = strrep(C,['-',D{2}],'');
end


% Remove space if present
C = strrep(C,' ','');

% Get the numbers
numIdx = regexp(C,'\d');

if isempty(numIdx) == 0
    
    % remove leading zeros
    if strcmp(C(numIdx(1)),'0') == 1
        C(numIdx(1)) = [];
    end
end

% Special thing to make ECG->EKG
%C = strrep(C,'ECG','EKG');

% Final channel name
chName = C;
newStr = chName;

end