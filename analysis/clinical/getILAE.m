function ilae = getILAE(outcomes)


%{
This function takes outcome scores (generated from getOutcome) and returns
ILAE scores, using the definitions from getOutcome
%}
ilae = zeros(size(outcomes));

for i = 1:length(outcomes)
    if outcomes(i) == 1
        ilae(i) = 1;
    elseif outcomes(i) == 1.25
        ilae(i) = 2;
    elseif outcomes(i) == 2.25
        ilae(i) = 3;
    elseif outcomes(i) == 3
        ilae(i) = 4;
    elseif outcomes(i) == 4
        ilae(i) = 5;
    elseif outcomes(i) == 4.5
        ilae(i) = 6;
    else
        error('weird outcome\n');
    end
end

end