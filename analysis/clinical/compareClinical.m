function compareClinical

outcomes = [...
    1; %HUP68
    1.3; %HUP70
    2; %HUP78
    1.7; %HUP86
    1; %HUP107
    1.3; %HUP111
    1; %HUP116
    4; %Study22
    3]; %Study28

vary = [...
    0; %HUP68
    1; %HUP70
    1; %HUP78
    1; %HUP86
    1; %HUP107
    1; %HUP111
    0; %HUP116
    1; %Study22
    1]; %Study 28

preIcDiff = [...
    0; %HUP68
    1; %HUP70
    1; %HUP78
    1; %HUP86
    0; %HUP107
    1; %HUP111
    0; %HUP116
    0; %Study22
    0]; %Study 28

bar([mean(outcomes(vary==0)),mean(outcomes(vary==1))])
    
[p,h,stats] = ranksum(outcomes(vary==0),outcomes(vary==1))

bar([mean(outcomes(preIcDiff==0)),mean(outcomes(preIcDiff==1))])
    
[p,h,stats] = ranksum(outcomes(preIcDiff==0),outcomes(preIcDiff==1))



end