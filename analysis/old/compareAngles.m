angles = getAngles(pt(1).sz(1).seq_matrix,pt(1).electrodeData);
ictal = icOrInteric(pt(1).sz(1).seq_matrix,[pt(1).sz(1).onset,pt(1).sz(1).offset]);

ic_yes = find(ictal == 1);
ic_no = find(ictal == 0);

%{
scatter(ic_yes,angles(ic_yes),'r');
hold on
scatter(ic_no,angles(ic_no),'b');
%}

angles = angles*pi/180;

[pval, table] = circ_wwtest(angles(ic_yes), angles(ic_no))