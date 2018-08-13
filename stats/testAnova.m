

spike_freq =  [1 2 5 7 10 12 2 3];
%spike_freq = [1 2 5 7 10 12 22 23];
patient = {'Steve','Steve','Steve','Steve','Erin','Erin','Erin','Erin'};
ic = {'interictal','interictal','ictal','ictal','interictal','interictal',...
    'ictal','ictal'};
whichSz = [1 2 1 2 3 4 3 4];

p1 = anovan(spike_freq,{patient,ic});
p2 = anovan(spike_freq,{patient,ic},'model','interaction','varnames',...
    {'patient','ic'});