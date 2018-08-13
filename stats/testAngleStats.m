angles = [90 50 8 5 -90 90 -135 -135 -90 -45];
ic = {'ic','ic','ic','ic','ic','interic','interic','interic','interic','interic'};

[h,p,ci,stats] = ttest(angles(1:5),angles(6:10));


%alpha = [pi/4, pi/5, 3*pi/4, 3.5*pi/4, pi/4, pi/6, 3*pi/4, 3.4*pi/4];
alpha = [pi/4, pi/5, 3*pi/4, 3.5*pi/4, pi/2, 1.1*pi/2, pi/4, 1.2*pi/4];
idp = [1 1 1 1 2 2 2 2];
idq = [1 1 2 2 1 1 2 2];
inter = 1;
fn = {'patient','ictal'};