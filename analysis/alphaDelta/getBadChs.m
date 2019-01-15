function getBadChs(power,whichPt)

alpha = power(whichPt).alpha;
delta = power(whichPt).delta;
ad_rat = alpha./delta;
nchs = size(alpha,1);


scatter(1:nchs,mean(alpha,2))
hold on
plot(get(gca,'xlim'),3*[median(mean(alpha,2)) median(mean(alpha,2))])


badChs = (mean(alpha,2) > 3*median(mean(alpha,2)));

figure
plot(mean(ad_rat,1));
hold on
plot(mean(ad_rat(~badChs,:),1));


end