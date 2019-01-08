function plotAlphaDelta(power,whichPt)

figure
plot(power(whichPt).times/3600,mean(power(whichPt).ad_rat,1));


end