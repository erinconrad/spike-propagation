fakePoints = randi(1000,1000,1);
fakeSOZ = [200 300];
fakeDist = zeros(1000,1);
for i = 1:size(fakePoints,1)
    fakeDist(i) = min(abs(fakePoints(i)-fakeSOZ));
end

nboot = 1e3;
dist_diff = zeros(nboot,1);
for ib = 1:nboot
small = fakeDist(randperm(1000,900));
small = randperm(1000,900);
all = 1:1000;
large = setdiff(all,small);
sdist = fakeDist(small);
ldist = fakeDist(large);
dist_diff(ib) = mean(sdist)-mean(ldist);
end

hist(dist_diff)
hold on
plot([0 0],ylim)