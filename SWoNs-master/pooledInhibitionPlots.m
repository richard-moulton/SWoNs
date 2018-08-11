%%
rep = 81;
maxstep = max(decision(rep,:, 2));
maxstep = maxstep;
if maxstep == 0
    maxstep = maxsteps;
end
figure;
subplot(1,2,1)
plot(1:maxstep, r{rep}(:, 1:maxstep))
subplot(1,2,2)
plot(1:maxstep, lams{rep}(:, 1:maxstep))

decision(decision==0) = NaN;
figure;
hist(decision(:,:,2),50);

col = brewermap(1000,'reds');
ind = sort(decision(:,:,2));
figure;
for i = 1:size(r,2)

    subplot(1,2,1)
    plot(1:length(r{i}(:)), r{i}(:),'Color',col(i,:))
    hold on
    subplot(1,2,2)
    plot(1:length(r{i}(:)), lams{i}(:),'Color',col(i,:))
    hold on
end
