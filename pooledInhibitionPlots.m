%%
rep = 801;
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
