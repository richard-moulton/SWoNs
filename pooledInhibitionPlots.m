%%
maxstep = max(decision(:, 2))
if maxstep == 0
    maxstep = maxsteps
end
figure;
subplot(1,2,1)
plot(1:maxstep, r(:, 1:maxstep))
subplot(1,2,2)
plot(1:maxstep, lams(:, 1:maxstep))
