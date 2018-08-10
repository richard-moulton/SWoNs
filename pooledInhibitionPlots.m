%%
rep = 5;
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

%% RT distributions
nBins = 150;
rtNonZero = rt;
rtNonZero(rtNonZero == 0) = nan; %maxsteps+1;
figure('Color', 'w', 'Position', [800 100 1000 500]);
subplot(1,2,1)
[rtCounts, rtCenters] = hist(rtNonZero, nBins);
totalCount = sum(cumsum(rtCounts), 2);
xCutoff = find(totalCount > totalCount(end) - 1, 1);
bar(rtCenters, rtCounts / totalCount(end), 1)
xlim([0 xCutoff])
ylabel('p(t_{D} = t)')
xlabel('Timestep')
title('PDF');
subplot(1,2,2)
cdf1 = cdfplot(rtNonZero(:, 1));
set(cdf1, 'LineWidth', 2);
hold on
cdf2 = cdfplot(rtNonZero(:, 2));
%xlim([0 maxsteps])
set(cdf2, 'LineWidth', 2);
set(gca,'XScale','log')
xlabel('Timestep')
ylabel('p(t_{D} ? t)')
title('CDF');
legend('L', 'R', 'Location', 'southeast')

%% r traces (not for presentation)
figure; 
subplot(1,2,1); hold on; 
subplot(1,2,2); hold on; 
for q = 1:250
    rs = r{q}';
    if anyDecision(q)
        subplot(1,2,1)
        plot(rs(:, 1), 'r')
    	plot(rs(:, 2), 'b')
    else
        subplot(1,2,2)
        plot(rs(:, 1), 'r')
        plot(rs(:, 2), 'b')
    end

end
