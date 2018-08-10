%%
rep = 5;
maxstep = max(decision(rep, :, 2));
maxstep = maxstep;
if maxstep == 0
    maxstep = maxsteps;
end
figure;
subplot(1,2,1)
plot(1:maxstep, r(rep, :, 1:maxstep))
subplot(1,2,2)
plot(1:maxstep, lams(rep, :, 1:maxstep))

%% DT distributions (Habiba)
set(0, 'defaultAxesFontSize', 14);
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
winFracs = wins / sum(wins);
undecidedFrac = 1 - sum(wins) / nReps;
annotTxt = sprintf(' f_{L} = %0.3f\n f_{R} = %0.3f\n f_{?} = %0.3f', winFracs(1), winFracs(2), undecidedFrac);
annot = annotation('textbox', [.35 .6 .3 .3], 'String', annotTxt, 'FitBoxToText', 'on');
annot.FontSize = 12;
subplot(1,2,2)
cdf1 = cdfplot(rtNonZero(:, 1));
set(cdf1, 'LineWidth', 2);
hold on
cdf2 = cdfplot(rtNonZero(:, 2));
%xlim([0 maxsteps])
set(cdf2, 'LineWidth', 2);
set(gca,'XScale','log')
xlabel('Timestep')
ylabel('p(t_{D} \geq t)')
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


%% sequential decision analysis
% decision = zeros(seqTrials, nNets, nSeq, 2);
figure;
q = 1;  % which repeat to plot
decisions = squeeze(decision(q, :, :, 1));  
subplot(2,1,1); hold on;
plot(1:seqTrials, trueSeq{q}, '*k');
plot(1:seqTrials, decisions(1, :), '.b')
plot(1:seqTrials, decisions(2, :), '.r')
ylim([1 2])
subplot(2,1,2);
plot(1:seqTrials, lams{q}(:, :, 1))