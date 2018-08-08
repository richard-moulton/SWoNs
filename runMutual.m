%% Simulate Mutual Inhibition
function runMutual (K,B)

dt = .1;
time = 0:dt:1000;
[a1,a2] = deal (zeros(1,length(time)));
dR = randn(length(time),1)./1000;

hold on
xlim([0 1000]);
ylim([-.05 .05]);

for j = 1:length (time)
        % Need here an update for difference of synchrony b/w networks
        [dA1,dA2] = updateInhib (a1(j),a2(j),dR(j),B,K,eps,dt);
        a1(j+1) = a1(j) + dA1;
        a2(j+1) = a2(j) + dA2;
end
p = plot(a1);
p2 = plot(a2);


title('Mutual Inhibition');
xlabel('Time')
ylabel('Activity');
legend([p, p2], 'Network 1','Network 2')
end


