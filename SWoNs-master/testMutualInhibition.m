%% Test Mutual Inhibition
K = .001;
B = 0:.001:5;
figure
for i = 1:length (B)
    hold on
    runMutual(K,B(i));
end