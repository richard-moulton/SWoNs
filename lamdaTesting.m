%% Test Lambda Parameter for Time to Synchrony
clear;
N = 50;
K = 4;
q = 0;
Lam = .08;
displayFlag = false;
testNet = createNetwork (N, K, q,displayFlag);
Edges = testNet.Edges;
iter = 2500;
numNodes = numel(unique(Edges.EndNodes(:,1)));

omega = randn(1,numNodes)*.1;  %initialize nodes with random intrinsic frequency
theta = [(2*pi*rand(1,numNodes))',zeros(numNodes,iter-1)];
f = waitbar (0,'Kuramoto Small World Testing');

for i = 1:length (Lam)
    waitbar(i/length(Lam));
    [theta2(:,:,i), r(:,i), psi(:, i)] = kuramNetwork (testNet,Edges,N,Lam(i),omega,theta,numNodes);
end
% set(0,'defaultfigurevisible','on');
% close (f);
% figure;
% imagesc(r')
% yticklabels(linspace(0,4,10))
% xlabel('time')
% ylabel('Coupling Constant')
% 
% col = brewermap (50,'*Reds');
% figure; 
% subplot(1,2,1)
% for i = 1:50
%     hold on
%     plot(r(:,i),'Color',col(i,:))
% end
% ylim([0.75 1.1]);
% subplot(1,2,2)
% imagesc(r');
