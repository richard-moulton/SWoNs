%% Test Lambda Parameter for Time to Synchrony
clear;
N = 50;
K = 1;
q = 0.1;
Lam = linspace(0,4,50);
displayFlag = false;
testNet = createNetwork (N, K, q,displayFlag);
Edges = testNet.Edges;
iter = 500;
numNodes = numel(unique(Edges.EndNodes(:,1)));

omega = randn(1,numNodes)*.1;  %initialize nodes with random intrinsic frequency
theta = [(2*pi*rand(1,numNodes))',zeros(numNodes,iter-1)];
set(0,'defaultfigurevisible','off');
f = waitbar (0,'Kuramoto Small World Testing');
for i = 1:length (Lam)
    waitbar(i/length(Lam));
    [theta2(:,:,i), r(:,i), psi(:, i)] = kuramNetwork (Edges,N,Lam(i),omega,theta,numNodes);
end
set(0,'defaultfigurevisible','on');

figure;
imagesc(r')
yticklabels(linspace(0,4,10))
xlabel('time')
ylabel('Coupling Constant')


stds = squeeze(std(theta2,[], 1))';
figure;
imagesc(stds);

figure;
imagesc(log(r' ./ stds(:, 1:end-1)))
% figure;
% for i = 1:500
%     
%         x=cos(theta2(:,i,4)' - psi(i, :));
%         y=sin(theta2(:,i,4)' - psi(i, :));
%         
%         s=linspace(0,2*pi,100);
%         cx=cos(s);
%         cy=sin(s);
%         p = plot(x,y,'.',cx,cy);
%         set(p,'MarkerSize',20);
%         axis([-1 1 -1 1])
%         axis square
%         drawnow
% end
% 


