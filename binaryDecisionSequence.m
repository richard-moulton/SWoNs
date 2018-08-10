clear;

% TODO
% - plot confidence bounds
% - RT distributions
% - save and load large networks
% - 

% network parameters
thres = 0.95;
N = [100 100];
K = [10 10];
q = [0.1 0.1];
nNets = size(N, 2);

% simulation parameters
h = 1e-4;
maxsteps = 1000;
nReps = 1000;

% decision model paramters
alpha = 0.8;  % 0 = perfect memory, 1 = no memory
beta = 0.2;
eps = 0.05;

% lambda update weight
delta = 0.1;

% networks initialization
%for b = 1:nNets
%    nets{b} = createNetwork(N(b), K(b), q(b), false, false);
%end
tic
net = createNetwork(N(1), K(1), q(1), false, false);
toc
disp('Network created!')
nets{1} = net; nets{2} = net;

for q = 1:nReps
    % initial oscillator conditions
    omegas{1} = rand(1, N(1));
    omegas{2} = rand(1, N(2));
    thetas0{1} = 2 * pi * rand(1, N(1));
    thetas0{2} = 2 * pi * rand(1, N(2));
    lams0 = [0, 0];
    if ~mod(q, 100)
        disp(q)
    end
    [trialThetas{q}, r(q, :, :, :), ~, lams(q, :, :, :), decision(q, :, :, :)] = pooledInhibBinaryDecision(nets, trueSeq, thres, alpha, beta, delta, eps, omegas, thetas0, lams0, maxsteps, h);
end

wins = sum(decision(:,:,1), 1);
rt = decision(:, :, 2);
anyDecision = logical(sum(decision(:,:,1), 2));
fracNoDecision = sum(~anyDecision) / nReps;