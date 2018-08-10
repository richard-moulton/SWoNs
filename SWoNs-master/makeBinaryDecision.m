clear;

% network parameters
thres = 0.99;
N = [200 200];
K = [10 10];
q = [0.1 0.1];

% simulation parameters
h = 1e-4;
maxsteps = 1000;
nReps = 1000;

% decision model paramters
alpha = 0.9;
beta = 0.2;
eps = 0.05;

% networks initialization
for b = 1:length(N)
    nets{b} = createNetwork(N(b), K(b), q(b), false, false);
end

for q = 1:nReps
    % initial oscillator conditions
    omegas{1} = rand(1, N(1));
    omegas{2} = rand(1, N(2));
    thetas0{1} = 2 * pi * rand(1, N(1));
    thetas0{2} = 2 * pi * rand(1, N(2));
    lams0 = [0, 0];

    [thetas{q}, r{q}, lams{q}, decision(q, :, :)] = pooledInhibBinaryDecision(nets, thres, alpha, beta, eps, omegas, thetas0, lams0, maxsteps, h);
end

ohboy = sum(decision(:,:,1), 1);