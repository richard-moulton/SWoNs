clear;

% TODO
% - plot confidence bounds
% - RT distributions
% - save and load large networks
% - 

% network parameters
thres = 0.95;
N = [25 25];
K = [6 6];
q = [0.1 0.1];
nNets = size(N, 2);

% simulation parameters
h = 1e-4;
maxsteps = 1000;
nReps = 5;

% decision model paramters
alpha = 0.8;  % 0 = perfect memory, 1 = no memory
beta = 0.2;
eps = 0.05;

% decision sequence length and lambda update weight
seqTrials = [100 100 100]';
pL = [0.5 0.7 0.3]';  % probability of left (1) on any trial in sequence

% networks initialization
%for b = 1:nNets
%    nets{b} = createNetwork(N(b), K(b), q(b), false, false);
%end
net = createNetwork(N(1), K(1), q(1), false, false);
nets{1} = net; nets{2} = net;

optFcn = @(params) correctBinSeqFromDelta(params, nets, seqTrials, pL, alpha, beta, eps, nReps, maxsteps, h, thres);

%deltas = logspace(-5,0,20);
%for d = 1:length(deltas)
%    disp(deltas(d))
%    fracCorrect(d) = optFcn(deltas(d));
%end

LB = [0 0 0];
UB = [5 1 10];
PLB = [0.01 0.01 1];
PUB = [1 1 4];
options = optimset('PlotFcns',@optimplotfval);
[bestDelta, minFracIncorrect] = bads(optFcn, [0.15 0.6 2], LB, UB, PLB, PUB, options);
%[bestDelta, minFracIncorrect] = fminsearch(optFcn, [0.15 0.6 3], options);

%plot(deltas, fracCorrect);
