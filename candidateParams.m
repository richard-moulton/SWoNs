%% first smooth decision variables. quick divergence. flank wiggles. 
% network parameters
N = [50 50];
K = [25 25];
q = [0 0];

% simulation parameters
h = 1e-2;
maxsteps = 200;
nReps = 25;

% decision model paramters
thres = 0.99;
alpha = 0.9;
beta = 0.2;
eps = 0;

%% fast sigmoid rise, decreased K
% network parameters
N = [50 50];
K = [10 10];
q = [0 0];

% simulation parameters
h = 1e-4;
maxsteps = 200;
nReps = 10;

% decision model paramters
thres = 0.99;
alpha = 0.9;
beta = 0.2;
eps = 0;

%% wins = 1540, 891; fracNoDecision = 0.0316
% network parameters
thres = 0.95;
N = [100 100];
K = [10 10];
q = [0.1 0];

% simulation parameters
h = 1e-4;
maxsteps = 1000;
nReps = 500;

% decision model paramters
alpha = 0.8;  % 0 = perfect memory, 1 = no memory
beta = 0.2;
eps = 0.05;