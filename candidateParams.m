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