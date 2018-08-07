function [theta, z] = kuramNetwork(net, lam, omega, theta0, steps, h)
% This script  runs a kuramoto oscillator through a pre-generated network
% conditioned on local connections
% Numerical integration through fourth-order Runge-Kutta Method
% BC/ML/SWoNS/2018
% Code adapted from appmath.wordpress.com, courtesy
% Jeongho Kim, Mathematical Sciences, Seoul National University.

% network parameters
N = net.numnodes;
A = full(adjacency(net));
w = net.Edges.Weight;

% pre-allocate arrays and set initial conditions
r = zeros(length(steps-1), 1);
theta = zeros(N, steps);
theta(:, 1) = theta0;

% partial function for evaluating the numerical integral
kuramotoPartial = @(thetaDiffs) kuramoto(thetaDiffs, lam, N, omega);

for iter = 1:steps-1
    % pairwise inter-node phase differences
    thetaDiffs = theta(:, iter)*ones(1, N) - (ones(N, 1)*theta(:, iter)'); 
    thetaDiffs(~A) = 0;  % exclude non-edges
    
    % numerical integration step
    theta(:, iter+1) = theta(:, iter) + rk4StepNonZero(kuramotoPartial, thetaDiffs, h, A)';
    
    % boundary conditions to keep theta in [0, 2*pi)
    indOver = theta(:, iter+1) >= 2*pi;
    indUnder = theta(:, iter+1) < 0;
    theta(indOver, iter+1) = theta(indOver, iter+1) - 2*pi;
    theta(indUnder, iter+1) = theta(indUnder, iter+1) + 2*pi;

    % network mean phase vector
    z = sum(exp(1i*sin(theta(:, iter+1)))) / N;  
end

function dTh = kuramoto(x, lam, N, omega)
    % Kuramoto oscillator differential equation 
    dTh = omega + (lam / N) * sum(sin(x));
end

function res = rk4StepNonZero(dxdt, x, h, posInds)
    % 4-th order Runge-Kutta method on non-zero values in x
    f1 = dxdt(x);
    x2 = broadcastAddToNonZero(x, 0.5*h*f1, posInds);
    f2 = dxdt(x2);
    x3 = broadcastAddToNonZero(x, 0.5*h*f2, posInds);
    f3 = dxdt(x3);
    x4 = broadcastAddToNonZero(x, h*f3, posInds);
    f4 = dxdt(x4);        
    res = (h/6)*f1 + 2*f2 + 2*f3 + f4;
end

function s = broadcastAddToNonZero(x, y, posInds)
    % broadcast addition of vector y to matrix x
    % but only to non-zero values in x, given by logical posInds
    yBroadcast = y + zeros(length(y));
    s = x;
    inds = logical(posInds);
    s(inds) = s(inds) + yBroadcast(inds);
end
end