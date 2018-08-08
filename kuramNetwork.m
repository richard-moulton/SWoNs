%% Run Kuramoto oscillator through network
function [theta, z] = kuramNetwork (net, lam, omega, theta0, steps, h)

% This script  runs a kuramoto oscillator through a pre-generated network
% conditioned on local connections
% Numerical integration through fourth-order Runge-Kutta Method
% BC/ML/SWoNS/2018
% Code adapted from appmath.wordpress.com, courtesy
% Jeongho Kim, Mathematical Sciences, Seoul National University.
% Someone please cleanv

% network parameters
N = net.numnodes;

% allocate and set initial conditions
theta = zeros(N, steps);
theta(:, 1) = theta0;
z = zeros(steps, 1);
z(1) = sum(exp(1i*theta(:, 1))) / N;

for k = 1:steps-1    
    thetaPairwiseDiffs = theta(:,k) - theta(:,k)'; %Generates phase matrix\
    A = adjacency(net);
    thetaPairwiseDiffs(~A) = 0; %no connection == 0

    f1 = kuramoto (thetaPairwiseDiffs,lam,N,omega);
    f2=kuramoto(thetaPairwiseDiffs+0.5*h*f1(:,1),lam,N,omega);
    f3=kuramoto(thetaPairwiseDiffs+0.5*h*f2(:,1),lam,N,omega);           %4-th order Runge-Kutta method.
    f4=kuramoto(thetaPairwiseDiffs+h*f3(:,1),lam,N,omega);        
    theta(:,k+1) = theta(:,k)+(h/6)*(f1(:,1))+2*f2(:,1)+2*f3(:,1)+f4(:,1);

    theta(:, k+1) = wrapTo2Pi(theta(:, k+1));

    z(k+1) = sum(exp(1i*theta(:,k+1)))/N;
end

function f=kuramoto(x,lam,N,omega)
 
f=omega+(lam/N)*sum(sin(x))'; %Take out rand for noise
 
end

end