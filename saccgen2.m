%% Simulate Eye
clear
dt = 0.001;
t = 0:dt:.5; % time array
deltaE = zeros(size(t));
deltaE(200:end) = 20; % position error (for model)
deltaE0 = zeros(size(t));
deltaE0(100:end) = 20; % position error (for plots) - used to implement saccade latency

T1 = 0.175; % visco-elastic time constant
T2 = 0.013; % inertial time constant
gPS = 600; % pulse generator gain


y(1:2) = 0; z(1:2) = 0; Eint(1:2) = 0; x(1:2) = 0; dEstar(1:2) = 0; Edot(1:2) = 0; Eint(1:2) = 0;
dEstar = zeros(length(t)); 
gRI = 1; 
for i = 1:length(t)-1
    dEstar(i+1) = dEstar(i) + dt*gRI*Edot(i); % resettable integrator
    Edot(i+1) = gPS; 
    Eint(i+1) = Eint(i) + dt*Edot(i); % neural integrator
    x(i) = (T1*Edot(i) + Eint(i)); % pulse-step addition
    z(i+1) = z(i) + dt/T2*(-z(i) + x(i)); % second ocular time constant (inertia)
    y(i+1) = y(i) + dt/T1*(-y(i) + z(i)); % output after second ocular time constant (visco-elasticity)
end

%% plot results
figure
hold on;
plot(t,deltaE0,'k:')
plot(t(1:end-1),y(1:end-1),'r')
ylabel('angular position (deg)')
