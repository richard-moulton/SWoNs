clear all; 
steps = 50000;
N = 50;
K = 20;
q = 0.15;
lambda = 20; 

%h = 0.1;

k_range = linspace(1, floor(N/2), 100);

network = createNetwork(N, K, q, 0);
A = table2array(network.Edges(:, 1));
E = table2array(network.Edges(:, 2));
theta = floor(rand(N, 1) * 2 * pi);
omega = randn(N, 1) * 0.1;
omega = omega - min(omega);
figure;
for i = 1:steps
    z(i) = sum(exp(1i * theta(:, i))) / N;
    r(i) = real(z(i));
    psi(i) = imag(z(i)); 
    for n = 1:N
        adjs = find(A(:, 1)==n);
        %coupling = 0;
        %for j = 1:length(adjs)
        %    coupling = coupling + E(adjs(j)) * sin(theta(n, i) - theta(A(adjs(j), 2), i)) / N;
        %end
        % not completely sure this vectorization is correct...
        coupling = dot(E(adjs), sin(theta(n, i) - theta(A(adjs, 2), 1)));
        dot_theta(n, i) = omega(n) + lambda * coupling / N;
        %dotTheta = @(theta) dottheta(omega(n), lambda, r(i), psi(i));
        %k1=dotTheta(theta(n, i));
        %k2=dotTheta(theta(n, i)+0.5*h*k1);
        %k3=dotTheta(theta(n, i)+0.5*h*k2);    %4-th order Runge-Kutta method.
        %k4=dotTheta(theta(n, i)+h*k3);
    end
    for n = 1:N
        theta(n, i + 1) = theta(n, i) + dot_theta(n, i);
        if theta(n, i + 1) > 2 * pi
            theta(n, i + 1) = theta(n, i + 1) - 2 * pi;
        elseif theta(n, i + 1) < 0
            theta(n, i + 1) = theta(n, i + 1) + 2 * pi;
        end
    end
    s = linspace(0, 2 * pi, 100);
    x = cos(theta(:, i));
    y = sin(theta(:, i));
    cx = cos(s);
    cy = sin(s);
    plot(x, y, 'o', cx, cy)
    axis([-1 1 -1 1])
    axis square
    
    % some more information to plot
    txt_fmt = '%8.3f';
    r_txt = sprintf(txt_fmt, r(i));
    psi_txt = sprintf(txt_fmt, psi(i));
    mufreq_txt = sprintf(txt_fmt, mean(theta(:, i)));
    text(-0.4, 0.1, 'r', 'FontWeight', 'bold');
    text(0.4, 0.1, r_txt, 'HorizontalAlignment', 'right')
    text(-0.4, 0, 'psi', 'FontWeight', 'bold');
    text(0.4, 0, psi_txt, 'HorizontalAlignment', 'right')
    text(-0.4, -0.1, 'mean freq', 'FontWeight', 'bold');
    text(0.4, -0.1, mufreq_txt, 'HorizontalAlignment', 'right')
    drawnow
    %pause(0.1)
end

