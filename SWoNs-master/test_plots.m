r = abs(z);
psi = angle(z);

%%
stds = squeeze(std(theta2,[], 1))';
figure;
imagesc(stds);

figure;
imagesc(log(r' ./ stds(:, 1:end-1)))

%% synch animation
theta_i = theta{end, 20};

figure;
for i = 1:size(theta_i, 2)
    
        x=cos(theta_i(:, i)');% - psi(i, :));
        y=sin(theta_i(:, i)');% - psi(i, :));
        
        s=linspace(0,2*pi,100);
        cx=cos(s);
        cy=sin(s);
        p = plot(x,y,'.',cx,cy);
        set(p,'MarkerSize',20);
        axis([-1 1 -1 1])
        axis square
        drawnow
        pause(0.1);
end
%% r vs. step, lambda
figure;
lam = cell2mat(oscParamCombs);
r = abs(squeeze(z(1, 1, 1, :, :)));
imagesc(1:steps, lam, r)
%set(gca,'YScale','log')
xlabel('Time')
ylabel('Coupling Constant')

%% r vs. step, K (lambda also variable between trials)
figure;
nLam = 15;
r = abs(squeeze(z(:,nLam,:)));
imagesc(r)

%% sync step (r threshold) vs. K, lam
thres=1;
r = abs(z);
lam = cell2mat(oscParamCombs);
netParams = cell2mat(netParamCombs);
K = netParams(:, 2);
for j = 1:size(r, 1)  % each K
for i = 1:size(r, 2)  % each lambda
    passesThres = find(r(j, i, :) > thres, 1);
    if isempty(passesThres)
        sync_step(j, i) = steps;
    else
        sync_step(j, i) = passesThres;
    end
end
end
figure;
imagesc(K*N/2, lam, sync_step');
xlabel('K')
ylabel('Coupling Constant')
title(sprintf('Step when r exceeds %0.2f', thres));
colorbar;

%% trying to animate synch along with synch vector
theta_p = theta{1, end};
z_p = squeeze(z(1, end, :));
r=abs(z_p);
psi=angle(z_p);

s=linspace(0,2*pi,100);
cx=cos(s);
cy=sin(s);

figure;
%p = plot(x,y,'.',cx,cy);
p=animatedline('Marker', '.', 'MarkerSize',20);
%p_z = plot([0 zx],[0 zy]); 
p_z=animatedline('LineWidth', 2);
%set(p_z, 'LineWidth', 2);
%set(p,'MarkerSize',20);
axis([-1 1 -1 1])
axis square

for i = 1:size(theta_p, 2)
        x=cos(theta_p(:, i)');% - psi(i, :));
        y=sin(theta_p(:, i)');% - psi(i, :));
        %set(p, 'XData', x);
        %set(p, 'YData', y);
        clearpoints(p);
        addpoints(p, x, y);
        
        
        zx = r(i) * cos(psi(i));
        zy = r(i) * sin(psi(i));
        %set(p_z, 'XData', zx);
        %set(p_z, 'YData', zy);  
        clearpoints(p_z);
        addpoints(p_z, zx, zy)
        drawnow 
        pause(0.1);
end

