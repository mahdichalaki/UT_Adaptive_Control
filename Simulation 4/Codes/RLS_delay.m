clear;
close all;
clc;

%% Transfer Function

z = tf('z');
Ts = 0.02;
sysd = (z-0.4)/((z-0.22)*(z-0.73)*z^3);
[numd,dend]=tfdata(sysd,'v');

%%

theta = [dend(2:end) numd]';
N = 10000;
theta_hat_init = zeros(length(theta),N);
P = 100000 * eye(length(theta));

%% Y - noise included

rng(2)
noise_var = 0.001;
noise = sqrt(noise_var) * randn(N,1);
noise = noise - mean(noise);

rng(3)
u_noise = sqrt(noise_var) * randn(N,1);
u_noise = u_noise - mean(u_noise);

u = sqrt(noise_var) * randn(N,1);
% u = zeros(N,1);

%% MDPP

t = N;

for i=1:t
    
y_star(i,:,1) = 1 + u_noise(i);
y_star(i,:,2) = (-1).^ceil(i/2000);
% u(i) = (-1).^ceil(i/2000) + u_noise(i);

Disturbance(i,:,1) = 0;
Disturbance(i,:,2) = 0;

end

i=2;

%% Solver

t = N;
y = zeros(N,1);

theta_hat = theta_hat_init;
y_hat = zeros(N,1);
Error = zeros(N,1);
Phi = zeros(N,length(theta));
    
for k = 6:N
    
%     if mod(k,4e3) == 0
%             P = 100000 * eye(size(theta,1));
%     end


    y(k)=[-y(k-1:-1:k-(length(numd)-1))', u(k:-1:k-(length(dend)-1))'] * theta + noise(k);
    Phi(k,:)=[-y(k-1:-1:k-(length(numd)-1))', u(k:-1:k-(length(dend)-1))'];
    P = P - P*Phi(k,:)'*inv(1+Phi(k,:)*P*Phi(k,:)')*Phi(k,:)*P;
    K = P * Phi(k,:)';
    theta_hat(:,k) = theta_hat(:,k-1) + K*(y(k)-Phi(k,:)*theta_hat(:,k-1));
    y_hat(k) = Phi(k,:)*theta_hat(:,k);
    Error(k) = y(k) - y_hat(k);
    
end

%% 
dend_hat = theta_hat(1:length(dend)-1,end);
threshold = 0.01;
degree = 0;

for i = 1:length(dend)-1
    
    if abs(dend_hat(end+1-i)) < threshold
        degree = degree+1;
    else
        break
    end 
end

disp (['Delay = ', num2str(degree)])

%% plotting

x = (linspace(1,N,N))';
limx = N; %for plot xlim

figure()
for i = 1:length(numd)-1
    subplot(2,3,i);
    title_text = "a_%d";
    plot(x,theta(i)*ones(1,N)) 
    hold on;
    plot(x,theta_hat(i,:))
    title(sprintf(title_text, i));
    xlim([0 limx])
    xlabel('Samples')
    legend('Real','Estimated');
    
end

figure()
for i = length(numd):length(numd)+length(dend)-1
    subplot(2,3,i-5);
    title_text = "b_%d";
    plot(x,theta(i)*ones(1,N)) 
    hold on;
    plot(x,theta_hat(i,:))
    title(sprintf(title_text, i-6));
    xlim([0 limx])
    xlabel('Samples')
    legend('Real','Estimated');
    
end
