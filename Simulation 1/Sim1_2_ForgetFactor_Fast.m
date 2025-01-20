close all
clear
clc

%% Parameters
m1 = 0.1;
m2 = m1;
k1 = 2/3;
b2 = k1;
k2 = 0.1;
b1 = k2;
k3 = 2*k1/3;

%% Transfer Function

s=tf('s');
G = k3/((m1*s^2+b1*s+k1+k3)*(m2*s^2+b2*s+k2+k3) - k3^2);
bw = bandwidth(G);
Ts = (2*pi)/(20*bw);
sysd = c2d(G, Ts, 'zoh');


%% Parameter extraction for LS estimation

[num,den]=tfdata(sysd,'v');
n = 8;
N = 100000;
% lambda = linspace(0.95,0.999,N);
lambda = 0.9999*ones(1,N);
theta = zeros(n,N);
theta(:,1:N/5) = repmat([den(2:5) num(2:5)]',1,N/5);
theta(:,1*N/5+1:2*N/5) = 1.05*repmat([den(2:5) num(2:5)]',1,N/5);
theta(:,2*N/5+1:3*N/5) = 0.95*repmat([den(2:5) num(2:5)]',1,N/5);
theta(:,3*N/5+1:4*N/5) = 1.1*repmat([den(2:5) num(2:5)]',1,N/5);
theta(:,4*N/5+1:N) = 1*repmat([den(2:5) num(2:5)]',1,N/5);

theta_hat_init = zeros(size(theta,1),N);
P = 1e9 * eye(size(theta,1));

%% Y - noise included

noise_var = 0.01;
noise = sqrt(noise_var) * randn(N,1);
noise = noise - mean(noise);

u = sqrt(noise_var) * randn(N,1);
[theta_hat_noisy,cost_noisy] = RLS(u,theta,theta_hat_init,noise,N,P,lambda);

%% Function

function [theta_hat, Cost] = RLS(u,theta,theta_hat_init,noise,N,P,lambda)

    y = zeros(N,1);
    y_zero = 0;
    y_minus1 = 0;
    y_minus2 = 0;

    u_zero = 0;
    u_minus1 = 0;
    u_minus2 = 0;
    
    y(1) = 0;
    y(2) = [-y(1) -y_zero -y_minus1 -y_minus2 u(1) u_zero u_minus1 u_minus2] * theta(:,2) + noise(2);
    y(3) = [-y(2) -y(1) -y_zero -y_minus1 u(2) u(1) u_zero u_minus1] * theta(:,3) + noise(3);
    y(4) = [-y(3) -y(2) -y(1) -y_zero u(3) u(2) u(1) u_zero] * theta(:,4) + noise(4);
    
    Phi = zeros(N,size(theta,1));
    Phi(2,:) = [-y(1) -y_zero -y_minus1 -y_minus2 u(1) u_zero u_minus1 u_minus2];
    Phi(3,:) = [-y(2) -y(1) -y_zero -y_minus1 u(2) u(1) u_zero u_minus1];
    Phi(4,:) = [-y(3) -y(2) -y(1) -y_zero u(3) u(2) u(1) u_zero];

    theta_hat = theta_hat_init;
    y_hat = zeros(N,1);
    Error = zeros(N,1);
    
    for k = 5:N 
        
        y(k)=[-y(k-1) -y(k-2) -y(k-3) -y(k-4) u(k-1) u(k-2) u(k-3) u(k-4)] * theta(:,k) + noise(k);
        Phi(k,:)=[-y(k-1) -y(k-2) -y(k-3) -y(k-4) u(k-1) u(k-2) u(k-3) u(k-4)];
        K = P * Phi(k,:)'*inv(lambda(k) + Phi(k,:)*P*Phi(k,:)');
        P = (P - K * Phi(k,:) * P)/lambda(k);
        theta_hat(:,k) = theta_hat(:,k-1) + K*(y(k)-Phi(k,:)*theta_hat(:,k-1));
        y_hat(k) = Phi(k,:)*theta_hat(:,k);
        Error(k) = y(k) - y_hat(k);
        
    end
    
    %plotting
    x = (linspace(1,N,N))';
    Cost = 1/2 * (Error' * Error);
    limx = 100000; %for plot xlim
    
    figure()
    plot(x,y)
    hold on
    plot(x,y_hat)
    xlabel('sample number')
    ylabel('Output')
    title('Output of Actual and Estimated systems')
    legend('Actual system','Estimated system')
    xlim([0 2000])
    
    figure()
    xlim([0 200])
    subplot(2,2,1);
    plot(x,theta(1,:))
    hold on
    plot(x,theta_hat(1,:))
    title('a_1')
    xlim([0 limx])
    
    subplot(2,2,2);
    plot(x,theta(2,:))
    hold on
    plot(x,theta_hat(2,:))
    title('a_2')
    xlim([0 limx])
    
    subplot(2,2,3);
    plot(x,theta(3,:))
    hold on
    plot(x,theta_hat(3,:))
    title('a_3')
    xlim([0 limx])
    
    subplot(2,2,4);
    plot(x,theta(4,:))
    hold on
    plot(x,theta_hat(4,:))
    title('a_4')
    xlim([0 limx])
    
    
    figure()
    subplot(2,2,1);
    plot(x,theta(5,:))
    hold on
    plot(x,theta_hat(5,:))
    title('b_1')
    xlim([0 limx])
    
    subplot(2,2,2);
    plot(x,theta(6,:))
    hold on
    plot(x,theta_hat(6,:))
    title('b_2')
    xlim([0 limx])
    
    subplot(2,2,3);
    plot(x,theta(7,:))
    hold on
    plot(x,theta_hat(7,:))
    title('b_3')
    xlim([0 limx])
    
    subplot(2,2,4);
    plot(x,theta(8,:))
    hold on
    plot(x,theta_hat(8,:))
    title('b_4')
    xlim([0 limx])
    
end