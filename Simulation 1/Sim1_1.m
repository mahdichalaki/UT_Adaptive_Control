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

%% Stability Checking
figure()
subplot(2,1,1); 
rlocus(G)
subplot(2,1,2);
rlocus(sysd)

%% Parameter extraction for LS estimation

[num,den]=tfdata(sysd,'v');
n = 8;
theta = [den(2:5) num(2:5)]';
N = 1000000;

noise_var = 0.01;
noise = sqrt(noise_var) * randn(N,1);
noise = noise - mean(noise);

%% Y - noise included

u = sqrt(noise_var) * randn(N,1);
[theta_hat_noisy,cost_noisy] = ls(u,theta,noise,N)

%% Y - Pulse

u = zeros(N,1);
u(1) = 1;
[theta_hat_pulse,cost_pulse] = ls(u,theta,noise,N)

%% Y - Step

u = ones(N,1);
[theta_hat_step,cost_step] = ls(u,theta,noise,N)

%% Y - Sinsusoid

u = zeros(N,1);
omega = 0.01;
for i=1:N
   u(i,1) = sin(omega*i*Ts);
end 
[theta_hat_sinusoid,cost_sinusoid] = ls(u,theta,noise,N)

%% Y - Ramp

u = (linspace(1,N,N))';
[theta_hat_ramp,cost_ramp] = ls(u,theta,noise,N)

%% Y - Under Parameter

u = sqrt(noise_var) * randn(N,1);
[theta_hat_under,cost_under] = ls_under(u,theta,noise,N)

%% Y - Over Parameter

u = sqrt(noise_var) * randn(N,1);
[theta_hat_over,cost_over] = ls_over(u,theta,noise,N)

%% LS Method

function [theta_hat, Cost] = ls(u,theta,noise,N)
    y = zeros(N,1);
    y_zero = 0;
    y_minus1 = 0;
    y_minus2 = 0;

    u_zero = 0;
    u_minus1 = 0;
    u_minus2 = 0;
    
    y(1) = 0;
    y(2) = [-y(1) -y_zero -y_minus1 -y_minus2 u(1) u_zero u_minus1 u_minus2] * theta + noise(2);
    y(3) = [-y(2) -y(1) -y_zero -y_minus1 u(2) u(1) u_zero u_minus1] * theta + noise(3);
    y(4) = [-y(3) -y(2) -y(1) -y_zero u(3) u(2) u(1) u_zero] * theta + noise(4);

    for k = 5:N
        y(k)=[-y(k-1) -y(k-2) -y(k-3) -y(k-4) u(k-1) u(k-2) u(k-3) u(k-4)] * theta + noise(k);
    end
    
    Phi = zeros(N,length(theta));
    Phi(2,:) = [-y(1) -y_zero -y_minus1 -y_minus2 u(1) u_zero u_minus1 u_minus2];
    Phi(3,:) = [-y(2) -y(1) -y_zero -y_minus1 u(2) u(1) u_zero u_minus1];
    Phi(4,:) = [-y(3) -y(2) -y(1) -y_zero u(3) u(2) u(1) u_zero];

    for k = 5:N
        Phi(k,:)=[-y(k-1) -y(k-2) -y(k-3) -y(k-4) u(k-1) u(k-2) u(k-3) u(k-4)];
    end

    theta_hat = inv(Phi'*Phi)*Phi'*y;
    
    %plotting
    x = (linspace(1,N,N))';
    y_hat = Phi*theta_hat;
    Error = y-y_hat;
    Cost = 1/2 * (Error' * Error);
    
    figure()
    plot(x,y)
    hold on
    plot(x,y_hat)
    xlabel('sample number')
    ylabel('Output')
    title('Output of Actual and Estimated systems')
    legend('Actual system','Estimated system')
    xlim([0 100])
end 

%% LS - Under Parameter

function [theta_hat, Cost] = ls_under(u,theta,noise,N)
    y = zeros(N,1);
    y_zero = 0;
    y_minus1 = 0;
    y_minus2 = 0;

    u_zero = 0;
    u_minus1 = 0;
    u_minus2 = 0;
    
    y(1) = 0;
    y(2) = [-y(1) -y_zero -y_minus1 -y_minus2 u(1) u_zero u_minus1 u_minus2] * theta + noise(2);
    y(3) = [-y(2) -y(1) -y_zero -y_minus1 u(2) u(1) u_zero u_minus1] * theta + noise(3);
    y(4) = [-y(3) -y(2) -y(1) -y_zero u(3) u(2) u(1) u_zero] * theta + noise(4);

    for k = 5:N
        y(k)=[-y(k-1) -y(k-2) -y(k-3) -y(k-4) u(k-1) u(k-2) u(k-3) u(k-4)] * theta + noise(k);
    end
    
    Phi = zeros(N,6);
    Phi(2,:) = [-y(1) -y_zero -y_minus1 u(1) u_zero u_minus1];
    Phi(3,:) = [-y(2) -y(1) -y_zero u(2) u(1) u_zero];

    for k = 4:N
        Phi(k,:)=[-y(k-1) -y(k-2) -y(k-3) u(k-1) u(k-2) u(k-3)];
    end

    theta_hat = inv(Phi'*Phi)*Phi'*y;
%     RMSE = sqrt(mean((theta - theta_hat).^2));
    
    %plotting
    x = (linspace(1,N,N))';
    y_hat = Phi*theta_hat;
    Error = y-y_hat;
    Cost = 1/2 * (Error' * Error);
    
    figure()
    plot(x,y)
    hold on
    plot(x,y_hat)
    xlabel('sample number')
    ylabel('Output')
    title('Output of Actual and Estimated systems')
    legend('Actual system','Estimated system')
    xlim([0 100])
end 

%% LS - Over Parameter

function [theta_hat, Cost] = ls_over(u,theta,noise,N)
    y = zeros(N,1);
    y_zero = 0;
    y_minus1 = 0;
    y_minus2 = 0;
    y_minus3 = 0;

    u_zero = 0;
    u_minus1 = 0;
    u_minus2 = 0;
    u_minus3 = 0;
    
    y(1) = 0;
    y(2) = [-y(1) -y_zero -y_minus1 -y_minus2 u(1) u_zero u_minus1 u_minus2] * theta + noise(2);
    y(3) = [-y(2) -y(1) -y_zero -y_minus1 u(2) u(1) u_zero u_minus1] * theta + noise(3);
    y(4) = [-y(3) -y(2) -y(1) -y_zero u(3) u(2) u(1) u_zero] * theta + noise(4);

    for k = 5:N
        y(k)=[-y(k-1) -y(k-2) -y(k-3) -y(k-4) u(k-1) u(k-2) u(k-3) u(k-4)] * theta + noise(k);
    end
    
    Phi = zeros(N,10);
    Phi(2,:) = [-y(1) -y_zero -y_minus1 -y_minus2 -y_minus3 u(1) u_zero u_minus1 u_minus2 u_minus3];
    Phi(3,:) = [-y(2) -y(1) -y_zero -y_minus1 -y_minus2 u(2) u(1) u_zero u_minus1 u_minus2];
    Phi(4,:) = [-y(3) -y(2) -y(1) -y_zero -y_minus1 u(3) u(2) u(1) u_zero u_minus1];
    Phi(5,:) = [-y(4) -y(3) -y(2) -y(1) -y_zero u(4) u(3) u(2) u(1) u_zero];

    for k = 6:N
        Phi(k,:)=[-y(k-1) -y(k-2) -y(k-3) -y(k-4) -y(k-5) u(k-1) u(k-2) u(k-3) u(k-4) u(k-5)];
    end

    theta_hat = inv(Phi'*Phi)*Phi'*y;
    
    %plotting
    x = (linspace(1,N,N))';
    y_hat = Phi*theta_hat;
    Error = y-y_hat;
    Cost = 1/2 * (Error' * Error);
    
    figure()
    plot(x,y)
    hold on
    plot(x,y_hat)
    xlabel('sample number')
    ylabel('Output')
    title('Output of Actual and Estimated systems')
    legend('Actual system','Estimated system')
    xlim([0 100])
end 