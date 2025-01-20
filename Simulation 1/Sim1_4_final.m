close all
clear
clc

%% Parameters
z1 = -5;
z2 = -2;
p1 = -1+4i;
p2 = -1-4i;
p3 = -3;

%% Transfer Function

s = tf('s');
G = (s-z1)*(s-z2)/((s-p1)*(s-p2)*(s-p3));
bw = bandwidth(G);
Ts = (2*pi)/(20*bw);
sysd = c2d(G, Ts, 'zoh');
step(G)

[C,~] = pidtune(sysd,'PI');

H = feedback(sysd*C,1);
step(H)

%% Parameter extraction for LS estimation
[num,den]=tfdata(sysd,'v');
n = 6;
theta = [den(2:4) num(2:4)]';
I = eye(length(theta));
R = 0.01;
theta_hat_init = zeros(length(theta),N-3);

%% From Simulink 

u = out.u_noisy;
y = out.y_noisy;
N=length(y);

%% Kalman Filter

x_k_k = zeros(length(theta),1);
P_k_k = 1e6 * I;

theta_hat = theta_hat_init;

% Sampling

 for k = 4:N
     
     Phi = [-y(k-1) -y(k-2) -y(k-3) u(k-1) u(k-2) u(k-3)]';
     x_k_1_k = I * x_k_k;
     P_k_1_k = I * P_k_k * I';
     K = P_k_1_k * Phi * inv(Phi'*P_k_1_k*Phi+R);
     x_k_k = x_k_1_k + K*(y(k) - Phi'*x_k_1_k);
     P_k_k = (I - K*Phi') * P_k_1_k * (I - K*Phi')' + K*R*K';
     theta_hat(: , k-3) = x_k_k;
     
 end


Y = y(4:end);


figure
subplot(2,2,1)
samples = 1:size(theta_hat , 2);
theta = [den(2:end),num(2:end)];
plot(samples , theta_hat(1,:))
hold on
plot(samples , ones(length(samples))*theta(1),'Color','red')
xlabel('sample')
title('a_1')

subplot(2,2,2)
plot(samples , theta_hat(2,:))
hold on
plot(samples , ones(length(samples))*theta(2),'Color','red')
xlabel('sample')
title('a_2')

subplot(2,2,3)
plot(samples , theta_hat(3,:))
hold on
plot(samples , ones(length(samples))*theta(3),'Color','red')
xlabel('sample')
title('a_3')


figure()
subplot(2,2,1)
plot(samples , theta_hat(4,:))
hold on
plot(samples , ones(length(samples))*theta(4),'Color','red')
xlabel('sample')
title('b_1')

subplot(2,2,2)
plot(samples , theta_hat(5 , :))
hold on
plot(samples , ones(length(samples))*theta(5),'Color','red')
xlabel('sample')
title('b_2')

subplot(2,2,3)
plot(samples , theta_hat(6 , :))
hold on
plot(samples , ones(length(samples))*theta(6),'Color','red')
xlabel('sample')
title('b_3')