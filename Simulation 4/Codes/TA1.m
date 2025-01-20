clear;
close all;
clc;

%% Transfer Function

z = tf('z');
Ts = 0.02;
sysd = (z-0.4)/((z-0.22)*(z-0.73)*z^3);
[numd,dend]=tfdata(sysd,'v');
numd = numd(end-1:end);

B = numd;
A = dend;

N = 10000;

%% MDPP with zero cancellation

% Desired Model
wn = 2.59;
zeta = 0.5912;
z1 = -20;
p3 = -25;
p4 = -27;
p5 = -30;
k1 = wn^2;
k2 = (p3*p4*p5)/(z1);

G1 = tf([1 2*zeta*wn wn^2], 1);
G2 = zpk (z1,[p3, p4, p5], k1*k2);

Csys = tf(G2/(G1));
Dsys = c2d(Csys,Ts,'zoh');

[numD,denD] = tfdata(Dsys,'v');

Am = denD;
Bm = numD;
Ao = [1, 0, 0, 0]; %observer

Bplus = B/B(1);
Bminus = B(1);
BmPrime = Bm/B(1);

[R,S,T] =  solve_diophantin(Am, A, Ao, Bminus, Bplus, BmPrime);
Ac = conv(A,R) + [zeros(1, length(Ao)), conv(B,S)]; % depends on degree

%% Y - noise included

rng(2)
noise_var = 0.001;
noise = sqrt(noise_var) * randn(N,1);
noise = noise - mean(noise);

u_noise = sqrt(noise_var) * randn(N,1);
u_noise = u_noise - mean(u_noise);

%% MDPP

t = N;

for i=1:t
    
u(i,:,1) = 1 + u_noise(i);
u(i,:,2) = (-1).^ceil(i/2000);

Disturbance(i,:,1) = 0;
Disturbance(i,:,2) = 0;

end

i=2;

%% Solver

y = filter(conv(B,T),Ac,u(:,:,2));
Control_input = filter(conv(A,T),Ac,u(:,:,2));

%% Plots

figure()
plot(Control_input, 'b', 'linewidth',1)
hold on
plot (u(:,:,i),'--r', 'linewidth',1)
xlabel ('Sample')
ylabel('Amplitude')
title('Control Effort')
grid on

figure()
hold on
plot (y, 'b', 'linewidth',1)
plot (u (:,1,i),'--r', 'linewidth',1)
legend ('Output', 'Reference Input')
xlabel ('Sample')
ylabel ('Amplitude')
title('Response')
grid on
