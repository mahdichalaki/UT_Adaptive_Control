clear;
close all;
clc;

%% Transfer Function

z = tf('z');
Ts = 0.02;
sysd = (z-1/0.65)/((z-0.35)*(z-0.2));
[numd,dend]=tfdata(sysd,'v');
numd = numd(2:end);

B = numd;
A = dend;


%% Parameter extraction for LS estimation

n = 6;
theta = [dend(2:end) numd]';
N = 10000;
theta_hat_init = 0.1*ones(length(theta)+2,N);
theta_hat_init(1:2,:) = ones(2,N);
P = 100000 * eye(length(theta)+2);

%% MDPP with zero cancellation

% Desired Model
wn = 2.59;
zeta = 0.5912;
z1 = -20;
k1 = wn^2;
k2 = -1/(z1);

G1 = tf([1 2*zeta*wn wn^2], 1);
G2 = zpk (z1,[], k1*k2);

Csys = tf(G2/(G1));
Dsys = c2d(Csys,Ts,'zoh');

[numD,denD] = tfdata(Dsys,'v');
numD = numD(2:end);

Am = denD;
Bm = numD;
Ao = 1; %observer

%% Y - noise included

C=(z-0.25)*(z-0.3);
[numc,denc]=tfdata(C,'v');
% numc = numc(2:end);


noise_var = 0.001;
noise = sqrt(noise_var) * randn(N,1);
noise = noise - mean(noise);

colored_noise = zeros(N,1);
% colored_noise(2:N) = 0.9*noise(2:N)+ 0.1 * noise(1:N-1);
colored_noise(3:N) = numc(1)*noise(3:N)+ numc(2)*noise(2:N-1) + numc(3)*noise(1:N-2);

u_noise = sqrt(noise_var) * randn(N,1);
u_noise = u_noise - mean(u_noise);

%% MDPP

t = N;

for i=1:t
    
u(i,:,1) = 1 + u_noise(i);
% u(i,:,1) = 1;
u(i,:,2) = (-1)^ceil(i/3000) + u_noise(i);
% u(i,:,2) = (-1)^ceil(i/3000);

Disturbance(i,:,1) = 0;
Disturbance(i,:,2) = 0;

end

%%
i=2;
%%
% [theta_hat, Cost] = RLS(i,u,theta,theta_hat_init,noise,N,P,Disturbance,A,B,Am,Ao);

%% SR

a0=A(1); a1=A(2); a2=A(3);
b0=B(1); b1=B(2);
c0=numc(1); c1=numc(2); c2=numc(3);



V_out=zeros(N,1);
V_control=zeros(N,1);
% [nums,dens]=tfdata(S,'v');
% Scoef=cell2mat(tfdata(S));
% s1=Scoef(1);
% s0=Scoef(2);
% [numr,denr]=tfdata(R,'v');
% numr = numr(2:end);
% Rcoef=cell2mat(tfdata(R));
% r0=Rcoef(2);


%% RLS

% function [theta_hat, Cost] = RLS(i,u,theta,theta_hat_init,noise,N,P,Disturbance,A,B,Am,Ao)
    
t = N;
y = zeros(N,1);
y_zero = 0;
y_minus1 = 0;

u_zero = 0;
u_minus1 = 0;

y(1) = 0;
y(2) = [-y(1) -y_zero u(1,:,i) u_zero] * theta + colored_noise(2);
y(3) = [-y(2) -y(1) u(2,:,i) u(1,:,i)] * theta + colored_noise(3);


y_hat = zeros(N,1);
Error = zeros(N,1);

ContrINPUT = zeros(N,1);
OUTPUT = zeros(N,1);
epsilon = zeros(N,2);
EE = zeros(N,2);

SS = zeros(N,2);
RR = zeros(N,2);

for k = 4:N
    
    y(k)=[-y(k-1) -y(k-2) ContrINPUT(k-1) ContrINPUT(k-2)] * theta + colored_noise(k);
    ContrINPUT(k) = S * [-y(k) -y(k-1)]' - R(2) * [ContrINPUT(k-1)]';
    OUTPUT(k) = y(k);
    
    V_out(k) = V_out(k-1) + OUTPUT(k)^2;
    V_control(k) = V_control(k-1) + ContrINPUT(k)^2;

end

%% Plotter
figure()
plot(ContrINPUT(:,1,:), 'b', 'linewidth',1)
hold on
plot (u(:,:,i),'--r', 'linewidth',1)
xlabel ('Sample')
ylabel('Amplitude')
title('Control Effort-Reference Input is Square Wave')
grid on

figure()
hold on
plot (OUTPUT (:,1,:), 'b', 'linewidth',1)
plot (u (:,1,i),'--r', 'linewidth',1)
legend ('Output', 'Reference Input')
xlabel ('Sample')
ylabel ('Amplitude')
title('Square Wave Response')
grid on

%% Variance

x = (linspace(1,N,N))';

figure()
plot(x,V_out)
xlabel('sample number')
ylabel('Variance')
title('Output_var')

figure()
plot(x,V_control)
xlabel('sample number')
ylabel('Variance')
title('input_var')

figure()
plot(x,OUTPUT)
xlabel('sample number')
ylabel('Magnitude')
title('Output')

figure()
plot(x,ContrINPUT)
xlabel('sample number')
ylabel('Magnitude')
title('Control input')


%% Hats
%plotting

limx = N; %for plot xlim
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
plot(x,theta(1)*ones(1,N))
hold on
plot(x,theta_hat(1,:))
title('a_1')
xlim([0 limx])

subplot(2,2,2);
plot(x,theta(2)*ones(1,N))
hold on
plot(x,theta_hat(2,:))
title('a_2')
xlim([0 limx])


figure()
subplot(2,2,1);
plot(x,theta(3)*ones(1,N))
hold on
plot(x,theta_hat(3,:))
title('b_1')
xlim([0 limx])


subplot(2,2,2);
plot(x,theta(4)*ones(1,N))
hold on
plot(x,theta_hat(4,:))
title('b_2')
xlim([0 limx])



% end