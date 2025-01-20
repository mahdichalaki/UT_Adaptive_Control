clear;
close all;
clc;

%% Transfer Function

z = tf('z');
Ts = 0.02;
sysd = (z-0.4)/((z-0.22)*(z-0.73)*z^3);
[numd,dend]=tfdata(sysd,'v');

%%
numd = numd(end-1:end);
% dend = dend(2:3);

B = numd;
A = dend;

n = length(A) - 1;
n1 = length(B) -1;


%% Parameter extraction for LS estimation

theta = [dend(2:end) numd]';
N = 10000;
theta_hat_init = 0.1*ones(length(theta),N);
theta_hat_init(1:3,:) = ones(3,N);
P = 100000 * eye(length(theta));

%% Diophantine

d = 5;
a1 = [A(1:3), zeros(1,d-1)]; 
na1 = length(a1)-1;
b1 = [1];
D = [1, zeros(1, na1+d-1)];


[F, G] = Diophantine(a1, b1, D);

F = F(2:end);

%% Alpha Beta

alpha = G;
beta = conv(F, B);
D1 = conv(beta, [1, zeros(1, d-3)]);

R_star = D1(1:d-3+1);
R_bar = D1(d-3+2:end);

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
    
y_star(i,:,1) = 1 + u_noise(i);
y_star(i,:,2) = (-1).^ceil(i/2000);


end

Dist = [zeros(1,1000), zeros(1,N-1000)];

i=2;

%% Solver

t = N;
y = zeros(N,1);
y_zero = 0;
y_minus1 = 0;
y_minus2 = 0;

u_zero = 0;
u_minus1 = 0;

y(1) = 0;
y(2) = 0;
y(3) = 0;

Phi = zeros(N,length(theta));


theta_hat = theta_hat_init;
y_hat = zeros(N,1);
Error = zeros(N,1);

ContrINPUT = zeros(N,1);
OUTPUT = zeros(N,1);
sumR = sum(R_star);

    
    for k = 6:N-d
        

        y(k)=[-y(k-1) -y(k-2) -y(k-3) -y(k-4) -y(k-5) ContrINPUT(k-3)+Dist(k-3) ContrINPUT(k-4)+Dist(k-4)] * theta;
        Phi(k,:)=[-y(k-1) -y(k-2) -y(k-3) -y(k-4) -y(k-5) ContrINPUT(k-3) ContrINPUT(k-4)];
                  
        ContrINPUT(k) = (y_star(k+d,:,2) - G(1:2) * [y(k) y(k-1)]'...
                    - R_bar(1:3) * [ContrINPUT(k-1) ContrINPUT(k-2) ContrINPUT(k-3)]')/sumR;

         OUTPUT(k) = y(k);
         
    end

%% Plots

figure()
plot(ContrINPUT, 'b', 'linewidth',1)
hold on
plot (y_star(:,:,i),'--r', 'linewidth',1)
xlabel ('Sample')
ylabel('Amplitude')
title('Control Effort')
grid on

figure()
hold on
plot (OUTPUT, 'b', 'linewidth',1)
plot (y_star(:,1,i),'--r', 'linewidth',1)
legend ('Output', 'Reference Input')
xlabel ('Sample')
ylabel ('Amplitude')
title('Response')
grid on
