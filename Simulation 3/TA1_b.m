clear;
close all;
clc;

%% Transfer Function

z = tf('z');
Ts = 0.02;
sysd = (z-0.65)/((z-0.35)*(z-0.2));
[numd,dend]=tfdata(sysd,'v');
numd = numd(2:end);

B = numd;
A = dend;


%% Parameter extraction for LS estimation

n = 4;
theta = [dend(2:end) numd]';
N = 10000;
theta_hat_init = 0.1*ones(length(theta),N);
theta_hat_init(1:2,:) = 0.85 * ones(2,N);
P = 100000 * eye(length(theta));

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

rng(50)
noise_var = 0.001;
noise = sqrt(noise_var) * randn(N,1);
noise = noise - mean(noise);

u_noise = sqrt(noise_var) * randn(N,1);
u_noise = u_noise - mean(u_noise);

%% MDPP

t = N;

for i=1:t
    
u(i,:,1) = 1 + u_noise(i);
u(i,:,2) = 0;

Disturbance(i,:,1) = 0;
Disturbance(i,:,2) = 0;

end

i=2;

%% Solver

uf = zeros(N,1);
yf = zeros(N,1);

theta_hat = theta_hat_init;
ContrINPUT = zeros(N,1);
OUTPUT = zeros(N,1);

Phi = zeros(N,length(theta));
Phi_f = zeros(N,length(theta));

SS = zeros(N,2);
RR = zeros(N,2);
TT = zeros(N,2);

y = zeros(N,1);
y_zero = 0;
y_minus1 = 0;

u_zero = 0;
u_minus1 = 0;

y(1) = 0;
y(2) = [-y(1) -y_zero u(1,:,i) u_zero] * theta + noise(2);
y(3) = [-y(2) -y(1) u(2,:,i) u(1,:,i)] * theta + noise(3);

Phi(2,:) = [-y(1) -y_zero u(1,:,i) u_zero];
Phi(3,:) = [-y(2) -y(1) u(2,:,i) u(1,:,i)];

T = [sum(Am) 0.1];
R = 1*ones(1,2);
S = 1*ones(1,2);

V_out=zeros(N,1);
V_control=zeros(N,1);

    for k = 4:N
        
        y(k)=[-y(k-1) -y(k-2) ContrINPUT(k-1) ContrINPUT(k-2)] * theta + noise(k);
        Phi(k,:)=[-y(k-1) -y(k-2) ContrINPUT(k-1) ContrINPUT(k-2)];

        ContrINPUT(k) = T * [u(k,:,i) u(k-1,:,i)]'...
                      + S * [-y(k) -y(k-1)]'...
                      - R(2:end) * [ContrINPUT(k-1)]';
        
        ContrINPUT(k) = ContrINPUT(k)/R(1);
        
        uf(k) = ContrINPUT(k) - [uf(k-1) uf(k-2)]*Am(2:end)';
        yf(k) = y(k) - [yf(k-1) yf(k-2)]*Am(2:end)';
        
        Phi_f(k,:) = [uf(k-1) uf(k-2) yf(k-1) yf(k-2)];
       
        P = P - P*Phi_f(k,:)'*((1+Phi_f(k,:)*P*Phi_f(k,:)')\Phi_f(k,:))*P;
        K = P * Phi_f(k,:)';
        theta_hat(:,k) = theta_hat(:,k-1) + K*(y(k)-Phi_f(k,:)*theta_hat(:,k-1));

        
        R = theta_hat(1:2,k)';
        S = theta_hat(3:end,k)';
        T = [sum(Am) 0];
        
        OUTPUT(k) = y(k);
        
        V_out(k) = V_out(k-1) + OUTPUT(k)^2;
        V_control(k) = V_control(k-1) + ContrINPUT(k)^2;
     
        SS(k,:) = S;
        TT(k,:) = T;
        RR(k,:) = R;
        
    end

var_y=var(OUTPUT);
mean_y=mean(OUTPUT);
var_u=var(ContrINPUT);
mean_u=mean(ContrINPUT);

sprintf('output variance is %d',var_y)
sprintf('output mean is %d',mean_y)
sprintf('control signal variance is %d',var_u)
sprintf('signal control mean is %d',mean_u)

%% Plots   

figure()
plot(ContrINPUT(:,1,:), 'b', 'linewidth',1)
hold on
plot (u(:,:,i),'--r', 'linewidth',1)
xlabel ('Sample')
ylabel('Amplitude')
title('Control Effort')
grid on

figure()
hold on
plot (OUTPUT (:,1,:), 'b', 'linewidth',1)
plot (u (:,1,i),'--r', 'linewidth',1)
legend ('Output', 'Reference Input')
xlabel ('Sample')
ylabel ('Amplitude')
title('Response')
grid on

%% Loss

x = (linspace(1,N,N))';
limx = N; %for plot xlim

figure()
subplot(1,2,1);
plot(x,V_out)
xlabel('sample number')
ylabel('loss')
title('Acc loss - Output')

subplot(1,2,2);
plot(x,V_control)
xlabel('sample number')
ylabel('loss')
title('Acc loss - Input')

%% S

figure()
subplot(2,1,1)
plot (SS (:,1), 'b', 'linewidth',1)
hold on
plot (-1.388*ones(N,1),'--r', 'linewidth',1)
legend ('Estimated', 'Actual')
xlabel ('Sample')
ylabel ('Amplitude')
title('S_0')
grid on

subplot(2,1,2)
plot (SS (:,2), 'b', 'linewidth',1)
hold on
plot (0.8708*ones(N,1),'--r', 'linewidth',1)
legend ('Estimated', 'Actual')
xlabel ('Sample')
ylabel ('Amplitude')
title('S_1')
grid on
    
%% R
    
figure()
subplot(2,1,1)
plot (RR (:,1), 'b', 'linewidth',1)
hold on
plot (1*ones(N,1),'--r', 'linewidth',1)
legend ('Estimated', 'Actual')
xlabel ('Sample')
ylabel ('Amplitude')
title('R_0')
grid on

subplot(2,1,2)
plot (RR (:,2), 'b', 'linewidth',1)
hold on
plot (-0.6508*ones(N,1),'--r', 'linewidth',1)
legend ('Estimated', 'Actual')
xlabel ('Sample')
ylabel ('Amplitude')
title('R_1')
grid on

    
%% T

figure()
subplot(2,1,1)
plot (TT (:,1), 'b', 'linewidth',1)
hold on
plot (0.0078*ones(N,1),'--r', 'linewidth',1)
legend ('Estimated', 'Actual')
xlabel ('Sample')
ylabel ('Amplitude')
title('T_0')
grid on

subplot(2,1,2)
plot (TT (:,2), 'b', 'linewidth',1)
hold on
plot (-0.0052*ones(N,1),'--r', 'linewidth',1)
legend ('Estimated', 'Actual')
xlabel ('Sample')
ylabel ('Amplitude')
title('T_1')
grid on
