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
theta_hat_init(1:2,:) = 0.85 * ones(2,N);
P = 100000 * eye(length(theta)+2);

%% Y - noise included

C=(z-0.25)*(z-0.3);
[numc,denc]=tfdata(C,'v');

rng(2)
noise_var = 0.001;
noise = sqrt(noise_var) * randn(N,1);
noise = noise - mean(noise);

colored_noise = zeros(N,1);
colored_noise(3:N) = numc(1)*noise(3:N)+ numc(2)*noise(2:N-1) + numc(3)*noise(1:N-2);

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

%% SR

a0=A(1); a1=A(2); a2=A(3);
b0=B(1); b1=B(2);
c0=numc(1); c1=numc(2); c2=numc(3);

V_out=zeros(N,1);
V_control=zeros(N,1);

r1 = (c2 + b0*c1/b1 - a2 - b0^2*c2/(b1^2) - 1 - b1*(c1-a1)/b0 ) / (a1 - a2*b0/b1 - b1/b0);
s0 = 1/b1 + (c1 - r1 -a1)/b0;
s1 = b0*c2/b1^2 - a2*r1/b1;

S = [s0 s1];
R = [1 r1];

%% Solver
    
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

var_y=var(OUTPUT);
mean_y=mean(OUTPUT);
var_u=var(ContrINPUT);
mean_u=mean(ContrINPUT);

sprintf('output variance is %d',var_y)
sprintf('output mean is %d',mean_y)
sprintf('control signal variance is %d',var_u)
sprintf('signal control mean is %d',mean_u)

%% Plotter

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

%% Variance

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