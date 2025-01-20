clear;
close all;
clc;

%% Transfer Function

s = tf('s');
sysc = (0.45*(s+1.9)*(s+3.4))/((0.55*s+4)*(0.2*s+1.4)*(s-1.35));
bw = bandwidth(sysc);
Ts = (2*pi)/(20*bw);

sysd = c2d(sysc, Ts,'zoh');
[numd,dend]=tfdata(sysd,'v');
numd = numd(2:end);

B = numd;
A = dend;


%% Parameter extraction for LS estimation

n = 6;
theta = [dend(2:end) numd]';
N = 10000;
theta_hat_init = 0.1*ones(length(theta),N);
theta_hat_init(1:3,:) = ones(3,N);
P = 100000 * eye(length(theta));

%% MDPP with zero cancellation

% Desired Model
wn = 2.59;
zeta = 0.5912;
z1 = -10;
z2 = -11;
p1 = -15;
k1 = wn^2;
k2 = -p1/(z1*z2);

G1 = tf([1 2*zeta*wn wn^2], 1);
G2 = zpk (p1,[],1);
G3 = zpk ([z1 z2],[], k1*k2);

Csys = tf(G3/(G1*G2));
Dsys = c2d(Csys,Ts,'zoh');

[numD,denD] = tfdata(Dsys,'v');
numD = numD(2:end);

Am = denD;
Bm = numD;
Ao = 1; %observer

%% Y - noise included

noise_var = 0.001;
noise = sqrt(noise_var) * randn(N,1);
noise = noise - mean(noise);

u_noise = sqrt(noise_var) * randn(N,1);

%% MDPP

t = N;

for i=1:t
    
u(i,:,1) = 1 + u_noise(i);
% u(i,:,1) = 1;
u(i,:,2) = (-1)^ceil(i/3000) + u_noise(i);
% u(i,:,2) = (-1)^ceil(i/3000);

end

Disturbance = zeros(N,1);
Disturbance(40:end) = 0.1;

%%
i=2;
%%
% [theta_hat, Cost] = RLS(i,u,theta,theta_hat_init,noise,N,P,Disturbance,A,B,Am,Ao);

%% RLS

% function [theta_hat, Cost] = RLS(i,u,theta,theta_hat_init,noise,N,P,Disturbance,A,B,Am,Ao)
    
    t = N;
    y = zeros(N,1);
    y_zero = 0;
    y_minus1 = 0;

    u_zero = 0;
    u_minus1 = 0;
    
    y(1) = 0;
    y(2) = [-y(1) -y_zero -y_minus1 u(1,:,i) u_zero u_minus1] * theta + noise(2);
    y(3) = [-y(2) -y(1) -y_zero u(2,:,i) u(1,:,i) u_zero] * theta + noise(3);
    
    Phi = zeros(N,length(theta));
    Phi(2,:) = [-y(1) -y_zero -y_minus1 u(1,:,i) u_zero u_minus1];
    Phi(3,:) = [-y(2) -y(1) -y_zero u(2,:,i) u(1,:,i) u_zero];

    theta_hat = theta_hat_init;
    y_hat = zeros(N,1);
    Error = zeros(N,1);
    
    ContrINPUT = zeros(N,1);
    OUTPUT = zeros(N,1);
    
    SS = zeros(N,3);
    RR = zeros(N,3);
    TT = zeros(N,3);
    
    for k = 4:N

        y(k)=[-y(k-1) -y(k-2) -y(k-3) ContrINPUT(k-1) ContrINPUT(k-2) ContrINPUT(k-3)] * theta;
        Phi(k,:)=[-y(k-1) -y(k-2) -y(k-3) ContrINPUT(k-1) ContrINPUT(k-2) ContrINPUT(k-3)];
        P = P - P*Phi(k,:)'*((1+Phi(k,:)*P*Phi(k,:)')\Phi(k,:))*P;
        K = P * Phi(k,:)';
        theta_hat(:,k) = theta_hat(:,k-1) + K*(y(k)-Phi(k,:)*theta_hat(:,k-1));
        y_hat(k) = Phi(k,:)*theta_hat(:,k);
        Error(k) = y(k) - y_hat(k);
        
        A_hat = [1 theta_hat(1:3,k)'];
        B_hat = [theta_hat(4:end,k)]';
        [R,S,T] = WZC(Am,Bm,Ao,A_hat,B_hat);
        

        ContrINPUT(k) = T * [u(k,:,i) u(k-1,:,i) u(k-2,:,i)]'...
                      + S * [-y(k) -y(k-1) -y(k-2)]'...
                      - R(2:end) * [ContrINPUT(k-1) ContrINPUT(k-2)]'...
                      + Disturbance(k);


         OUTPUT(k) = y(k);
         SS(k,:) = S;
         TT(k,:) = T;
         RR(k,:) = R;
    end
    
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
    
    
    %plotting
    x = (linspace(1,N,N))';
    Cost = 1/2 * (Error' * Error);
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
    
    subplot(2,2,3);
    plot(x,theta(3)*ones(1,N))
    hold on
    plot(x,theta_hat(3,:))
    title('a_3')
    xlim([0 limx])
    
    figure()
    subplot(2,2,1);
    plot(x,theta(4)*ones(1,N))
    hold on
    plot(x,theta_hat(4,:))
    title('b_1')
    xlim([0 limx])
    
    
    subplot(2,2,2);
    plot(x,theta(5)*ones(1,N))
    hold on
    plot(x,theta_hat(5,:))
    title('b_2')
    xlim([0 limx])
    
    subplot(2,2,3);
    plot(x,theta(6)*ones(1,N))
    hold on
    plot(x,theta_hat(6,:))
    title('b_3')
    xlim([0 limx])
    
    
% end