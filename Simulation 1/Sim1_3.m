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
N = 1e6;
theta_hat_init = zeros(length(theta),N);
P = 100000 * eye(length(theta));

%% From Simulink 

u_noisy = out.u_noisy(2:N+1);
y_noisy = out.y_noisy(1:N);

%% Y - noisy input

[theta_hat_noisy,cost_noisy] = RLS(u_noisy,y_noisy,theta,theta_hat_init,N,P);


%% Function

function [theta_hat, Cost] = RLS(u,y,theta,theta_hat_init,N,P)

    y_zero = 0;
    y_minus1 = 0;
    y_minus2 = 0;

    u_zero = 0;
    u_minus1 = 0;
    u_minus2 = 0;
    
    
    Phi = zeros(N,length(theta));
    Phi(2,:) = [-y(1) -y_zero -y_minus1 u(1) u_zero u_minus1];
    Phi(3,:) = [-y(2) -y(1) -y_zero u(2) u(1) u_zero];
    Phi(4,:) = [-y(3) -y(2) -y(1) u(3) u(2) u(1)];

    theta_hat = theta_hat_init;
    y_hat = zeros(N,1);
    Error = zeros(N,1);
    
    for k = 5:N
        
        Phi(k,:)=[-y(k-1) -y(k-2) -y(k-3) u(k-1) u(k-2) u(k-3)];
        P = P - P*Phi(k,:)'*inv(1+Phi(k,:)*P*Phi(k,:)')*Phi(k,:)*P;
        K = P * Phi(k,:)';
        theta_hat(:,k) = theta_hat(:,k-1) + K*(y(k)-Phi(k,:)*theta_hat(:,k-1));
        y_hat(k) = Phi(k,:)*theta_hat(:,k);
        Error(k) = y(k) - y_hat(k);
        
    end
    
    
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
    
end

