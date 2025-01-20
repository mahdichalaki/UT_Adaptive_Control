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

t = N;
y = zeros(N,1);

Phi = zeros(N,length(theta));
theta_hat = theta_hat_init;
y_hat = zeros(N,1);

ContrINPUT = zeros(N,1);
OUTPUT = zeros(N,1);

%% Y - noise included

rng(2)
noise_var = 0.001;
noise = sqrt(noise_var) * randn(N,1);
noise = noise - mean(noise);

u_noise = sqrt(noise_var) * randn(N,1);
u_noise = u_noise - mean(u_noise);

%% MDPP

for i=1:t
    
y_star(i,:,1) = 1 + u_noise(i);
y_star(i,:,2) = (-1).^ceil(i/2000);

Disturbance(i,:,1) = 0;
Disturbance(i,:,2) = 0;

end

i=2;

%% Diophantine

d = 3;
% a1 = A; 
b1 = [1];
D = [1, zeros(1, n+d-1)];

%% Solver

    
    for k = 6:N-d
        
        y(k)=[-y(k-1) -y(k-2) -y(k-3) -y(k-4) -y(k-5) ContrINPUT(k-3) ContrINPUT(k-4)] * theta;
        Phi(k,:)=[-y(k-1) -y(k-2) -y(k-3) -y(k-4) -y(k-5) ContrINPUT(k-3) ContrINPUT(k-4)];
        
        P = P - P*Phi(k,:)'*inv(1+Phi(k,:)*P*Phi(k,:)')*Phi(k,:)*P;
        K = P * Phi(k,:)';
        theta_hat(:,k) = theta_hat(:,k-1) + K*(y(k)-Phi(k,:)*theta_hat(:,k-1));
        y_hat(k) = Phi(k,:)*theta_hat(:,k);
        
        a1_hat = [1, theta_hat(1:5,k)'] ;
        [F, G] = Diophantine(a1_hat, b1, D);
        F = F(end-2:end);
        alpha = G;
        B_hat = theta_hat(6:7,k)';
        beta = conv(F, B_hat);
        
        
                  
        ContrINPUT(k) = y_star(k+d,:,2) - alpha * [y(k) y(k-1) y(k-2) y(k-3) y(k-4)]'...
                    - beta(2:end) * [ContrINPUT(k-1) ContrINPUT(k-2) ContrINPUT(k-3)]';
                
        if ContrINPUT(k)>2
            ContrINPUT(k) = 2;
        elseif ContrINPUT(k)<-2
            ContrINPUT(k) = -2;
        end

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

%% plotting

x = (linspace(1,N,N))';
limx = N; %for plot xlim

figure()
for i = 1:length(dend)-1
    subplot(2,3,i);
    title_text = "a_%d";
    plot(x,theta(i)*ones(1,N)) 
    hold on;
    plot(x,theta_hat(i,:))
    title(sprintf(title_text, i));
    xlim([0 limx])
    xlabel('Samples')
    legend('Real','Estimated');
    
end

figure()
for i = length(dend):length(numd)+length(dend)-1
    subplot(2,3,i-5);
    title_text = "b_%d";
    plot(x,theta(i)*ones(1,N)) 
    hold on;
    plot(x,theta_hat(i,:))
    title(sprintf(title_text, i-6));
    xlim([0 limx])
    xlabel('Samples')
    legend('Real','Estimated');
    
end
