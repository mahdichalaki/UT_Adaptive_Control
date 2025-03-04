%% discretization of the model
clear;
sysC  = zpk([-30, -0.6], [-1.3,-2.5,1.1], 0.2);
BW = bandwidth(sysC); 

disc_ratio = 20;
fs = BW * disc_ratio/(2*pi);
Ts = 1/fs;

sysD = c2d(sysC, Ts, 'zoh');
[numD, denD] = tfdata(sysD, 'v');

% bodeplot(sysC);
% hold on;
% bodeplot(sysD);

B = numD;
A = denD;
%%
if B(1) == 0
    B = B(2:end);
end
%% desired system
overshoot = 10;
settling_time = 3;

zeta = cos(atan2(pi,-log(overshoot/100)));
wn = 4/(zeta*settling_time); 

zeta=0.999;
wn = 4/(zeta*settling_time); 

% sys_2d_desired  = tf([wn^2], [1, 2*zeta*wn, wn^2]);
% damp(sys_desired)
z1 = -30;
z2 = -0.6;
p3 = -55;

k2 = -p3;

G1 = tf([wn^2],[1, 2*zeta*wn, wn^2]);
G2 = zpk([z1,z2],[p3],k2);

sys_desC = G1*G2;
sys_desD = c2d(sys_desC, Ts, 'zoh');
[num_desD, den_desD] = tfdata(sys_desD, 'v');

Bm_prim2 = num_desD; 
Am = den_desD;

%% without zero cancellation, with noise
main_title = "2-1-with-noise";

syms q;
% input properties
num_samples = 500;
step_mag = 1;
t = 0:num_samples-1;
uc = (-1).^-ceil(t/200);
uc(1) = -1;

% must be changed for every problem
deg_desA = 4;
deg_desB = 3;
A_estimated = [1,0,0,0]; % initial conditions
B_estimated = [0.01, 0.035,0.2];

assert(length(A_estimated)==deg_desA && length(B_estimated)==deg_desB, "initial condistions are not true")

% initial B
Bplus = [1];
Bminus = B_estimated;
Bm = conv(Bplus, Bminus);
BmPrime = (sum(Am)/sum(Bm)); 


% initial R S T and Ao
R_solved = [1,-0.1,-0.1];
S_solved = [1,0.1,0.1];
T = [1,0,0];
Ao = [1,0,0];

% R S and T which were calculated in the last part
R_real = [1,-0.7386,-0.5508];
S_real = [6.7933,-7.7807,1.9653];
T_real = [1.5652,0,0];

% initial conditions for y
skip_instances = max(deg_desA, deg_desB);
total_parameters = deg_desA - 1 + deg_desB;
y = [];
y(1:skip_instances) = 0;
u(1:skip_instances) = 0;


theta_real = [A(2:end).'; B.'];

% noise and its degree
noise_poly = [1,0.1,0.1].';
deg_noise = length(noise_poly);
noise_variance = 0.001; % system noise variance
noise = sqrt(noise_variance) * randn(1, num_samples);

R_solved_toPlot = zeros([num_samples, length(R_solved)]);
S_solved_toPlot = zeros([num_samples, length(S_solved)]);
T_toPlot = zeros([num_samples, length(T)]);
theta_hat_toPlot = zeros([num_samples, total_parameters]);

%
theta_epsilon_zero = [1,0,0].';
els_solver = ELSClass(100 * eye(total_parameters + length(theta_epsilon_zero)), 0.1 * ones([total_parameters,1]), theta_epsilon_zero);

for i = skip_instances:length(uc)
    phi_t = [-y(i-1:-1:i-(deg_desA - 1)), u(i-1:-1:i-deg_desB)].';   

    noise_t = [noise(i:-1:i-(deg_noise-1))] * noise_poly; %only for simulation. not involved in controller calculations directly

    y(i) = phi_t.' * theta_real + noise_t;

    u(i) = T * [uc(i:-1:i-(length(T)-1))].' + S_solved * [-y(i:-1:i-(length(S_solved)-1))].' - R_solved(2:end) * [u(i-1:-1:i-(length(R_solved)-1))].';

    theta_hat_new = els_solver.update_ELS(y(i), phi_t);


    A_estimated = theta_hat_new(1:(deg_desA - 1)).';
    B_estimated = theta_hat_new(deg_desA:total_parameters).';
    
    Bplus = [1];
    Bminus = B_estimated;
    Bm = conv(Bplus, Bminus);
    BmPrime = (sum(Am)/sum(Bm)); 
    
    [R_solved,S_solved,T] = solve_diophantin(Am, [1,A_estimated], Ao, Bminus, Bplus, BmPrime);

    theta_hat_toPlot(i, :) = theta_hat_new(1:total_parameters).';
    R_solved_toPlot(i, :) = R_solved;
    S_solved_toPlot(i, :) = S_solved;
    T_toPlot(i, :) = T;
%     disp(R_solved)
end
%% plotters

figure()
plot(t,y,'b', 'linewidth',1)
hold on;
plot(t, uc, "--r",'linewidth',1)
xlabel("Sample");
ylabel("Amplitude")
title("Square Wave Response");
legend('Output', 'Reference Input');

figure()
plot(t, u,'b', 'linewidth',1)
hold on;
plot(t, uc, "--r",'linewidth',1)
xlabel("Sample");
ylabel("Amplitude")
title("Control Effort-Reference Input is Square Wave");


%% R 
total_plot_rows = ceil(length(R_solved)/2) + ceil(length(S_solved)/2)+ ceil(length(T)/2) - 1;
f2 = figure();
for i = 2:(length(R_solved)) % first parameter of R is 1
    subplot(1,2,i-1);
    title_text = "R_%d";
    plot(1:num_samples, R_solved_toPlot(:,i),'b', 'linewidth',1) 
%     hold on;
%     plot(1:num_samples, ones([num_samples,1]) * R_real(i), 'DisplayName','Real') 
    title(sprintf(title_text, i-1));
    legend('Estimated');
    xlabel("Sample")
    ylabel("Amplitude")
    xlim([0 500])
    
end
end_of_r_plot = i-1;

%% S
figure()
for i = 1:(length(S_solved)) % first parameter of R is 1
    subplot(2,2,i);
    title_text = "S_%d";
    plot(1:num_samples, S_solved_toPlot(:,i), 'b', 'linewidth',1) 
%     hold on;
%     plot(1:num_samples, ones([num_samples,1]) * S_real(i), 'DisplayName','Real') 
title(sprintf(title_text, i-1));
    legend('Estimated');
    xlabel("Sample")
    ylabel("Amplitude")
    xlim([0 500])
end

%% T
for i = 1:(length(T)) % first parameter of R is 1
    subplot(2,2,i);
    title_text = "T_%d";
    plot(1:num_samples, T_toPlot(:,i), 'b', 'linewidth',1) 
%     hold on;
%     plot(1:num_samples, ones([num_samples,1]) * T_real(i), 'DisplayName','Real') 
    title(sprintf(title_text, i));
    legend('Estimated');
    xlabel("Sample")
    ylabel("Amplitude")
    xlim([0 500])
end


%% Theta
% f5 = figure();
% for i = 1:length(theta_hat_toPlot(1,:))
%     title_text = "θ_%d";
%     subplot(ceil(length(theta_hat_toPlot(1,:)))/2,2,i);
%     plot(1:num_samples, theta_hat_toPlot(:,i))
%     hold on;
%     plot(1:num_samples, ones([num_samples,1]) * theta_real(i))
%     title(sprintf(title_text, i));
%     legend('Estimated','Actual');
%     xlabel("Sample")
%     ylabel("Amplitude")
%     xlim([0 500])
% end
% Theta
f5 = figure();
for i = 1:3
    title_text = "a_%d";
    subplot(2,2,i);
    plot(1:num_samples, theta_hat_toPlot(:,i), 'DisplayName','Estimated')
    hold on;
    plot(1:num_samples, ones([num_samples,1]) * theta_real(i) , 'DisplayName','Actual')
    title(sprintf(title_text, i));
    legend()
    xlim([0 500])
    xlabel("Sample")
%     legend('Location','best');
%     xlabel("sample number")
end
figure()
for i = 4:6
    title_text = "b_%d";
    subplot(2,2,i-3);
    plot(1:num_samples, theta_hat_toPlot(:,i-3), 'DisplayName','Estimated')
    hold on;
    plot(1:num_samples, ones([num_samples,1]) * theta_real(i-3) , 'DisplayName','Actual')
    title(sprintf(title_text, i-3));
    legend()
    xlim([0 500])
    xlabel("Sample")
%     legend('Location','best');
%     xlabel("sample number")
end




% sys_new = tf([conv(B,T)],conv(A,T));
