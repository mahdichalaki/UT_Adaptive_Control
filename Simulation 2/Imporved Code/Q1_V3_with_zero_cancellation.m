%% discretization of the model

sysC  = zpk([-30, -0.6], [-1.3,-2.5,1.1], 0.2);
BW = bandwidth(sysC); 

disc_ratio = 10;
fs = BW * disc_ratio/(2*pi);
Ts = 1/fs;

sysD = c2d(sysC, Ts, 'zoh');
[numD, denD] = tfdata(sysD, 'v');

% bodeplot(sysC);
% hold on;
% bodeplot(sysD);

B = numD;
A = denD;

%% desired system
overshoot = 10;
settling_time = 3;

zeta = cos(atan2(pi,-log(overshoot/100)));
wn = 4/(zeta*settling_time); 

% sys_2d_desired  = tf([wn^2], [1, 2*zeta*wn, wn^2]);
% damp(sys_desired)
z1 = -15;
z2 = -25;
p3 = -20;

k2 = -p3/(z1*z2);

G1 = tf([wn^2],[1, 2*zeta*wn, wn^2]);
G2 = zpk([z1,z2],[p3],k2);

sys_desC = G1*G2;
sys_desD = c2d(sys_desC, Ts, 'zoh');
[num_desD, den_desD] = tfdata(sys_desD, 'v');

BmPrime = num_desD; 
Am = den_desD;
%% without pole-zero cancellation

syms q;

% Bm = nonzeros(Bm);
% Am = nonzeros(Am);
if B(1) == 0
    B = B(2:end);
end
% A = nonzeros(A);

Bplus = [B/B(1)];
Bminus = B(1);
BmPrime = BmPrime/B(1);
Ao = [1];

[R_solved,S_solved,T] =  solve_diophantin(Am, A, Ao, Bminus, Bplus, BmPrime);
Ac_solved = conv(A,R_solved) + [0,conv(B,S_solved)];
%% step
% num_samples = 100;
% step_mag = 1;
% t = 0:num_samples-1;
% uc = step_mag * t>=0;
% 
% y = filter(conv(B,T),Ac_solved,uc);
% u = filter(conv(A,T),Ac_solved,uc);
% 
% figure()
% plot(t,y)
% figure()
% plot(t, u)

%%
num_samples = 1000;
step_mag = 1;
t = 0:num_samples-1;
uc = (-1).^ceil(t/200);

y = filter(conv(B,T),Ac_solved,uc);
u = filter(conv(A,T),Ac_solved,uc);

figure()
plot(t,y)
hold on;
plot(t, uc, "--r")
title("output")
figure()
plot(t, u)
title("control input")

% sys_new = tf([conv(B,T)],conv(A,T));
