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

zeta=0.999;
wn = 4/(zeta*settling_time); 

% sys_2d_desired  = tf([wn^2], [1, 2*zeta*wn, wn^2]);
% damp(sys_desired)

p3 = -55;

k2 = -p3;

G1 = tf([wn^2],[1, 2*zeta*wn, wn^2]);
G2 = zpk([z1,z2],[p3],k2);

sys_desC = G1*G2;
sys_desD = c2d(sys_desC, Ts, 'zoh');
[num_desD, den_desD] = tfdata(sys_desD, 'v');

Bm_prim2 = num_desD; 
Am = den_desD;

%% without pole-zero cancellation

syms q;

% Bm = nonzeros(Bm);
% Am = nonzeros(Am);
% B = nonzeros(B);
% A = nonzeros(A);

Bplus = [1];
Bminus = B;
Ao = [1,0,0];
Bm = conv(Bplus, Bminus);
BmPrime = (sum(Am)/sum(Bm)); 

[R_solved,S_solved,T] = solve_diophantin(Am, A, Ao, Bminus, Bplus, BmPrime);
Ac_solved = conv(A,R_solved) + conv(B,S_solved);

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
