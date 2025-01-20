clear;
close all;
clc;
%% System
a1 = 3.5;
a2 = 1.5;
b = 3;
s = tf('s');
sysc = 3/((s+0.5)*(s+3));

%% Desired System
wn = 6;
zeta = 0.65;
am1 = 2*zeta*wn;
am2 = wn^2;
bm = wn^2;
sys = tf(bm, [1 am1 am2]);
[y,t] = step(sys,5);

%% Continuous STR
a0=30;
P0=1e6*eye(3);
tet0=[2;0;0.5];
theta1 = [a1; a2; b];