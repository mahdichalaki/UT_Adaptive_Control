clear;
close all;
clc;

%% Transfer Function

a = 9;
b = 0;

c1 = (b+2)/10;
c2 = (a+2)/10;
R1 = (10*a+b)*0.01;
R2 = (10*a+b+20)*0.01;

tau1 = 63.85;
tau2 = 1048.2575;
A = 1;

b = R2/(tau1*tau2);
a0 = 1/(tau1*tau2);
a1 = (tau1+tau2+A*R2)/(tau1*tau2);

s = tf('s');
G = b/(s^2 + a1*s + a0);

bm = 10;
a0m = 10;
a1m = 14;

Gm = bm/(s^2 + a1m*s + a0m);

gamma = 3000;
alpha = 0.001;



