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

figure
step (sysc,10,'b')
hold on
step (sysd,10,'r')
grid on
legend ('Continuos', 'Discrete')

%% MDPP with Zero Cancellation
% Desired Model
wn = 2.59;
zeta = 0.5912;
z1 = -20;
z2 = -25;
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

figure
step(Csys,'b')
hold on
step(Dsys,'r')
grid on
legend('Continuos','Discrete')
stepinfo (Dsys)

% MDPP
Ao=1;
[R,S,T]=WZC(Am,Bm,Ao,A,B);

t=600;
for i=1:t
    
INPUT (i,:,1)=1;
INPUT (i,:,2)=(-1)^ceil(i/200);
Disturbance (i,:,1)=0;
Disturbance (i,:,2)=0;

end

OUTPUT1 = Sim(t, INPUT, [conv(A,R) + [0 conv(B,S)]], conv(B,T));
OUTPUT2 = Sim(t, Disturbance, [conv(A,R) + [0 conv(B,S)]], conv(B,R));
OUTPUT = OUTPUT1 + OUTPUT2;

ContrINPUT1 = Sim(t, INPUT, [conv(A,R) + [0 conv(B,S)]], conv(A,T));
ContrINPUT2 = Sim(t, Disturbance, [conv(A,R) + [0 conv(B,S)]], conv(B,S));
ContrINPUT = ContrINPUT1 - ContrINPUT2;


figure()
hold on
plot (OUTPUT (:,1,1), 'b', 'linewidth',1)
plot (INPUT (:,1,1),'--r', 'linewidth',1)
legend ('Output', 'Reference Input')
xlabel('Sample')
ylabel('Amplitude')
title('Unit Step Response')
grid on

figure()
plot(INPUT(:,1,1)-OUTPUT(:,1,1),'b','linewidth',1)
xlabel('Sample')
ylabel('Amplitude')
title('Unit Step Response Error')
grid on

figure()
plot(ContrINPUT(:,1,2), 'b', 'linewidth',1)
hold on
plot (INPUT(:,:,2),'--r', 'linewidth',1)
xlabel ('Sample')
ylabel('Amplitude')
title('Control Effort-Reference Input is Square Wave')
grid on

figure
hold on
plot (OUTPUT (:,1,2), 'b', 'linewidth',1)
plot (INPUT (:,1,2),'--r', 'linewidth',1)
legend ('Output', 'Reference Input')
xlabel ('Sample')
ylabel ('Amplitude')
title('Square Wave Response')
grid on

%% MDPP without zero cancellation

% Desired Model
wn = 2.59;
zeta = 0.5912;
z1 = -20;
z2 = -25;
p1 = -15;
kl = wn^2;
k2 = -p1/(z1*z2);

G1=tf([1 2*zeta*wn wn^2],1);
G2=zpk(p1,[],1);
G3=zpk([z1 z2],[],kl*k2);

Csys = tf(G3/(G1*G2));
Dsys = c2d(Csys,Ts,'zoh');

[numD,denD] = tfdata(Dsys,'v');
numD = numD(2:end);

beta = sum(denD)/sum(numd);
Dsys = tf(beta*numd,denD,Ts);
% Csys=d2c (Dsys);

Am = denD;
Bm = beta*numd;

% MDPP
Ao=[1 0 0]; %observer
[R,S,T] = WOZC(Am,Ao,A,B);

t=900;
for i=1:t
    
INPUT (i, :,1)=1;
INPUT (i,:,2)=(-1)^ceil(i/300);
Disturbance (i,:,1)=0;
Disturbance (i,:,2)=0;

end

OUTPUT1 = Sim(t, INPUT, [conv(A,R)+[0 conv(B,S)]],conv(B,T)); 
OUTPUT2 = Sim(t, Disturbance, [conv(A,R) + [0 conv(B,S)]], conv(B,R));
OUTPUT = OUTPUT1 + OUTPUT2;

ContrINPUT1 = Sim(t, INPUT, [conv(A,R) + [0 conv(B,S)]], conv(A,T));
ContrINPUT2 = Sim(t, Disturbance, [conv(A,R) + [0 conv(B,S)]], conv(B,S));
ContrINPUT = ContrINPUT1 - ContrINPUT2;

figure
hold on
plot (OUTPUT (:,1,1), 'b', 'linewidth',1)
plot (INPUT (:,1,1),'--r', 'linewidth',1)
legend ('Output', 'Reference Input')
xlabel ('Sample')
ylabel ('Amplitude')
title('Unit Step Response')
grid on

figure
plot (INPUT (:,1,1) -OUTPUT (:,1,1), 'linewidth', 1)
xlabel('Sample')
ylabel ('Amplitude')
title('Unit Step Response Error')
grid on

figure()
plot(ContrINPUT(:,1,2), 'b', 'linewidth',1)
hold on
plot (INPUT(:,:,2),'--r', 'linewidth',1)
xlabel ('Sample')
ylabel('Amplitude')
title('Control Effort-Reference Input is Square Wave')
grid on

figure
hold on
plot (OUTPUT (:,1,2), 'b', 'linewidth',1)
plot (INPUT (:,1,2),'--r', 'linewidth',1)
legend ('Output', 'Reference Input')
xlabel ('Sample')
ylabel ('Amplitude')
title('Square Wave Response')
grid on











