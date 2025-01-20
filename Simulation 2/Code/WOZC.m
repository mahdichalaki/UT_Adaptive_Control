function [R, S, T] = WOZC(Am, Ao, A, B)

a1 = A(2);
a2 = A(3);
a3 = A(4);

b0 = B(1);
b1 = B(2);
b2 = B(3);

am1 = Am(2);
am2 = Am(3);
am3 = Am(4);

ao1 = Ao(2);
ao2 = Ao(3);

R=[];
T=[];
S=[];

H = [1 0 b0 0 0;...
    a1 1 b1 b0 0;...
    a2 a1 b2 b1 b0;...
    a3 a2 0 b2 b1;...
    0 a3 0 0 b2];

h = [am1+ao1-a1;...
    am2+ao2+ao1*am1-a2;...
    am3+ao1*am2+ao2*am1-a3;...
    am3*ao1+am2*ao2;...
    am3*ao2];

RS = (inv(H)*h)';
R = [1 RS(1:2)];
S = RS(3:5);

beta=(1+am1+am2+am3)/(b0+b1+b2);
T = beta*Ao;

end



