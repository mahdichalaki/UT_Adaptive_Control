function [R,S,T] = WZC(Am, Bm, Ao, A, B)

a1 = A(2);
a2 = A(3);
a3 = A(4);

b0 = B(1);
b1 = B(2);
b2 = B(3);

am1 = Am(2);
am2 = Am(3);
am3 = Am(4);

R=[];
T=[];
S=[];

Bplus = B/b0;
Bminus = b0;

Rprime = 1;
R = conv(Rprime, Bplus);
s0 = (am1-a1)/b0;
sl = (am2-a2)/b0;
s2 = (am3-a3)/b0;
S = [s0 sl s2];

% Bmprime = b0^(-1)*Bm;
Bmprime = Bm/b0;
T = conv(Ao,Bmprime);

end


