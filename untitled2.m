clear
clc
close all

r=85;
Te = 1;
d=2;

delta1=0;tt1=0;delta2=0;tt2=0;
w1 =1;
w2=w1;
zeta1=0.9;
zeta2= zeta1;
A = [1 -0.7];
B= [0 0.3];


% Hp = c2d(Hp,Te,'zoh');
Hp = tf(B,A,Te,'Variable','z');
Hp  = d2c(Hp,'zoh');
[R,S,T,Bm,Am,B,A]=rst(Hp,Te,d,delta1,tt1,delta2,tt2,zeta1,zeta2,w1,w2);