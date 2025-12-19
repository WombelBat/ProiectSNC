clc
close all
clear

y_st= 15.2;
r =173;
c = 173;
amp = y_st/c;


% load("arx_model.mat")
load("bj_model_data.mat")

amp2 = amp/dcgain(m_bj);

nb = 2;
nA =2;
nk = 1;
d=nk;



A = m_bj.F;
% A = m_bj.A;
% B= m_bj.B * amp2*c;
B= m_bj.B * amp2;
Te = m_bj.Ts;

tt2 = 360;

tt1 = tt2-60;
delta1=1.5;
zeta1 = tseta_fun(delta1)
% w1 = calc_omega(tt1,zeta1)  


delta2 = 0.15;
zeta2 = tseta_fun(delta2)
% tt2 = tt1;


Hp = tf(B,A,Te,'Variable','z^-1');
Hp  = d2c(Hp,'zoh');
[R,S,T,Bm,Am,B,A,M,P] = rst_int(Hp,Te,d,delta1,tt1,delta2,tt2);


% B = B/c

% B= B * amp2;





































% s=tf('s');
% Ho1 = w1^2/(s^2 + 2* w1 * zeta1* s + w1^2);
% 
% H0d1 = c2d(Ho1,5,'zoh');
% 
% H0d1 =tf(H0d1.Numerator, H0d1.Denominator, H0d1.Ts, 'Variable', 'z^-1');
% 
% P = H0d1.Denominator{1};
% 
% for i= 1:nk
%     B = [0 B];
%     A = [A 0];
% end
% 
% % for i = 0:nk
% %     zer = zeros(1,i)
% %     A
% %      M = [ M, zer;A']
% % 
% % end
% H= 0.4^2/(s^2 +2 * 0.4* 0.99 * s + 0.4^2);
% 
% w1 =0.8;
% zeta1 = 0.8;
% Ho1 = w1^2/(s^2 + 2* w1 * zeta1* s + w1^2);
% 
% Ho1 = 2/( ( 10*s +1 ) * ( 0.1*s +1));
% H0d1 = c2d(Ho1,0.1,'zoh');
% 
% 
% 
% 
% M = [1 -1.358 0.36 0 0;
%     0 1 -1.358 0.36 0;
%     0 0 1 -1.358 0.36;
%     0 0 0.007 0.005 0;
%     0 0 0 0.007 0.005]';
% 
% P = [1 -1.87 0.87 0 0];
% 
% x = M\P'
% 
% S = x(1:3)
% R = x(4:5)
% T = sum(P)/sum(B)
% 




