function [R,S,T,Bm,Am,B,A] = rst(Hp,Te,d,delta1,tt1,delta2,tt2,zeta1,zeta2,w1,w2)
    % aflare rst
    %% P1
    if nargin <7
        error("mai pune parametri");
    end
    if nargin <8
        zeta1 = tseta_fun(delta1);
        zeta2 = tseta_fun(delta2);
    end
    if nargin <10
        w1 = getOmega(zeta1,tt1);
        w2 = getOmega(zeta2,tt2);
    end

    %% P2 amfla Bm,Am din urmarie
    
    
    H01 = get_H02(zeta1,w1);
    H01 = c2d(H01,Te,'zoh');
    
    
    H01 =tf(H01.Numerator, H01.Denominator, H01.Ts, 'Variable', 'z^-1');
    [Bm,Am] = tfdata(H01,'v');

    %% P3 afla P din reglare
    
    H02 = get_H02(zeta2,w2);
    H02 = c2d(H02,Te,'zoh');
    
    H02 = tf(H02.Numerator,H02.Denominator,Te,'Variable','z^-1');
    P = H02.Denominator{1};
    [~,P] = tfdata(H02 ,'v');
    
    %% P4 aflam A,B din Hp
       
    Hp = c2d(Hp,Te,'zoh');
    Hp = tf(Hp.Numerator,Hp.Denominator,Te,'Variable','z^-1');
    nb = length(Hp.Numerator{1}) -1;
    na = length(Hp.Denominator{1}) -1;
    z = tf('z');
    Hp = Hp * (z^-d);
    
    [B,A] = tfdata(Hp,'v');
    A = A(1:na+1);
    
    %% P5 afla ns nR
    nS = nb+d-1;
    nR = na -1;
    
    nS = nS+1;% ca sa includa si termenul liber
    nR = nR+1; % ca sa includa si termenul liber
    
    
    %% P6 afla M
    
    % M= zeros(nb+na+d, na+nb+d);
    M=[];
    %partea cu A
    for i =0: nb + d-1
        z1 = zeros(i,1);
        z2 = zeros(nb+d-1 -i,1);
        M = [M [z1;A';z2] ];
    end
    
    %partea cu B
    for i = 0:na-1
        z1 = zeros(i,1);
        z2 = zeros(na-1 -i,1);
        M = [M [z1;B';z2] ];
    end
    
    %% P7 afla vector SR
    
    % refacem un pic P
    P =[P zeros(1,na+nb +d -length(P) )];
    
    SR = M\P';
    S = SR(1:nS);
    
    R = SR(nS+1:end);
    if length(S) == nS && length(R) == nR
        disp("good");
    else
        disp("refa functia de SR")
    end
    
    %% P8 afla T
    % B = B(2:end)
    T = sum(P)/sum(B);  
    % sau
    % T = tf(P,B,Te,'Variable','z^-1')

end
