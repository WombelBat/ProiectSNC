clc
close all
clear

load("raspunsIndicial.mat")

% ursuS1.plot
u0 = 68;
du = 15;

tmort = 3;
tc = 75;
tt=88;

Te = 5;
ti = 3000;

p =2;

N = ( tc/( p* Te ) );

N = ceil(N);
Lspab = 2^N -1;

tspab = p*Lspab * Te;

%vec_spab=SPAB_generator(N,p,Te,u0,du);

% figure
% plotFreq(vec_spab,Te,'s')

spab_data_load = load("Ursu_spab_60min.mat");
spab_data_temp = spab_data_load.Ursu_spab_60min;
% 
y_spab = spab_data_temp.simout.Data((2*N*p+1) +27 +22 : end-1 );
u_spab = spab_data_temp.comanda.Data((2*N*p+1) +27 +22 : end-1 );

% y_spab = spab_data_temp.simout.Data((2*N*p+1)  : end-1 );
% u_spab = spab_data_temp.comanda.Data((2*N*p+1) : end-1 );

Ts_spab = spab_data_load.Te;
spab_data = iddata(y_spab,u_spab,Ts_spab);


% figure;
% subplot(3,1,1)
% plot((y_spab))
% title("hopa")
% subplot(3,1,2)
% plot((spab_data_temp.simout.Data))
% subplot(3,1,3)
% plot(spab_data_temp.comanda.Data)

t_med = getTrend(spab_data,0);
spab_data_cent = detrend(spab_data,t_med);

% filtram
[b_butt,a_butt] =butter(1, 0.2);
% y_spab_filtrat = filter(b_butt,a_butt,y_spab);
y_spab_filtrat = filter(b_butt,a_butt,spab_data_cent.OutputData); % cu detrend

% se intampla pentru ca dupa filtrarre tot nu arata bine datele nu trebuie
% y_spab_filtrat = y_spab_filtrat(6+11: end );
% u_spab_filtrat = u_spab(6+11:end);

u_spab_filtrat = spab_data_cent.InputData;

% figure
% subplot(2,1,1)
% plot(y_spab)
data_spab_filt = iddata(y_spab_filtrat,u_spab_filtrat,Ts_spab);
% subplot(2,1,2)
% plot(y_spab_filtrat)

% separate data
eData = data_spab_filt(1:Lspab);
vData = data_spab_filt(Lspab+1 :end);
% timp mort
nk =min([tmort, delayest(data_spab_filt)]);

% find data na nb
M = struc(1:10,1:10,nk);

V = arxstruc(eData.InputData,eData.OutputData,vData.InputData,vData.OutputData,M);

% 

% selstruc(V,0)
% ordin = selstruc(V,'plot');

%  din model cele mai bune
na = 8;
nb  =5;

m_arx = arx(eData, [na nb nk]);
[ymod,fit,ic]=compare(m_arx,vData);
compare(m_arx,vData)
