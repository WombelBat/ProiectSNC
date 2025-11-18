clc
close all
clear

load("raspunsIndicial.mat")

% ursuS1.plot
u0 = 68;
du = 15;
c=173;
y_st= 15.62;

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
t_out_spab =spab_data_temp.tout;

t_out_spab =t_out_spab((2*N*p+1) +27 +22 :length(spab_data_temp.comanda.Data) -1) ;
y_spab = spab_data_temp.simout.Data((2*N*p+1) +27 +22 : end-1 );
u_spab = spab_data_temp.comanda.Data((2*N*p+1) +27 +22 : end-1 );

winpin = [y_spab u_spab ];
save('DateBXY.txt','-ascii','winpin')

% y_spab = spab_data_temp.simout.Data((2*N*p+1)  : end-1 );
% u_spab = spab_data_temp.comanda.Data((2*N*p+1) : end-1 );

Ts_spab = spab_data_load.Te;
spab_data = iddata(y_spab,u_spab,Ts_spab);



%% 
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

%% 

% filtram
[b_butt,a_butt] =butter(1, 0.2);
% y_spab_filtrat = filter(b_butt,a_butt,y_spab);
y_spab_filtrat = filter(b_butt,a_butt,spab_data_cent.OutputData); % cu detrend

% y_spab_filtrat =spab_data_cent.OutputData;

% se intampla pentru ca dupa filtrarre tot nu arata bine datele nu trebuie
% y_spab_filtrat = y_spab_filtrat(6+11: end );
% u_spab_filtrat = u_spab(6+11:end);

u_spab_filtrat = spab_data_cent.InputData;

plotFreq(y_spab,1/Te)
legend( "Spectrul Iesirii Nefiltrate (y_spab)")
title("Spectrul Iesirii Nefiltrate (y_spab)")

data_spab_filt = iddata(y_spab_filtrat,u_spab_filtrat,Ts_spab);

plotFreq(y_spab_filtrat,1/Te)
title("Spectrul Iesirii Filtrate")
legend("Spectrul Iesirii Filtrate")


plotFreq(data_spab_filt.OutputData,1/Te,legend = "Spectrul Iesirea filtrata")
title("Spectrul Iesirea filtrata")
legend("Spectrul Iesirea filtrata")

plotFreq(data_spab_filt.InputData, 1/Te)
title("Spectrul Intrarea filtrata")
legend("Spectrul Intrarea filtrata")

% separate data
eData = data_spab_filt(1:Lspab);
vData = data_spab_filt(Lspab+1 :end);

% test gets fit 72.22% but [9 10] parameters
t_med = getTrend(eData,0);
eData = detrend(eData,t_med);
t_med = getTrend(vData,0);
vData = detrend(vData,t_med);

% timp mort
nk =min([tmort, delayest(data_spab_filt)]);

%%
% arx model

% find data na nb
M = struc(1:10,1:10,nk);

V = arxstruc(eData.InputData,eData.OutputData,vData.InputData,vData.OutputData,M);

% 

% selstruc(V,0)
% ordin = selstruc(V,'plot');

%  din model cele mai bune
na = 9; %9
nb  =5; %5

m_arx = arx(eData, [na nb nk]);
% this just gives fittness also get the others
[~,fit_arx,~]=compare(m_arx,vData);
figure
sgtitle('Validare si Analiza Reziduuri Model ARX')
subplot(2,1,1)
compare(m_arx,vData)
[E_arx,R_arx] = resid(vData,m_arx);
subplot(2,1,2)
resid(vData,m_arx);

loss_arx = m_arx.Report.Fit.LossFcn;
fpe_arx = fpe(m_arx);
mse_arx = mean(E_arx.OutputData.^2);
figure('Name', 'Analiza Stabilitate - ARX');
subplot(2,2,1); step(m_arx); title('Raspuns la treapta (ARX)');
subplot(2,2,2); pzmap(m_arx); title('Diagrama Poli-Zerouri (ARX)');
subplot(2,2,3); nyquist(m_arx); title('Diagrama Nyquist (ARX)');
subplot(2,2,4); bode(m_arx); title('Diagrama Bode (ARX)');


% figure
% step(m_final)

amp = y_st/c;
amp2 = amp/dcgain(m_arx);

figure("Name","Grafic de comparatie a corectie amplitdini (ARX)")
subplot(2,1,1)
step(m_arx*amp2 *c);
title("raspuns corectat al modelului (ARX)")
subplot(2,1,2)
step(c * m_arx)
title("raspuns original modelului (ARX)")
%%
% armax model

% fit_armax= 0;
% nc_limit = 10;
% for i =1 :nc_limit 
%     t_armax = armax(eData,[na,nb,i,nk]);
%     [~,t_fit,~]=compare(t_armax,vData);
% 
%     if t_fit > fit_armax
%         fit_armax = t_fit;
%         nc = i;
%         m_armax = t_armax;
%     end
% end
load("armax_model.mat")

figure
sgtitle('Validare si Analiza Reziduuri Model ARMAX')
subplot(2,1,1)
compare(m_armax,vData)
[E_armax,R_armax] = resid(vData,m_armax);
subplot(2,1,2)
resid(vData,m_armax);


loss_armax = m_armax.Report.Fit.LossFcn;
fpe_armax = fpe(m_armax);
mse_armax = mean(E_armax.OutputData.^2);
figure('Name', 'Analiza Stabilitate - ARMAX');
    subplot(2,2,1); step(m_armax); title('Raspuns la treapta (ARMAX)');
    subplot(2,2,2); pzmap(m_armax); title('Diagrama Poli-Zerouri (ARMAX)');
    subplot(2,2,3); nyquist(m_armax); title('Diagrama Nyquist (ARMAX)');
    subplot(2,2,4); bode(m_armax); title('Diagrama Bode (ARMAX)');

% figure
% step(m_final)

amp = y_st/c;
amp2 = amp/dcgain(m_armax);

figure("Name","Grafic de comparatie a corectie amplitdini ARMAX")
subplot(2,1,1)
step(m_armax*amp2 *c);
title("raspuns corectat al modelului ARMAX")
subplot(2,1,2)
step(c * m_armax)
title("raspuns original modelului ARMAX")
 
%%
% bj model

% fit_bj= 0;
% nb_limit = 4;
% nc_limit  = 4;
% nd_limit  = 4;
% nf_limit  = 4;
% for i = 1 : nb_limit
%     for j = 1 : nc_limit
%         for k = 1 : nd_limit
%             for l = 1 : nf_limit
% 
%                 t_bj = bj(eData, [i j k l nk]);
%                 [~,t_fit,~]=compare(t_bj,vData);
% 
%                 if fit_bj < t_fit
%                     fit_bj = t_fit;
%                     nb_bj = i;
%                     nc_bj = j;
%                     nd = k;
%                     nf = l;
%                     m_bj = t_bj;
% 
%                 end
%             end
%         end
%     end
% end  

load("bj_model_data.mat")

figure
sgtitle('Validare si Analiza Reziduuri Model BJ')
subplot(2,1,1)
compare(m_bj,vData);

[E_bj,R_bj] = resid(vData,m_bj);
subplot(2,1,2)
resid(vData,m_bj);

loss_bj = m_bj.Report.Fit.LossFcn;
fpe_bj = fpe(m_bj);
mse_bj = mean(E_bj.OutputData.^2);
figure('Name', 'Analiza Stabilitate - BJ');
    subplot(2,2,1); step(m_bj); title('Raspuns la treapta (BJ)');
    subplot(2,2,2); pzmap(m_bj); title('Diagrama Poli-Zerouri (BJ)');
    subplot(2,2,3); nyquist(m_bj); title('Diagrama Nyquist (BJ)');
    subplot(2,2,4); bode(m_bj); title('Diagrama Bode (BJ)');

% figure
% step(m_final)

amp = y_st/c;
amp2 = amp/dcgain(m_bj);

figure("Name","Grafic de comparatie a corectie amplitdini BJ")
subplot(2,1,1)
step(m_bj*amp2 *c);
title("raspuns corectat al modelului BJ")
subplot(2,1,2)
step(c * m_bj)
title("raspuns original modelului BJ")


  %%
% oe model 

% fit_oe= 0;
% nb_limit = 6;
% nf_limit  = 6;
% 
% eData.InterSample = 'foh';
% vData.InterSample = 'foh';
% for i = 1 : nb_limit
%     for j = 1 : nf_limit
%                 t_oe = oe(eData, [i j nk]);
%                 [~,t_fit,~]=compare(t_oe,vData);
% 
%                 if fit_oe < t_fit
%                     fit_oe = t_fit;
%                     nb_oe = i;         
%                     nf_oe = j;
%                     m_oe = t_oe;
%                 end
%     end
% end  

load("oe_model_data.mat")

figure
sgtitle('Validare si Analiza Reziduuri Model OE')
subplot(2,1,1)
compare(m_oe,vData)

[E_oe,R_oe] = resid(vData,m_oe);
subplot(2,1,2)
resid(vData,m_oe);

loss_oe = m_oe.Report.Fit.LossFcn;
fpe_oe = fpe(m_oe);
mse_oe = mean(E_oe.OutputData.^2);
figure('Name', 'Analiza Stabilitate - OE');
    subplot(2,2,1); step(m_oe); title('Raspuns la treapta (OE)');
    subplot(2,2,2); pzmap(m_oe); title('Diagrama Poli-Zerouri (OE)');
    subplot(2,2,3); nyquist(m_oe); title('Diagrama Nyquist (OE)');
    subplot(2,2,4); bode(m_oe); title('Diagrama Bode (OE)');

% figure
% step(m_final)

amp = y_st/c;
amp2 = amp/dcgain(m_oe);

figure("Name","Grafic de comparatie a corectie amplitdini OE")
subplot(2,1,1)
step(m_oe*amp2 *c);
title("raspuns corectat al modelului OE")
subplot(2,1,2)
step(c * m_oe)
title("raspuns original modelului OE")
   %%
% comparatie

Modele = ["ARX"; "ARMAX"; "BJ"; "OE"];
FIT_val = [fit_arx; fit_armax; fit_bj; fit_oe];
Loss_val = [loss_arx; loss_armax; loss_bj; loss_oe]*1e4;
FPE_val = [fpe_arx; fpe_armax; fpe_bj; fpe_oe]*1e4;
MSE_val = [mse_arx; mse_armax; mse_bj; mse_oe]*1e4;

% Creare si afisare tabel
T = table(Modele, FIT_val, Loss_val, FPE_val, MSE_val, ...
    'VariableNames', {'Model', 'FIT (%)', 'Loss Function', 'FPE', 'MSE'});

disp(T);

%%
% aplificare model
m_final = m_bj;

% figure
% step(m_final)

amp = y_st/c;
amp2 = amp/dcgain(m_final);

figure("Name","Grafic de comparatie a corectie amplitdini DIN MODELUL ALES")
subplot(2,1,1)
step(m_final*amp2 *c);
title("raspuns corectat al modelului final")
subplot(2,1,2)
step(c * m_final)
title("raspuns original modelului final ")
%%


