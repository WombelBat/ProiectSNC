clc
close all
clear

load("raspunsIndicial.mat")

ursuS1.plot
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

vec_spab=SPAB_generator(N,p,Te,u0,du);

figure
plotFreq(vec_spab,Te,'s')