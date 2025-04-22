% define parameter values (rates, activation thresholds, baseline
% diffusivity), run simulations for 4 conditions, plot results

addpath('./functions');
addpath(pwd); %add current folder as path
close all;

% define Nodal parameters:
%k_N = 1.25*1e-6;

% reduction factor<1 for k_N in order to increase the non-relay lambda, and
% time scale
factor = 1;%0.4;

k_N = factor*1e-2;%1e-3;
D_N = 1.95; 
L = 300;
dx = 1;
s0 = 10;
xs = dx/2; %size of tiny source (discretised point source)
lambda = sqrt(D_N/k_N);
c0 = s0*xs/(k_N*lambda*(1-exp(-L/lambda)))

N0 = 0.95*c0;
sigma_crit = k_N*N0;
sigma_N =2.5*sigma_crit;%2*sigma_crit;%1;%0.01;%0.28;   1.75
N_a = N0 %0.01*31.6228; 

%with relay, the peak magnitude is roughly: sigmaN/kN + c0
Nrelay = sigma_N/k_N + c0;

% define Lefty parameters:
%k_L = 7.5*1e-7;%5*
N_L = 0.9*Nrelay %0.95*Nrelay
%2.5*N0;%117;%50;%150;%200;
k_L = (1e-4)*k_N;
alpha_L = 10*k_L;
L_N = 2;%2*15.8114;
k_NL = (1e+8)*k_N;%1e+4;%0.1;%1;%0.1

k_E = 0.1*(1e-1)*k_N; %keep k_E to original value when changing k_N (1/factor)*0.25*
N_E = N_a %0.5*N_L
%(N_a+N_L)/2
%alpha_L = 3000*1e-3;

Movie = 'On';

foldername = strcat('Combined_kNL=',num2str(k_NL),'_kN=',num2str(k_N),'_kL=',num2str(k_L),'_kE=',num2str(k_E),'_LN=',num2str(L_N),'_NL=',num2str(N_L),'_Na=',num2str(N_a),'_NE=',num2str(N_E),'_alphaL=',num2str(alpha_L),'_sigmaN',num2str(sigma_N));%,'_higher2xAlphaL'  ,'_withDcutoff'
mkdir(foldername)
cd(foldername)

% run simulations for WT, Wnt11, Lefty, LeftyNoPack
Mode = "Lefty";
CompactingTissue_1D(k_NL,k_N,k_L,L_N,N_L,N_a,N_E,k_E,sigma_N,0,Movie,Mode);

Mode = "WT";
CompactingTissue_1D(k_NL,k_N,k_L,L_N,N_L,N_a,N_E,k_E,sigma_N,alpha_L,Movie,Mode);

Mode = "Wnt11";
CompactingTissue_1D(k_NL,k_N,k_L,L_N,N_L,N_a,N_E,k_E,sigma_N,alpha_L,Movie,Mode);

% Mode = "Lefty";
% CompactingTissue_1D(k_NL,k_N,k_L,L_N,N_L,ca,N_E,k_E,Movie,Mode);
% Mode = "LeftyNoPack";
% CompactingTissue_1D(k_NL,k_N,k_L,L_N,N_L,ca,N_E,k_E,Movie,Mode);

% run plotting file: inside function, make the 'Figures' subfolder

mkdir('Figures_png')
mkdir('Figures_fig')

plotting(k_NL,k_N,k_L,L_N,N_L,N_a,N_E,k_E);


cd ..