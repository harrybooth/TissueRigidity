% define parameter values (rates, activation thresholds, baseline
% diffusivity), run simulations for 4 conditions, plot results

% parameters as in: January2025/Results_kNL=0.1_kN=1.25e-06_kL=7.5e-07_kE=0.0001_LN=15.8114_NL=200_Na=31.6228_NE=23.7171_alphaL=0.012_withDcutoff”

addpath('./functions');
addpath(pwd); %add current folder as path
close all;

% define parameters
k_NL = 0.1;%1;%0.1
k_N = 1.25*1e-6;
k_L = 7.5*1e-7;%5*
k_E = 5*1e-4;
L_N = 15.8114;
N_L = 200;%150;%200;
ca = 31.6228; 
sigma_N = 1e-2;
N_E = 0.75*ca; %this is actually lower than N_a!
alpha_L = 12*1e-3;
Movie = 'Off';

foldername = strcat('Results_kNL=',num2str(k_NL),'_kN=',num2str(k_N),'_kL=',num2str(k_L),'_kE=',num2str(k_E),'_LN=',num2str(L_N),'_NL=',num2str(N_L),'_Na=',num2str(ca),'_NE=',num2str(N_E),'_alphaL=',num2str(alpha_L),'_Dcutoff0_newphi');%,'_higher2xAlphaL'  ,'_withDcutoff'
mkdir(foldername)
cd(foldername)

% run simulations for WT, Wnt11, Lefty, LeftyNoPack

Mode = "WT";
CompactingTissue_1D(k_NL,k_N,k_L,L_N,N_L,ca,N_E,k_E,sigma_N,alpha_L,Movie,Mode);
Mode = "Wnt11";
CompactingTissue_1D(k_NL,k_N,k_L,L_N,N_L,ca,N_E,k_E,sigma_N,alpha_L,Movie,Mode);

% Mode = "Lefty";
% CompactingTissue_1D(k_NL,k_N,k_L,L_N,N_L,ca,N_E,k_E,Movie,Mode);
% Mode = "LeftyNoPack";
% CompactingTissue_1D(k_NL,k_N,k_L,L_N,N_L,ca,N_E,k_E,Movie,Mode);

% run plotting file: inside function, make the 'Figures' subfolder
mkdir('Figures_png')
mkdir('Figures_fig')

plotting(k_NL,k_N,k_L,L_N,N_L,ca,N_E,k_E);


cd ..