function CompactingTissue_1D(k_NL,k_N,k_L,L_N,N_L,ca,N_E,k_E,sigma_N,alpha_L,Movie,Mode)

addpath('./functions');
close all

%clear all;
set(0,'defaultfigurecolor',[1 1 1])
%Movie = 'On';
Plotting = 'Off';

%define the parameters:
L=300; %microns
dx = 1;
D_N = 1.95; %microns^2/s; value of eff. diffusion constant at psi=0.8
D_L = 15;%1.25*
D_N_min = 0.*D_N;
D_L_min = 0.*D_L;
% to change the scaling Deff(phi), we need to adjust the functions:
% Dscaling() and Dscaling_inv()

dt = 0.9*(dx^2)/(2*max([D_N D_L])); %time step according to Neumann criterion
tmax = 9000;
t=0;
nsteps = 0;
nPlot = 10000;
tvec = [];
tvec_l = [];%for lengthscale measurement
dcdtvec = [];
maxvelo =1;
maxvelo_threshold = 1e-3;
maxvelo_threshold2 = 1e-3;

% parameterise alpha angle parameter:
xgrid = 0:dx:L;
% for calculation of IF:
lim1 = 50;%50;
ind1 = find(xgrid<lim1);
ind1 = ind1(end);
length1 = 50;%50;

N = length(xgrid);
cN = zeros(N,1);
cL = zeros(N,1);

cE = zeros(N,1);
alphapar_crit = 0.8662;
alpha_init = 0.89;%25;
alphapar_E = 1/alpha_init;%1/0.9;

max_cNcurrent = 0;

alphapar = 1./(cE+alphapar_E);
psi_t = 1-packfracBernat(alphapar);
psi0 = mean(psi_t)

beta = 0.01; %absorption rate at maximal packing fraction
s0 = 10;%3;%10; %initial production at boundary
%psi0 = 0.9; %initial packing fraction
%k_N = 1e-4;%psi*beta;
%k_L = 0.6*1e-4;
%k_E = 1e-4;%*k_N; % for now
k_N = k_N/(psi0/(1-psi0));%/4; %adjust for packfrac
k_L = k_L/(psi0/(1-psi0));%/4;
D_N_0 = Dscaling_inv(1-psi0,D_N_min,D_N);%D_N_min + (D_N-D_N_min)/(1-psi0);% D_N/(1-psi0);
D_L_0 = Dscaling_inv(1-psi0,D_L_min,D_L);%D_L_min + (D_L-D_L_min)/(1-psi0);%D_L/(1-psi0);


c_first = [];
c_second = [];
IF_first = [];
IF_second = [];
lengthscale50 = [];
lengthscale10 = [];
lengthscale20 = [];

% profiles c(x,t):
Nodal_dyn = []; 
Lefty_dyn = [];
Nprod_dyn = [];
Lprod_dyn = [];
NLprod_dyn = [];
Ldeg_dyn = [];
Ndeg_dyn = [];
keffN_dyn =[];
keffL_dyn =[];
DeffN_dyn =[];
DeffL_dyn =[];
packfrac_dyn = [];

alpha_dyn = [];
ECad_dyn = [];

%initialise with influx at boundary
cN(1) = 0;
prod = zeros(N,1);
prod(1) = s0; %constant influx; total influx is dx*s0

% parametrise Nodal relay activation:
alpha = sigma_N;%1e-2;
alpha = alpha/psi0;
%ca = sqrt(1e+3);
ma = 2;

% parametrise Nodal-Lefty inhibition:
%k_NL = 1e-2;
%L_N = sqrt(1e+3);
m_NL = 2;

% parametrise Lefty activation:
if strcmp(Mode,'Lefty')||strcmp(Mode,'LeftyNoPack')
    alpha_L = 0;
else
    %alpha_L = 3*1e-3;
end
% rescale by porosity:
alpha_L = alpha_L/psi0;
%N_L = sqrt(1e+3);
m_L = 8;%2;

% parameterise E-Cad activation:
if strcmp(Mode,'Wnt11')||strcmp(Mode,'LeftyNoPack')
    alpha_E = 0;
else
    alpha_E = 0.5*k_E;%10*k_E;%1e-2; 
end
%k_E = alpha_E;
%N_E = 0.75*ca;%0.1*ca;
m_E = 1;

% parametrise packing fraction:
cp = 0.005; % 5 %50 higher than relay
mp = 50;%4

% % parameterise alpha angle parameter:
% alphapar_crit = 0.8662;
% alpha_init = 0.925;
% alphapar_E = 1/alpha_init;%1/0.9;
% 
% max_cNcurrent = 0;
% 
% alphapar = 1./(cE+alphapar_E);

%% read in steady state:
%data = importdata('steadystate.dat');
%c = data(:,2);

%% folder name:
foldername = Mode;%strcat('movie_kNL=',num2str(k_NL),'_kN=',num2str(k_N),'_kL=',num2str(k_L),'_kE=',num2str(k_E),'_LN=',num2str(L_N),'_NL=',num2str(N_L),'_Na=',num2str(ca),'_Dn=',num2str(D_N),'_NE=',num2str(N_E),'_cp=',num2str(cp),'_newECad_noMech_psi0=0.9');%_sharperPFlowerTh');
mkdir(foldername)
cd(foldername)

if strcmp(Movie,'On')
    vidObj1 = VideoWriter('dynamics.mp4','MPEG-4');
    vidObj1.Quality = 100;
    vidObj1.FrameRate = 20;
    open(vidObj1);
end

%% Construct matrices for time update:
%define matrix M for Nodal:
uvec = ones(N,1);
%D_Neff = D_N*(1-psi_t);
Mdiff = spdiags([uvec -2*uvec uvec],-1:1,N,N);
% no diffusive flux boundary
Mdiff(1,2) = 2; 
Mdiff(N,N-1) = 2;
%Mdiff = (D_Neff/dx^2).*Mdiff;
% add degradation:
MN = Mdiff + spdiags(-k_N*uvec,0,N,N);

%define matrix M for Lefty:
uvec = ones(N,1);
Mdiff = spdiags([uvec -2*uvec uvec],-1:1,N,N);
% no diffusive flux boundary
Mdiff(1,2) = 2; 
Mdiff(N,N-1) = 2;
Mdiff = (D_L/dx^2)*Mdiff;
% add degradation:
ML = Mdiff + spdiags(-k_L*uvec,0,N,N);

%% Analytical solution: (from old version)
xs = dx/2; %size of tiny source (discretised point source)
lambda = sqrt(D_N/k_N);
c0 = s0*xs/(k_N*lambda*(1-exp(-L/lambda)));
%c_analyt = (s0*lambda/D)*exp(-xgrid/lambda); this is with given
%diffusive influx
c_analyt = c0*exp(-xgrid/lambda);
c_lambda = interp1(xgrid,c_analyt,lambda);

x_anal = 0:(xs/10):L; % this is the same as above, actually
c_anal = zeros(size(x_anal));
cl = -(s0/(2*k_N))*(exp(xs/lambda)-exp((2*L-xs)/lambda))/(1-exp(2*L/lambda));
cr = (s0/(2*k_N))*(exp((2*L-xs)/lambda))*(-1+exp(2*xs/lambda))/(-1+exp(2*L/lambda));
c_anal(x_anal<=xs) = cl*(exp(x_anal(x_anal<=xs)/lambda)+exp(-x_anal(x_anal<=xs)/lambda))+s0/k_N;
c_anal(x_anal>xs) = cr*exp(-2*L/lambda)*exp(x_anal(x_anal>xs)/lambda)+cr*exp(-x_anal(x_anal>xs)/lambda);

%plot(x_anal, c_anal)
%psi_t = packfrac(cN,cp,mp);

% %% Time evolution: no plotting, just obtain max values 
% max_cN = 0;
% max_terms = 0;
% while (t<tmax)&&(maxvelo>maxvelo_threshold2)
%     % time-update of concentrations
%     maxvelo = max(abs(MN*cN + prod + relay(cN,alpha,ca,ma) - inhibition_NL(cN,cL,k_NL,L_N,m_NL)));
% 
%     cN = cN + dt*(MN*cN + prod + relay(cN,alpha,ca,ma) - inhibition_NL(cN,cL,k_NL,L_N,m_NL));% +prod;% + s(c));
%     cL = cL + dt*(ML*cL + relay(cN,alpha_L,N_L,m_L));
%     t = t + dt;
%     nsteps = nsteps+1;
%     max_cN = max([max_cN max(cN)]);
%     max_terms = max([max_terms max(relay(cN,alpha,ca,ma)) max(relay(cN,alpha_L,N_L,m_L)) max(inhibition_NL(cN,cL,k_NL,L_N,m_NL))]);
% end
%% Time evolution:
% without mechanical feedback
maxvelo = 1;
t = 0;
cN= zeros(N,1);
cL= zeros(N,1);

max_cN = 176.0878; %this is from the noMech simulation (Wnt11 mutant)
%psi_t = packfrac(cN,cp,mp);
psi_t = 1-packfracBernat(alphapar);
%psi_t = packfrac(alphapar,alphapar_crit,mp);

% first measurement point:
IF = 1-psi_t;%packfrac(c,cp,mp);
IF_first = [IF_first trapz(xgrid(1:ind1),IF(1:ind1))/length1];
tvec = 0;
lengthscale10 = 0;
lengthscale20 = 0;
lengthscale50 = 0;

while (t<tmax)&&(maxvelo>maxvelo_threshold)

%define matrix M for Nodal:
uvec = ones(N,1);
%D_Neff = D_N*(1-psi_t);
D_Neff = Dscaling(1-psi_t,D_N_min,D_N_0);% D_N_min+(D_N_0-D_N_min)*(1-psi_t);
dxD_Neff = gradient(D_Neff,dx);
Mdiff = spdiags([uvec -2*uvec uvec],-1:1,N,N);
% no diffusive flux boundary
Mdiff(1,2) = 2; 
Mdiff(N,N-1) = 2;
Mdiff = (D_Neff./dx^2).*Mdiff;

% diffusivity-gradient:
Mdiff2 = spdiags([-uvec 0*uvec uvec],-1:1,N,N); 
% boundary conditions:
Mdiff2(1,2) = 0; 
Mdiff2(N,N-1) = 0;

Mdiff_N2 = (dxD_Neff/(2*dx)).*Mdiff2;

% add degradation:
k_Nvec = k_N*(psi_t./(1-psi_t));
MN = Mdiff + spdiags(-k_Nvec.*uvec,0,N,N) +Mdiff_N2;

%define matrix M for Lefty:
uvec = ones(N,1);
%D_Leff = D_L*(1-psi_t);
D_Leff = Dscaling(1-psi_t,D_L_min,D_L_0);%D_L_min+(D_L_0-D_L_min)*(1-psi_t);
dxD_Leff = gradient(D_Leff,dx);
Mdiff = spdiags([uvec -2*uvec uvec],-1:1,N,N);
% no diffusive flux boundary
Mdiff(1,2) = 2; 
Mdiff(N,N-1) = 2;
Mdiff = (D_Leff./dx^2).*Mdiff;

Mdiff_L2 = (dxD_Leff/(2*dx)).*Mdiff2;
% add degradation:
k_Lvec = k_L*(psi_t./(1-psi_t));
ML = Mdiff+ spdiags(-k_Lvec.*uvec,0,N,N)+ Mdiff_L2 ;

ME = spdiags(-k_E*uvec,0,N,N);


%     % concentration-dependent diffusion:
%     D0 = D/0.2;
%     Deff = D0*(1-packfrac(c,cp));
%     dxDeff = gradient(Deff,dx);
% 
%     Mdiff = spdiags([uvec -2*uvec uvec],-1:1,N,N);
%     % no flux boundary
%     Mdiff(1,2) = 2; 
%     Mdiff(N,N-1) = 2;
%     Mdiff = (Deff/dx^2).*Mdiff;
% 
%     Mdiff2 = spdiags([-uvec 0*uvec uvec],-1:1,N,N); 
%     % boundary conditions:
%     Mdiff2(1,2) = 0; 
%     Mdiff2(N,N-1) = 0;
%     Mdiff2 = (dxDeff/(2*dx)).*Mdiff2;
% 
%     % concentration-dependent degradation:
%     k = beta*packfrac(c,cp);
%     M = Mdiff + Mdiff2 + spdiags(-k.*uvec,0,N,N);
% 
%     % update in-flux term:
%     %prod(1) = 2*Deff(1)*s0/dx;% - s0*dxDeff(1)/Deff(1);
%     prod(1) = 2*Deff(1)*s0/dx; % - s0*dxDeff(1);

    % time-update of concentrations
    maxvelo = max(abs(MN*cN + prod + relay(cN,alpha,ca,ma) - inhibition_NL(cN,cL,k_NL,L_N,m_NL)))
    alphavec = alpha.*psi_t;
    alpha_Lvec = alpha_L*psi_t;

    %cN = cN + dt*(MN*cN + prod + relay(cN,alpha,ca,ma) - inhibition_NL(cN,cL,k_NL,L_N,m_NL));% +prod;% + s(c));
    %cL = cL + dt*(ML*cL + relay(cN,alpha_L,N_L,m_L));
    
    cN = cN + dt*(MN*cN + prod + relay(cN,alphavec,ca,ma) - inhibition_NL(cN,cL,k_NL,L_N,m_NL));% +prod;% + s(c));
    cL = cL + dt*(ML*cL + relay(cN,alpha_Lvec,N_L,m_L));
    cE = cE + dt*(ME*cE + relay(cN,alpha_E,N_E,m_E));

    alphapar = 1./(cE+alphapar_E);

    max_cNcurrent = max([max_cNcurrent max(cN)]);

    t = t + dt;
    nsteps = nsteps+1;

    %psi_tplus = packfrac_irrev(psi_t,cN,cp,mp);
    %psi_tplus = packfrac_irrev(psi_t,cE,cp,mp);
    psi_tplus = packfrac_irrev(psi_t,alphapar,alphapar_crit,mp);
    psi_t = psi_tplus;

    % measure avg. concentration and IF fraction
    % integrate betweeen 0..50 
    
    %c_first = [c_first trapz(xgrid(1:ind1),cN(1:ind1))/length1];
    
    % % integrate between 150..200
    % lim2 = 150;
    % ind2 = find(xgrid<lim2);
    % ind2 = ind2(end);
    % lim3 = 200;
    % ind3 = find(xgrid<lim3);
    % ind3 = ind3(end);
    % length2 = lim3-lim2;
    % c_second = [c_second trapz(xgrid(ind2:ind3),c(ind2:ind3))/length2];
    % IF_second = [IF_second trapz(xgrid(ind2:ind3),IF(ind2:ind3))/length2];
    % position of half-maximum concentration (at that time point)

  

    % % integrate up to that length scale:
    % %lim2 = 150;
    % %ind2 = find(xgrid<lim2);
    % if xgrid(ind_l)>0
    %     ind2 = 1;%ind2(end);
    %     %lim3 = 200;
    %     %ind3 = find(xgrid<lim3);
    %     ind3 = ind_l;%ind3(end);
    %     length2 = xgrid(ind_l);%lim3-lim2;
    %     c_second = [c_second trapz(xgrid(ind2:ind3),cN(ind2:ind3))/length2];
    %     IF_second = [IF_second trapz(xgrid(ind2:ind3),IF(ind2:ind3))/length2];
    %     tvec_l = [tvec_l t];
    % end

    if (mod(nsteps,nPlot)==0)||(t==dt)
        tvec = [tvec t];
        dcdtvec = [dcdtvec maxvelo];
        maxvelo

          % position of half-concentration rel to max over whole run:
          current_max = max_cN;%max(cN);
          ind_x = find(cN>current_max/2);
          if ~isempty(ind_x)
              ind_l = ind_x(end);
              lengthscale50 = [lengthscale50 xgrid(ind_l)];
          else
              lengthscale50 = [lengthscale50 0];
          end

          ind_x = find(cN>0.1*current_max);
          if ~isempty(ind_x)
              ind_l = ind_x(end);
              lengthscale10 = [lengthscale10 xgrid(ind_l)];
          else
              lengthscale10 = [lengthscale10 0];
          end

          ind_x = find(cN>0.2*current_max);
          if ~isempty(ind_x)
              ind_l = ind_x(end);
              lengthscale20 = [lengthscale20 xgrid(ind_l)];
          else
              lengthscale20 = [lengthscale20 0];
          end

        IF = 1-psi_t;%packfrac(c,cp,mp);
        IF_first = [IF_first trapz(xgrid(1:ind1),IF(1:ind1))/length1];

        % save spatial profiles
        Nodal_dyn = [Nodal_dyn; cN']; 
        Lefty_dyn = [Lefty_dyn; cL'];
        Nprod_dyn = [Nprod_dyn; relay(cN,alphavec,ca,ma)'];
        Lprod_dyn = [Lprod_dyn; relay(cN,alpha_Lvec,N_L,m_L)'];
        NLprod_dyn = [NLprod_dyn; inhibition_NL(cN,cL,k_NL,L_N,m_NL)'];
        Ldeg_dyn = [Ldeg_dyn; (k_Nvec.*cN)'];
        Ndeg_dyn = [Ndeg_dyn; (k_Lvec.*cL)'];
        keffN_dyn =[keffN_dyn; k_Nvec'];
        keffL_dyn =[keffL_dyn; k_Lvec'];
        DeffN_dyn =[DeffN_dyn; D_Neff'];
        DeffL_dyn =[DeffL_dyn; D_Leff'];
        packfrac_dyn = [packfrac_dyn; psi_t'];
        ECad_dyn = [ECad_dyn; cE'];
        alpha_dyn = [alpha_dyn; alphapar'];

        if strcmp(Plotting,'On')
            figure1 = figure(1);
            figure1.Position = [10 10 1400 400];
            subplot(1,4,1)
            plot(xgrid,cN,'LineWidth',2);hold on
            plot(xgrid,cL,'LineWidth',2);hold off
            %plot(xgrid,c_analyt,'--','LineWidth',2);hold off
            %legend({strcat('\alpha=',num2str(alpha)),'\alpha=0'});
            legend('c_N','c_L')
            ylim([0 max_cN]);
            xlim([0 L]);
            title(strcat('t=',num2str(t/60),' min'));
            xlabel("position in tissue [µm]")
            ylabel("concentration [mol/µm^2]")
            pbaspect([1 0.75 1])
            %axis square
            subplot(1,4,2)
            semilogy(xgrid,cN,'LineWidth',2);hold on
            semilogy(xgrid,cL,'LineWidth',2);
            semilogy(xgrid,cE,'LineWidth',2);
            semilogy(xgrid,alphapar,'LineWidth',2);
            hold off
            %ylim([0 max_cN])
            xlim([0 L]);
            xlabel("position in tissue [µm]")
            ylabel("concentration [mol/µm^2]")
            legend('c_N','c_L','c_E','\alpha')
            title('log scale')
            pbaspect([1 0.75 1])
            subplot(1,4,4)
            semilogy(xgrid,relay(cN,alphavec,ca,ma),'LineWidth',2);hold on;
            semilogy(xgrid,relay(cN,alpha_Lvec,N_L,m_L),'LineWidth',2);
            semilogy(xgrid,inhibition_NL(cN,cL,k_NL,L_N,m_NL),'LineWidth',2);hold off;
            %ylim([0 max_terms])
            legend('N prod.','L prod.','N-L inhib.')
            title('magnitude of prod. and inhib. terms');
            xlabel("position in tissue [µm]")
            xlim([0 L]);
            %ylabel("production [mol/(µm^2 s)]")
            pbaspect([1 0.75 1])
            subplot(1,4,3)
            plot(xgrid,psi_t,'LineWidth',2)
            %semilogy(xgrid,cN,'LineWidth',2);hold on
            %semilogy(xgrid,cL,'LineWidth',2);hold off
            %plot(xgrid,c_analyt,'--','LineWidth',2);hold off
            %legend({strcat('\alpha=',num2str(alpha)),'\alpha=0'});
            %ylim([0 max_cN])
            ylim([0 1])
            xlim([0 L]);
            %legend('c_N','c_L')
            %title('concentrations in log-scale');
            xlabel("position in tissue [µm]")
            %ylabel("concentration [mol/µm^2]")
            % plot(xgrid,packfrac(cN,cp,mp),'LineWidth',2);
            % xlabel("position in tissue [µm]")
            ylabel("packing fraction \psi")
            % ylim([0 1])
            pbaspect([1 0.75 1])
            % subplot(1,4,4)
            % loglog(tvec, dcdtvec,'LineWidth',2);
            % title('numerical convergence')
            % xlabel('time t')
            % ylabel('max|\partial_t c|');
            % pbaspect([1 0.75 1])
            %saveas(figure(1),strcat('Relay_t=',num2str(t),'.fig'));
            %saveas(figure(1),strcat('Relay_t=',num2str(t),'.png'));
            set(findall(gcf,'-property','FontSize'),'FontSize',14);
        end

        if strcmp(Movie,'On')
            writeVideo(vidObj1, getframe(gcf));
        end

        % %  figure(2)
        % % subplot(1,3,1)
        % % yyaxis left
        % % plot(tvec,IF_first,'LineWidth',2);
        % % xlabel('time t')
        % % ylabel('avg. IF fraction')
        % % ylim([0 max(IF_first)])
        % % yyaxis right
        % % plot(tvec,c_first,'LineWidth',2);
        % % ylabel('avg. concentration [mol/µm^2]')
        % % title({'region $[0\mu m, 50\mu m]$'},'Interpreter','Latex')
        % % axis square
        % % subplot(1,3,2)
        % % if xgrid(ind_l)>0
        % %     yyaxis left
        % %     plot(tvec_l,IF_second,'LineWidth',2);
        % %     xlabel('time t')
        % %     ylabel('avg. IF fraction')
        % %     ylim([0 max(IF_second)])
        % %     yyaxis right
        % %     plot(tvec_l,c_second,'LineWidth',2);
        % %     ylabel('avg. concentration [mol/µm^2]')
        % %     title({'region $[150\mu m, 200\mu m]$'},'Interpreter','Latex')
        % %     title({'region $[0, l]$'},'Interpreter','Latex')
        % % end
        % % axis square
        % % subplot(1,3,3)
        % % yyaxis left
        % % plot(tvec,IF_first,'LineWidth',2);
        % % ylabel('avg. IF fraction in region [0,50]')
        % % ylim([0 max(IF_first)])
        % % yyaxis right
        % % plot(tvec,lengthscale,'LineWidth',2);
        % % ylabel('position of half-max. Nodal conc. [µm]')
        % % xlabel('time t')
        % % pbaspect([1 0.75 1])
        % % title({'Nodal patterning length scale vs. time'},'Interpreter','Latex')
        % % axis square
    end
end
csteady = cN; 
prodsteady = relay(cN,alpha,ca,ma);
%c = zeros(N,1);
max_cNcurrent

if strcmp(Movie,'On')
    close(vidObj1);
end

figure1 = figure(1);
figure1.Position = [10 10 1400 400];
subplot(1,4,1)
    plot(xgrid,cN,'LineWidth',2);hold on
    plot(xgrid,cL,'LineWidth',2);hold off
    legend('c_N','c_L')
    ylim([0 max_cN]);
    xlim([0 L]);
    title(strcat('t=',num2str(t/60),' min'));
    xlabel("position in tissue [µm]")
    ylabel("concentration [mol/µm^2]")
    pbaspect([1 0.75 1])
subplot(1,4,2)
    semilogy(xgrid,cN,'LineWidth',2);hold on
    semilogy(xgrid,cL,'LineWidth',2);
    semilogy(xgrid,cE,'LineWidth',2);
    semilogy(xgrid,alphapar,'LineWidth',2);
    hold off
    xlim([0 L]);
    xlabel("position in tissue [µm]")
    ylabel("concentration [mol/µm^2]")
    legend('c_N','c_L','c_E','\alpha')
    title('log scale')
    pbaspect([1 0.75 1])
subplot(1,4,4)
    semilogy(xgrid,relay(cN,alphavec,ca,ma),'LineWidth',2);hold on;
    semilogy(xgrid,relay(cN,alpha_Lvec,N_L,m_L),'LineWidth',2);
    semilogy(xgrid,inhibition_NL(cN,cL,k_NL,L_N,m_NL),'LineWidth',2);hold off;
    legend('N prod.','L prod.','N-L inhib.')
    title('magnitude of prod. and inhib. terms');
    xlabel("position in tissue [µm]")
    xlim([0 L]);
    pbaspect([1 0.75 1])
subplot(1,4,3)
    plot(xgrid,psi_t,'LineWidth',2)
    ylim([0 1])
    xlim([0 L]);
    xlabel("position in tissue [µm]")
    ylabel("packing fraction \psi")
    pbaspect([1 0.75 1])
set(findall(gcf,'-property','FontSize'),'FontSize',14);

figure(2)
    % yyaxis left
    % plot(tvec,IF_first,'LineWidth',2);
    % ylabel('avg. IF fraction in region [0,50]')
    % ylim([0 max(IF_first)])
    % yyaxis right
    plot(tvec/60,lengthscale10,'LineWidth',2);hold on;
    plot(tvec/60,lengthscale20,'LineWidth',2);
    plot(tvec/60,lengthscale50,'LineWidth',2);
    ylabel('position of x% max. Nodal conc. [µm]')
    xlabel('time t [min]')
    xlim([0 200])
    legend({'10%','20%','50%'})
    pbaspect([1 0.75 1])
    title({'Nodal patterning length scale vs. time'},'Interpreter','Latex')
    set(findall(gcf,'-property','FontSize'),'FontSize',14);

figure(3)
    % yyaxis left
    plot(tvec/60,IF_first,'LineWidth',2);
    ylabel('avg. IF fraction in region [0,50\mu m]')
    xlabel('time t [min]')
    xlim([0 200])
    ylim([0 0.2])
    pbaspect([1 0.75 1])
    set(findall(gcf,'-property','FontSize'),'FontSize',14);
    % ylim([0 max(IF_first)])
    % yyaxis right    

figure(4)
    yyaxis left
    plot(tvec/60,IF_first,'LineWidth',2);
    xlabel('time t [min]')
    ylabel('avg. IF fraction in region [0,50µm]')
    ylim([0 0.15])
    xlim([0 100])
    yyaxis right
    plot(tvec/60,lengthscale20,'LineWidth',2);
    ylabel('position of 20% max. conc. [µm]')
    title({'Nodal length scale and IFF'},'Interpreter','Latex')
    xlim([0 100])
    pbaspect([1 0.75 1])
    set(findall(gcf,'-property','FontSize'),'FontSize',14);

saveas(figure(1),'Relay_steadystate.fig');
saveas(figure(1),'Relay_steadystate.png');
saveas(figure(2),'Observables.fig');
saveas(figure(2),'Observables.png');
saveas(figure(3),'IFfraction.fig');
saveas(figure(3),'IFfraction.png');
saveas(figure(4),'Combined.fig');
saveas(figure(4),'Combined.png');

maxvelo = 1;

writematrix(max_cN,'Nodal_maximum.dat');
writematrix([tvec; IF_first],'IFfraction.dat');
writematrix([xgrid; cN'; cL'],'NodalLefty_concentration_final.dat');
writematrix([tvec; lengthscale10; lengthscale50; lengthscale20],'Nodal_lengthscale.dat');
writematrix(tvec,'tvec.dat');
writematrix(xgrid,'xgrid.dat');

%write all dynamic profiles:
writematrix(Nodal_dyn, 'Nodal_dynamics.dat');
writematrix(Lefty_dyn, 'Lefty_dynamics.dat');
writematrix(Nprod_dyn, 'Nprod_dynamics.dat');
writematrix(Lprod_dyn, 'Lprod_dynamics.dat');
writematrix(NLprod_dyn, 'NLprod_dynamics.dat');
writematrix(Ldeg_dyn, 'Ldeg_dynamics.dat');
writematrix(Ndeg_dyn, 'Ndeg_dynamics.dat');
writematrix(keffN_dyn, 'keffN_dyn_dynamics.dat');
writematrix(keffL_dyn, 'keffL_dyn_dynamics.dat');
writematrix(DeffN_dyn, 'DeffN_dyn_dynamics.dat'); 
writematrix(DeffL_dyn, 'DeffL_dyn_dynamics.dat');
writematrix(packfrac_dyn, 'Packfrac_dynamics.dat');
writematrix(ECad_dyn, 'ECad_dynamics.dat');
writematrix(alpha_dyn, 'Alpha_dynamics.dat');


% dt = dt;
% psi_t = packfrac(cN,cp,mp);
% % turn on mechanical feedback
% while (t<tmax)&&(maxvelo>maxvelo_threshold)
% 
%     % concentration-dependent diffusion:
%     D0 = D_N/0.2;
%     Deff = D0*(1-psi_t);
%     % what if D does not depend on psi:
%     %Deff = D*ones(size(xgrid));
%     dxDeff = gradient(Deff,dx);
% 
%     Mdiff = spdiags([uvec -2*uvec uvec],-1:1,N,N);
%     % no flux boundary
%     Mdiff(1,2) = 2; 
%     Mdiff(N,N-1) = 2;
%     Mdiff = (Deff/dx^2).*Mdiff;
% 
%     Mdiff2 = spdiags([-uvec 0*uvec uvec],-1:1,N,N); 
%     % boundary conditions:
%     Mdiff2(1,2) = 0; 
%     Mdiff2(N,N-1) = 0;
%     Mdiff2 = (dxDeff/(2*dx)).*Mdiff2;
% 
%     % concentration-dependent degradation:
%     k_N = beta*psi_t;
%     M = Mdiff + Mdiff2 + spdiags(-k_N.*uvec,0,N,N);
% 
%     % update in-flux term:
%     %prod(1) = 2*Deff(1)*s0/dx - s0*dxDeff(1)/Deff(1);
%     %prod(1) = 2*Deff(1)*s0/dx; % - s0*dxDeff(1);
%     prod(1) = s0; %corrected for constant influx
% 
%     % time-update of concentrations
%     maxvelo = max(abs(M*cN +prod + relay(cN,alpha,ca,ma)))
%     tvec = [tvec t];
%     dcdtvec = [dcdtvec maxvelo];
% 
%     cN = cN + dt*(M*cN +prod + relay(cN,alpha,ca,ma));% +prod;% + s(c));
%     t = t + dt;
%     nsteps = nsteps+1;
%     % update psi only if it increases
%     psi_tplus = packfrac_irrev(psi_t,cN,cp,mp);
%     psi_t = psi_tplus;
% 
%     % measure avg. concentration and IF fraction
%     % integrate betweeen 0..50 
%     lim1 = 50;
%     ind1 = find(xgrid<lim1);
%     ind1 = ind1(end);
%     length1 = 50;
%     c_first = [c_first trapz(xgrid(1:ind1),cN(1:ind1))/length1];
%     IF = 1-psi_t;%packfrac(c,cp,mp);
%     IF_first = [IF_first trapz(xgrid(1:ind1),IF(1:ind1))/length1];
%     % % integrate between 150..200
%     % lim2 = 150;
%     % ind2 = find(xgrid<lim2);
%     % ind2 = ind2(end);
%     % lim3 = 200;
%     % ind3 = find(xgrid<lim3);
%     % ind3 = ind3(end);
%     % length2 = lim3-lim2;
%     % c_second = [c_second trapz(xgrid(ind2:ind3),c(ind2:ind3))/length2];
%     % IF_second = [IF_second trapz(xgrid(ind2:ind3),IF(ind2:ind3))/length2];
%     % position of half-maximum concentration (at that time point)
%     current_max = max(cN);
%     ind_x = find(cN>current_max/2);
%     ind_l = ind_x(end);
%     lengthscale50 = [lengthscale50 xgrid(ind_l)];
% 
%     % integrate up to that length scale:
%     %lim2 = 150;
%     %ind2 = find(xgrid<lim2);
%     if xgrid(ind_l)>0
%         ind2 = 1;%ind2(end);
%         %lim3 = 200;
%         %ind3 = find(xgrid<lim3);
%         ind3 = ind_l;%ind3(end);
%         length2 = xgrid(ind_l);%lim3-lim2;
%         c_second = [c_second trapz(xgrid(ind2:ind3),cN(ind2:ind3))/length2];
%         IF_second = [IF_second trapz(xgrid(ind2:ind3),IF(ind2:ind3))/length2];
%         tvec_l = [tvec_l t];
%     end
% 
%     if mod(nsteps,10*nPlot)==0
%         figure1 = figure(1);  
%         figure1.Position = [10 10 1000 400]; 
%         subplot(1,3,1)
%         plot(xgrid,csteady,'k:','LineWidth',2);hold on;
%         plot(xgrid,cN,'LineWidth',2);hold off;
%         %plot(xgrid,(s0*sqrt(D/beta)/D)*exp(-xgrid/sqrt(D/beta)),'--');hold off
%         %legend({strcat('\alpha=',num2str(alpha)),'\alpha=0'});
%         legend({'relay st.st.','with mech. feedback'});
%         %legend('concentration','relay production')
%         title(strcat('t=',num2str(t)));
%         xlabel("position in tissue [µm]")
%         ylabel("concentration [mol/µm^2]")
%         axis square
%         subplot(1,3,2)
%         plot(xgrid,prodsteady,'k:','LineWidth',2);hold on;
%         plot(xgrid,relay(cN,alpha,ca,ma),'LineWidth',2);hold off;
%         %legend('concentration','relay production')
%         title(strcat('t=',num2str(t)));
%         xlabel("position in tissue [µm]")
%         ylabel("production [mol/(µm^2 s)]")
%         axis square
%         % if xgrid(ind_l)>0
%         %     yyaxis left
%         %     plot(tvec_l,IF_second,'LineWidth',2);
%         %     xlabel('time t')
%         %     ylabel('avg. IF fraction')
%         %     ylim([0 max(IF_second)])
%         %     yyaxis right
%         %     plot(tvec_l,c_second,'LineWidth',2);
%         %     ylabel('avg. concentration [mol/µm^2]')
%         %     %title({'region $[150\mu m, 200\mu m]$'},'Interpreter','Latex')
%         %     title({'region $[0, l]$'},'Interpreter','Latex')
%         %     axis square
%         % end
%         subplot(1,3,3)
%         plot(xgrid,psi_t,'LineWidth',2);%packfrac(c,cp,mp)
%         xlabel("position in tissue [µm]")
%         ylabel("packing fraction \psi")
%         ylim([0 1])
%         axis square
% %         subplot(1,4,4)
% %         loglog(tvec, dcdtvec,'LineWidth',2);
% %         title('numerical convergence')
% %         xlabel('time t')
% %         ylabel('max|\partial_t c|');
% %         axis square
%         %saveas(figure1,strcat('RelayMechFeedback_t=',num2str(t),'.fig'));
%         %saveas(figure1,strcat('RelayMechFeedback_t=',num2str(t),'.png'));
% 
%         writeVideo(vidObj1, getframe(gcf));
% 
%         figure(2)
%         subplot(1,3,1)
%         yyaxis left
%         plot(tvec,IF_first,'LineWidth',2);
%         xlabel('time t')
%         ylabel('avg. IF fraction')
%         ylim([0 max(IF_first)])
%         yyaxis right
%         plot(tvec,c_first,'LineWidth',2);
%         ylabel('avg. concentration [mol/µm^2]')
%         title({'region $[0\mu m, 50\mu m]$'},'Interpreter','Latex')
%         axis square
%         subplot(1,3,2)
%         yyaxis left
%         plot(tvec_l,IF_second,'LineWidth',2);
%         xlabel('time t')
%         ylabel('avg. IF fraction')
%         ylim([0 max(IF_second)])
%         yyaxis right
%         plot(tvec_l,c_second,'LineWidth',2);
%         ylabel('avg. concentration [mol/µm^2]')
%         title({'region $[150\mu m, 200\mu m]$'},'Interpreter','Latex')
%         axis square
%         subplot(1,3,3)
%         yyaxis left
%         plot(tvec,IF_first,'LineWidth',2);
%         xlabel('time t')
%         ylabel('avg. IF fraction in region [0,50]')
%         ylim([0 max(IF_first)])
%         yyaxis right
%         plot(tvec,lengthscale50,'LineWidth',2);
%         ylabel('position of half-max. conc. [µm]')
%         title({'patterning length scale vs. IF fraction'},'Interpreter','Latex')
%         axis square
%     end
% end
% saveas(figure1,'RelayMechFeedback_steadystate.fig');
% saveas(figure1,'RelayMechFeedback_steadystate.png');
% saveas(figure(2),'Observables.fig');
% saveas(figure(2),'Observables.png');
% close(vidObj1)

cd ..
end
