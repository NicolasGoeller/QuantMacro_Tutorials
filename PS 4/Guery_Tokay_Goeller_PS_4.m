% ==========================
% Problem set solving the RBC model different
% ways
% ==========================
clear
close all
%addpath('C:\Users\tbroer\Dropbox\Teaching\PSE\2021 Quantitative Macro\Problem sets\PS RBC')
%addpath('C:\Users\tbroer\Dropbox\Teaching\PSE\2021 Quantitative Macro\Problem sets')

%%
% ============
% Calibration targets
% ============

Kshare=1/3;
Rbar=1.01;
Lbar=1/3;
InvKrat=0.05;

% ============
% Exogenous parameters
% ============
gamma = 2; % CRRA
rho = 0.95;   % persistence of TFP shock
sigmaepsilon = 0.007; % volatility of TFP shock
psi=1; %inverse Frisch elasticity

% ============
% options and convergence criteria
% ============
criter_V = 1e-6; % conv criterion for value function
N=80; % number of grid points
linear=1; % grid linear or not
M=3; % number of support points for the shock
T=100; % periods for transition
N_sim=100; % number of simulations

% ==============
% Problem I: Calibration
% ==============
delta= InvKrat;
alpha=Kshare;
beta=1/Rbar;
margprod=1/beta-1+delta;
KNrat=((alpha*beta)/(1-beta*(1-delta)))^(1/(1-alpha));
wbar=(1-alpha)*KNrat^alpha;
cKrat=(margprod/alpha)-delta;
theta=wbar/(Lbar^(gamma+psi)*((margprod/alpha)-delta)^gamma*(margprod/alpha)^(gamma/(alpha-1)));
kbar=KNrat*Lbar;
% ==============
% Grids, transition probabilities, etc
% ==============
kmax=1.1*kbar;
kmin=kbar*0.9;
if delta==1
    if linear==1
        kgrid=linspace(kmin,kmax,N);
    else
        temp=linspace(0,0.5,N).^5/0.5^5*(kmax-kmin);
        kgrid=kmin+temp;
    end
else
    kgrid=linspace(kmin ,kmax,N);
end

%%
% ============
% Markov chain
% ============
%[Z, P] = tauchen(M,0,rho,sigmaepsilon,1)
[Z,P] = rouwenhorst(M,0,rho,sigmaepsilon);

%Need function, to iterate MC to z_t+1 based on z_t and P
% Function needs to map into itself

%for discrete

% simulate discrete-time, discrete-state Markov Process
z = dismc_sim(P,3,T,N_sim);
%Below has T on cols, N_sim on rows?
% for j=1:N_sim
%     z(j,1)=3; %% We are starting in the high state here?
%     for t=2:T
%         z(j,t)=;
%     end
% end

% simulate discrete-time, continuous-state Markov Process
z_cont = conmc_sim(0, sigmaepsilon, rho, 0, T,N_sim);
% for j=1:N_sim
%     z_cont(j,1)=1;
%     for t=2:T
%         %z_cont(j,t)=[Fill this in if you do linearisation or deterministic sequence];
%     end
% end


Z_cont_lev = exp(z_cont);
Z_lev=exp(Z);
Z_sim=Z_lev(z);


%%

% =====================
% Option III:   Log-linearization
% =====================

% some ratios as function of parameters
ybar=kbar^alpha*Lbar^(1-alpha);
cbar=ybar-delta*kbar;
ckrat=cbar/kbar;
Rbar=alpha*(ybar/kbar);

%%

% a. write system as A_1 E[x_t+1]=A_2  x_t+B_1 epsilon_t
% order k,a,c,l
a=[     kbar/ybar       ,       0       ,       0       ,       0           ;
    beta*Rbar*(alpha-1) ,   beta*Rbar   ,   -gamma      , -beta*Rbar*(alpha-1);
            0           ,       0       ,       0       ,                   0;
            0           ,       1       ,       0       ,                   0];



b=[alpha + (1-delta)*kbar/ybar  ,       1     ,        -cbar/ybar   ,   1-alpha ; 
            0                   ,       0     ,         -gamma      ,       0   ; 
            -alpha              ,       -1    ,         gamma       , 1 + alpha ; 
            0                   ,       rho   ,       0             ,     0     ];

%%

% re-order to have state variables first
% order k,a,c,n
nk=2;
[f,p] = solab(a,b,nk);

%%
% extract cons and lab policies
c_polfunc=f(1,:);
n_polfunc=f(2,:);
LOM=p(1,:);

%%
% policy functions are in log deviations. so need to pre-multiply by cbar
% and add cbar to get levels
% for j=1:M
%     c_pol_lin(:,j)=cbar*exp(c_polfunc(1)*log(kgrid/kbar)+c_polfunc(2)*Z(j));
%     %c_pol_lin(:,j)=cbar*(c_polfunc(1)*) 
%     n_pol_lin(:,j)=[Fill this in];
%     k_pol_lin(:,j)=[Fill this in];
% end

%%


k_sim_lin(:,1)=kbar*ones(N_sim,1);
for j=1:N_sim
    for t=2:T
        c_sim_lin(j,t-1)=(c_polfunc(1)*(k_sim_lin(j,t-1)-kbar)/kbar+c_polfunc(2)*(Z_cont_lev(j,t-1)-1))* cbar +cbar;
        L_sim_lin(j,t-1)=(n_polfunc(1)*(k_sim_lin(j,t-1)-kbar)/kbar+n_polfunc(2)*(Z_cont_lev(j,t-1)-1))* Lbar + Lbar;
        k_sim_lin(j,t)=(LOM(1)*(k_sim_lin(j,t-1)-kbar)/kbar+LOM(2)*(Z_cont_lev(j,t-1)-1))*kbar + kbar;
        inv_sim_lin(j,t-1)=k_sim_lin(j,t)-(1-delta)*k_sim_lin(j,t-1);
        y_sim_lin(j,t-1)=Z_cont_lev(j,t-1)*(k_sim_lin(j,t-1)^alpha)*L_sim_lin(j,t-1)^(1-alpha);
    end
    
    c_sim_lin(j,T)=(c_polfunc(1)*(k_sim_lin(j,T)-kbar)/kbar+c_polfunc(2)*(Z_cont_lev(j,T)-1))* cbar +cbar;
    L_sim_lin(j,T)=(n_polfunc(1)*(k_sim_lin(j,T)-kbar)/kbar+n_polfunc(2)*(Z_cont_lev(j,T)-1))* Lbar + Lbar;
    k_sim_lin(j,T+1)= (LOM(1)*(k_sim_lin(j,T)-kbar)/kbar+LOM(2)*(Z_cont_lev(j,T)-1))*kbar + kbar;
    inv_sim_lin(j,T)=k_sim_lin(j,T+1)-(1-delta)*k_sim_lin(j,T);
    y_sim_lin(j,T)=Z_cont_lev(j,T)*(k_sim_lin(j,T)^alpha)*L_sim_lin(j,T)^(1-alpha);
end



%%


% ==============
% Figures
% ==============

plotz=3;
figure('Name','Policy functions'); % this is not asked
subplot(3,1,1);
title('Consumption','Interpreter','Latex','fontsize',13);
hold on;
% plot([c_disc(:,plotz)])
% plot([c_interp(:,plotz)])
plot([c_sim_lin(:,plotz)]);
xlabel('Capital','Interpreter','Latex','fontsize',13);
subplot(3,1,2);
title('Labor supply','Interpreter','Latex','fontsize',13);
hold on;
% plot([L_disc(:,plotz)])
% plot([L_interp(:,plotz)])
plot([L_sim_lin(:,plotz)]);
xlabel('Capital','Interpreter','Latex','fontsize',13);
subplot(3,1,3);
title('K prime','Interpreter','Latex','fontsize',13);
hold on;
% plot([kprime_disc(:,plotz)])
% plot([kprime_interp(:,plotz)])
plot([k_sim_lin(:,plotz)]);
xlabel('Capital','Interpreter','Latex','fontsize',13);
h = legend('Linearised','Location', 'best','Orientation','Vertical');
set(h,'fontsize',13,'Interpreter','Latex');%'Orientation', 'horizontal'

%%

figure('Name','Simulated time series')

subplot(3,1,1)
hold on

title('Simulation of the capital stock')
%for sim=1:N_sim
sim=1
% plot(k_disc_sim(sim,1:T),'k','Linewidth',1)
% plot(k_interp_sim(sim,1:T),'k:','Linewidth',1)
plot(k_sim_lin(sim,1:T),'k--','Linewidth',1)
%plot(k_det_sim(sim,1:T),'k-+','Linewidth',1)
%end
h = legend('Linear','Location', 'best','Orientation','Vertical');
set(h,'fontsize',12,'Interpreter','Latex');%'Orientation', 'horizontal'

subplot(3,1,2)
hold on

title('Simulation of consumption')
%for sim=1:N_sim
sim=1
% plot(c_disc_sim(sim,1:T),'k','Linewidth',1)
% plot(c_interp_sim(sim,1:T),'k:','Linewidth',1)
plot(c_sim_lin(sim,1:T),'k--','Linewidth',1)
%plot(c_det_sim(sim,1:T),'k-+','Linewidth',1)
h = legend('Linear','Location', 'best','Orientation','Vertical');
set(h,'fontsize',12,'Interpreter','Latex');%'Orientation', 'horizontal'

subplot(3,1,3)
hold on

title('Simulation of labour supply')
%for sim=1:N_sim
sim=1;
% plot(L_disc_sim(sim,1:T),'k','Linewidth',1)
% plot(L_interp_sim(sim,1:T),'k:','Linewidth',1)
plot(L_sim_lin(sim,1:T),'k--','Linewidth',1)
%plot(L_det_sim(sim,1:T),'k-+','Linewidth',1)
h = legend('Linear','Location', 'best','Orientation','Vertical');
set(h,'fontsize',12,'Interpreter','Latex');%'Orientation', 'horizontal'

%%

disp('Standard devations of Z discrete, continuous, theoretical')
disp([mean(std(Z_sim')),mean(std(z_cont')),(sigmaepsilon^2/(1-rho^2))^0.5])

%%

for i=1:N_sim
    rho_emp(i)=corr(Z_sim(i,1:end-1)',Z_sim(i,2:end)');
    rho_emp_cont(i)=corr((z_cont(i,1:end-1))',(z_cont(i,2:end))');
    corr_Cy_lin(i)=corr(y_sim_lin(i,:)',c_sim_lin(i,:)');
    corr_invy_lin(i)=corr(y_sim_lin(i,:)',inv_sim_lin(i,:)');
    corr_Ly_lin(i)=corr(y_sim_lin(i,:)',L_sim_lin(i,:)');
    std_c(i) = std(c_sim_lin(i,:)');
    std_inv(i) = std(inv_sim_lin(i,:)');
    std_L(i) = std(L_sim_lin(i,:)');
    std_y(i) = std(y_sim_lin(i,:)');
end

mean_std_c = mean(std_c);
mean_std_inv = mean(std_inv);
mean_std_L = mean(std_L);
mean_std_y = mean(std_y);

%%

disp('Autocorrelation of Z discrete, continuous');
disp(mean(rho_emp_cont));
disp('Standard devations of Z discrete, continuous');
disp(std(rho_emp_cont));

%%

disp('Correlation with y of C,I,L')
% disp('VFI discrete')
% disp([Fill this in])
% disp('VFI interpolation')
% disp([Fill this in])
disp('Linear')
disp([mean(corr_Cy_lin), mean(corr_invy_lin), mean(corr_Ly_lin)]);
% disp('IRF')
% disp([Fill this in])

%%

disp('Standard deviations of C,I,L relative to y')
% disp('VFI discrete')
% disp([Fill this in])
% disp('VFI interpolation')
% disp([Fill this in])
disp('Linear')
disp([mean_std_c/mean_std_y, mean_std_inv/mean_std_y, mean_std_L/mean_std_y])
% disp('IRF')
% disp([Fill this in])
