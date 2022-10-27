% ==========================
% Problem set solving the stochastic neoclassical growth model different
% ways
% ==========================
% clear the workspace
clear
% close all figures
close all

% ============
% parameters  - you may have to change this according to instructions
% ============

params.alpha=0.4; % capital share was named theat before
params.beta = 0.99; % disfcount factor
params.sigma = 1.0001; % CRRA coefficient (for 1 equals log, but need to replace the function, so set close to 1)
params.delta=1;
params.phipi = 1.5 ; %parameter of the monetary policy feedback rule
params.varphi = 1; %intratemporal first order condition for labor supply parameter
params.theta = 0.666667; %probability of a firm not resetting its price
params.rhoa = 0.95 ; %AR1 parameter for labor productivity
params.rhonu = 0.5 ; %AR1 parameter for monetary policy shocks
params.lambda = 0.3 ; %share of hands-to-mouth households
params.taud = 0.1; %firm profit redistribution share
params.vartheta = 0.666667; %Frisch elasticity of labor supply

% =============
% Shock parameters
% =============

params.epsilon_nu = 0.005; %monetary policy shock
params.epsilon_a = 0.005; %productivity shock
params.maxiter = 100; %maxiter for broyden

% ============
% options and convergence criteria
% ============

params.criter_V = 1e-10; % conv criterion for value function
params.T=100; % periods for transition
T=params.T;

% ===========
% calculated parameters at steady state
% ===========



%  ===========================
%% Blueprint for solution
%  ===========================

% Define equilibrium error function
% - TANK model
% - RANK model

% Set initial guess vectors to eq-ss values
x_tank_init=zeros(3*params.T,1);
x_rank_init=zeros(3*params.T,1);


% Set perturbance vectors to assigned t=0 shock
epsi_nu_shock = [params.epsilon_nu; zeros(params.T-1,1)];
epsi_a_no_shock = zeros(params.T,1);
epsi_nu_no_shock = zeros(params.T,1);
epsi_a_shock = [params.epsilon_a; zeros(params.T-1,1)];

nu_val_nu_shock = zeros(T,1); %vector of values of nu given its persistence and the shock at t=0
nu_val_nu_no_shock = zeros(T,1); %vector of values of nu without shock (always 0)
a_val_a_shock = zeros(T,1); %vector of values of a given its persistence and the shock at t=0
a_val_a_no_shock = zeros(T,1); %vector of values of a without shock (always 0)

nu_val_nu_shock(1)=epsi_nu_shock(1); %initializing nu value at t=0 (shock of epsi)
a_val_a_shock(1)=epsi_a_shock(1); %initalizing a value at t=0 (shock of epsi)

for t=1:T-1
    nu_val_nu_shock(t+1)=params.rhonu*nu_val_nu_shock(t)+epsi_nu_shock(t+1); %constructing value given AR1 process
    a_val_a_shock(t+1)=params.rhoa*a_val_a_shock(t)+epsi_a_shock(t+1); %AR1 process of a
end

z_nu_shock = [epsi_a_no_shock; nu_val_nu_shock];
z_a_shock = [a_val_a_shock; epsi_nu_no_shock];

%  ==============================
%% Option I : Broydens Method for TANK model
%  ==============================

% Broydens method with monetary shock epsilon_nu at t=0
% T1 = 10;
% jacob_tank = tank_jacob(x_tank_init,z_nu_shock, params);
% trans_tank_nu_shock = TANK_broyden(x_tank_init, z_nu_shock, jacob_tank, params.maxiter, params);  
% 
% 
% %investment_tank_nu_shock= params.phipi*trans_tank_nu_shock(T+1:2T,1) + nu_val_nu_shock(2:end);
% 
% subplot(3,1,1);
% plot(trans_tank_nu_shock(1:T,1));
% subplot(3,1,2);
% plot(trans_tank_nu_shock(T+1:2*T,1));
% subplot(3,1,3);
% plot(trans_tank_nu_shock(2*T+1:3*T,1));

output1 = trans_tank_nu_shock1(1:T);
inflation1 = trans_tank_nu_shock1(T+1:2*T);
consumption1 = trans_tank_nu_shock1(2*T+1:end);
interest1 = params.phipi*inflation1 + z_nu_shock(T+1:end);
labor1 = output1 - z_nu_shock(1:T);
shock1 = z_nu_shock(T+1:end);
%%

% Broydens method with production shock epsilon_a at t=0
jacob_tank2 = tank_jacob(x_tank_init,z_a_shock, params);
trans_tank_nu_shock2 = TANK_broyden(x_tank_init, z_a_shock, jacob_tank2, params.maxiter, params);  

output2 = trans_tank_nu_shock2(1:T);
inflation2 = trans_tank_nu_shock2(T+1:2*T);
consumption2 = trans_tank_nu_shock2(2*T+1:end);
interest2 = params.phipi*inflation2 + z_nu_shock(T+1:end);
labor2 = output2 - z_nu_shock(1:T);
shock2 = z_nu_shock(1:T);
%% Monetary shock
x_fsolve_nushock= fsolve(@(x) tank_error(x,z_nu_shock,params),x_tank_init);

cH_nushock=(x_fsolve_nushock(1:T)-(1-params.lambda)*x_fsolve_nushock(2*T+1:3*T))/params.lambda;
interest_nushock=params.phipi*x_fsolve_nushock(T+1:2*T)+nu_val_nu_shock;
n_nushock = x_fsolve_nushock(1:T) - z_nu_shock(1:T);
%% Productivity shock
x_fsolve_ashock= fsolve(@(x) tank_error(x,z_a_shock,params),x_tank_init);

cH_ashock=(x_fsolve_ashock(1:T)-(1-params.lambda)*x_fsolve_ashock(2*T+1:3*T))/params.lambda;
interest_ashock=params.phipi*x_fsolve_ashock(T+1:2*T)+nu_val_nu_no_shock;
n_ashock = x_fsolve_ashock(1:T) - a_val_a_shock;

%%
% Input dictionary
% x_fsolve_nushock(1:T) output (y)
% x_fsolve_nushock(T+1:2T) inflation (pi)
% x_fsolve_nushock(2T+1:3T) consumption of smoothers (cS)


%for large T
subplot(3,1,1);
plot(x_fsolve_nushock(1:T/5,1));
subplot(3,1,2);
plot(x_fsolve_nushock(T+1:1.2*T,1));
subplot(3,1,3);
plot(x_fsolve_nushock(2*T+1:2.2*T,1));


% subplot(3,1,1);
% plot(x_fsolve_nushock(1:T,1));
% subplot(3,1,2);
% plot(x_fsolve_nushock(T+1:2*T,1));
% subplot(3,1,3);
% plot(x_fsolve_nushock(2*T+1:3*T,1));




%  ==============================
%% Option II : Log-lin method for TANK models (eq system is already linear)
%  ==============================

%once the system is reduced, we have a set of 5 eq and 5 variables with 2
%AR1 process, we can then compute a log linear solution

kappa= (params.sigma + params.vartheta)*(1-params.theta)*(1-params.beta*params.theta)/params.beta;

% a. write system as A_1 E[x_t+1]=A_2  x_t+B_1 epsilon_t
% order nu, a, y, pi, cS
a=[ 0 ,       0       ,       0       , 1/params.sigma,       0           ;
    0 , 0, 1, 1/params.sigma,0;
    0,0,0,params.beta,0;
    0,1,0,0,0;
    1,0,0,0,0];



b=[1/params.sigma,0,0, params.phipi/params.sigma,1;
    1/params.sigma,0,1,params.phipi/params.sigma,0;
    0, kappa, -kappa, 1, 0;
    0,params.rhoa,0,0,0;
    params.rhonu, 0, 0, 0, 0];

%%

% re-order to have state variables first
% order nu, a, y, pi, cS
nk=2;
[f,p] = solab(a,b,nk);

%%
% extract cons and lab policies
y_polfunc=f(1,:);
pi_polfunc=f(2,:);
cS_polfunc=f(3,:);
LOM_nu=p(1,:);
LOM_a =p(2,:);

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

y_loglin=zeros(T,1);
pi_loglin=zeros(T,1);
cS_loglin=zeros(T,1);

for t=2:T
    y_loglin(t)=y_polfunc
end


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






%% Plotting error in equations
error_test_tank= tank_error(trans_tank_nu_shock,z_nu_shock,params);
subplot(2,2,1);
plot(error_test_tank(1:T,1));
subplot(2,2,2);
plot(error_test_tank(T+1:2*T,1));
subplot(2,2,3);
plot(error_test_tank(2*T+1:3*T,1));
subplot(2,2,4);
plot(error_test_tank(3*T+1:4*T,1));


%% IRF plotting TANK with monetary shock

figure(1);
title('Monetary Shock - Impulse Response Functions'); % this is not asked

subplot(4,2,1);
title('Smoother Consumption','Interpreter','Latex','fontsize',13);
hold on;
% plot([c_interp(:,plotz)])
plot([x_fsolve_nushock(2*T+1:2.2*T,1)]);
%ylabel('SConsumption','Interpreter','Latex','fontsize',13);

subplot(4,2,2);
title('Hand-Mouth Consumption','Interpreter','Latex','fontsize',13);
hold on;
% plot([c_interp(:,plotz)])
plot([cH_nushock(1:T/5)]);
%ylabel('H Consumption','Interpreter','Latex','fontsize',13);

subplot(4,2,3);
title('Labor supply','Interpreter','Latex','fontsize',13);
hold on;
% plot([L_disc(:,plotz)])
% plot([L_interp(:,plotz)])
plot([n_nushock(1:T/5)]);
%ylabel('Labor supply','Interpreter','Latex','fontsize',13);

subplot(4,2,4);
title('Inflation','Interpreter','Latex','fontsize',13);
hold on;
% plot([kprime_disc(:,plotz)])
% plot([kprime_interp(:,plotz)])
plot([x_fsolve_nushock(T+1:1.2*T,1)]);
%ylabel('Inflation','Interpreter','Latex','fontsize',13);

subplot(4,2,5);
title('Interest rate','Interpreter','Latex','fontsize',13);
hold on;
% plot([kprime_disc(:,plotz)])
% plot([kprime_interp(:,plotz)])
plot([interest_nushock(1:T/5)]);
%ylabel('Interest rate','Interpreter','Latex','fontsize',13);

subplot(4,2,6);
title('Output','Interpreter','Latex','fontsize',13);
hold on;
% plot([kprime_disc(:,plotz)])
% plot([kprime_interp(:,plotz)])
plot([x_fsolve_nushock(1:T/5,1)]);
%ylabel('Output','Interpreter','Latex','fontsize',13);

subplot(4,2,7);
title('Shock','Interpreter','Latex','fontsize',13);
hold on;
% plot([kprime_disc(:,plotz)])
% plot([kprime_interp(:,plotz)])
plot([z_nu_shock(T+1:1.2*T)]);
%ylabel('Shock','Interpreter','Latex','fontsize',13);

% h = legend('Linearised','Location', 'best','Orientation','Vertical');
% set(h,'fontsize',13,'Interpreter','Latex');%'Orientation', 'horizontal'

%% IRF plotting TANK with productivity shock

figure(2);
title('Productivity Shock - Impulse Response Functions'); % this is not asked

subplot(4,2,1);
title('Smoother Consumption','Interpreter','Latex','fontsize',13);
hold on;
% plot([c_interp(:,plotz)])
plot([x_fsolve_ashock(2*T+1:2.2*T,1)]);
%ylabel('SConsumption','Interpreter','Latex','fontsize',13);

subplot(4,2,2);
title('Hand-Mouth Consumption','Interpreter','Latex','fontsize',13);
hold on;
% plot([c_interp(:,plotz)])
plot([cH_ashock(1:T/5)]);
%ylabel('H Consumption','Interpreter','Latex','fontsize',13);

subplot(4,2,3);
title('Labor supply','Interpreter','Latex','fontsize',13);
hold on;
% plot([L_disc(:,plotz)])
% plot([L_interp(:,plotz)])
plot([n_ashock(1:T/5)]);
%ylabel('Labor supply','Interpreter','Latex','fontsize',13);

subplot(4,2,4);
title('Inflation','Interpreter','Latex','fontsize',13);
hold on;
% plot([kprime_disc(:,plotz)])
% plot([kprime_interp(:,plotz)])
plot([x_fsolve_ashock(T+1:1.2*T,1)]);
%ylabel('Inflation','Interpreter','Latex','fontsize',13);

subplot(4,2,5);
title('Interest rate','Interpreter','Latex','fontsize',13);
hold on;
% plot([kprime_disc(:,plotz)])
% plot([kprime_interp(:,plotz)])
plot([interest_ashock(1:T/5)]);
%ylabel('Interest rate','Interpreter','Latex','fontsize',13);

subplot(4,2,6);
title('Output','Interpreter','Latex','fontsize',13);
hold on;
% plot([kprime_disc(:,plotz)])
% plot([kprime_interp(:,plotz)])
plot([x_fsolve_ashock(1:T/5,1)]);
%ylabel('Output','Interpreter','Latex','fontsize',13);

subplot(4,2,7);
title('Shock','Interpreter','Latex','fontsize',13);
hold on;
% plot([kprime_disc(:,plotz)])
% plot([kprime_interp(:,plotz)])
plot([z_a_shock(1:0.2*T)]);
%ylabel('Shock','Interpreter','Latex','fontsize',13);

% h = legend('Linearised','Location', 'best','Orientation','Vertical');
% set(h,'fontsize',13,'Interpreter','Latex');%'Orientation', 'horizontal'

%%
% Broydens method with monetary shock epsilon_a at t=0 for RANK model

rank_params = params;
rank_params.taud = 0;
rank_params.lambda = 0.000001;
%jacob_tank3 = tank_jacob(x_tank_init,z_nu_shock, params);
%trans_tank_nu_shock3 = TANK_broyden(x_tank_init, z_nu_shock, jacob_tank3, params.maxiter, params);  

x_fsolve_rank= fsolve(@(x) tank_error(x,z_nu_shock,params),x_tank_init);

cH_rank=(x_fsolve_rank(1:T)-(1-params.lambda)*x_fsolve_rank(2*T+1:3*T))/params.lambda;
interest_rank=params.phipi*x_fsolve_rank(T+1:2*T)+nu_val_nu_shock;
n_rank = x_fsolve_rank(1:T) - a_val_a_no_shock;

%%

figure(3)
title('Monetary Shock in TANK and RANK','Impulse Response Functions'); % this is not asked

subplot(4,2,1);
title('Smoother Consumption','Interpreter','Latex','fontsize',13);
hold on;
plot([x_fsolve_nushock(2*T+1:2.2*T,1)])
plot([x_fsolve_rank(2*T+1:2.2*T,1)]);
%ylabel('SConsumption','Interpreter','Latex','fontsize',13);

subplot(4,2,2);
title('Hand-Mouth Consumption','Interpreter','Latex','fontsize',13);
hold on;
plot([cH_nushock(1:T/5)]);
plot([cH_rank(1:T/5)]);
%ylabel('H Consumption','Interpreter','Latex','fontsize',13);

subplot(4,2,3);
title('Labor supply','Interpreter','Latex','fontsize',13);
hold on;
plot([n_nushock(1:T/5)]);
plot([n_rank(1:T/5)]);
%ylabel('Labor supply','Interpreter','Latex','fontsize',13);

subplot(4,2,4);
title('Inflation','Interpreter','Latex','fontsize',13);
hold on;
plot([x_fsolve_nushock(T+1:1.2*T,1)]);
plot([x_fsolve_rank(T+1:1.2*T,1)]);
%ylabel('Inflation','Interpreter','Latex','fontsize',13);

subplot(4,2,5);
title('Interest rate','Interpreter','Latex','fontsize',13);
hold on;
plot([interest_nushock(1:T/5)]);
plot([interest_rank(1:T/5)]);
%ylabel('Interest rate','Interpreter','Latex','fontsize',13);

subplot(4,2,6);
title('Output','Interpreter','Latex','fontsize',13);
hold on;
plot([x_fsolve_nushock(1:T/5,1)]);
plot([x_fsolve_rank(1:T/5,1)]);
%ylabel('Output','Interpreter','Latex','fontsize',13);

subplot(4,2,7);
title('Shock','Interpreter','Latex','fontsize',13);
hold on;
plot([z_nu_shock(T+1:1.2*T)]);
%ylabel('Shock','Interpreter','Latex','fontsize',13);

% h = legend('Linearised','Location', 'best','Orientation','Vertical');
% set(h,'fontsize',13,'Interpreter','Latex');%'Orientation', 'horizontal'
