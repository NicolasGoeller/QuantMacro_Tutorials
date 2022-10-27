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
params.maxiter = 1000; %maxiter for broyden

% ============
% options and convergence criteria
% ============

params.criter_V = 1e-10; % conv criterion for value function
params.T=20; % periods for transition
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

% start Broyden Function (including initial jacobian)

%h= tank_error(x_tank_init, z_nu, params);

%  ==============================
%% Option I : Broydens Method for TANK model
%  ==============================

% Broydens method with monetary shock epsilon_nu at t=0
% T1 = 10;
% jacob_tank = tank_jacob(x_tank_init,z_nu_shock, params);
% trans_tank_nu_shock = TANK_broyden(x_tank_init, z_nu_shock, jacob_tank, params.maxiter, params);  

%investment_tank_nu_shock= params.phipi*trans_tank_nu_shock(T+1:2T,1) + nu_val_nu_shock(2:end);

x_fsolve= fsolve(@(x) tank_error(x,z_nu_shock,params),x_tank_init);

% Input dictionary
% trans_tank_nu_shock(1:T) output (y)
% trans_tank_nu_shock(T+1:2T) inflation (pi)
% trans_tank_nu_shock(2T+1:3T) consumption of smoothers (cS)
% trans_tank_nu_shock(3T+1:end) consumption of Hand-to-mouth (cH)

%for large T
% subplot(3,1,1);
% plot(x_fsolve(1:T/5,1));
% subplot(3,1,2);
% plot(x_fsolve(T+1:1.2*T,1));
% subplot(3,1,3);
% plot(x_fsolve(2*T+1:2.2*T,1));


subplot(3,1,1);
plot(x_fsolve(1:T,1));
subplot(3,1,2);
plot(x_fsolve(T+1:2*T,1));
subplot(3,1,3);
plot(x_fsolve(2*T+1:3*T,1));




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


%%

% Broydens method with productivity shock epsilon_a at t=0
% T1 = 10;
% jacob_tank = tank_jacob(x_tank_init,z_a, params);
% trans1 = ncgm_broyden(x1, jacob1, 1e-6, T1, params);



%% IRF plotting TANK with monetary shock


figure('Monetary Shock','Impulse Response functions'); % this is not asked

subplot(3,2,1);
title('Consumption','Interpreter','Latex','fontsize',13);
hold on;
% plot([c_disc(:,plotz)])
% plot([c_interp(:,plotz)])
plot([c_sim_lin(:,plotz)]);
ylabel('Consumption','Interpreter','Latex','fontsize',13);

subplot(3,2,2);
title('Labor supply','Interpreter','Latex','fontsize',13);
hold on;
% plot([L_disc(:,plotz)])
% plot([L_interp(:,plotz)])
plot([L_sim_lin(:,plotz)]);
ylabel('Labor supply','Interpreter','Latex','fontsize',13);

subplot(3,2,3);
title('Inflation','Interpreter','Latex','fontsize',13);
hold on;
% plot([kprime_disc(:,plotz)])
% plot([kprime_interp(:,plotz)])
plot([k_sim_lin(:,plotz)]);
ylabel('Inflation','Interpreter','Latex','fontsize',13);

subplot(3,2,4);
title('Interest rate','Interpreter','Latex','fontsize',13);
hold on;
% plot([kprime_disc(:,plotz)])
% plot([kprime_interp(:,plotz)])
plot([k_sim_lin(:,plotz)]);
ylabel('Interest rate','Interpreter','Latex','fontsize',13);

subplot(3,2,5);
title('Output','Interpreter','Latex','fontsize',13);
hold on;
% plot([kprime_disc(:,plotz)])
% plot([kprime_interp(:,plotz)])
plot([k_sim_lin(:,plotz)]);
ylabel('Output','Interpreter','Latex','fontsize',13);

subplot(3,2,6);
title('Shock','Interpreter','Latex','fontsize',13);
hold on;
% plot([kprime_disc(:,plotz)])
% plot([kprime_interp(:,plotz)])
plot([k_sim_lin(:,plotz)]);
ylabel('Shock','Interpreter','Latex','fontsize',13);

h = legend('Linearised','Location', 'best','Orientation','Vertical');
set(h,'fontsize',13,'Interpreter','Latex');%'Orientation', 'horizontal'





%%
% ==============
% 4. Figures
% ==============
% plot policy function
% plot policy function
figure(1)
% levels
subplot(2,1,1); subplot
title('Policy functions');
hold on;
plot(k_lin,k_lina); %% here is a graph of k'(k)
plot(k_lin_num',k_lina_num');
xlabel('k_{t}','FontSize',10);
ylabel('k_{t+1}','FontSize',10);
subplot (2,1,2);
title('Percentage of difference');
plot(k_lin,abs((k_lina-k_lina_num)/k_lina));
xlabel('k_{t}','FontSize',10);
ylabel('% of difference','FontSize',10);
%plot(k_lin_num,k_lina_num);
hold off

%%
% plot the transition
% Time variable
figure(2)
subplot(2,1,1)
T = 100;
crop_off = 70;
x = 1:1:crop_off;
figure(2)
title('Simulated deterministic transition - infinite time')

hold on
plot(x, k_analyt(1:crop_off)', x, trans2(1:crop_off)), xlabel('Time steps'), ylabel('Capital level');

subplot(2,1,2);
hold on
plot(x, consum_analyt(1:crop_off)', x, trans2(T+1:T+crop_off)), xlabel('Time steps'), ylabel('Consumption level');
hold off

if delta==1 && abs(sigma-1)<0.001
    h = legend('Analytical solution', 'Broydens method' ,'Location', 'bestoutside','Orientation','Vertical');
    h.Title.String = 'Analytically solvable paths';
else
    h = legend(['Analytical solution', 'Broydens method'] ,'Location', 'bestoutside','Orientation','Vertical');
    h.Title.String = 'Non-linearised paths';
end
set(h,'fontsize',12,'Interpreter','Latex');%'Orientation', 'horizontal'




