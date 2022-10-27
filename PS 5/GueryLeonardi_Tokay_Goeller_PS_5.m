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
x_tank_init=ones(3*params.T,1);
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
jacob_tank1 = tank_jacob(x_tank_init,z_nu_shock, params);
trans_tank_nu_shock1 = TANK_broyden(x_tank_init, z_nu_shock, jacob_tank1, params.maxiter, params);

output1 = trans_tank_nu_shock1(1:T);
inflation1 = trans_tank_nu_shock1(T+1:2*T);
consumption1 = trans_tank_nu_shock1(2*T+1:end);
interest1 = params.phipi*inflation1 + z_nu_shock(T+1:end);
labor1 = output1 - z_nu_shock(1:T);
shock1 = z_nu_shock(T+1:end);

% Broydens method with production shock epsilon_a at t=0
jacob_tank2 = tank_jacob(x_tank_init,z_a_shock, params);
trans_tank_nu_shock2 = TANK_broyden(x_tank_init, z_a_shock, jacob_tank2, params.maxiter, params);  

output2 = trans_tank_nu_shock2(1:T);
inflation2 = trans_tank_nu_shock2(T+1:2*T);
consumption2 = trans_tank_nu_shock2(2*T+1:end);
interest2 = params.phipi*inflation2 + z_nu_shock(T+1:end);
labor2 = output2 - z_nu_shock(1:T);
shock2 = z_nu_shock(1:T);
%%
% Input dictionary
% trans_tank_nu_shock(1:T) output (y)
% trans_tank_nu_shock(T+1:2T) inflation (pi)
% trans_tank_nu_shock(2T+1:3T) consumption of smoothers (cS)
% trans_tank_nu_shock(3T+1:end) consumption of Hand-to-mouth (cH)

subplot(2,2,1);
plot(trans_tank_nu_shock(1:T,1));
subplot(2,2,2);
plot(trans_tank_nu_shock(T+1:2*T,1));
subplot(2,2,3);
plot(trans_tank_nu_shock(2*T+1:3*T,1));
subplot(2,2,4);
plot(trans_tank_nu_shock(3*T+1:4*T,1));


%  ==============================
%% Option II : Log-lin method for TANK models (eq system is already linear)
%  ==============================



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

subplot(3,2,1);
title('Consumption','Interpreter','Latex','fontsize',13);
hold on;
% plot([c_interp(:,plotz)])
plot([consumption1(1:T/5)]);
ylabel('Consumption','Interpreter','Latex','fontsize',13);

subplot(3,2,2);
title('Labor supply','Interpreter','Latex','fontsize',13);
hold on;
% plot([L_disc(:,plotz)])
% plot([L_interp(:,plotz)])
plot([labor1(1:T/5)]);
ylabel('Labor supply','Interpreter','Latex','fontsize',13);

subplot(3,2,3);
title('Inflation','Interpreter','Latex','fontsize',13);
hold on;
% plot([kprime_disc(:,plotz)])
% plot([kprime_interp(:,plotz)])
plot([inflation1(1:T/5)]);
ylabel('Inflation','Interpreter','Latex','fontsize',13);

subplot(3,2,4);
title('Interest rate','Interpreter','Latex','fontsize',13);
hold on;
% plot([kprime_disc(:,plotz)])
% plot([kprime_interp(:,plotz)])
plot([interest1(1:T/5)]);
ylabel('Interest rate','Interpreter','Latex','fontsize',13);

subplot(3,2,5);
title('Output','Interpreter','Latex','fontsize',13);
hold on;
% plot([kprime_disc(:,plotz)])
% plot([kprime_interp(:,plotz)])
plot([output1(1:T/5)]);
ylabel('Output','Interpreter','Latex','fontsize',13);

subplot(3,2,6);
title('Shock','Interpreter','Latex','fontsize',13);
hold on;
% plot([kprime_disc(:,plotz)])
% plot([kprime_interp(:,plotz)])
plot([shock1(1:T/5)]);
ylabel('Shock','Interpreter','Latex','fontsize',13);

h = legend('Linearised','Location', 'best','Orientation','Vertical');
set(h,'fontsize',13,'Interpreter','Latex');%'Orientation', 'horizontal'

%% IRF plotting TANK with productivity shock

figure('Productivity Shock','Impulse Response Functions'); % this is not asked

subplot(3,2,1);
title('Consumption','Interpreter','Latex','fontsize',13);
hold on;
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
plot(z_nu_shock(T+1:end))
% Broydens method with monetary shock epsilon_a at t=0 for RANK model

rank_params = params;
rank_params.taud = 0;
rank_params.lambda = 0.000001;
jacob_tank3 = tank_jacob(x_tank_init,z_nu_shock, params);
trans_tank_nu_shock3 = TANK_broyden(x_tank_init, z_nu_shock, jacob_tank3, params.maxiter, params);  

output3 = trans_tank_nu_shock3(1:T);
inflation3 = trans_tank_nu_shock3(T+1:2*T);
consumption3 = trans_tank_nu_shock3(2*T+1:end);
interest3 = params.phipi*inflation3 + z_nu_shock(T+1:end);
labor3 = output3 - z_nu_shock(1:T);
shock3 = z_nu_shock(1:T);

figure('Monetary Shock in TANK and RANK','Impulse Response Functions'); % this is not asked

subplot(3,2,1);
title('Consumption','Interpreter','Latex','fontsize',13);
hold on;
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
