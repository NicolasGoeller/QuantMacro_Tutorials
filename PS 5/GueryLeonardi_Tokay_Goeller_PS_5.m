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
params.vartheta = %Frisch elasticity of labor supply

% =============
% Shock parameters
% =============




% ============
% options and convergence criteria
% ============

params.criter_V = 1e-7; % conv criterion for value function
params.T=100; % periods for transition

% ===========
% calculated parameters at steady state
% ===========



% ==============
% 2. Solving deterministic equation systems
% ==============

% a&b. analytical policies, finite and infinite horizon; delta = 1

%% Analytical part and sequence
alphabeta = zeros(1,T); %creating alpahbeta matrix with alphabeta(i)=(alpha*beta)^i
alphabeta(1)= 1;
for t=2:T
    alphabeta(t)=alphabeta(t-1)*alpha*beta;
end

if delta==1 % only makes sense for delta=1
    k_analyt=zeros(1,T);
    consum_analyt=zeros(1,T); %vectors with consumption at t
    consum_analyt_finite=zeros(1,T);
    
    k_analyt(1)=k_0;
    k_analyt_finite(1)=k_0;
    for t=2:T
        % use the recursive expression of k calculated in Part1.3
        k_analyt(t)= (1-1/(sum(alphabeta(1:T-t+2))))*k_analyt(t-1)^alpha ; %sum to T-t+2 because we begin at 1 and not 0
        k_analyt_finite(t)= (1-(1/sum(alphabeta(1:T-t+2))))*k_analyt_finite(t-1)^alpha;
        % k_analyt_finite(t)=[xxxx fill this in xxxx];
        consum_analyt(t-1)=k_analyt(t-1)^alpha - k_analyt(t);
        consum_analyt_finite(t-1)=k_analyt_finite(t-1)^alpha - k_analyt_finite(t);
    end
    consum_analyt_finite(T)=k_analyt_finite(T);
    consum_analyt(T)=k_analyt(T);
end

%%

% c. numerical solution algorithms to above problems

%set guesses
kbar=((1/params.beta-1+params.delta)/(params.alpha*params.beta))^(1/(params.alpha-1));
cbar = kbar^params.alpha-params.delta*kbar;

params.k0 = 0.75*kbar;
params.cbar = cbar;

% Compile inputs
%kt = [ones(T,1)*0.75*kbar; 0];
%ct = [ones(T,1)*0.8*cbar; 0];

%% Blueprint for solution

% Define equilibrium error function
% - TANK model
% - RANK model

% Set initial guess vectors to eq-ss values

% Set perturbance vectors to assigned t=0 shock
epsi_ = [0.005; zeros(T,1)];
epsi_a = zeros(T+1,1);

% start Broyden Function (including initial jacobian)



%% Broydens Method

% Broydens method Finite time T=10
T1 = 10;
kt1 = ones(T1,1)*0.75*kbar;
ct1 = ones(T1,1)*0.8*cbar;
x1 = [kt1; ct1];
jacob1 = ncgm_jacob(x1, params);
trans1 = ncgm_broyden(x1, jacob1, 1e-6, T1, params);



%% IRF plotting


figure('Name','Impulse Response functions'); % this is not asked

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




