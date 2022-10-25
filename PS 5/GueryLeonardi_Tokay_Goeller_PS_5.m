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

% ============
% options and convergence criteria
% ============

params.criter_V = 1e-7; % conv criterion for value function
params.T=100; % periods for transition

% ===========
% calculated parameters in 1
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

% param values
params.alpha = 0.4;
params.beta = 0.99;
params.sigma = 1.000001;
params.delta = 1;
params.kterm = 0;


%set guesses
kbar=((1/params.beta-1+params.delta)/(params.alpha*params.beta))^(1/(params.alpha-1));
cbar = kbar^params.alpha-params.delta*kbar;

params.k0 = 0.75*kbar;
params.cbar = cbar;

% Compile inputs
%kt = [ones(T,1)*0.75*kbar; 0];
%ct = [ones(T,1)*0.8*cbar; 0];

%% Broydens Method

% Broydens method Finite time T=10
T1 = 10;
kt1 = ones(T1,1)*0.75*kbar;
ct1 = ones(T1,1)*0.8*cbar;
x1 = [kt1; ct1];
jacob1 = ncgm_jacob(x1, params);
trans1 = ncgm_broyden(x1, jacob1, 1e-6, T1, params);

% Broydens method Finite time T=100
T2 = 100;
kt2 = ones(T2,1)*0.75*kbar;
ct2 = ones(T2,1)*0.8*cbar;
x2 = [kt2; ct2];
jacob2 = ncgm_jacob(x2, params);
trans2 = ncgm_broyden(x2, jacob2, 1e-6, T2, params);

% Broydens method Finite time T=200
T3 = 200;
kt3 = ones(T3,1)*0.75*kbar;
ct3 = ones(T3,1)*0.8*cbar;
x3 = [kt3; ct3];
jacob3 = ncgm_jacob(x3, params);
trans3 = ncgm_broyden(x3, jacob3, 1e-6, T3, params);

% Broydens method Infinite time - set T=?
T4 = 15;
kt4 = ones(T4,1)*0.75*kbar;
ct4 = ones(T4,1)*0.8*cbar;
x4 = [kt4; ct4];
jacob4 = ncgm_jacob(x4, params);
trans4 = ncgm_broyden(x4, jacob4, 1e-6, T4, params);

%% Multiple Shooting
T=20;
criter_V=1e-8;
c0_shots_T10 =ncgm_multiple_shooting(alpha, beta, delta, sigma, T, criter_V);

% T=100;
% criter_V=1e-8;
% c0_shots_T100=ncgm_multiple_shooting(alpha, beta, delta, sigma, T, criter_V);
% 
% T=200;
% criter_V=1e-8;
% c0_shots_T200=ncgm_multiple_shooting(alpha, beta, delta, sigma, T, criter_V);

%%
% =====================
% 3. Log-linearization
% =====================

% some ratios as function of parameters
ybar=kbar^alpha;
cbar=ybar-delta*kbar;
% need for matrix form of loglin
ckrat=cbar/kbar;

% a. write system as A E[y_t+1]+B y_t=0



% order c,k,z
A=[-sigma, beta *(alpha-1)*alpha*kbar^(alpha-1) ; 0 , 1 ];

B= [-sigma , 0 ; -ckrat,(1/beta)];

D = inv(A)*B;

% note that these are right-hand eigenvectors
[ ev lambda]=eig(D); %%%give the egein vectors and eigen values of D
aaa=inv(ev); %%% give the invert of the eigen vector matrix

% find eigenvalues equal or larger than one, and check that they equal the
% number of jump variables - in that case, set BKcond to 1

BKcond = 1;
eigen = [lambda(1,1) lambda(2,2)];
if all(eigen > 1) || all(eigen < 1)
    BKcond = 0;
end

%If the Blanchard Kahn condition is satisfied, we can find the expression
%of "alpha1"
if BKcond~=1
    disp('BK conditions not satisfied')
else
    bkev =find(abs(diag(lambda))>1);
    invP=aaa(bkev,:);%%Select the element of the invert of the vector matrix needed to compute the policy function
    polfunc= -invP(1,2)/invP(1);% you need to find the policy for consumption here : derived analytically 
end

% policy functions are in log deviations. so need to pre-multiply by cbar
% and add cbar to get levels


% calculate the deterministic transition  using the linearised policy
% functions, and law of motion
%%
% Broydens method Infinite time - set T=?
T4 = 51;
T=50;
kt4 = ones(T4,1)*0.75*kbar;
ct4 = ones(T4,1)*0.8*cbar;
x4 = [kt4; ct4];
jacob4 = ncgm_jacob(x4, params);
k_lina_num2 = ncgm_broyden(x4, jacob4, 1e-6, T4, params);
k_lina_num = k_lina_num2(1:T4);


k_lina_num=ncgm_broyden(x4, jacob4, 1e-6, 50, params); %values goes from k1 to kT (k(t) vector)
k_lin(1)=k_0;
for t=2:T
    c_lin(t-1)=polfunc*((k_lin(t-1)-kbar)/kbar)*cbar + cbar;
    k_lin(t)=k_lin(t-1)^alpha +(1-delta)*k_lin(t-1) - c_lin(t-1);
end

%Here is generating a sequence of k lagged by 1 unit of time 
k_lina(1)=k_lin(1,2);

k_lin_num = k_lina_num; 
k_lin_num(end)=[];% values goes from k0 to kT-1 (k(t-1))
k_lina_num(1)=[];

for t=2:(T)
    c_lina(t-1)=polfunc*((k_lina(t-1)-kbar)/kbar)*cbar + cbar;
    k_lina(t)=k_lina(t-1)^alpha +(1-delta)*k_lina(t-1) - c_lina(t-1);
end

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




