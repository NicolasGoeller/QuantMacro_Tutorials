% ==========================
% Problem set solving the stochastic neoclassical growth model different
% ways
% ==========================
% clear the workspace
clear
% close all figures
close all
addpath([xxx put your path here in '' xxx])

% ============
% parameters  - you may have to change this according to instructions
% ============
alpha=0.4; % capital share was named theat before
beta = 0.99; % discount factor
%rho = 0.95;   % persistence of TFP shock this doesnt eist in the PS
sigma = 1.00001; % CRRA coefficient (for 1 equals log, but need to replace the function, so set close to 1)
delta=1;

% ============
% options and convergence criteria
% ============

criter_V = 1e-7; % conv criterion for value function
T=150; % periods for transition

%mean of capital non-stochastic steady state
kbar=((1/beta-1+delta)/(alpha*beta))^(1/(alpha-1)); %% here we inserted a *beta
% initial level of capital in the transition
k_0=kbar*0.75; % you may have to change this from the problem set instructions


% ==============
% 1. analytical case delta=1, finite T and infinite T
% ==============

% a. analytical policies, finite and infinite horizon
% use fsolve

if delta==1 % only makes sense for delta=1
    k_analyt=zeros(1,T);
    
    k_analyt(1)=k_0;
    k_analyt_finite(1)=k_0;
    for t=2:T
        k_analyt(t)= [xxxx fill this in xxxx] ;
        k_analyt_finite(t)=[xxxx fill this in xxxx];
    end
end
%% This is nico's generator for system of equations

x0 = [0,0];
x = fsolve(@root2d,x0)

%%

% =====================
% 4. Log-linearization
% =====================

% some ratios as function of parameters
ybar=kbar^alpha;
cbar=ybar-delta*kbar;
% Dont know what these should be used for
%ckrat=cbar/kbar;
%R=1/beta;
%margprod=R-1+delta;

% a. write system as A E[y_t+1]+B y_t=0



% order c,k,z
A=[[xxxx fill this in xxxx] ];

B= [[xxxx fill this in xxxx] ];

D = [xxxx fill this in xxxx] ;

% note that these are right-hand eigenvectors
[ ev lambda]=eig(D);
aaa=inv(ev);

% find eigenvalues equal or larger than one, and check that they equal the
% number of jump variables - in that case, set BKcond to 1

[xxxx fill this in xxxx] 

if BKcond~=1
    disp('BK conditions not satisfied')
else
    [xxxx fill this in xxxx] 
    polfunc= [xxxx fill this in xxxx] % you need to find the policy for consumption here
end

% policy functions are in log deviations. so need to pre-multiply by cbar
% and add cbar to get levels

clev=cbar+cbar*polfunc(1)*(kgrid-kbar)/kbar;
kprimelev=kgrid'.^alpha+(1-delta)*kgrid'-clev(:,i);

% calculate the deterministic transition  using the linearised policy
% functions, and law of motion
k_lin(1)=k_0;
for t=2:T
    c_lin(t-1)=[xxxx fill this in xxxx] 
    k_lin(t)=[xxxx fill this in xxxx] 
end


% ==============
% 4. Solve deterministic sequence
% ==============

% as usual we need an initial guess for the sequences of k and c - here we
% just use steady state
x0=[kbar*ones(T,1);cbar*ones(T,1)];

zpath=ones(T,1);
zpath(1)=1.01;
for i=2:T-1
zpath(i)=exp(log(zpath(i-1)));
end

% now we need a function that returns errors of the equation system,
% constraining k(0)=k_0, and c(T)=cbar

% here we use Broyden's method

%Initial guess for jacobian - use finite difference at steady state
%sequences of k and c

clear dx J
for i=1:2*T
    dx = zeros(2*T,1);
    dx(i)=x0(i)*0.001;
    J(:,i)=[xxxx fill this in xxxx] ;
end

crit=1e-10;
x=x0;
f=rbc_obj(x,alpha,beta,sigma,delta,k_0,cbar,zpath);

while max(abs(f))>crit

dx = [xxxx fill this in xxxx] 
x=x+dx;
f = rbc_obj(x,alpha,beta,gamma,delta,k_0,cbar,ones(T,1));
J = [xxxx fill this in xxxx] ;

end

k_trans_br=x(1:T,1);
c_trans_br=x(T+1:2*T,1);

% or we just let matlab solve it

[x,fsol]=fsolve([xxxx fill this in xxxx] );

k_trans=x(1:T,1);
c_trans=x(T+1:2*T,1);



% ==============
% Figures
% ==============
% plot policy function
figure(1)
% levels
subplot(2,1,1)
title('Policy functions')
hold on
[xxxx fill this in xxxx] 

% plot the transition
figure(2)
title('Simulated transition - deterministic')
hold on
[xxxx fill this in xxxx] 
if delta==1 && abs(sigma-1)<0.001
    h = legend([xxxx fill this in xxxx] ,'Location', 'best','Orientation','Vertical');
else
    h = legend([xxxx fill this in xxxx] ,'Location', 'best','Orientation','Vertical');
end
set(h,'fontsize',12,'Interpreter','Latex');%'Orientation', 'horizontal'
%%
params.alpha = 0.4;
params.beta = 0.99;
params.sigma = 1;
params.delta = 1;
params.kterm = 0;
kbar=((1/params.beta-1+params.delta)/(params.alpha*params.beta))^(1/(params.alpha-1));
cbar = kbar^params.alpha-params.delta*kbar;
kt = ones(9,1)*0.75*kbar;
ct = ones(10,1)*0.8*cbar;
x = [kt;ct];
%x = [0.75*kbar;0.1;0.1;0.1;0.1; 1;1;1;1;1;1];

fsolve(ncgw_seq(),x)

%% Broydens method

%symbolise all unknowns
T=10
syms x [1 2*T+1]
x
%create numerical Jacobian

%%
function F = ncgw_seq()
T = 9;
par.alpha = 0.4;
par.beta = 0.99;
par.sigma = 1;
par.delta = 1;
par.kterm = 0;

%k = x(1:T,1);
%c = x(T+1:end,1);
syms x [1 2*T+1];

for t = 1:T
    %t
    funcit = [2*(t-1)+1, 2*(t-1)+2];
    if t == T
        F(2*T+1) = x(T+t) + par.kterm - x(t)^par.alpha - (1-par.delta)*x(t);
        break
    end
    F(funcit(1)) = x(T+t) + x(t+1) -x(t)^par.alpha - (1-par.delta)*x(t);
    F(funcit(2)) = x(T+t)^(par.sigma*(-1)) - par.beta*(par.alpha*x(t+1)^(par.alpha - 1) + ...
        1 - par.delta)*x(T+t+1)^(par.sigma*(-1));
end
end

