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
alpha=0.4; % capital share was named theat before
beta = 0.99; % discount factor
%rho = 0.95;   % persistence of TFP shock this doesnt eist in the PS
sigma = 1.0001; % CRRA coefficient (for 1 equals log, but need to replace the function, so set close to 1)
delta=1;

% ============
% options and convergence criteria
% ============

criter_V = 1e-7; % conv criterion for value function
T=10; % periods for transition

%mean of capital non-stochastic steady state
kbar=((1/beta-1+delta)/(alpha))^(1/(alpha-1)); %% here we inserted a *beta
cbar= kbar^alpha-delta*kbar
% initial level of capital in the transition
k_0=kbar*0.75; % you may have to change this from the problem set instructions


%% ==============
% 1. analytical case delta=1, finite T and infinite T
% ==============

% a. analytical policies, finite and infinite horizon
% use fsolve
%% Result of the analytical part

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

%% ==============
% 4. Solve deterministic sequence
% ==============

% as usual we need an initial guess for the sequences of k and c - here we
% just use steady state
x0=[kbar*ones(T,1);cbar*ones(T,1)];
%%
% zpath=ones(T,1);
% zpath(1)=1.01;
% for i=2:T-1
%     zpath(i)=exp(log(zpath(i-1)));
% end

% now we need a function that returns errors of the equation system,
% constraining k(0)=k_0, and c(T)=cbar

% here we use Broyden's method

%Initial guess for jacobian - use finite difference at steady state
%sequences of k and c
%rbc_obj(x,alpha,beta,sigma,delta,kbar,cbar)
%%
clear dx J
for i=1:2*T
    dx = zeros(2*T,1);
    dx(i)=x0(i)*0.001;
    diffF = rbc_obj_start(x0+dx,alpha,beta,sigma,delta,kbar,cbar)-rbc_obj_start(x0,alpha,beta,sigma,delta,kbar,cbar);
    diffX= dx(i);
    J(:,i)=[diffF/diffX] ;
end
clear x
%%
crit=1e-4;
x=x0;
f0=rbc_obj_start(x,alpha,beta,sigma,delta,k_0,cbar);
f=rbc_obj_start(x,alpha,beta,sigma,delta,k_0,cbar);
%%
while max(abs(f))>crit

%dx = [xxxx fill this in xxxx] 
%     dx = zeros(2*T,1);
%     for i=1:2*T
%         for j=1:2*T
%             dx(i)= dx(i)-J(i,j)^(-1)*f(i);
%         end
%     end
    f0 = rbc_obj_start(x,alpha,beta,sigma,delta,k_0,cbar);
    dx = -J\f; 
    %dx = dx*(x'*x)/(10*(dx'*dx));
    %while x+dx<1e-6
      %  dx=dx/2;
    %end
    xn1=max(ones(2*T,1)*1e-8,x+dx);
    dx=xn1-x;
    f = rbc_obj_start(xn1,alpha,beta,sigma,delta,k_0,cbar);
    %J = [] ;
    J = J + ((f-f0)-J*dx)*dx'/(dx'*dx);
    %B = B + ((y' - B*s)*s')/(s'*s);
    x=xn1;


end


%%
%k_trans_br=x(1:T,1);
%c_trans_br=x(T+1:2*T,1);

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
% param values
params.alpha = 0.4;
params.beta = 0.99;
params.sigma = 1.000001;
params.delta = 1;
params.kterm = 0;
params.cterm = 0;

%set guesses
kbar=((1/params.beta-1+params.delta)/(params.alpha*params.beta))^(1/(params.alpha-1));
cbar = kbar^params.alpha-params.delta*kbar;
T = 10;

params.k0 = 0.75*kbar;

% Compile inputs
kt = [ones(T,1)*0.75*kbar; 0];
ct = [ones(T,1)*0.8*cbar; 0];
%kt = ones(T,1)*0.75*kbar;
%ct = ones(T,1)*0.8*cbar;
x = [kt; ct];
jacob = eye(2*T); %Jacobian guess as identity matrix
%jacob
%a = ncgm_seq(x, params)
%%

a= ncgm_seq(x, params)
%fsolve(@ncgm_seq, x)
%%
ncgm_broyden(x, jacob, 1e-4, 50, params)

%%

0.3^(params.sigma*(-1)) - params.beta*(params.alpha*0^(params.alpha - 1) + 1 - params.delta)*0^(params.sigma*(-1))

