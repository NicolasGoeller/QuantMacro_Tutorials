%%Multiple shooting code
% first need of defining steps foir the try of k1

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

criter_V = 1e-6; % conv criterion for value function
T=10; % periods for transition

kbar=((1/beta-1+delta)/(alpha))^(1/(alpha-1)); %% here we inserted a *beta
cbar= kbar^alpha-delta*kbar
% initial level of capital in the transition
k0=kbar*0.75; % you may have to change this from the problem set instructions
c0=cbar*0.8

%pas=kbar/10; %pas relativement petit
%we solve kT for k0 and c0 known recursively
%k(t+1)=-c(t)+k(t)^alphe+(1-delta)k(t)
%c(t+1)=c(t)*(beta*(1-delta+alpha*k(t+1)^(alpha-1)))^(-sigma)

%point=1000;

c = [c0 ; zeros(T,1)];
k = [k0 ; zeros(T,1)];
i =1;


%c0_shots = linspace(cbar*0.25,1.1*cbar,point)';
%%kT_results = zeros(1,1);
%%kT_results_error = zeros(1,1);

%c0Low=0.8*cbar;
%c0Up=1.2*cbar;
c0=0.5*cbar;
% c02=0.51*cbar;
% 
% %initialization pour first shot
% c = [c0 ; zeros(T,1)];
% %cUp = [c0Up; zeros(T,1)];
% k = [k0 ; zeros(T,1)];
% 
%  for t=1:T
%      k(t+1)=-c(t)+k(t)^alpha+(1-delta)*k(t);
%      c(t+1)=c(t)*(beta*(1-delta+alpha*k(t+1)^(alpha-1)))^(-sigma);
%  end
%  kT_results(1)=k(T+1); 
%  kT_results_error(1)=kbar-k(T+1); %first result
% 
%  %initialisation second shot
% 
% c = [c02 ; zeros(T,1)];
% %cUp = [c0Up; zeros(T,1)];
% k = [k0 ; zeros(T,1)];
% 
%  for t=1:T
%      k(t+1)=-c(t)+k(t)^alpha+(1-delta)*k(t);
%      c(t+1)=c(t)*(beta*(1-delta+alpha*k(t+1)^(alpha-1)))^(-sigma);
%  end
%  kT_results =[kT_results; k(T+1)]; 
%  kT_results_error =[kT_results_error; kbar-k(T+1)]; %second result
% 
% %%
% 
% i=3;
epsilon = 1e-5;
i=1;

%%

kT_results = [];
k1T_results = [];
kT_results_error = []; %kT_results_error(i) give the error value for the ith shot
k1T_results_error = [];

err = abs(k(T+1)-kbar);

%%
while err>criter_V

    
    c = [c0 ; zeros(T,1)];
    c1 = [c0*(1+epsilon); zeros(T,1)]; %built to estimate the derivative
    %cUp = [c0Up; zeros(T,1)];
    k = [k0 ; zeros(T,1)];
    k1 = [k0 ; zeros(T,1)];

    for t=1:T
        k(t+1)=-c(t)+k(t)^alpha+(1-delta)*k(t);
        c(t+1)=c(t)*(beta*(1-delta+alpha*k(t+1)^(alpha-1)))^(-sigma);
        k1(t+1)=-c1(t)+k1(t)^alpha+(1-delta)*k1(t);
        c1(t+1)=c1(t)*(beta*(1-delta+alpha*k1(t+1)^(alpha-1)))^(-sigma);
    end

    kT_results = [kT_results ; k(T+1)];
    k1T_results = [k1T_results ; k1(T+1)];
    kT_results_error = [kT_results_error ;k(T+1)-kbar]; %kT_results_error(i) give the error value for the ith shot
    k1T_results_error = [k1T_results_error ;k1(T+1)-kbar];
    deriv_estimate = (k1T_results_error(i)-kT_results_error(i))/(c0*epsilon);
    c0new = c0 - kT_results_error(i)/deriv_estimate;
    c0=c0new;
    err=abs(k(T+1)-kbar);
    %c0new = c0-(k(T+1)-kbar)*(c0-c02)/(kT_results_error(i-2)-kT_results_error(i-1));
    i=i+1;
end


