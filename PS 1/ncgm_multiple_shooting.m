function out = ngcm_multiple_shooting(alpha, beta, delta, sigma, T, criter_V)
%function to comput c0 by imposing convergence to kT=kbar

kbar=((1/beta-1+delta)/(alpha))^(1/(alpha-1)); %compute the steady state kbar
cbar= kbar^alpha-delta*kbar; %compute the steady state cbar

% initial level of capital in the transition
k0=kbar*0.75; % you may have to change this from the problem set instructions
c0=cbar*0.5; %arbitrary input the first shot for c0

c = [c0 ; zeros(T,1)]; %vector storing the value of ct for initial value of c0
k = [k0 ; zeros(T,1)]; %vector storing the value of kt for initial value of c0
epsilon = 1e-5; %epsilon use to estimate the derivative of f(c0^j)=k(T+1)^j-kbar
%k(T+1)^j = f(c0^j)+kbar
i=1;

kT_results = []; %kT estimate for initial value of c0=c0^j : kT_results(j)=kT^j=k(T+1)^j = f(c0^j)+kbar
k1T_results = []; %estimate useful to approximate the derivative
kT_results_error = []; %kT_results_error(i) give the error value for the ith shot (this is the objective function)
k1T_results_error = []; %useful to approx the derivative

err = abs(k(T+1)-kbar);

%%
while err>criter_V %convergence criteria

    
    c = [c0 ; zeros(T,1)];
    c1 = [c0*(1+epsilon); zeros(T,1)]; %built to estimate the derivative
    k = [k0 ; zeros(T,1)];
    k1 = [k0 ; zeros(T,1)]; %built to estimate derivative

    for t=1:T %we define the function f(c0)=kT recursively
        k(t+1)=-c(t)+k(t)^alpha+(1-delta)*k(t);
        c(t+1)=c(t)*(beta*(1-delta+alpha*k(t+1)^(alpha-1)))^(-sigma);
        k1(t+1)=-c1(t)+k1(t)^alpha+(1-delta)*k1(t);
        c1(t+1)=c1(t)*(beta*(1-delta+alpha*k1(t+1)^(alpha-1)))^(-sigma);
    end

    kT_results = [kT_results ; k(T+1)];
    k1T_results = [k1T_results ; k1(T+1)];
    
    kT_results_error = [kT_results_error ;k(T+1)-kbar]; %kT_results_error(i) give the error value for the ith shot
    k1T_results_error = [k1T_results_error ;k1(T+1)-kbar];
    
    deriv_estimate = (k1T_results_error(i)-kT_results_error(i))/(c0*epsilon);%estimation fo the derivative of f at c0j
    
    c0new = c0 - kT_results_error(i)/deriv_estimate;
    c0=c0new;
    err=abs(k(T+1)-kbar);
    i=i+1;
end
out = c0;
