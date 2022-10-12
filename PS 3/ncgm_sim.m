function [kpath, cpath, zpath] = ncgm_sim(T, M, N, n_sim, par)
%Function ncgm_sim
%
%Purpose:    Run a discretized, stochastic NCGM for the indicated
%            parameters, grid size, shock structure and numerb of
%            simulations.
%
%Format:     {kpath, cpath, zpath} = ncgm_sim(T, M, N, Nsim, par)
%
%Input:      T        scalar, number of time steps
%            M        scalar, grid size of capital k
%            N        scalar, grid size of shock z
%            n_sim    scalar, number of simulations to run
%            par      structure, numerical parameters for model
%
%Output:     kpath       T+1*n_sim matrix, optimal paths for kprime
%            cpath       T*n_sim matrix, optimal paths for consumption 
%            zpath       T*n_sim matrix, derived stochastic paths of shocks

% Preallocate output formats
kpath = zeros(T+1,n_sim);
cpath = zeros(T,n_sim);
zpath = zeros(T,n_sim);
%kgrid = zeros(1,M);

% Generate kgrid
kbar=((1/par.beta-1+par.delta)/(par.alpha))^(1/(par.alpha-1));
% the grid
if par.delta==1
    if par.linear==1
        % linear sequence kbar -2kbar in N steps
        kgrid=linspace(kbar*0.85,1.15*kbar,M);
    else
        % linear sequence 0-0.5 in N steps, divided by 0.5 all to the power
        % of 5; times 1.5kbar
        temp=linspace(0,0.5,M).^5/0.5^5*(0.85*kbar-1.15*kbar);
        % 0.5kbar + temp
        kgrid=0.85*kbar+temp;
    end
else
    % if delta ~= 1, linear sequence 0.25kbar-2kbar in N steps
    kgrid=linspace(0.85*kbar ,1.15*kbar,M);
end

% Generate shock matrix
rng(4);
[Z_tauchen, P_tauchen] = tauchen(N,0,par.rho,par.sigma,2);
p = dtmc(P_tauchen);
X0 = [0 0 N_sim 0 0]; %simulate N_sim starting at initial value of Z=0 
X = simulate(p,T,"X0",X0); %X0 define the number of simulation to be made starting with a given initial condition
for i=1:N 
      zpath(X==i)=Z_tauchen(i,1);
end

% Generate kprime matrix
[V, kprime, index] = discrete_search();

kpath(1,:) = par.k0;

for i = 1:n_sim
    
    for t = 1:T
        % Get index of k in kgrid
        index = find(kgrid == k);
        % Get kprime optimal for cuurent z and k
        kpath(t+1,i) = kprime(index,zpath(t,i));
        % Compute optimal c for k, kprime and z
        cpath(t,i) = exp(zpath(t,i))*kpath(t,i)^par.alpha - kpath(t,i) + (1-par.delta)*kpath(t,i);
        k = kpath(t,i);
    end
end


