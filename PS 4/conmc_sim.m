function conmc = conmc_sim(mu, sigma, rho, start, T, N_sim)
% Code to simulate a discrete-time, continuous state markov chain
%Purpose:    Simulates the evolution of continuous markov chains according to
%            a given starting value
%
%Format:     dismc = dismc_sim(mu, sigma, rho, start, T, N_sim)
%
%Input:      mu         scalar, mean of disturbance distribution
%            sigma      scalar, std. dev of disturbance distribution
%            rho        scalar, persistence parameter
%            start      scalar, starting value of all simulations 
%            T          scalar, length of markov chain series
%            N_sim      scalar, number of simulations
%
%Output:     conmc       N_sim*T matrix, of markov chain value series

% Pre-allocate for speed
conmc = zeros(N_sim,T);

% Run loop over each simulation, then over all time periods
for j=1:N_sim
    conmc(j,1)=start; % Set start values for all simulations
    for t=2:T
        % Get a random value for disturbance
        epsilon = normrand(mu, sigma);
        % compute new shock value
        conmc(j,t) = rho*conmc(j,t-1)+ epsilon;
    end
end
