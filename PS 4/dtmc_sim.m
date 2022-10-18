function dtmc = dtmc_sim(trans, start, T, N_sim)
% Code to simulate a discrete-time, discrete state markov chain
%Purpose:    Simulates the evolution of discrete markov chains according to
%            a given starting value
%
%Format:     dtmc = dtmc_sim(trans, start, T, N_sim)
%
%Input:      trans      scalar, transition matrix
%            start      scalar, starting value of all simulations 
%            T          scalar, length of markov chain series
%            N_sim      scalar, number of simulations
%
%Output:     dtmc       T*N_sim matrix, of markov chain value series

% Pre-allocate for speed
dtmc = zeros(T,N_sim);

% Run loop over each simulation, then over all time periods
for j=1:N_sim
    dtmc(j,1)=start; % Set start values for all simulations
    for t=2:T
        % Get out conditional transition probabilities
        thresholds = cumsum(trans(dtm(j,t-1),:)); 
        % get a random value
        r = rand();
        % Find position of value that exceeds random and set as new value
        dtmc(j,t) = find(thresholds> r,1);
    end
end
