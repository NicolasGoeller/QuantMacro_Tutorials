function out = ncgm_broyden(guess, jacob_guess, convergence_crit, maxiter, ncgm_par)
%function for executing broydens method of rootfinding
% jacob_guess: initial value of Jacobian matrix
% convergence_crit: criterion for accepting the convergence of the equations
% maxiter: maximim number of iterations before breaking
error = ncgm_seq(guess, ncgm_par);
B = jacob_guess;
iter = 0;

%check if ssum squared of dev from 0 in equation system
while ncgm_error(error) > convergence_crit
    iter = iter + 1;
    %Compute difference of new guess
    s = -inv(B)*ncgm_seq(guess, ncgm_par); 
    size(s)
    size(guess)
    
    %compute new guess 
    guess_new = guess + s; 
    %Compute function value difference for two guesses
    y = ncgm_seq(guess_new, ncgm_par) - ncgm_seq(guess, ncgm_par);
    %y = ncgm_seq(guess_new) - ncgm_seq(guess);
   
    guess = guess_new;

    error = ncgm_seq(guess, ncgm_par); %evaluate error of new guess

    if abs(s'*s) > 1e-2 % evaluate if norm of x difference is not 0
        B = B + ((y' - B*s)*s')/(s'*s); %Compute updated Jacobian
    end
    if iter > maxiter % if iteration reaches very large value terminate loop 
        break
    end

    out = guess;
end

