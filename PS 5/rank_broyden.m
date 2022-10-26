function out = rank_broyden(guess, shocks, jacob_guess, convergence_crit, maxiter, params)
%function for executing broydens method of rootfinding
% guess: sequence of endogeneous variable values to use as initial guess
% shocks: Sequence of exogenous variable values
% jacob_guess: initial value of Jacobian matrix
% convergence_crit: criterion for accepting the convergence of the equations
% maxiter: maximim number of iterations before breaking
% params: parameter values to be fed into the function

error = rank_error(guess, shocks, params);

B = jacob_guess;
iter = 0;

%check if sum squared of dev from 0 in equation system
while max(abs(error)) > convergence_crit
    iter = iter + 1;

    %Compute difference of new guess
    s = -inv(B)*rank_error(guess, shocks, params);

    %compute new guess 
    guess_new = guess + s; 
    %Impose non-negativity constraint
    guess_new = max(guess_new, 1e-8);

    %Compute function value difference for two guesses
    error_new = rank_error(guess, shocks, params);
    y = error_new - error;
    
    % Recycle values for next round
    guess = guess_new;
    error = error_new;

    if abs(s'*s) > 1e-2 % evaluate if norm of x difference is not 0
        B = B + ((y - B*s)*s')/(s'*s); %Compute updated Jacobian
    end
    if iter >= maxiter % if iteration reaches very large value terminate loop 
        break
    end

    out = guess;
end

