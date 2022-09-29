function out = ncgm_broyden(guess, jacob_guess, convergence_crit, maxiter, ncgm_par)
%function for executing broydens method of rootfinding
% guess: sequence of values to use as initial guess
% jacob_guess: initial value of Jacobian matrix
% convergence_crit: criterion for accepting the convergence of the equations
% maxiter: maximim number of iterations before breaking

error = rbc_obj_start(guess, ncgm_par);

B = jacob_guess;
iter = 0;

%check if ssum squared of dev from 0 in equation system
while ncgm_error(error) > convergence_crit
    iter = iter + 1;
    iter
    %Compute difference of new guess
    s = -inv(B)*rbc_obj_start(guess, ncgm_par); 
    %s
    %compute new guess 
    guess_new = guess + s; 
    %Impose non-negativity constraint
    guess_new = max(guess_new, 1e-8);
    %guess_new 

    %Compute function value difference for two guesses
    y = rbc_obj_start(guess_new, ncgm_par) - rbc_obj_start(guess, ncgm_par);
    %y = ncgm_seq(guess_new) - ncgm_seq(guess);
    %y
    guess = guess_new;
    
    error = rbc_obj_start(guess, ncgm_par); %evaluate error of new guess

    if abs(s'*s) > 1e-2 % evaluate if norm of x difference is not 0
      
        B = B + ((y - B*s)*s')/(s'*s); %Compute updated Jacobian
        %B
    end
    if iter >= maxiter % if iteration reaches very large value terminate loop 
        break
    end

    out = guess;
end

