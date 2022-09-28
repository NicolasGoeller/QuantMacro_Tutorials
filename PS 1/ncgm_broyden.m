function out = ncgm_broyden(guess, jacob_guess, convergence_crit, ncgm_par)

T = (length(guess)-1)/2;
%guess = [x(1:T,1); x(T+1:end,2)];
error = ncgm_seq(guess, ncgm_par);
B = jacob_guess;

while ncgm_error(error) > convergence_crit
    s = -inv(B)*ncgm_seq(guess, ncgm_par);  
    guess_new = guess + s;
    y = ncgm_seq(guess_new, ncgm_par) - ncgm_seq(guess, ncgm_par);
    guess_new[]
    guess = guess_new;

    error = ncgm_seq(guess, ncgm_par);

    if abs(s'*s) > 1e-2
        B = B + ((y-B*s)*s')/(s'*s);
    end
    out = guess;
end

