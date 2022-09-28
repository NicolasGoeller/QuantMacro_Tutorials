function out = ncgm_broyden(guess, jacob_guess, convergence_crit, ncgm_par)

error = ncgm_seq(guess, ncgm_par);
B = jacob_guess;

while ncgm_error(error) > convergence_crit
    s = -inv(B)*ncgm_seq(guess, ncgm_par);
    size(guess)
    size(s)
    guess_new = guess + s;
    y = ncgm_seq(guess_new, ncgm_par) - ncgm_seq(guess, ncgm_par);
    guess = guess_new;

    error = ncgm_seq(guess, ncgm_par);

    if abs(s'*s) > 1e-2
        B = B + ((y-B*s)*s')/(s'*s);
    end
    out = guess;
end

