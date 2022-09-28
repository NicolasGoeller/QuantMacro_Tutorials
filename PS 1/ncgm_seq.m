function F = ncgm_seq(x, par)

T = length(x);
k = x(1:T,1);
c = x(1:end,2);

F = zeros(2*T+1,1);

for t = 1:T
    %t
    funcit = [2*(t-1)+1, 2*(t-1)+2];
    if t == T
        F(2*T+1) = c(t) + par.kterm - k(t)^par.alpha - (1-par.delta)*k(t);
        break
    end 
    F(funcit(1)) = c(t) + k(t+1) -x(t)^par.alpha - (1-par.delta)*k(t);
    F(funcit(2)) = c(t)^(par.sigma*(-1)) - par.beta*(par.alpha*k(t+1)^(par.alpha - 1) + ...
        1 - par.delta)*c(t+1)^(par.sigma*(-1));
end
end