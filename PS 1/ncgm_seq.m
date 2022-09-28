function F = ncgm_seq(x, par)
%F = zeros(1, 2);
%F(1) = exp(-exp(-(x(1)+x(2)))) - x(2)*(1+x(1)^2);
%F(2) = x(1)*cos(x(2)) + x(2)*sin(x(1)) - 0.5;
%Input of one vector with 2T+2 length

%par.alpha = 0.4;
%par.beta = 0.99;
%par.sigma = 1.000001;
%par.delta = 1;
%par.kterm = 0;
%par.cterm = 0;

T = (length(x)-2)/2;
k = [x(1:T+1,1); par.kterm];
c = [x(T+2:end,1); par.cterm];

F = zeros(2*T, 1);

for t = 1:T
    F(t) = c(t) + k(t+1) -k(t)^par.alpha - (1-par.delta)*k(t);
    F(T+t) = c(t)^(par.sigma*(-1)) - par.beta*(par.alpha*k(t+1)^(par.alpha - 1) + ...
        1 - par.delta)*c(t+1)^(par.sigma*(-1));

    %funcit = [2*(t-1)+1, 2*(t-1)+2]'; %index values of equations
    %if t == T+1
    %    F(2*T+1) = c(t) + par.kterm - k(t)^par.alpha - (1-par.delta)*k(t);
    %    break
    %end 
    %F(funcit(1)) = c(t) + k(t+1) -k(t)^par.alpha - (1-par.delta)*k(t);
    %F(funcit(2)) = c(t)^(par.sigma*(-1)) - par.beta*(par.alpha*k(t+1)^(par.alpha - 1) + 1 - par.delta)*c(t+1)^(par.sigma*(-1));
%end
end