function error = rank_error(x, z, params)

% Input dictionary
% x(1:T) output
% x(T+1:2T) inflation
% x(2T+1:end) interest rate

% Shock dictionary
% z(1:T) productivity shock a
% z(T+1,end) monetary shock nu

% Preallocate error vector for speed
error = zeros(length(x),1);

kappa = (params.sigma + params.vartheta)*(1-params.theta)*(1-params.beta*params.theta)/params.beta;

% Here we assume bc of RE that E(x_{t+1}) = x_{+1}
for t=1:T
    % DIS equation erros
    error(t,1)=  x(t,1) + params.sigma^(-1)*(x(2*T+t,1) - x(T+t+1,1)) - x(t+1,1);
    % NKPC equation error
    error(T+t,1) = x(T+t,1) - params.beta*x(T+t+1,1) - kappa*(x(t,1) - params.varphi*z(t,1));
    % taylor rule error
    error(2*T+t,1) = x(2*T+t,1) - params.phipi*x(T+t,1) + z(T+t,1);
end
