function error = tank_error(x, z, params)

% Input dictionary
% x(1:T) output (y)
% x(T+1:2T) inflation (pi)
% x(2T+1:3T) consumption of smoothers (cS)
% x(3T+1:end) consumption of Hand-to-mouth (cH)

% Shock dictionary
% z(1:T) productivity shock a
% z(T+1,end) monetary shock nu

% Preallocate error vector for speed
error = zeros(length(x),1);

T= params.T;

kappa = (params.sigma + params.vartheta)*(1-params.theta)*(1-params.beta*params.theta)/params.beta;

% calculate errors in the following equation that define equilibrium in
% TANK model :
% y(t)-lambda*cH(t) - (1-lambda)*cS(t) =0
% cS(t+1)-cS(t)-(1/sigma)*(phi_pi*pi(t)+nu(t)-pi(t+1)) =0
% y(t+1)-y(t)-(1/sigma)*(phi_pi*pi(t)+nu(t)-pi(t+1))=0
% pi(t)-beta*pi(t+1)-kappa*(y(t)-a(t))=0
for t=1:T
    % y(t)-lambda*cH(t) - (1-lambda)*cS(t) equation error
    error(t,1)=  x(t)-params.lambda*x(3*T+t)-(1-params.lambda)*x(2*T+t);
    % cS(t+1)-cS(t)-(1/sigma)*(phi_pi*pi(t)+nu(t)-pi(t+1)) equation error
    error(T+t,1) = x(2*T+t+1)-x(2*T+t)-(1/params.sigma)*(params.phipi*x(T+t)+z(T+t)-x(T+t+1));
    % y(t+1)-y(t)-(1/sigma)*(phi_pi*pi(t)+nu(t)-pi(t+1)) equation error
    error(2*T+t,1) = x(t+1)-x(t)-(1/params.sigma)*(params.phipi*x(T+t)+z(T+t)-x(T+t+1));
    % pi(t)-beta*pi(t+1)-kappa*(y(t)-a(t)) equation error
    error(3*T+t,1) = x(T+t)-params.beta*x(T+t+1)-kappa+(x(t)-z(t));
end