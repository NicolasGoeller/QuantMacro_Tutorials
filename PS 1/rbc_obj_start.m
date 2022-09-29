    function error=rbc_obj_start(x,par)%alpha,beta,sigma,delta,kbar,cbar)

%x(1:T) is capital
%x(T+1:2T) is consumption
alpha = par.alpha;
beta = par.beta;
sigma = par.sigma;
delta = par.delta;
cbar = par.cbar;
k0 = par.k0;

%par.alpha = 0.4;
%par.beta = 0.99;
%par.sigma = 1.000001;
%par.delta = 1;
%par.cbar = params.cbar;
%par.k0 = params.k0;

% Euler equation
% -sigma*chat_t=(-sigma)*chat_t+1
% +beta*margprod*khat_t+1+beta*margprod*zbar*zhat_t+1

% res constraint
% khat_t+1 - ckrat*chat_t - (margprod+(1-delta)) khat_t -kbar^alpha*zbar*zhat_t=0

T=length(x)/2;

%% are these really k0 and cbar???? or steady state values
error(1,1)=x(1,1)-k0; 
error(2*T,1)=x(2*T,1)-cbar;

% calculate errors in budget constraint - capital k
for t=2:T
   error(t,1)= x(t,1) - x(t-1,1)*(1-delta) - x(t-1,1)^alpha + x(T+t-1,1);
end

% calculate errors in EE - consumption c
for t=1:T-1
   error(T+t,1) = x(T+t,1)^(-sigma) - beta*(alpha*x(t+1,1)^(alpha-1)+1-delta)*x(T+1+t,1)^(-sigma);
   %error(2*T-t,1) = x(2*T-t,1)^(sigma*(-1)) - beta*(alpha*x(T-t,1)^(alpha-1) + 1 - delta)*x(2*T+1-t,1)^(sigma*(-1));
end

end