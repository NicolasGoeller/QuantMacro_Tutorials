    function error=rbc_obj_start(x,par)%alpha,beta,sigma,delta,kbar,cbar)

%x(1:T) is capital
%x(T+1:2T) is consumption

%Set paramater values
alpha = par.alpha;
beta = par.beta;
sigma = par.sigma;
delta = par.delta;
cbar = par.cbar;
k0 = par.k0;

T=length(x)/2;

% Set error values for known terms
error(1,1)=x(1,1)-k0; 
error(2*T,1)=x(2*T,1)-cbar;

% calculate errors in budget constraint - capital k
for t=2:T
   error(t,1)= x(t,1) - x(t-1,1)*(1-delta) - x(t-1,1)^alpha + x(T+t-1,1);
end

% calculate errors in EE - consumption c
for t=1:T-1
   error(T+t,1) = x(T+t,1)^(-sigma) - beta*(alpha*x(t+1,1)^(alpha-1)+1-delta)*x(T+1+t,1)^(-sigma);
end

end