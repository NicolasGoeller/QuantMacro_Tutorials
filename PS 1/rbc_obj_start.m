    function error=rbc_obj(x,alpha,beta,sigma,delta,kbar,cbar,z_trans)

%x(1:T) is capital
%x(T+1:2T) is consumption

% Euler equation
% -sigma*chat_t=(-sigma)*chat_t+1
% +beta*margprod*khat_t+1+beta*margprod*zbar*zhat_t+1

% res constraint
% khat_t+1 - ckrat*chat_t - (margprod+(1-delta)) khat_t -kbar^alpha*zbar*zhat_t=0

T=length(x);

%% are these really k0 and cbar???? or steady state values
error(1,1)=x(1,1)-kbar; %% this was initally k0
error(2*T,1)=x(2*T,1)-cbar;

% calculate errors in budget constraint
for t=2:T
   error(t,1)= x(t,1) - kbar;
end

% calculate errors in EE
for t=1:T-1
   error(T+t,1)= x(T+t,1) - cbar;
end

end