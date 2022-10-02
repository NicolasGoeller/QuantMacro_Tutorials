% ==============
% 4. Solve deterministic sequence
% ==============

% as usual we need an initial guess for the sequences of k and c - here we
% just use steady state
x0=[kbar*ones(T,1);cbar*ones(T+1,1)];

finite = @(x)rbc_obj_finite(x,theta,beta,gamma,delta,k_0,T);
%Uncomment to the infinite case
%infinite = @(x) rbc_obj_infinite(x,theta,beta,gamma,delta,k_0,kbar,T);

%Start easy
tic
[x_sol,fsol]=fsolve(@(x)rbc_obj_finite(x,theta,beta,gamma,delta,k_0,T),x0);
toc

%Then Broyden's method

J = zeros(2*T+1);
h = 0.01;

for i=1:2*T+1
    dx = zeros(2*T+1,1);
    dx(i)=h;
    J(:,i)=(finite(x0+dx)-finite(x0))./h ;
end

crit=1e-10;
x=x0;
f = finite(x);

while max(abs(f))>crit
    dx = -J\f;
    x=x+dx;
    f = finite(x);
    J = J +f*dx'/(dx'*dx) ;
    max(abs(f))
end

k_trans_br=x(1:T,1);
c_trans_br=x(T+1:2*T+1,1);