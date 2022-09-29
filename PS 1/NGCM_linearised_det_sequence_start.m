% ==========================
% Problem set solving the stochastic neoclassical growth model different
% ways
% ==========================
% clear the workspace
clear
% close all figures
close all

% ============
% parameters  - you may have to change this according to instructions
% ============
alpha=0.4; % capital share was named theat before
beta = 0.99; % disfcount factor
%rho = 0.95;   % persistence of TFP shock this doesnt eist in the PS
sigma = 1.0001; % CRRA coefficient (for 1 equals log, but need to replace the function, so set close to 1)
delta=1;

% ============
% options and convergence criteria
% ============

criter_V = 1e-7; % conv criterion for value function
T=10; % periods for transition

%mean of capital non-stochastic steady state
kbar=((1/beta-1+delta)/(alpha))^(1/(alpha-1)); %% here we inserted a *beta
cbar= kbar^alpha-delta*kbar
% initial level of capital in the transition
k_0=kbar*0.75; % you may have to change this from the problem set instructions


% ==============
% Problem 2 analytical case delta=1, finite T and infinite T
% ==============

% a&b. analytical policies, finite and infinite horizon

alphabeta = zeros(1,T); %creating alpahbeta matrix with alphabeta(i)=(alpha*beta)^i
alphabeta(1)= 1;
for t=2:T
    alphabeta(t)=alphabeta(t-1)*alpha*beta;
end

if delta==1 % only makes sense for delta=1
    k_analyt=zeros(1,T);
    consum_analyt=zeros(1,T); %vectors with consumption at t
    consum_analyt_finite=zeros(1,T);
    
    k_analyt(1)=k_0;
    k_analyt_finite(1)=k_0;
    for t=2:T
        % use the recursive expression of k calculated in Part1.3
        k_analyt(t)= (1-1/(sum(alphabeta(1:T-t+2))))*k_analyt(t-1)^alpha ; %sum to T-t+2 because we begin at 1 and not 0
        k_analyt_finite(t)= (1-(1/sum(alphabeta(1:T-t+2))))*k_analyt_finite(t-1)^alpha;
        % k_analyt_finite(t)=[xxxx fill this in xxxx];
        consum_analyt(t-1)=k_analyt(t-1)^alpha - k_analyt(t);
        consum_analyt_finite(t-1)=k_analyt_finite(t-1)^alpha - k_analyt_finite(t);
    end
    consum_analyt_finite(T)=k_analyt_finite(T);
    consum_analyt(T)=k_analyt(T);
end

% c. numerical solution algorithms to above problems

% param values
params.alpha = 0.4;
params.beta = 0.99;
params.sigma = 1.000001;
params.delta = 1;
params.kterm = 0;
params.cterm = 0;

%set guesses
kbar=((1/params.beta-1+params.delta)/(params.alpha*params.beta))^(1/(params.alpha-1));
cbar = kbar^params.alpha-params.delta*kbar;
params.k0 = 0.75*kbar;

% Compile inputs
%kt = [ones(T,1)*0.75*kbar; 0];
%ct = [ones(T,1)*0.8*cbar; 0];

%% Broydens Method

% Broydens method Finite time T=10
T1 = 10;
kt1 = ones(T1,1)*0.75*kbar;
ct1 = ones(T1,1)*0.8*cbar;
x1 = [kt1; ct1];
jacob1 = ncgm_jacob(x1, params);
trans1 = ncgm_broyden(x1, jacob1, 1e-6, T1, params);

% Broydens method Finite time T=100
T2 = 100;
kt2 = ones(T2,1)*0.75*kbar;
ct2 = ones(T2,1)*0.8*cbar;
x2 = [kt2; ct2];
jacob2 = ncgm_jacob(x2, params);
trans2 = ncgm_broyden(x2, jacob2, 1e-6, T2, params);

% Broydens method Finite time T=200
T3 = 200;
kt3 = ones(T3,1)*0.75*kbar;
ct3 = ones(T3,1)*0.8*cbar;
x3 = [kt3; ct3];
jacob3 = ncgm_jacob(x3, params);
trans3 = ncgm_broyden(x3, jacob3, 1e-6, T3, params);

% Broydens method Infinite time - set T=?
T4 = 15;
kt4 = ones(T4,1)*0.75*kbar;
ct4 = ones(T4,1)*0.8*cbar;
x4 = [kt4; ct4];
jacob4 = ncgm_jacob(x4, params);
trans4 = ncgm_broyden(x4, jacob4, 1e-6, T4, params);

x = 1:1:50; % Time variable
plot(x, , x, trans4, '.-'), legend('Analytical solution', 'Broydens method', 'Multiple Shooting');

%% Multiple Shooting


%%
% =====================
% 3. Log-linearization
% =====================

% some ratios as function of parameters
ybar=kbar^alpha;
cbar=ybar-delta*kbar;
% need for matrix form of loglin
ckrat=cbar/kbar;

% a. write system as A E[y_t+1]+B y_t=0



% order c,k,z
A=[-sigma, beta *(alpha-1)*alpha*kbar^(alpha-1) ; 0 , 1 ];

B= [-sigma , 0 ; -ckrat,(1/beta)];

D = inv(A)*B;

% note that these are right-hand eigenvectors
[ ev lambda]=eig(D); %%%give the egein vectors and eigen values of D
aaa=inv(ev); %%% give the invert of the eigen vector matrix

% find eigenvalues equal or larger than one, and check that they equal the
% number of jump variables - in that case, set BKcond to 1

BKcond = 1;
eigen = [lambda(1,1) lambda(2,2)];
if all(eigen > 1) || all(eigen < 1)
    BKcond = 0;
end

%If the Blanchard Kahn condition is satisfied, we can find the expression
%of "alpha1"
if BKcond~=1
    disp('BK conditions not satisfied')
else
    bkev =find(abs(diag(lambda))>1);
    invP=aaa(bkev,:);%%Select the element of the invert of the vector matrix needed to compute the policy function
    polfunc= -invP(1,2)/invP(1);% you need to find the policy for consumption here : derived analytically 
end

% policy functions are in log deviations. so need to pre-multiply by cbar
% and add cbar to get levels
%%

% calculate the deterministic transition  using the linearised policy
% functions, and law of motion
%Here is generating a sequence of k from the policy function
k_lin(1)=k_0;
for t=2:T
    c_lin(t-1)=polfunc*((k_lin(t-1)-kbar)/kbar)*cbar + cbar; %% compute the linear c from the formula derived in the latex part 3.
    k_lin(t)=k_lin(t-1)^alpha +(1-delta)*k_lin(t-1) - c_lin(t-1);%% same thing as above
end

%Here is generating a sequence of k lagged by 1 unit of time 
k_lina(1)=k_lin(1,2);
for t=2:(T)
    c_lina(t-1)=polfunc*((k_lina(t-1)-kbar)/kbar)*cbar + cbar;
    k_lina(t)=k_lina(t-1)^alpha +(1-delta)*k_lina(t-1) - c_lina(t-1);
end


%graphpart 4 q1
T= [10,100,200];
k_analyt_finite = zeros(3,T(3));
c_analyt_finite =zeros (3,T(3));
for i=1:length(T)
    k_analyt_finite(i,1)=k_0;
    for t=2:T(i)+1
        k_analyt_finite(i,t) = alpha*beta*((1-(alpha*beta)^(T(i)+1-t))/(1-(alpha*beta)^(T(i)+1-t+1)))*k_analyt_finite(i,t-1)^alpha;
        c_analyt_finite(i,1) = k_0^alpha +(1-delta)*k_analyt_finite(i,1) - k_analyt_finite(i,2);
        c_analyt_finite(i,t) = ((1-alpha*beta)/(1-(alpha*beta)^(T(i)-t+1)))*k_analyt_finite(i,t)^alpha;
    end
end



T=15;
k_analyt = zeros(1,T);
c_analyt = zeros(1,T);
k_analyt(1) = k_0;



for t=2:T+1
    k_analyt(t) = alpha*beta*k_analyt(t-1)^alpha;
    c_analyt(t) = (1-alpha*beta)*k_analyt(t)^alpha;
end


hold on 
plot(1:16, k_analyt);
plot(1:16, c_analyt);
xlabel('Time','FontSize',10);
ylabel('Value','FontSize',10);
legend('Capital','Consumption');
hold off




% ==============
% 4. Solve deterministic sequence
% ==============

% as usual we need an initial guess for the sequences of k and c - here we
% just use steady state
x0=[kbar*ones(T,1);cbar*ones(T,1)];
%%
% zpath=ones(T,1);
% zpath(1)=1.01;
% for i=2:T-1
%     zpath(i)=exp(log(zpath(i-1)));
% end

% now we need a function that returns errors of the equation system,
% constraining k(0)=k_0, and c(T)=cbar

% here we use Broyden's method

%Initial guess for jacobian - use finite difference at steady state
%sequences of k and c
%rbc_obj(x,alpha,beta,sigma,delta,kbar,cbar)
%%
clear dx J
for i=1:2*T
    dx = zeros(2*T,1);
    dx(i)=x0(i)*0.001;
    diffF = rbc_obj_start(x0+dx,alpha,beta,sigma,delta,kbar,cbar)-rbc_obj_start(x0,alpha,beta,sigma,delta,kbar,cbar);
    diffX= dx(i);
    J(:,i)=[diffF/diffX] ;
end
clear x
%%
crit=1e-4;
x=x0;
f0=rbc_obj_start(x,alpha,beta,sigma,delta,k_0,cbar);
f=rbc_obj_start(x,alpha,beta,sigma,delta,k_0,cbar);
%%
while max(abs(f))>crit

%dx = [xxxx fill this in xxxx] 
%     dx = zeros(2*T,1);
%     for i=1:2*T
%         for j=1:2*T
%             dx(i)= dx(i)-J(i,j)^(-1)*f(i);
%         end
%     end
    f0 = rbc_obj_start(x,alpha,beta,sigma,delta,k_0,cbar);
    dx = -J\f; 
    %dx = dx*(x'*x)/(10*(dx'*dx));
    %while x+dx<1e-6
      %  dx=dx/2;
    %end
    xn1=max(ones(2*T,1)*1e-8,x+dx);
    dx=xn1-x;
    f = rbc_obj_start(xn1,alpha,beta,sigma,delta,k_0,cbar);
    %J = [] ;
    J = J + ((f-f0)-J*dx)*dx'/(dx'*dx);
    %B = B + ((y' - B*s)*s')/(s'*s);
    x=xn1;


end


%%
% ==============
% 4. Figures
% ==============
% plot policy function
figure(1)
% levels
subplot(2,1,1)
title('Policy functions')
hold on
plot(k_lin,k_lina) %% here is a graph of k'(k)
xlabel('k_{t}','FontSize',10);
ylabel('k_{t+1}','FontSize',10);
hold off

%%
% plot the transition
x = 1:1:50;
figure(2)
title('Simulated transition - deterministic')
hold on
% Time variable
plot(x, k_analyt, x, trans4, x, ), legend('Analytical solution', 'Broydens method', 'Multiple Shooting'),
xlabel('Time steps'), ylabel('Capital level');
hold on
plot(x, consum_analyt, x, trans4, x, ), legend('Analytical solution', 'Broydens method', 'Multiple Shooting'),
xlabel('Time steps'), ylabel('Consumption level');

if delta==1 && abs(sigma-1)<0.001
    h = legend([xxxx fill this in xxxx] ,'Location', 'best','Orientation','Vertical');
else
    h = legend([xxxx fill this in xxxx] ,'Location', 'best','Orientation','Vertical');
end
set(h,'fontsize',12,'Interpreter','Latex');%'Orientation', 'horizontal'




