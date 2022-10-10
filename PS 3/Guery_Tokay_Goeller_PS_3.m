% ==========================
% Problem set solving the stochastic neoclassical growth model different
% ways
% ==========================
% clear the workspace
clear
% close all figures
close all
%addpath('C:\Users\tbroer\Dropbox\Teaching\PSE\2021 Quantitative Macro\Problem sets\PS RBC')
%addpath('C:\Users\tbroer\Dropbox\Teaching\PSE\2021 Quantitative Macro\Problem sets')

% ============
% parameters
% ============
alpha=0.4; % capital share - this is alpha
beta = 0.9; % discount factor
rho = 0.95;   % persistence of TFP shock
sigma = 2.00000001; % CRRA coefficient (for 1 equals log, but need to replace the function, so set close to 1)
delta=0.1;

% ============
% options and convergence criteria
% ============

Howard =1; % set to 1 if you want to do policy fct iteration / Howard improvement
criter_V = 1e-6; % conv criterion for value function
N=50; % number of grid points
linear=1; % grid linear or not
N_sim = 100; % nbr of simulation
T = 150;%period of transition
%mean of capital non-stochastic steady state
kbar=((1/beta-1+delta)/(alpha))^(1/(alpha-1));
k_0 = 0.75*kbar;









% ==============
% 0. Grids,  etc
% ==============

% center the grid around mean of capital non-stochastic steady state
kbar=((1/beta-1+delta)/(alpha))^(1/(alpha-1));
% the grid
if delta==1
    if linear==1
        % linear sequence kbar -2kbar in N steps
        kgrid=linspace(kbar/2,2*kbar,N);
    else
        % linear sequence 0-0.5 in N steps, divided by 0.5 all to the power
        % of 5; times 1.5kbar
        temp=linspace(0,0.5,N).^5/0.5^5*(2*kbar-kbar/2);
        % 0.5kbar + temp
        kgrid=kbar/2+temp;
    end
else
    % if delta ~= 1, linear sequence 0.25kbar-2kbar in N steps
    kgrid=linspace(kbar/4 ,2*kbar,N);
end


% Problem 1 - discretize income process and simulate
% ==============
[Z_tauchen, P_tauchen] = tauchen(5,0,0.95,0.007,2);
p = dtmc(P_tauchen);
X0 = Z_tauchen';
X = simulate(p,150,"A", Z_tauchen);
graphplot(p,'ColorEdges',true);

figure;
simplot(p,X);
A = zeros(1,lenght(Z_tauchen));
A=A'


disp('Standard devations of Z discrete, continuous')
disp([mean(std(log(Z_sim)')),(sigmaepsilon^2/(1-rho^2))^0.5])
for i=1:N_sim
    rho_emp(i)=corr(log(Z_sim(i,1:end-1))',log(Z_sim(i,2:end))');
end
disp('Autocorrelation of Z discrete, continuous')

disp([mean(rho_emp),rho])





% ==============
% 1. analytical case delta=1, finite T and infinite T
% ==============

alpha=0.4; % capital share - this is alpha
beta = 0.99; % discount factor
rho = 0.95;   % persistence of TFP shock
sigma = 2.00000001; % CRRA coefficient (for 1 equals log, but need to replace the function, so set close to 1)
delta=0.1;

% a. analytical policies, finite and infinite horizon


if delta==1 % only makes sense for delta=1
    k_analyt=zeros(1,T);
    k_analyt(1)=k_0;
    k_analyt_finite(1)=k_0;
    for t=2:T
        k_analyt(t)=(alpha*beta)*k_analyt(t-1)^alpha;
        temp=0;
        if t<T
            for i=t:T-1
                temp=alpha*beta*(1+temp);
            end
        end
        coef(t)=temp/(1+temp);
        k_analyt_finite(t)=coef(t)*k_analyt_finite(t-1)^alpha;
    end
    for i=1:N
        kprime_analyt(i)=(alpha*beta)*kgrid(i)^alpha;
    end
end

% ==============
% 2. Value function iteration (solve with and w/o Howard / policy function iteration), infinite T
% ==============

%pre-allocate matrices c, u for speed
u = zeros(N,N);
c = zeros(N,N);

% one period return
for i=1:N
    for j=1:N %Think if kgrid(j) is correct or should maybe be kprime for the - kgrid(j)
        c(i,j)= kgrid(i)^alpha - kgrid(j) + (1-delta)*kgrid(i);
        if c(i,j)>0
            u(i,j)=(c(i,j)^(1-sigma) -1)/(1-sigma);
        else
            %But I think the above should be correct, bc of same structure
            %here
            u(i,j)=-1e50*((kgrid(i)^alpha+(1-delta)*kgrid(i)-kgrid(j))<=0);
            %we penalize the negative value of c by giving them an utility
            %of minus infinity
        end
    end
end

% initial guesses and preallocations
dV=1;

%V= k_0^alpha - kgrid + (1-delta)*k_0;
%V= (V.^(1-sigma) -1)/(1-sigma);
% Alternative initial value definition
V = zeros(1,N);
VV = zeros(N,N);
Vnew = zeros(1,N);
kprime = zeros(1,N);
Vnewhow = zeros(1,N);
index = zeros(1,N); %stores the index of the optimal policy kprime(kgrid(i))
%kgrid(index(i))=kprime(i)
iter=0;
tic

% main loop
while dV>criter_V
    iter=iter+1;
    for i=1:N % loop over capital today
        for j=1:N %... and over capital tomorrow
            VV(i,j)= u(i,j) + beta*V(j); %this is without the interpolation
        end
        % take maximum over capital tomorrow
        [Vnew(i), maxI]= max(VV(i,:));
        index(i) = maxI;
        % record policy function - doesnt make sense
        kprime(i) = kgrid(maxI);
    end
    % Howard - doesn't help much here
    if Howard==1 && iter>3 
        dVV=1;
        while dVV>criter_V
        for i=1:N 
            temp = kgrid(i)^alpha - kprime(i) + (1-delta)*kgrid(i);
            temp = (temp^(1-sigma) -1)/(1-sigma);
            Vnewhow(i)=temp + beta*Vnew(index(i)); %Vnew(index(i))=V evaluated in kprime(i)
            clear temp
        end
        dVV=max(abs(Vnewhow-Vnew));
        Vnew=Vnewhow;
        disp("dVHoward");
        disp(dV);
         

        end
    end 
    % calculate convergence criterion
    % but basically, we stop if the old-new difference in values is small
    dV=max(abs(Vnew-V)); %works also without doubling max
    % updated value function
    V=Vnew;
    disp('dV')
    disp(dV)
    iter=iter+1;
end

V_disc_VFI=V;
kprime_VFI = kprime;
toc

%% Build plots for policy functions - Plot policy function k'(k)
hold on
title("K' policy function plot")
plot(kgrid, kprime), xlabel('Capital values at t'), ylabel('Capital level at t+1');
hold off
%%
hold on
title("Value function plot by capital values")
plot(kgrid, V), xlabel('Capital values at t'), ylabel('Value function result at t');
hold off
%%

% Euler equation errors in percent

% consumption vector today
c1= kgrid.^alpha + (1- delta)*kgrid - kprime;
% consumption vector at choice kprime tomorrow
c2= interp1(kgrid, c1, kprime,'linear','extrap');
% marginal productivity
margprod=alpha.*kprime.^(alpha-1) + 1 - delta;

EEerror_disc=(c1 - beta.*margprod.^(-1/sigma).*c2)./c1;
maxEEerror_disc=max(abs(EEerror_disc));

hold on
title("Euler equation error with corresponding values")
plot(kgrid, EEerror_disc), xlabel('Capital values at t'), ylabel('Euler equation error');
hold off
%%
% ==============
% 3. Value function iteration with interpolation
% ==============
% set options for fminsearch
options=optimset('MaxIter',5000,'MaxFunEval',5000,'TolFun',1e-12);
% initial guesses

dV=1;
V=V_disc_VFI;%zeros(N,1);
iter=0;
tic
kprime_VFI_cont=kprime_VFI; %initial guess


tic

% fminsearch
while dV>criter_V
    iter=iter+1;
    kprimelow=min(kgrid); kprimehigh=1.3*min(kgrid);
    for i=1:N % loop over capital today
        % find maximum over capital tomorrow - now using interpolation
        [kprime_VFI_cont(i),Vnew(i)]=fminsearch(@(x) Valuefun(x,kgrid,kgrid(i),alpha,sigma,V,delta,beta),kprime_VFI_cont(i),options);
    end
    Vnew=-Vnew; % take negative as Valuefun is for minimisation
    % calculate convergence criterion
    dV=max(max(abs(Vnew-V)));
    % updated value function
    V=Vnew;
    %disp('dV')
    %disp(dV)
end
t_fmins=toc;
%%
% with Golden Search
dV=1;
ctemp=kgrid.^alpha+(1-delta)*kgrid;
V=V_disc_VFI;%(ctemp.^(1-sigma)-1)/(1-sigma)
kprime_VFI_contG = zeros(1,N);
VnewG = zeros(1,N);
iter=0;        
alpha1 = (3-sqrt(5))/2;
alpha2 = (sqrt(5)-1)/2;
tic
while dV>criter_V
    iter=iter+1;

    for i=1:N % loop over capital today
        
        if i==1
            kprimelow=min(kgrid); 
            kprimehigh=max(kgrid);
        else
            kprimelow=max(kprimelow,min(kgrid)); 
            kprimehigh=1.2*kprimelow;
        end
        b1=kprimelow+alpha1*(kprimehigh-kprimelow); %lower bound candidate
        b2=kprimelow+alpha2*(kprimehigh-kprimelow); %upper bound candidate
        %Starting value function results
        Vlow=Valuefunpos(kprimelow,kgrid,kgrid(i),alpha,sigma,V,delta,beta);
        Vhigh=Valuefunpos(kprimehigh,kgrid,kgrid(i),alpha,sigma,V,delta,beta);
        %Candidate value function results
        % Lower bound cand result
        Vb1=Valuefunpos(b1,kgrid,kgrid(i),alpha,sigma,V,delta,beta);
        % Upper bound cand result
        Vb2=Valuefunpos(b2,kgrid,kgrid(i),alpha,sigma,V,delta,beta);
        dk=1;
        criter_k=1e-12;
        while dk>criter_k
        % use golden search
        if Vb2>Vb1                                            
            kprimelow=b1; 
            Vlow=Valuefunpos(kprimelow,kgrid,kgrid(i),alpha,sigma,V,delta,beta);

            b1=kprimelow+alpha1*(kprimehigh-kprimelow);
            Vb1=Valuefunpos(b1,kgrid,kgrid(i),alpha,sigma,V,delta,beta);

            b2=kprimelow+alpha2*(kprimehigh-kprimelow);
            Vb2=Valuefunpos(b2,kgrid,kgrid(i),alpha,sigma,V,delta,beta);
        else
            kprimehigh=b2;
            Vhigh=Valuefunpos(kprimehigh,kgrid,kgrid(i),alpha,sigma,V,delta,beta);

            b2=kprimelow+alpha2*(kprimehigh-kprimelow);
            Vb2=Valuefunpos(b2,kgrid,kgrid(i),alpha,sigma,V,delta,beta);

            b1=kprimelow+alpha1*(kprimehigh-kprimelow);
            Vb1=Valuefunpos(b1,kgrid,kgrid(i),alpha,sigma,V,delta,beta);
        end
        dk=abs(Vhigh-Vlow);
        end
        kprime_VFI_contG(i)=1/2*(kprimelow+kprimehigh);
        VnewG(i)=1/2*(Vlow+Vhigh);

        % disp([Vnew(i,1),kprime_VFI_cont(i)])
    end
    
    % calculate convergence criterion
    dV=max(max(abs(VnewG-V)));
    % updated value function
    V=VnewG;    
    disp('dV')
    disp(dV)
end
t_G=toc;

%%

% Euler equation errors in percent

% consumption vector today
c1_cont= kgrid.^alpha + (1- delta)*kgrid - kprime_VFI_contG;
% consumption vector at choice kprime tomorrow
c2_cont= interp1(kgrid, c1_cont, kprime_VFI_contG,'linear','extrap');
% marginal productivity
margprod_cont=alpha.*kprime_VFI_contG.^(alpha-1) + 1 - delta;

EEerror_cont=(c1_cont - beta.*margprod_cont.^(-1/sigma).*c2_cont)./c1_cont;
maxEEerror_cont=max(abs(EEerror_cont));

hold on
title("Euler equation error with corresponding values - Golden Search")
plot(kgrid, EEerror_cont), xlabel('Capital values at t'), ylabel('Euler equation error');
hold off

%%
% ==============
% Figures
% ==============
% plot policy function
figure(1)
% levels
subplot(2,1,1)
title("K' policy function plot with Golden Search");
hold on
if delta==1 && abs(sigma-1)<0.001
    kprime_analyt_pol=alpha*beta*kgrid.^alpha;
    plot(kgrid,kprime_analyt_pol,'b-','Linewidth',1)
end

%plot(kgrid,kprime_VFI,'k-','Linewidth',1)
%plot(kgrid,kprime_VFI_cont,'k--','Linewidth',1)
%plot(kgrid,kprime_VFI_contG,'k--','Linewidth',1)
plot(kgrid, kprime_VFI, kgrid, kprime_VFI_contG)
xlabel('Capital values at t'), ylabel('Capital values at t+1');

h = legend('VFI path', 'Golden Search path' ,'Location', 'best','Orientation','Vertical');
h.Title.String = 'Search methods';

set(h,'fontsize',12,'Interpreter','Latex')
