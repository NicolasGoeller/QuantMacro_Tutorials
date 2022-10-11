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

%% Problem Set 3


% ============
% parameters
% ============
alpha=0.4; % capital share - this is alpha
beta = 0.987; % discount factor
rho = 0.95;   % persistence of TFP shock
gamma_c = 2.00000001; % CRRA coefficient (for 1 equals log, but need to replace the function, so set close to 1)
delta=0.1;
sigma = 0.007;


% ============
% options and convergence criteria
% ============

Howard =1; % set to 1 if you want to do policy fct iteration / Howard improvement
criter_V = 1e-6; % conv criterion for value function
M=50; % number of grid points
N=5; % grid for z
linear=1; % grid linear or not
N_sim = 100; % nbr of simulation
N=5;
T = 150;%period of transition
%mean of capital non-stochastic steady state
kbar=((1/beta-1+delta)/(alpha))^(1/(alpha-1));
k_0 = kbar;



% ==============
% 0. Grids,  etc
% ==============

% center the grid around mean of capital non-stochastic steady state
kbar=((1/beta-1+delta)/(alpha))^(1/(alpha-1));
% the grid
if delta==1
    if linear==1
        % linear sequence kbar -2kbar in N steps
        kgrid=linspace(kbar*0.85,1.15*kbar,M);
    else
        % linear sequence 0-0.5 in N steps, divided by 0.5 all to the power
        % of 5; times 1.5kbar
        temp=linspace(0,0.5,M).^5/0.5^5*(0.85*kbar-1.15*kbar);
        % 0.5kbar + temp
        kgrid=0.85*kbar+temp;
    end
else
    % if delta ~= 1, linear sequence 0.25kbar-2kbar in N steps
    kgrid=linspace(0.85*kbar ,1.15*kbar,M);
end


%% Problem 1 - discretize income process and simulate
% ==============

[Z_tauchen, P_tauchen] = tauchen(N,0,rho,sigma,2);
p = dtmc(P_tauchen);
X0 = [0 0 N_sim 0 0]; %simulate N_sim starting at initial value of Z=0 
X = simulate(p,T,"X0",X0); %X0 define the number of simulation to be made starting with a given initial condition
for i=1:N 
      X(X==i)=Z_tauchen(i,1);
end
Mean_X= mean(X,1);
a = mean(Mean_X);
Std_X = std(X);
b = mean(Std_X);
[Acf_x,lag] = autocorr(X(:,4));
graphplot(p,'ColorEdges',true);

figure;
simplot(p,X);
A = zeros(1,length(Z_tauchen));
A=A';


% disp('Standard devations of Z discrete, continuous')
% disp([mean(std(log(Z_sim)')),(sigmaepsilon^2/(1-rho^2))^0.5])
% for i=1:N_sim
%     rho_emp(i)=corr(log(Z_sim(i,1:end-1))',log(Z_sim(i,2:end))');
% end
% disp('Autocorrelation of Z discrete, continuous')
% 
% disp([mean(rho_emp),rho])


%% Problem 2 : Discrete grid value function iteration 

%pre-allocate matrices c, u for speed
u = zeros(M,M,N);
c = zeros(M,M,N);

% one period return
for i=1:M
    for l=1:M %Think if kgrid(j) is correct or should maybe be kprime for the - kgrid(j)
        for j=1:N %loop on value of z
            c(i,l,j)= exp(Z_tauchen(j))*kgrid(i)^alpha - kgrid(l) + (1-delta)*kgrid(i);
            if c(i,l,j)>0
                u(i,l,j)=(c(i,l,j)^(1-gamma_c) -1)/(1-gamma_c);
            else
                %But I think the above should be correct, bc of same structure
                %here
                u(i,l,j)=-1e50*((exp(Z_tauchen(j))*kgrid(i)^alpha+(1-delta)*kgrid(i)-kgrid(l))<=0);
                %we penalize the negative value of c by giving them an utility
                %of minus infinity
            end
        end
    end



%V= k_0^alpha - kgrid + (1-delta)*k_0;
%V= (V.^(1-sigma) -1)/(1-sigma);
% Alternative initial value definition
V = zeros(M,N);
VV = zeros(M,M,N);
Vnew = zeros(M,N);
kprime = zeros(M,N);
%Vnewhow = zeros(M,N);
index = zeros(M,N); %stores the index of the optimal policy kprime(kgrid(i))
%kgrid(index(i))=kprime(i)
iter=0;
tic

% main loop
while dV>criter_V
    iter=iter+1;
    for i=1:M % loop over capital today
        for j=1:N %... and over value of Z
            for l=1:M %... and over capital tomorrow
                VV(i,l,j)= u(i,l,j) + beta*V(l,j); %this is without the interpolatio
            end
             % take maximum over capital tomorrow
            [Vnew(i,j), maxI]= max(VV(i,:,j));
            index(i,j) = maxI;
            % record policy function - doesnt make sense
            kprime(i,j) = kgrid(maxI);
        end
    end
    % Howard - doesn't help much here
%     if Howard==1 && iter>3 
%         dVV=1;
%         while dVV>criter_V
%         for i=1:M 
%             temp = kgrid(i)^alpha - kprime(i) + (1-delta)*kgrid(i);
%             temp = (temp^(1-gamma_c) -1)/(1-gamma_c);
%             Vnewhow(i)=temp + beta*Vnew(index(i)); %Vnew(index(i))=V evaluated in kprime(i)
%             clear temp
%         end
%         dVV=max(abs(Vnewhow-Vnew));
%         Vnew=Vnewhow;
%         disp("dVHoward");
%         disp(dV);
%          
% 
%         end
%     end 
    % calculate convergence criterion
    % but basically, we stop if the old-new difference in values is small
    dV=max(max(abs(Vnew-V))); 
    % updated value function
    V=Vnew;
    disp('dV')
    disp(dV)
    iter=iter+1;
end

V_disc_VFI=V;
kprime_VFI = kprime;
toc

% Build plots for policy functions - Plot policy function k'(k)
hold on
title("K' policy function plot")
plot(kgrid, kprime), xlabel('Capital values at t'), ylabel('Capital level at t+1');
h = legend('log(z)=-0.0448','log(z)=-0.0224','log(z)=0','log(z)=0.0224','log(z)=0.0448', ...
    'Location', 'best','Orientation','Vertical');
h.Title.String = 'log(z) values';
set(h,'fontsize',12,'Interpreter','Latex')
hold off

%
hold on
title("Value function plot by capital values")
plot(kgrid, V), xlabel('Capital values at t'), ylabel('Value function result at t');
h = legend('log(z)=-0.0448','log(z)=-0.0224','log(z)=0','log(z)=0.0224','log(z)=0.0448', ...
    'Location', 'best','Orientation','Vertical');
h.Title.String = 'log(z) values';
set(h,'fontsize',12,'Interpreter','Latex')
hold off

%% Euler equation errors in percent

c1 = zeros(M,N);
c2 = zeros(M,N);
margprod = zeros(M,N);
EEerror_disc = zeros(M,N);
%maxEEerror_disc = 


% consumption vector today
for j=1:N
    c1(:,j)=exp(Z_tauchen(j))*kgrid'.^alpha + (1- delta)*kgrid' - kprime(:,j);
    % consumption vector at choice kprime tomorrow
    c2(:,j)= interp1(exp(Z_tauchen(j))*kgrid', c1(:,j), kprime(:,j),'linear','extrap');
    % marginal productivity
    margprod(:,j)=alpha.*kprime(:,j).^(alpha-1) + 1 - delta;
    EEerror_disc(:,j)=abs((c1(:,j) - beta.*margprod(:,j).^(-1/gamma_c).*c2(:,j))./c1(:,j));
    %maxEEerror_disc(:,j)=max(abs(EEerror_disc(:,j)));
end 

%Graph of the EE error 

hold on
title("Euler equation error with corresponding values")
plot(kgrid, EEerror_disc), xlabel('Capital values at t'), ylabel('Euler equation error');
h = legend('log(z)=-0.0448','log(z)=-0.0224','log(z)=0','log(z)=0.0224','log(z)=0.0448', ...
    'Location', 'best','Orientation','Vertical');
h.Title.String = 'log(z) values';
set(h,'fontsize',12,'Interpreter','Latex')
hold off

%% 3. Log-linearization
% =====================

% some ratios as function of parameters
ybar=kbar^alpha;
cbar=ybar-delta*kbar;
% need for matrix form of loglin
ckrat=cbar/kbar;

% a. write system as A E[y_t+1]+B y_t=0



% order c,k,z
A=[-gamma_c, beta *(alpha-1)*alpha*kbar^(alpha-1) ; 0 , 1 ];

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

% initial guesses and preallocations
dV=1;
%If the Blanchard Kahn condition is satisfied, we can find the expression
%of "alpha1"
if BKcond~=1
    disp('BK conditions not satisfied')
else
    bkev =find(abs(diag(lambda))>1);
    invP=aaa(bkev,:);%%Select the element of the invert of the vector matrix needed to compute the policy function
    polfunc= -invP(1,2)/invP(1);% you need to find the policy for consumption here : derived analytically 
end



%% Problem SET 2



% ==============
% 1. analytical case delta=1, finite T and infinite T
% ==============

alpha=0.4; % capital share - this is alpha
beta = 0.99; % discount factor
rho = 0.95;   % persistence of TFP shock
gamma_c = 2.00000001; % CRRA coefficient (for 1 equals log, but need to replace the function, so set close to 1)
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
    for i=1:M
        kprime_analyt(i)=(alpha*beta)*kgrid(i)^alpha;
    end
end

% ==============
% 2. Value function iteration (solve with and w/o Howard / policy function iteration), infinite T
% ==============

%pre-allocate matrices c, u for speed
u = zeros(M,M);
c = zeros(M,M);

% one period return
for i=1:M
    for l=1:M %Think if kgrid(j) is correct or should maybe be kprime for the - kgrid(j)
        c(i,l)= kgrid(i)^alpha - kgrid(l) + (1-delta)*kgrid(i);
        if c(i,l)>0
            u(i,l)=(c(i,l)^(1-gamma_c) -1)/(1-gamma_c);
        else
            %But I think the above should be correct, bc of same structure
            %here
            u(i,l)=-1e50*((kgrid(i)^alpha+(1-delta)*kgrid(i)-kgrid(l))<=0);
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
V = zeros(1,M);
VV = zeros(M,M);
Vnew = zeros(1,M);
kprime = zeros(1,M);
Vnewhow = zeros(1,M);
index = zeros(1,M); %stores the index of the optimal policy kprime(kgrid(i))
%kgrid(index(i))=kprime(i)
iter=0;
tic

% main loop
while dV>criter_V
    iter=iter+1;
    for i=1:M % loop over capital today
        for l=1:M %... and over capital tomorrow
            VV(i,l)= u(i,l) + beta*V(l); %this is without the interpolation
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
        for i=1:M 
            temp = kgrid(i)^alpha - kprime(i) + (1-delta)*kgrid(i);
            temp = (temp^(1-gamma_c) -1)/(1-gamma_c);
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

EEerror_disc=(c1 - beta.*margprod.^(-1/gamma_c).*c2)./c1;
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
    for i=1:M % loop over capital today
        % find maximum over capital tomorrow - now using interpolation
        [kprime_VFI_cont(i),Vnew(i)]=fminsearch(@(x) Valuefun(x,kgrid,kgrid(i),alpha,gamma_c,V,delta,beta),kprime_VFI_cont(i),options);
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
kprime_VFI_contG = zeros(1,M);
VnewG = zeros(1,M);
iter=0;        
alpha1 = (3-sqrt(5))/2;
alpha2 = (sqrt(5)-1)/2;
tic
while dV>criter_V
    iter=iter+1;

    for i=1:M % loop over capital today
        
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
        Vlow=Valuefunpos(kprimelow,kgrid,kgrid(i),alpha,gamma_c,V,delta,beta);
        Vhigh=Valuefunpos(kprimehigh,kgrid,kgrid(i),alpha,gamma_c,V,delta,beta);
        %Candidate value function results
        % Lower bound cand result
        Vb1=Valuefunpos(b1,kgrid,kgrid(i),alpha,gamma_c,V,delta,beta);
        % Upper bound cand result
        Vb2=Valuefunpos(b2,kgrid,kgrid(i),alpha,gamma_c,V,delta,beta);
        dk=1;
        criter_k=1e-12;
        while dk>criter_k
        % use golden search
        if Vb2>Vb1                                            
            kprimelow=b1; 
            Vlow=Valuefunpos(kprimelow,kgrid,kgrid(i),alpha,gamma_c,V,delta,beta);

            b1=kprimelow+alpha1*(kprimehigh-kprimelow);
            Vb1=Valuefunpos(b1,kgrid,kgrid(i),alpha,gamma_c,V,delta,beta);

            b2=kprimelow+alpha2*(kprimehigh-kprimelow);
            Vb2=Valuefunpos(b2,kgrid,kgrid(i),alpha,gamma_c,V,delta,beta);
        else
            kprimehigh=b2;
            Vhigh=Valuefunpos(kprimehigh,kgrid,kgrid(i),alpha,gamma_c,V,delta,beta);

            b2=kprimelow+alpha2*(kprimehigh-kprimelow);
            Vb2=Valuefunpos(b2,kgrid,kgrid(i),alpha,gamma_c,V,delta,beta);

            b1=kprimelow+alpha1*(kprimehigh-kprimelow);
            Vb1=Valuefunpos(b1,kgrid,kgrid(i),alpha,gamma_c,V,delta,beta);
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

EEerror_cont=(c1_cont - beta.*margprod_cont.^(-1/gamma_c).*c2_cont)./c1_cont;
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
if delta==1 && abs(gamma_c-1)<0.001
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

%%


