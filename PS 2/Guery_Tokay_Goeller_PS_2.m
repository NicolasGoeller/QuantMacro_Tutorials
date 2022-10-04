% ==========================
% Problem set solving the stochastic neoclassical growth model different
% ways
% ==========================
% clear the workspace
clear
% close all figures
close all
addpath('C:\Users\tbroer\Dropbox\Teaching\PSE\2021 Quantitative Macro\Problem sets\PS RBC')
addpath('C:\Users\tbroer\Dropbox\Teaching\PSE\2021 Quantitative Macro\Problem sets')

% ============
% parameters
% ============
alpha=0.4; % capital share - this is alpha
beta = 0.99; % discount factor
rho = 0.95;   % persistence of TFP shock
sigma = 2.00000001; % CRRA coefficient (for 1 equals log, but need to replace the function, so set close to 1)
delta=1;

% ============
% options and convergence criteria
% ============

Howard =0; % set to 1 if you want to do policy fct iteration / Howard improvement
criter_V = 1e-7; % conv criterion for value function
N=50; % number of grid points
linear=1; % grid linear or not

%mean of capital non-stochastic steady state
kbar=((1/beta-1+delta)/(alpha))^(1/(alpha-1));


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
        c(i,j)= kgrid(i)^alpha*pshock(i) - kgrid(j) + (1-delta)*kgrid(i);
        if c(i,j)>0
            u(i,j)=(c(i,j)^(1-sigma) -1)/(1-sigma);
        else
            %But I think the above should be correct, bc of same structure
            %here
            u(i,j)=-1e50*((kgrid(i)^alpha+(1-delta)*kgrid(i)-kgrid(j))<=0);
        end
    end
end

% initial guesses
dV=1;
V=[ xxxx FILL THIS IN xxxx ]
iter=0;
tic

% main loop
while dV>criter_V
    iter=iter+1;
    for i=1:N % loop over capital today
        for j=1:N %... and over capital tomorrow
            VV(i,j)=value; % calculate value for each i,j
        end
        % take maximum over capital tomorrow
        [ xxxx FILL THIS IN xxxx ]
        % record policy function
        [ xxxx FILL THIS IN xxxx ]
    end
    % Howard - doesn't help much here
    if Howard==1 && iter>3
        dVV=1;
        while dVV>criter_V
            for i=1:N
                Vnewhow(i)=[ xxxx FILL THIS IN xxxx ]
                clear temp
            end
            dVV=max(max(abs(Vnewhow-Vnew)))
            %disp(dVV)
            Vnew=Vnewhow;
        end
    end
    
    % calculate convergence criterion - why the double max() here?
    % but basically, we stop if the old-new difference in values is small
    dV=max(max(abs(Vnew-V)));
    % updated value function
    V=Vnew;
    disp('dV')
    disp(dV)
end

V_disc_VFI=V;
toc

% Euler equation errors in percent

% consumption vector today
c1=[ xxxx FILL THIS IN xxxx ]
% consumption vector at choice kprime tomorrow
c2=[ xxxx FILL THIS IN xxxx ]
% marginal productivity
margprod=[ xxxx FILL THIS IN xxxx ]

EEerror_disc=[ xxxx FILL THIS IN xxxx ]
maxEEerror_disc=max(abs(EEerror_disc))

% ==============
% 3. Value function iteration with interpolation
% ==============
% set options for fminsearch
options=optimset('MaxIter',5000,'MaxFunEval',5000,'TolFun',1e-12);
% initial guesses
dV=1;
V=V_disc_VFI%zeros(N,1);
iter=0;
tic
kprime_VFI_cont=kprime_VFI; %initial guess

        alpha1 = (3-sqrt(5))/2;
        alpha2 = (sqrt(5)-1)/2;
tic

% fminsearch
while dV>criter_V
    iter=iter+1;
    kprimelow=min(kgrid); kprimehigh=1.3*min(kgrid);
    for i=1:N % loop over capital today
        % find maximum over capital tomorrow - now using interpolation
        [kprime_VFI_cont(i),Vnew(i,1)]=fminsearch(@(x) Valuefun(x,kgrid,kgrid(i),alpha,sigma,V,delta,beta),kprime_VFI_cont(i),options);
    end
    Vnew=-Vnew; % take negative as Valuefun is for minimisation
    % calculate convergence criterion
    dV=max(max(abs(Vnew-V)));
    % updated value function
    V=Vnew;
    %disp('dV')
    %disp(dV)%
end
t_fmins=toc;

% with Golden Search
dV=1;
ctemp=kgrid.^alpha+(1-delta)*kgrid;
V=V_disc_VFI;%(ctemp.^(1-sigma)-1)/(1-sigma)
iter=0;
tic
while dV>criter_V
    iter=iter+1;

    for i=1:N % loop over capital today
        
        if i==1
            kprimelow=min(kgrid); kprimehigh=max(kgrid);
        else
            kprimelow=max(kprimelow,min(kgrid)); kprimehigh=1.2*kprimelow;
        end
        b1=kprimelow+alpha1*(kprimehigh-kprimelow);
        b2=kprimelow+alpha2*(kprimehigh-kprimelow);
        Vlow=Valuefunpos(kprimelow,kgrid,kgrid(i),alpha,sigma,V,delta,beta);
        Vhigh=Valuefunpos(kprimehigh,kgrid,kgrid(i),alpha,sigma,V,delta,beta);
        Vb1=Valuefunpos(b1,kgrid,kgrid(i),alpha,sigma,V,delta,beta);
        Vb2=Valuefunpos(b2,kgrid,kgrid(i),alpha,sigma,V,delta,beta);
        dk=1;
        criter_k=1e-12;
        while dk>criter_k
        % use golden search
        if Vb2>Vb1                                            
        kprimelow=[ xxxx FILL THIS IN xxxx ]; Vlow=[ xxxx FILL THIS IN xxxx ];
        b1=[ xxxx FILL THIS IN xxxx ];Vb1=[ xxxx FILL THIS IN xxxx ]2;
        b2=[ xxxx FILL THIS IN xxxx ];
        Vb2=[ xxxx FILL THIS IN xxxx ];
        else
        kprimehigh=[ xxxx FILL THIS IN xxxx ];Vhigh=[ xxxx FILL THIS IN xxxx ];
        b2=[ xxxx FILL THIS IN xxxx ];Vb2=[ xxxx FILL THIS IN xxxx ];
        b1=[ xxxx FILL THIS IN xxxx ];
        Vb1=[ xxxx FILL THIS IN xxxx ];
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



% Euler equation errors in percent
% consumption vector today
c1=[ xxxx FILL THIS IN xxxx ]
% consumption vector at choice kprime tomorrow
c2=[ xxxx FILL THIS IN xxxx ]
% marginal productivity
margprod=[ xxxx FILL THIS IN xxxx ]

EEerror_cont=[ xxxx FILL THIS IN xxxx ]
maxEEerror_cont=max(abs(EEerror_cont));


% ==============
% Figures
% ==============
% plot policy function
figure(1)
% levels
subplot(2,1,1)
hold on
if delta==1 && abs(sigma-1)<0.001
    kprime_analyt_pol=alpha*beta*kgrid.^alpha;
    plot(kgrid,kprime_analyt_pol,'b-','Linewidth',1)
end
plot(kgrid,kprime_VFI,'k-','Linewidth',1)
plot(kgrid,kprime_VFI_cont,'k--','Linewidth',1)
plot(kgrid,kprime_VFI_contG,'k--','Linewidth',1)
