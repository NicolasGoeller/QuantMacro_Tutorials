function [V,kprime,index] = discrete_search(alpha, delta, gamma_c, beta, criter_V, M,N,Z_tauchen,kgrid)

%This function return the grid of optimal value of the value function, of
%kprime given k and z

%index is the matrix of index of kprime in kgrid defined as follows kgrid(index(i,j))=kprime(i,j)


%M is the size of the grid for value of z
%N is the size of the grid for value of k
%Z_tauchen is the grid of value of z
%kgrid is the grid of value of k


u = zeros(M,M,N);
c = zeros(M,M,N);
dV=1;

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
end

V = zeros(M,N);
VV = zeros(M,M,N);
Vnew = zeros(M,N);
kprime = zeros(M,N);
index = zeros(M,N); %stores the index of the optimal policy kprime(kgrid(i))
%we have the relation kgrid(index(i,j))=kprime(i,j)
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