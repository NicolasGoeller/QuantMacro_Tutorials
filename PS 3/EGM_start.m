% grid for kprime
kprimegrid=linspace(0.9*kbar,1.1*kbar,100);
% cash-on-hand next period
xprime=[fill this in]; % first dimension is prod, 2nd kprime
% guess for policy
kprime=0.01*kbar+0.99*kprimegrid;
% initialise 
dkprime=1;
iter=0;
tic
% loop over policy function
while dkprime>0.00000001
    iter=iter+1;
    % consumption next period given xprime, kprime
          CPRIME=max(0.0000000001,[fill this in]);
          % Marg utility
            UPRIME=[fill this in];
            % marg utility as function of kprime, zprime
        MPROD=t[fill this in];
        % expected marg utility given Z(k) today
        for k=1:M
            EMU=[fill this in];
            c_EGM(k,:)=([fill this in]);
        end

        % cash-on-hand implied by new c() and kprime

x=[fill this in];
% new kprime as interpolation of kprimegrid (x) at xprime
kprimenew=interpn([fill this in])

dkprime=max(max(abs(kprimenew-kprime)));
%update policy
kprime=0.9*kprime+0.1*kprimenew;

end
time_EGM=toc;