function xx=util(x,K,V,kgrid,z,gamma,theta,psi,alpha,beta,delta,p)
c=x(1);
L=x(2);
if c<=0 || L<0 || L>1 || c>exp(z)*K^alpha*L^(1-alpha)+(1-delta)*K
    xx=-100000000000;
    else
        kprime=exp(z)*K^alpha*L^(1-alpha)+(1-delta)*K-c;
        %keyboard
        xx=[Fill this in to calculate utility including the continuation value V at kprime],'linear')';
    end
end
% 
% function x=util(x,K,V,kgrid,z,gamma,theta,psi,alpha,beta,delta,p)
% c=x(1);
% L=x(2);
% if c<=0 || L<0 || L>1 || c>exp(z)*K^alpha*L^(1-alpha)+(1-delta)*K
%     x=-100000000000;
% else
%     kprime=exp(z)*K^alpha*L^(1-alpha)+(1-delta)*K-c;
%     %keyboard
%     x=(c^(1-gamma)-1)/(1-gamma)+theta*(1-L)^(1+psi)/(1+psi)+beta*p*interp1(kgrid,V,kprime,'linear','extrap')';
% end
% end