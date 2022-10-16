function error=FOC(x,K,z,gamma,theta,psi,alpha,beta,delta,kprime)
%c=x(1);
L=x(1);
if L<0 || L>1
    error(1)=100000;
    %error(2)=100000;
else
%error(1)=;%+1000000*(c<0)+1000000*(L<0)+1000000*(L>0);
error(1)=[Fill this in];
end
end