

    % function that gives rhs of bellman equation
    function xx=Valuefunpos(kprime,kgrid,k,alpha,sigma,V,delta,beta)
    c=k.^alpha+(1-delta)*k-kprime;
    if c<0
        xx=-10000000;
    else
    u=(c^(1-sigma)-1)/(1-sigma);
    Vprime=interp1(kgrid,V,kprime,'linear','extrap');
    xx=u+beta*Vprime;
  
    end
    % This here is the only difference to the other valuefunpos
      xx=xx;
    end