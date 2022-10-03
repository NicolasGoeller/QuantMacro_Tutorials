

    % function that gives rhs of bellman equation
        function xx=Valuefunpos(kprime,kgrid,k,theta,gamma,V,delta,beta)
    c=k.^theta+(1-delta)*k-kprime;
    if c<0
        xx=-10000000;
    else
    u=(c^(1-gamma)-1)/(1-gamma);
    Vprime=interp1(kgrid,V,kprime,'linear','extrap');
    xx=u+beta*Vprime;
  
    end
      xx=xx;
    end