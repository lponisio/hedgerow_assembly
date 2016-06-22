function x=uniform_rand(mu,fvar,m,n)

x=mu*(1+fvar*(2*rand(m,n)-1));

end