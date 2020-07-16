function p=FindMin(X0,spec,slope,intercept);

alpha=X0(1);
mu=X0(2);
sigma1=X0(3);
sigma2=X0(4);
pwr=X0(5);
p=Chi2(spec,alpha,mu,sigma1,sigma2,pwr,slope,intercept);
