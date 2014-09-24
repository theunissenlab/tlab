function [alpha_f,mu_f,sigma_f1,sigma_f2,pwr_f,d]=GetMotoSpec(spec,alpha,mu,sigma1,sigma2,pwr)
slope=-0.6025;
intercept=-0.5979;

x0=[alpha,mu,sigma1,sigma2,pwr];
options = optimset('MaxIter',50,'Display','off');

f=@(x)FindMin(x,spec,slope,intercept);
xf=fminsearch(f,x0,options);
alpha_f=xf(1);
mu_f=xf(2);
sigma_f1=xf(3);
sigma_f2=xf(4);
pwr_f=xf(5);%sum(spec);
d=Chi2(spec,alpha_f,mu_f,sigma_f1,sigma_f2,pwr_f,slope,intercept);

