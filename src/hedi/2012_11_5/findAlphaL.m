function mz=findAlphaL(f,spec,X0,slope,intercept,alpha0);
% give best alpha and beta using brute force methods



[pmax,imax]=max(spec);
spec=spec/pmax;
fmax=f(imax);
fband=f(2)-f(1);
fmin=500; % 500 Hz minimal frequency 
mz=[];

%alpha0=a;%%X0(1);
beta0=slope*alpha0+intercept;

a=X0(1);
L=X0(2);

x0=[a,L];
g=@(x)CHI2VC(x,spec,slope,intercept,alpha0);
options = optimset('MaxIter',100,'Display','off');
xf=fminsearch(g,x0,options);
mz=[xf g(xf)];