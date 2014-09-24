function mz=findAlphaBetaX0(f,spec,X0);
% give best alpha and beta using brute force methods



[pmax,imax]=max(spec);
spec=spec/pmax;
fmax=f(imax);
fband=f(2)-f(1);
fmin=500; % 500 Hz minimal frequency 
mz=[];

alpha0=X0(1);
beta0=X0(2);
fcenter=X0(3);
fwidth=X0(4);
x0=[alpha0,beta0,fcenter,fwidth];
g=@(x)CHI2spec(x,spec);
options = optimset('MaxIter',200);
xf=fminsearch(g,x0,options);
mz=[xf g(xf)];