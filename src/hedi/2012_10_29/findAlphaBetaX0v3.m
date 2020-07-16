function mz=findAlphaBetaX0v3(f,spec,X0);
% give best alpha and beta using brute force methods



[pmax,imax]=max(spec);
spec=spec/pmax;
fmax=f(imax);
fband=f(2)-f(1);
fmin=500; % 500 Hz minimal frequency 
mz=[];

alpha0=X0(1);
beta0=X0(2);
x1=X0(3);
y1=X0(4);
x2=X0(5);
y2=X0(6);
x3=X0(7);
y3=X0(8);

x0=[alpha0,beta0,x1,y1,x2,y2,x3,y3];
g=@(x)CHI2specV3(x,spec);
options = optimset('MaxIter',200);
xf=fminsearch(g,x0,options);
mz=[xf g(xf)];