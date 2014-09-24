function mz=findAlphaBetaX0v2(f,spec,X0);
% give best alpha and beta using brute force methods



[pmax,imax]=max(spec);
spec=spec/pmax;
fmax=f(imax);
fband=f(2)-f(1);
fmin=500; % 500 Hz minimal frequency 
mz=[];

alpha0=X0(1);
beta0=X0(2);
fcenter1=X0(3);
fwidth1=X0(4);
fcenter2=X0(5);
fwidth2=X0(6);
fcenter3=X0(7);
fwidth3=X0(8);

x0=[alpha0,beta0,fcenter1,fwidth1,fcenter2,fwidth2,fcenter3,fwidth3];
g=@(x)CHI2specV2(x,spec);
options = optimset('MaxIter',1000);
xf=fminsearch(g,x0,options);
mz=[xf g(xf)];