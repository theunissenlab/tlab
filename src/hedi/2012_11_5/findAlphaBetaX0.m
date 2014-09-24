function mz=findAlphaBetaX0(f,spec,X0,slope,intercept);
% give best alpha and beta using brute force methods



[pmax,imax]=max(spec);
spec=spec/pmax;
fmax=f(imax);
fband=f(2)-f(1);
fmin=500; % 500 Hz minimal frequency 
mz=[];

alpha0=X0(1);
beta0=slope*alpha0+intercept;

x1=X0(2);
y1=X0(3);
x2=X0(4);
y2=X0(5);
x3=X0(6);
y3=X0(7);

x0=[alpha0,x1,y1,x2,y2,x3,y3];
g=@(x)CHI2spec(x,spec,slope,intercept);
options = optimset('MaxIter',50,'Display','off');
xf=fminsearch(g,x0,options);
mz=[xf g(xf)];