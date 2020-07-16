function y=lowP(a,ttime,tau);
input=[ttime,a];
tspan=ttime;
x0=a(1);
g=@(t,x)dlowP(t,x,input,tau);

[t,y]=ode45(g,tspan,x0);
