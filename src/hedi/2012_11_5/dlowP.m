function dy=dlowP(t,x,input,tau)
in=interp1(input(:,1),input(:,2),t);
y=x(1);
dy=[(in - y)/tau];
