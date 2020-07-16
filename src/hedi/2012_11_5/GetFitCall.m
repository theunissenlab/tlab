function sxpa=GetFitCall(hashname,spec,t,f,x0,slope,intercept,thr);
% try to fit a whole call spec with time vector t and frequency f

tmax=t(end);
n=length(t);

xna=[];
for tt=fix(n/2):1:n;
    sp=spec(:,tt);
    xpa=fitProc(sp,x0,xna,f,tt,slope,intercept);
    xna=[xna;tt xpa];
    
end;
for tt=fix(n/2)-1:-1:1;
    sp=spec(:,tt);
    xpa=fitProc(sp,x0,xna,f,tt,slope,intercept);
    xna=[xna;tt xpa];
end;

sxpa=sortrows(xna,1);

save(hashname,'sxpa');
