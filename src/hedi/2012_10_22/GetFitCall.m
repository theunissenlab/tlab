function sxpa=GetFitCall(hashname,spec,t,f,X0);
% try to fit a whole call spec with time vector t and frequency f

tmax=t(end);
n=length(t);

xpa=[];


X=X0;
for tt=fix(n/2):1:n;
    sp=spec(:,tt);
    [pmax,imax]=max(sp);
    sp=sp/pmax;   
    mz=findAlphaBetaX0v2(f,sp,X);
    xpa=[xpa;tt mz];
    X=mz(1:8);
    fprintf(1,'%f %f %f %f %f %f %f %f %f\n',tt,mz(1),mz(2),mz(3),mz(4),mz(5),mz(6),mz(7),mz(8));
end;

X=mz(1,1:8); %fix(n/2) value

for tt=fix(n/2)-1:-1:1;      
    sp=spec(:,tt);
    [pmax,imax]=max(sp);
    sp=sp/pmax;   
    mz=findAlphaBetaX0v2(f,sp,X);
    xpa=[xpa;tt mz];
    X=mz(1:8);
    fprintf(1,'%f %f %f %f %f %f %f %f %f\n',tt,mz(1),mz(2),mz(3),mz(4),mz(5),mz(6),mz(7),mz(8));
end;

sxpa=sortrows(xpa,1);

save(hashname,'sxpa');
