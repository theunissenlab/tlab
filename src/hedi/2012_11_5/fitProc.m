function xpa=fitProc(sp,x0,xna,f,t,slope,intercept)

[pmax,imax]=max(sp);
sp=sp/pmax;
ffalpha=@(x)findAlphaBetaX0(f,sp,x,slope,intercept);

m0=ffalpha(x0);
if ~isempty(xna);
    X=xna(end,2:(end-1));
    mz=ffalpha(X);
    [mu,imu]=min(xna(:,end));
    xm=xna(imu,2:(end-1));
    mM=ffalpha(xm);
else
    mM=m0;
    mz=m0;
end;

en=[mz;m0;mM];
[k,ik]=min(en(:,end));
xpa=en(ik,:);
fprintf('%f %f ',t,xpa(end));
for k=1:length(xpa(1:(end-1)));
    fprintf(1,'%f ',xpa(k));
end;
fprintf('\n');
