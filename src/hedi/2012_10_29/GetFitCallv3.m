function sxpa=GetFitCallv3(hashname,spec,t,f,X0);
% try to fit a whole call spec with time vector t and frequency f

tmax=t(end);
n=length(t);

xpa=[];


X=X0;
% sp=spec(:,fix(n/2));
% [pmax,imax]=max(sp);
% X(5)=f(imax);
% X(3)=f(imax)/10;
% X(7)=f(imax)*2;
for tt=fix(n/2):1:n;
    sp=spec(:,tt);
    [pmax,imax]=max(sp);
    sp=sp/pmax;   
    mz=findAlphaBetaX0v3(f,sp,X);
    mz2=findAlphaBetaX0v3(f,sp,X0);   
    [mu,imu]=min(mz(:,9));
    mz3=findAlphaBetaX0v3(f,sp,mz(imu,1:8));
   if mz(:,9)<mz2(:,9)
       if mz(:,9)<mz3(:,9);           
           xpa=[xpa;tt mz];
           X=mz(1:8);
           fprintf(1,'%f %f %f %f %f %f %f %f %f %f\n',mz(9),tt,mz(1),mz(2),mz(3),mz(4),mz(5),mz(6),mz(7),mz(8));
       else
           xpa=[xpa;tt mz3];
           X=mz3(1:8);
           fprintf(1,'%f %f %f %f %f %f %f %f %f %f\n',mz3(9),tt,mz3(1),mz3(2),mz3(3),mz3(4),mz3(5),mz3(6),mz3(7),mz3(8));
       end
   else
       if mz2(:,9)<mz3(:,9);
           xpa=[xpa;tt mz2];
           X=mz2(1:8);
           fprintf(1,'%f %f %f %f %f %f %f %f %f %f\n',mz2(9),tt,mz2(1),mz2(2),mz2(3),mz2(4),mz2(5),mz2(6),mz2(7),mz2(8));
       else
           xpa=[xpa;tt mz3];
           X=mz3(1:8);
           fprintf(1,'%f %f %f %f %f %f %f %f %f %f\n',mz3(9),tt,mz3(1),mz3(2),mz3(3),mz3(4),mz3(5),mz3(6),mz3(7),mz3(8));
       end;
   end
   
end;

X=mz(1,1:8); %fix(n/2) value

for tt=fix(n/2)-1:-1:1;      
    sp=spec(:,tt);
    [pmax,imax]=max(sp);
    sp=sp/pmax;   
    mz=findAlphaBetaX0v3(f,sp,X);
    mz2=findAlphaBetaX0v3(f,sp,X0);   
    [mu,imu]=min(mz(:,9));
    mz3=findAlphaBetaX0v3(f,sp,mz(imu,1:8));
   if mz(:,9)<mz2(:,9)
       if mz(:,9)<mz3(:,9);           
           xpa=[xpa;tt mz];
           X=mz(1:8);
           fprintf(1,'%f %f %f %f %f %f %f %f %f %f\n',mz(9),tt,mz(1),mz(2),mz(3),mz(4),mz(5),mz(6),mz(7),mz(8));
       else
           xpa=[xpa;tt mz3];
           X=mz3(1:8);
           fprintf(1,'%f %f %f %f %f %f %f %f %f %f\n',mz3(9),tt,mz3(1),mz3(2),mz3(3),mz3(4),mz3(5),mz3(6),mz3(7),mz3(8));
       end
   else
       if mz2(:,9)<mz3(:,9);
           xpa=[xpa;tt mz2];
           X=mz2(1:8);
           fprintf(1,'%f %f %f %f %f %f %f %f %f %f\n',mz2(9),tt,mz2(1),mz2(2),mz2(3),mz2(4),mz2(5),mz2(6),mz2(7),mz2(8));
       else
           xpa=[xpa;tt mz3];
           X=mz3(1:8);
           fprintf(1,'%f %f %f %f %f %f %f %f %f %f\n',mz3(9),tt,mz3(1),mz3(2),mz3(3),mz3(4),mz3(5),mz3(6),mz3(7),mz3(8));
       end;
   end
end;

sxpa=sortrows(xpa,1);

save(hashname,'sxpa');
