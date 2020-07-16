function [fi,gi,freq]=getImageFilter(ttime,filters)

n=length(ttime);
tau=5;
freq=0:50:10000;
fi=[];
gi=[];
fx1=lowP(filters(:,1),ttime,tau);
fp1=lowP(filters(:,2),ttime,tau);
fx2=lowP(filters(:,3),ttime,tau);
fp2=lowP(filters(:,4),ttime,tau);
fx3=lowP(filters(:,5),ttime,tau);
fp3=lowP(filters(:,6),ttime,tau);

for i=1:n;
    x1=fx1(i);y1=fp1(i);x2=fx2(i);y2=fp2(i);x3=fx3(i);y3=fp3(i);
    filter=abs(interp1([-10 x1 x2 x3 10000],[0 y1 y2 y3 0],freq,'cubic'));
    me=sum(filter.*freq)/sum(filter);
    vs=sum(filter.*(freq-me).^2)/sum(filter);
    fi=[fi;i,me,sqrt(vs)];
    gi=[gi;filter];
end;

    