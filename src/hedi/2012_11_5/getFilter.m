function [fi,gi]=getFilter(ttime,filters)

n=length(ttime);

freq=200:8000;
fi=[];
gi=[];
for i=1:n;
    f=filters(i,:);
    x1=f(1);y1=f(2);x2=f(3);y2=f(4);x3=f(5);y3=f(6);
    filter=abs(interp1([-10 x1 x2 x3 10000],[0 y1 y2 y3 0],freq,'cubic'));
    me=sum(filter.*freq)/sum(filter);
    vs=sum(filter.*(freq-me).^2)/sum(filter);
    fi=[fi;me,sqrt(vs)];
    gi=[gi;filter];
end;
    