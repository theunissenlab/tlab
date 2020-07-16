function mf=findFun(x);
% find the fundamentals using bullshit correlation
tfrep=timefreq(x,44100,'stft');
z=sum(tfrep.spec,2);
z=z-mean(z);
n=length(tfrep.f);
My=-1e6;
mf=-1;
[m,fm]=max(z);
if fm==1;
    [m,fm]=max(z(3:end));
end;

for mu=1:n;
    f1=floor(fm/mu);
    y1=sum(z(f1:f1:n));
    f2=floor(fm*mu);
    y2=sum(z(f1:f1:n));
    if y1> My;
        My=y1;
        mf=f1/n*tfrep.f(end);
    end;
    if y2> My;
        My=y2;
        mf=f2/n*tfrep.f(end);
    end;
end;
    
% for f=3:n;
%     
%     if y>My;
%         My=y;
%         mf=f/n*tfrep.f(end);       
%     end;
% end;
% res=[];
% for f=3:60;
%     y=computeCOf(sum(z),f);
%     res=[res;f y];
% end
% [m,fm]=max(res(:,2));
