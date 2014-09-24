fband=25;
alpha0=-0.4;

out=[];
for freq=200:50:4000;
    f1=freq-fband;
    f2=freq+fband;
    sp=find(stf_fun(:,3)>f1 & stf_fun(:,3)<f2);
    x=stf_fun(sp,1);y=stf_fun(sp,2);  
    fit=polyfit(x,y,1);
    out=[out;freq,fit(1),fit(2)];
end;
n=length(200:50:4000);
lasta=1;
lastb=0;
lastc=-alpha0;
iso=[];
for i=1:n;
    a=out(i,2);
    c=out(i,3);
    
    %plot(-0.7:0.1:0.0,a*
    % compute intersection
    % beta - a*alpha -c=0
    % lastb*beta +lasta*alpha+lastc = 0
    %
    % lastb(a*alpha+c)+lasta*alpha + lastc=0
    % (lastb*a+lasta)*alpha+(lastb+lastc)*c=0
    
    % new_alpha=-(lastb*c+lastc)/(lastb*a+lasta);
    
    new_alpha =-(lastb*c+lastc)/(lastb*a+lasta);
    new_beta  =a*new_alpha+c;
    
    iso=[iso; new_alpha,new_beta]
    % take the orthogonal of (a,1) is (-1,a) 
    
    lastb=a
    lasta=1
    lastc=-lastb*new_beta-lasta*new_alpha
end;

fit2=polyfit(out(:,1),out(:,2),1);
