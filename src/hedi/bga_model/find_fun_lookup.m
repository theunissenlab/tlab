fband=25;
alpha=0.4;

out=[];
for freq=400:100:3500;
    f1=freq-fband;
    f2=freq+fband;
    sp=find(stf_fun(:,3)>f1 & stf_fun(:,3)<f2);
    x=stf_fun(sp,1);y=stf_fun(sp,2);  
    fit=polyfit(x,y,1);
    
    beta=fit(1)*alpha+fit(2);
    bprime=fit(1)*beta-alpha
    next_alpha=
    yfit=fit(1)*x+fit(2);
    yresid=y-yfit;
    SSresid=sum(yresid.^2);
    SStotal=(length(y)-1)*var(y);
    R2=1-SSresid/SStotal;
    out=[out;freq,fit(1),fit(2),R2]
end;