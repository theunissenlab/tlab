function [pred,b,r2s,psi]=RidgePlot(y,data,k,Stims,constraints)

r2s={};
% first do the ridge regression
b=ridge(y,data,k,0);   
pred=b(1)+data*b(2:end);
% compute its R2
r2=GetR2(y,pred);
r2s.all=r2;
psi={};
%create the subsets
nc=length(constraints);
%nrow for plotting
nrow=ceil(nc/2);

% plot all in first
subplot(nrow,2,1);
plot(y,pred,'+',min(y):0.1:max(y),min(y):0.1:max(y));
xlabel('psth');
ylabel('prediction');

% then the rest
for i=1:nc;
    fi=find(Stims==constraints(i));
    ys=y(fi);
    ds=data(fi,:);
    bsi=ridge(ys,ds,k,0);
    psi{i}={};
    psi{i}.fi=fi;
    psi{i}.ys=ys;
    psi{i}.ds=ds;
    psi{i}.bsi=bsi;
    
    ps=pred(fi);
    r2=GetR2(ys,ps);
    st=sprintf('d%dm',constraints(i));
    r2s.(st)=r2;
    subplot(nrow,2,i+1);
    plot(ys,ps,'+',min(y):0.1:max(y),min(y):0.1:max(y));
    xlabel('psth');
    ylabel('prediction');
end
  