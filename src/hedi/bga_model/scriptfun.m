alpha=unique(stf_fun(:,1));
beta=unique(stf_fun(:,2));
n=length(alpha);
m=length(beta);
fun=200+zeros(n,m);
har=zeros(n,m);

for i=1:n;
    alpha_row=stf_fun(stf_fun(:,1)==alpha(i),:);
    k=length(alpha_row(:,2));
    fun(i,(end-k+1):end)=alpha_row(:,3)';
    har(i,(end-k+1):end)=(sum(alpha_row(:,6:2:end),2)'-1)/10;
end
subplot(2,1,1)
imagesc(alpha,beta,fun);
colorbar
subplot(2,1,2)
imagesc(alpha,beta,har);
colorbar