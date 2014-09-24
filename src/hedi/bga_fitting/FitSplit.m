function [best_fit,archives]=FitSplit(split,freq_a,nsplit,tsplit)
ntime=length(split.t);
freqs=split.f;

guess=[];
best_fit=[];
archives=[];
stime=floor(ntime/2);
for n=stime:ntime;
    fprintf('fitting %d over %d time: %d - %.2f done\n',nsplit,tsplit,n,-(stime-n)/ntime*100);
    sp=split.spec(:,n);
    res=ProceedToFit(sp,freqs,freq_a,guess);
    [k,p]=size(res);
    [mi,im]=min(res(:,end));
    best_fit=[best_fit;n res(im,:)];
    guess=res(im,:);
    archives=[archives; [n*ones(k,1),res]];
    fprintf('best %.4f %.2f %.2f\n',res(im,1),res(im,2),res(im,end));
end
guess=best_fit(best_fit(:,1)==ntime/2,2:end);
for n=(stime-1):-1:1;
    fprintf('fitting %d over %d time: %d - %.2f done\n',nsplit,tsplit,n,50+(stime-n)/ntime*100);
    sp=split.spec(:,n);
    res=ProceedToFit(sp,freqs,freq_a,guess);
    [k,p]=size(res);
    [mi,im]=min(res(:,end));
    best_fit=[best_fit;n res(im,:)];
    guess=res(im,:);
    archives=[archives; [n*ones(k,1),sortrows(res,p)]];
    fprintf('best %.4f %.2f %.2f\n',res(im,1),res(im,2),res(im,end));
end;

best_fit=sortrows(best_fit,1);
archives=sortrows(archives,1);

