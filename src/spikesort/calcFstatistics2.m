function ev = calcFstatistics2(r)

repnum=size(r,2);
n=size(r,2);
K=size(r,1);
N=n*K;
% for ii=1:size(r,3)
%      progressdot(ii,1000,20000,size(r,3));
%      d=squeeze(r(:,:,ii));
%      dm=mean(d,2);
%      if isempty(dm)
% 	     continue;
%      end
%      dmm=mean(d(:));
%      bgv=n*sum((dm-dmm).^2)/(K-1);
%      wgv=bsxfun(@minus,d,dm).^2;
%      wgv=sum(wgv(:))/(N-K);
%      ev(ii)=bgv/wgv;
% end

sr=nanstd(r);
sr(sr<1e-3)=NaN; % cut out non-valid data on the moving edge
sr=nanmean(sr,2);

dm=nanmean(r,2);
dmm=nanmean(nanmean(r,1),2);
bgv=n*sum(bsxfun(@minus,dm,dmm).^2)/(K-1);
wgv=bsxfun(@minus,r,dm).^2;
wgv=nansum(nansum(wgv,1),2)/(N-K);
ev=bgv./wgv;
ev=ev(:);

ev(isnan(sr))=0;
