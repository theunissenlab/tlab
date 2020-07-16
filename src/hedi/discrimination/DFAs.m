function error=DFAs(z_data,group,nbird,percent,nPC,nDF)

[nc,np]=size(z_data);
ntrain=floor(nc*percent);
nvalid=nc-ntrain;
[Coeff, Score, latent] = princomp(z_data);

done=0;
igroup=unique(group);
while done==0
    shuffle=randperm(nc);
    itrain=shuffle(1:ntrain);
    ivalid=shuffle((ntrain+1):end);
    gtrain=group(itrain);
    gvalid=group(ivalid);
    if length(unique(gtrain))==length(unique(gtrain))
        done=1;
    end
end;

% train data set 
[d, p, stats] = manova1(Score(itrain, 1:nPC),gtrain);
[mean_grp, sem_grp, meanCI_grp, range_grp] = grpstats(stats.canon(:,1:nDF),gtrain, {'mean', 'sem', 'meanci', 'range'});

% valid dat set 
vc = Score(ivalid,1:nPC);
vCanon = vc*stats.eigenvec(:, 1:nDF);

Mah = zeros(nvalid,nbird); % 
for i = 1:nvalid
    for j = 1:nbird
        Mah(i,j) = sqrt((vCanon(i,:)-mean_grp(j,:))*(vCanon(i,:)-mean_grp(j,:))');
    end
end
%Mah
ConfMat=zeros(nbird,nbird);
 
for i = 1:nvalid;
    ix=ivalid(i);
    ic=find(igroup==group(ix));
    [row, col] = find(Mah(i,:)==min(Mah(i,:)));
    ConfMat(ic,col) = ConfMat(ic,col)+1; % classification incorrecte
end
%ConfMat
ntotal=sum(sum(ConfMat));
zMat=ConfMat;
zMat(1:(nbird+1):end)=0;
nout=sum(sum(zMat));
error=nout/ntotal*100;