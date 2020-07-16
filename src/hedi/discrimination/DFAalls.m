function error=DFAalls(z_data,group,nbird,nPC,nDF)

[nc,np]=size(z_data);

[Coeff, Score, latent] = princomp(z_data);


% train data  
[d, p, stats] = manova1(Score(:, 1:nPC),group);
[mean_grp, sem_grp, meanCI_grp, range_grp] = grpstats(stats.canon(:,1:nDF),group, {'mean', 'sem', 'meanci', 'range'});

% valid dat set 
vc = Score(:,1:nPC);
vCanon = vc*stats.eigenvec(:, 1:nDF);

Mah = zeros(nc,nbird); % 
for i = 1:nc
    for j = 1:nbird
        Mah(i,j) = sqrt((vCanon(i,:)-mean_grp(j,:))*(vCanon(i,:)-mean_grp(j,:))');
    end
end
%Mah
ConfMat=zeros(nbird,nbird);
 
for i = 1:nc;
    ic=group(i);
    [row, col] = find(Mah(i,:)==min(Mah(i,:)));
    ConfMat(ic,col) = ConfMat(ic,col)+1; % classification incorrecte
end
%ConfMat
ntotal=sum(sum(ConfMat));
zMat=ConfMat;
zMat(1:(nbird+1):end)=0;
nout=sum(sum(zMat));
error=nout/ntotal*100;