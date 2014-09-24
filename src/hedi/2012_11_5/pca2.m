% get all vocal stuff
close all;
clear all;
data_dir='../output';
nbird=11;

max_len=get_longest(nbird);
alpha_min=-0.40;
big_data2=[];
nbird=14;
ncall=16;
Group=zeros(nbird*ncall,1);
for bird=[1,3:9,11:16];
    
    ldir=sprintf('%s/res.*v4*mbird%dcall*.mat',data_dir,bird);
    d=dir(ldir);
    n=length(d);
    if n~=ncall;
        disp('problem\n');
        disp(bird)
    end;
    for i=1:n;
        fname=sprintf('%s/%s',data_dir,d(i).name);
        load(fname)       
        y=lowP(res.data(:,2),res.data(:,1),5);       
        if length(y)<max_len;
              y=[y;alpha_min*ones(max_len-length(y),1)];
        end;
        fi=getFilter(res.data(:,1),res.data(:,3:end));        
        yf=lowP(fi(:,1),res.data(:,1),20);
        if length(yf)<max_len;
              yf=[yf;,yf(end)*ones(max_len-length(yf),1)];
          end;
        big_data2=[big_data2;y(1:10:end)',yf(1:10:end)'];
        offset=(bird-1)*ncall + i;
        Group(offset)=bird;        
    end;
end;
    
nPC=10;
nDF=5;

[Coeff, Score, latent] = princomp(zscore(big_data2));
[d, p, stats] = manokva1(Score(:, 1:nPC),Group);
[mean_grp, sem_grp, meanCI_grp, range_grp] = grpstats(stats.canon(:,1:nDF),Group, {'mean', 'sem', 'meanci', 'range'});
%PC=Coeff(:, 1:nPC) * stats.eigenvec(:, 1:nDF);
