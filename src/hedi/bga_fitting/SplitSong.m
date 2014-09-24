function splits=SplitSong(tspec,dB)
spec=tspec.spec;
env=[sum(spec,1),0];
env=env/max(env)
dfs=zeros(1,length(tspec.t)+1);
spl=find(env>dB);
dfs(spl)=1;
ddfs=find(abs(diff(dfs))>0.5);
uddfs=diff(dfs);
splits={};
t=1;
i=1;
while t<length(ddfs);
    %t
    %ddfs(t)
    %uddfs(t)
    if uddfs(ddfs(t))>0;
        splits{i}={};
        start=ddfs(t);
        finish=ddfs(t+1);
        splits{i}.start=start;
        splits{i}.finish=finish;
        splits{i}.spec=spec(:,start:finish);
        splits{i}.t=tspec.t(start:finish);
        splits{i}.f=tspec.f;
        i=i+1;
        t=t+2;
    else
%       fprintf('not starting with calls - error\n');
   end
end
    
