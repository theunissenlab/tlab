function max_len=get_longest(nbird)
max_len=-1;
data_dir='../output';
for bird=[1:nbird];
    ldir=sprintf('%s/res.*v4*mbird%dcall*.mat',data_dir,bird);
    d=dir(ldir);
    n=length(d);
    for i=1:n;
        fname=sprintf('%s/%s',data_dir,d(i).name);
        load(fname)
        if max_len<length(res.data(:,1));
            max_len=length(res.data(:,1));           
        end;
    end;
end;
