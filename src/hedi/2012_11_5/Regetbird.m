function Regetbird(sex,individual,a);

data_path='../data';
addpath('../sound');
for i=1:16;
    fname=sprintf('%s/%cbird%dcall%d.wav',data_path,sex,individual,i);
    if exist(fname,'file');
        ReProcessFileFit(fname,a,20);
    end;
end;

