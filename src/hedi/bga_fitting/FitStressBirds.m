function FitStressBirds(sbird,data_dir,freq_a)

name=sbird;


for i=1:12;
    if i<10;
        c=sprintf('0%d',i);
    else
        c=sprintf('%d',i);
    end;
    wavname=sprintf('%s/%s_%s.wav',data_dir,name,c);
    fprintf('processing %s\n',wavname);
    spl=fitWave(wavname,.5,freq_a);
    savename=sprintf('fit.%s_%s.mat',name,c);
    save(savename,'spl');
end;
