function launch_stress_birds(sbird,data_dir)

name=sbird.name;
start=sbird.start;

for i=1:12;
    if i<10;
        c=sprintf('0%d',i);
    else
        c=sprintf('%d',i);
    end;
    wavname=sprintf('%s/%s_%s.wav',data_dir,name,c);
    fprintf('processing %s',wavname);
    ProcessFileStressFit(wavname,start);
    
end;