function fsave(sex,bird,call)

fname=sprintf('../output/pAv3.%cbird%dcall%d.mat',sex,bird,call);
load(fname);
outname=sprintf('../output_txt/%cbird%dcall%d',sex,bird,call);
fout=fopen(outname,'w');

[n,m]=size(sxpa)
for i=1:n;
    for j=1:m;
        fprintf(fout,'%f ',sxpa(i,j));
    end;
    fprintf(fout,'\n');
end;
fclose(fout);

