

data={}
data.sex='all';
msex=['m','f'];

len=[0,0,0];

n=length(birds);
for sex=1:n;
    bs=birds(sex).birds
    nbirds=length(bs);
    for bird=1:nbirds;
        nbirds
        bsc=bs(bird).calls
        ncall=length(bsc);
        for call=1:ncall;
            bsc(call)
            syr=birds(sex).birds(bird).calls(call).syr;
            if len(sex)<length(syr);
                len(sex)=length(syr);
            end;
            if len(3)<length(syr);
                len(3)=length(syr);
            end;
        end;
    end;
end;

big_data=[];
datam = {};
datam.sex='m';
datam.data=[];
datam.group=[];

dataf = {};
dataf.sex='f';
dataf.data=[];
dataf.group=[];

dataa = {};
dataa.sex='all';
dataa.data=[];
dataa.group=[];

big_data=[datam,dataf,dataa];

for sex=1:n;
    bs=birds(sex).birds
    nbirds=length(bs);
    for bird=1:nbirds;
        nbirds
        bsc=bs(bird).calls
        ncall=length(bsc);
        for call=1:ncall;
            bsc(call)
            syra=birds(sex).birds(bird).calls(call).syr;
            fila=birds(sex).birds(bird).calls(call).fil;
            syrs=birds(sex).birds(bird).calls(call).syr;
            fils=birds(sex).birds(bird).calls(call).fil;
            size(syrs)
            if length(syrs)<len(sex);
                syrs=[syrs',syrs(end)+zeros(1,len(sex)-length(syrs))];
                fils=[fils',fils(end)+zeros(1,len(sex)-length(fils))];
            else
                syrs=syrs';
                fils=fils';
            end;
            if length(syra)<len(3);
                syra=[syra',syra(end)+zeros(1,len(3)-length(syra))];
                fila=[fila',fila(end)+zeros(1,len(3)-length(fila))];
            else 
                syra=syra';
                fila=fila';
            end;           
            big_data(sex).data=[big_data(sex).data;[syrs,fils]];
            big_data(3).data=[big_data(3).data;[syra,fila]]; 
            big_data(sex).group=[big_data(sex).group;bird];
            big_data(3).group=[big_data(3).group;100*sex+bird];
        end;
    end;
end;