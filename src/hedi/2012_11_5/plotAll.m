

n=length(birds);
for sex=1:n;
    sex
    bs=birds(sex).birds
    nbirds=length(bs);
    for bird=1:nbirds;
        nbirds
        bsc=bs(bird).calls
        ncall=length(bsc);
        for call=1:ncall;
            call
            plotGuy(birds,sex,bird,call);
        end;
    end;
end;
