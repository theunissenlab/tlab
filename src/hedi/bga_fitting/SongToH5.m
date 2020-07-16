function SongToH5(song, f, groupName)

    h5 = h5utils();    
    %h5.create_group(f, groupName);
    
    slope=-0.6025;
    intercept=-0.5979;
    md5 = song.wavname;
    
    songGroup = sprintf('%s/%s', groupName, md5);

    nsplits = length(song.r_moto);
    
    h5.set_attr(f, songGroup, 'slope', slope);
    h5.set_attr(f, songGroup, 'intercept', intercept);
        
    for k = 1:nsplits
        
        splitGroup = sprintf('%s/split%d', songGroup, k);
        split = song.splits{k};
        mgram = song.r_moto{k}.b;
        
        
        h5.set_attr(f, splitGroup, 'start_time', split.start);
        h5.set_attr(f, splitGroup, 'end_time', split.finish);
        
        h5.set_ds(f, splitGroup, 'motogram', mgram(:, [1:5, 7]));                
        h5.set_ds(f, splitGroup, 'spectrogram', split.spec);
        h5.set_ds(f, splitGroup, 'frequencies', split.f);        
        
    end
    