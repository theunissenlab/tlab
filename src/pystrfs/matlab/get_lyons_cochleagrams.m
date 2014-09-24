function cgrams = get_lyons_cochleagrams()

    freqsFile = '~/netman/data/lyons_freqs.csv';    
    freqs = csvread(freqsFile);
    
    preprocFile = '/auto/k6/mschachter/pystrfs/preproc/lyons.agc_1.earQ_8.step_0.25.h5';
    
    
    h5 = h5utils();
    f = h5.open(preprocFile);
    md5s = h5.get_subgroups(f, '/');
    
    cgrams = cell(length(md5s), 1);
    
    for k = 1:length(md5s)
        
        s = struct();        
        s.spec = h5.get_ds(f, sprintf('%s/spectrogram', md5s{k}));
        s.freqs = freqs;
        s.sample_rate = 1000;
        s.wav_file = sprintf('/auto/k6/mschachter/pystrfs/all_stims/%s.wav', md5s{k});
        
        cgrams{k} = s;
    end
    
    h5.close(f);
    
    


