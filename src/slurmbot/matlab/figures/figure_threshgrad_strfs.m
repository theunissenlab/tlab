function figure_threshgrad_strfs(cellName)

    dataRootDir = '/auto/fdata/mschachter/data';
    cellDir = fullfile(dataRootDir, cellName);
    outputDir = fullfile(cellDir, 'conspecific', 'output');
    
   tvals = {'0.00', '0.25', '0.50', '0.75', '1.00'};
   
   tstrfs = cell(1, length(tvals));
   
   basename = 'strflab.stft.threshgrad.nlType_linear.thresh_%s.JN_%d.mat';
   
   njns = 3;
   
   for k = 1:length(tvals)
       strfSum = zeros(60, 75);
       for jn = 1:njns
           fname = sprintf(basename, tvals{k}, jn-1);
           fpath = fullfile(outputDir, fname);
           mvars = load(fpath);
           strfSum = strfSum + squeeze(mvars.modelParamsTrained.w1);
           clear mvars;
       end
       
       strfAvg = strfSum / 3;
       tstrfs{k} = strfAvg;       
   end
   
   
   for k = 1:length(tvals)
       strf = tstrfs{k};
       mval = max(strf(:));
       figure; hold on;
       imagesc(strf); axis tight;
       caxis([-mval mval]);
       title(sprintf('thresh=%s', tvals{k}));
   end
   
   