function figure_lars_strfs()

   dataRootDir = '/auto/fdata/mschachter/data';
   
   cells = {'yy1617_4_A','gg0116_1_A','oo0909_2_A', 'pupu1203_4_A', 'pipi1112_1_A',...
             'pupi0414_10','obla1305_7','gr0869_5', 'glb5202_9', 'gr1070_1',...
             'pupu2122_2_B','yy1617_5_B','ww1211_5_B', 'yy2728_4_B', 'blabla1515_2_B'};
    
   fname = 'strflab.stft.LARS.best.nlType_linear.mat';
   
   
   for k = 1:length(cells)
       cellName = cells{k};
       
       cellDir = fullfile(dataRootDir, cellName);
       outputDir = fullfile(cellDir, 'conspecific', 'output');
       fpath = fullfile(outputDir, fname);
       
       mvars = load(fpath);
       strf = squeeze(mvars.modelParamsTrained.w1);
       clear mvars;
       
       mval = max(strf(:));
       figure; hold on;
       imagesc(strf); axis tight;
       caxis([-mval mval]);
       title(sprintf('cell=%s', cellName));
   end
   
   