function figure_lin_comps(type)

    if nargin < 1
        type = 'linear';
    end

    dataDir = '/auto/fdata/mschachter/data';
    stimsDir = '/auto/fdata/mschachter/data/all_stims';


    cells = {'yy1617_4_A','gg0116_1_A','oo0909_2_A', 'pupu1203_4_A', 'pipi1112_1_A',...
             'pupi0414_10','obla1305_7','gr0869_5', 'glb5202_9', 'gr1070_1',...
             'pupu2122_2_B','yy1617_5_B','ww1211_5_B', 'yy2728_4_B', 'blabla1515_2_B'};
    
    clen = length(cells);
    
    types = {'df', 'tg', 'lars'};
    
    infoTrain = struct;
    infoBoundTrain = zeros(1, clen);
    infoValid = struct;
    infoBoundValid = zeros(1, clen);
    
    for k = 1:length(types)       
        tname = types{k};
        infoTrain.(tname) = zeros(1, clen);
        infoValid.(tname) = zeros(1, clen);
    end
    
    
    for k = 1:length(cells)        
        c = cells{k}
        mdata = get_method_data(dataDir, stimsDir, c);
        
        for m = 1:length(types)           
            tname = types{m};
            infoTrain.(tname)(k) = mdata.(tname).(type).infoTrain;
            infoValid.(tname)(k) = mdata.(tname).(type).infoValid;
            infoBoundTrain(k) = mdata.infobndTrain;
            infoBoundValid(k) = mdata.infobndValid;
        end        
    end
    
    
    dfPerfs = infoValid.df ./ infoBoundValid;
    tgPerfs = infoValid.tg ./ infoBoundValid;
    larsPerfs = infoValid.lars ./ infoBoundValid;
    
    h = figure; hold on;      
    plot(infoBoundValid, dfPerfs, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    plot(infoBoundValid, tgPerfs, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
    plot(infoBoundValid, larsPerfs, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
    axis([30 235 0 1]);
    title(sprintf('NL=%s', type));
    xlabel('Info Bound');
    ylabel('Performance Ratio');
    legend('DF', 'TG', 'LARS');
    
    
    dfMean = mean(dfPerfs);
    dfDev = std(dfPerfs);
    
    tgMean = mean(tgPerfs);
    tgDev = std(tgPerfs);
    
    larsMean = mean(larsPerfs);
    larsDev = std(larsPerfs);
    
    fprintf('  DF: avg perf=%f +/- %f\n', dfMean, dfDev);
    fprintf('  TG: avg perf=%f +/- %f\n', tgMean, tgDev);
    fprintf('LARS: avg perf=%f +/- %f\n', larsMean, larsDev);
    
    outputFile = '~/Desktop/sfn_poster/images/figure_lin_comps.svg';
    plot2svg_2d(outputFile, h);
    
    
    
    
    