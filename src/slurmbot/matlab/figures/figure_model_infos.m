function figure_model_infos()

    dataDir = '/auto/fdata/mschachter/data';
    stimsDir = '/auto/fdata/mschachter/data/all_stims';
    
    mldCells = {'yy1617_4_A','gg0116_1_A','oo0909_2_A', 'pupu1203_4_A', 'pipi1112_1_A'};

    ovCells = {'pupi0414_10','obla1305_7','gr0869_5', 'glb5202_9', 'gr1070_1'};
    
    lCells = {'pupu2122_2_B','yy1617_5_B','ww1211_5_B', 'yy2728_4_B', 'blabla1515_2_B'};
    
    
    cellsByRegion = {mldCells, ovCells, lCells};
    regNames = {'MLd', 'OV', 'Field L'};
    
    
    for m = 1:length(regNames)
    
        carr = cellsByRegion{m};

        infoBndVals = zeros(1, length(carr));
        dfLinInfos = zeros(1, length(carr));
        dfExpInfos = zeros(1, length(carr));
        tgLinInfos = zeros(1, length(carr));
        tgExpInfos = zeros(1, length(carr));
        larsLinInfos = zeros(1, length(carr));
        larsExpInfos = zeros(1, length(carr));

        for k = 1:length(carr)               
            cellName = carr{k};
            mdata = get_method_data(dataDir, stimsDir, cellName);

            infoBndVals(k) = mdata.infobndValid;

            dfLinInfos(k) = mdata.df.linear.infoValid;
            dfExpInfos(k) = mdata.df.exponential.infoValid;

            tgLinInfos(k) = mdata.tg.linear.infoValid;
            tgExpInfos(k) = mdata.tg.exponential.infoValid;

            larsLinInfos(k) = mdata.lars.linear.infoValid;
            larsExpInfos(k) = mdata.lars.exponential.infoValid;

        end

        hline = [0 infoBndVals];

        figure; hold on;

        plot(hline, hline, 'k-');

        plot(infoBndVals, dfLinInfos, 'kx');
        plot(infoBndVals, tgLinInfos, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
        plot(infoBndVals, larsLinInfos, 'k*', 'MarkerSize', 12);

        plot(infoBndVals, dfExpInfos, 'rx');
        plot(infoBndVals, tgExpInfos, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
        plot(infoBndVals, larsExpInfos, 'r*', 'MarkerSize', 12);
        legend('Bound', 'DF (lin)', 'TG (lin)', 'LARS (lin)', 'DF (exp)', 'TG (exp)', 'LARS (exp)');
        title(sprintf('Validation: %s', regNames{m}));
        xlabel('Info Upper Bound');
        ylabel('Info Captured');
        axis tight;
        
        for k = 1:length(carr)               
            
            xval = infoBndVals(k) - mean(infoBndVals)*2e-1;
            yval = infoBndVals(k) - mean(infoBndVals)*1e-1;
            text(xval, yval, strrep(carr{k}, '_', '\_'));
        end
        
    end
    
    