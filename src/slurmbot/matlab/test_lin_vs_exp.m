function test_lin_vs_exp()

    dataDir = '/auto/fdata/mschachter/data';
    
    mldCells = {'yy1617_4_A','gg0116_1_A','oo0909_2_A', 'pupu1203_4_A', 'pipi1112_1_A'};

    ovCells = {'pupi0414_10','obla1305_7','gr0869_5', 'glb5202_9', 'gr1070_1'};
    
    lCells = {'pupu2122_2_B','yy1617_5_B','ww1211_5_B', 'yy2728_4_B', 'blabla1515_2_B'};
    
    
    display = 1;
    
    mldLinScores = [];
    mldExpScores = [];
    for k = 1:length(mldCells)        
        cellName = mldCells{k};
        [linScores, expScores] = find_best_method(dataDir, cellName, display);
        
        mldLinScores = [mldLinScores linScores];
        mldExpScores = [mldExpScores expScores];                
    end
    
    
    ovLinScores = [];
    ovExpScores = [];
    for k = 1:length(ovCells)        
        cellName = ovCells{k};
        [linScores, expScores] = find_best_method(dataDir, cellName, display);
        
        ovLinScores = [ovLinScores linScores];
        ovExpScores = [ovExpScores expScores];                
    end
    
    
    lLinScores = [];
    lExpScores = [];
    for k = 1:length(lCells)        
        cellName = lCells{k};
        [linScores, expScores] = find_best_method(dataDir, cellName, display);
        
        lLinScores = [lLinScores linScores];
        lExpScores = [lExpScores expScores];                
    end
    
    alpha = 0.05;
    
    fprintf('t-testing MLd scores...\n');    
    [h, p, ci] = ttest(mldLinScores, mldExpScores, alpha, 'both')
    
    
    fprintf('t-testing OV scores...\n');    
    [h, p, ci] = ttest(ovLinScores, ovExpScores, alpha, 'both')
    
    
    fprintf('t-testing Field L scores...\n');    
    [h, p, ci] = ttest(lLinScores, lExpScores, alpha, 'both')
    
    
    


