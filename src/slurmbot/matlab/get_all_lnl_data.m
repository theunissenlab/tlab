function alldata = get_all_lnl_data()

    %allcells = {'yy1617_4_A','gg0116_1_A','oo0909_2_A','pupu1203_4_A','pipi1112_1_A','blublu0916_5_A','pipi1112_2_A','blabla1904_1_A','gg0318_2_A','lghp1616_9_A','pupi0414_10','obla1305_7','gr0869_5','glb5202_9','gr1070_1','glb5656_2','pupi0414_8','gr1070_4','obla17102_2','blahp4108_4','pupu2122_2_B','yy1617_5_B','ww1211_5_B','yy2728_4_B','blabla1515_2_B','yy0203_3_A','blabla1515_7_B','pupu2122_3_B','oo2015_5_B','ww1211_3_B'};
    allcells = {'yy1617_4_A','gg0116_1_A','oo0909_2_A','pupu1203_4_A','blublu0916_5_A','pipi1112_2_A','blabla1904_1_A','gg0318_2_A','lghp1616_9_A','pupi0414_10','obla1305_7','gr0869_5','glb5202_9','gr1070_1','glb5656_2','pupi0414_8','gr1070_4','obla17102_2','blahp4108_4','pupu2122_2_B','yy1617_5_B','ww1211_5_B','yy2728_4_B','blabla1515_2_B','yy0203_3_A','blabla1515_7_B','pupu2122_3_B','oo2015_5_B','ww1211_3_B'};
    
    
    alldata = cell(length(allcells));
    
    for k = 1:length(allcells)
        alldata{k} = lnl_comp_onecell(allcells{k}, {'stft', 'lyons', 'surprise'}, {'conspecific', 'flatrip'});                
    end
    
    