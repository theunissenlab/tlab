function figure_strf_responses()

   dataRootDir = '/auto/fdata/mschachter/data';
   stimsDir = fullfile(dataRootDir, 'all_stims');
   
   cells = {'yy1617_4_A','gg0116_1_A','oo0909_2_A', 'pupu1203_4_A', 'pipi1112_1_A',...
             'pupi0414_10','obla1305_7','gr0869_5', 'glb5202_9', 'gr1070_1',...
             'pupu2122_2_B','yy1617_5_B','ww1211_5_B', 'yy2728_4_B', 'blabla1515_2_B'};
    
   types = {'df', 'tg', 'lars'};      
         
   for k = 1:length(cells)
       cellName = cells{k};
       mdata = get_method_data(dataRootDir, stimsDir, cellName);
     
       vg = mdata.vGroups;
       vindx = findIdx(vg, mdata.groupIndex);
       respReal = mdata.resp(vindx);
       
       dfResp = [];
       dfStrf = [];
       
       figure('Name', cellName); hold on;
       
       for m = 1:length(types)
           
            t = types{m};
            
            spl = (m-1)*3 + 1;
            spr = spl + 1;
            
            s = mdata.(t).linear.strf;
            smax = max(abs(s(:)));
            
            l1 = sum(abs(s(:)));
            l2 = sum(abs(s(:).^2));
            vperf = mdata.(t).linear.perf;
            
            respValid = mdata.(t).linear.respValid;
            
            if m == 1
                dfResp = respValid;
                dfStrf = s;
            end
            
            subplot(3, 3, spl); hold on;
            imagesc(s); axis tight;
            t1 = sprintf('l1=%f | l2=%f | perf=%f', l1, l2, vperf);
            title(t1);
            caxis([-smax smax]);
            
            subplot(3, 3, spr); hold on;
            plot(respReal, 'k-');
            plot(respValid, 'r-');
            axis tight;
            
            if m > 1
                
               %rdiff = dfResp - respValid;
               sdiff = dfStrf - s;
               sdmax = max(abs(sdiff(:)));
               subplot(3, 3, spr+1);
               imagesc(flipud(sdiff)); axis tight;
               caxis([-sdmax sdmax]);
               nsdiff = norm(sdiff(:));
               title(sprintf('Diff | norm=%f', nsdiff));
                
            end
            
       end
       
   end
   