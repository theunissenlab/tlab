function [plotmin, plotmax] = findplotminmax(tmpv2)

plotokflag = 0;
cnt = 0;

plotmin = min(tmpv2(:));
plotmax = max(tmpv2(:));
nevents = size(tmpv2,1);

maxloop = min([100, fix(0.001*nevents)]);
stdEst = max(4*std(tmpv2,0,1));

fprintf(1, 'Finding min and max for plotting...\n');
while plotokflag == 0
    

    [maxidx, ~] = find(tmpv2 == plotmax);
    [minidx, ~] = find(tmpv2 == plotmin);

    if max([abs(plotmax), abs(plotmin)]) >  stdEst
        if ( abs(plotmax) > abs(plotmin) )
            tmpv2(maxidx, :) = 0;
        else
            tmpv2(minidx, :) = 0;
        end
        plotmin = min(tmpv2(:));
        plotmax = max(tmpv2(:));
    else
        plotokflag=1;
    end
    
    cnt = cnt+1;
    if cnt> maxloop
        plotokflag=1;
    end
end

fprintf(1, 'Done\n');
