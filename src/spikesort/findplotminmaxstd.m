function [plotmin, plotmax] = findplotminmax(tmpv2)

plotokflag = 0;
cnt = 0;

plotmin = -4*std(tmpv2(:));
plotmax = 4*std(tmpv2(:));



