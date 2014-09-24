function [b, c] = fit_gompertz_fmin_jacknife(xvals, yvals, groupIndex, a, b0, c0)

    groups = unique(groupIndex);
    
    nperms = round(length(groups) / 2);
    
    b = zeros(1, nperms);
    c = zeros(1, nperms);
    
    for k = 1:nperms
       
        h2 = k*2;
        h1 = h2-1;
        hg1 = groups(h1);
        hg2 = groups(h2);
        %holdoutGroups = groups(h1:h2);
        trainGroups = groups(groups ~= hg1 & groups ~= hg2);
       
        trainIndex = findIdx(trainGroups, groupIndex);
        
        [gb, gc] = fit_gompertz_fmin(xvals(trainIndex), yvals(trainIndex), a, b0, c0);

        b(k) = gb;
        c(k) = gc;        
        
    end
    