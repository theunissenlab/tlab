function [score] = avgCorrMetric(spiketrains)


for ii = 1:5
    cnt = 1;
    delta = inf;
    mclast = 0;

    while delta > 1e-7 || cnt < 1000

        idx = randperm(size(spiketrains,2));
        idx = idx(1:end-mod(size(idx,2),2));

        c(cnt) = corr(nanmean(spiketrains(:,idx(1:.5*length(idx))),2),nanmean(spiketrains(:,idx(.5*length(idx)+1:end)),2));

        score = nanmean(c);

        delta = abs(score - mclast);

        mclast = score;

        cnt = cnt+1;
    end
    
    s(ii) = score;
    
end

score = round(100*median(s));