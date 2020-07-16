function [bandwidth]=make_bandwidth(meanspike);
[y, fo, meanP, Pupper, Plower, stP]=mtpsd_JN([meanspike'-mean(meanspike) meanspike'-mean(meanspike)], 128, 1000);

%[Pxx,Pxxc,F] = pmtm(meanspike,[],128,1000);

cutindex=min(find(Plower==0));
if isempty(cutindex)
    cutindex=length(fo);
end
power=meanP(1:cutindex,1,2)./sum(meanP(1:cutindex,1,2));
fw=fo(1:cutindex);

bandwidth=sum((fw.^2).*power);
bandwidth=sqrt(bandwidth);