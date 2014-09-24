% Load all meanSnips and calculate the SNR...


% First store the information for the multi unit data
load meanSnipsNS.mat
nUnits = length(allNames);

sdMultiUnits = zeros(1, nUnits);
snrMinMaxMultiUnits = zeros(1, nUnits);
snrMaxMultiUnits = zeros(1, nUnits);
snrMeanMultiUnits = zeros(1, nUnits);


for i = 1:nUnits
    
    maxVal = max(allMeanSnips(i,:));
    maxInd = find(allMeanSnips(i,:) == maxVal);
    minVal = min(allMeanSnips(i,:));
    minInd = find(allMeanSnips(i,:) == minVal);
    
    sdMultiUnits(i) = mean(allSdSnips(i,:));
    snrMinMaxMultiUnits(i) = (maxVal-minVal)./sqrt((allSdSnips(i,maxInd).^2 + allSdSnips(i,minInd).^2)/2);
    
    if abs(maxVal) > abs(minVal)
        snrMaxMultiUnits(i) = abs(maxVal)./allSdSnips(i,maxInd);
    else
        snrMaxMultiUnits(i) = abs(minVal)./allSdSnips(i,minInd);
    end
    
    snrMeanMultiUnits(i) = mean(allMeanSnips(i,:).^2)./mean(allSdSnips(i,:).^2);
    
    
end



% Load the spike sorted data and find other variables data
load meanSnipsAll.mat
nUnits = length(allNames);

sdSingleUnits = zeros(1, nUnits);
snrMinMaxSingleUnits = zeros(1, nUnits);
snrMaxSingleUnits = zeros(1, nUnits);
snrMeanSingleUnits = zeros(1, nUnits);


for i = 1:nUnits
    
    maxVal = max(allMeanSnips(i,:));
    maxInd = find(allMeanSnips(i,:) == maxVal);
    minVal = min(allMeanSnips(i,:));
    minInd = find(allMeanSnips(i,:) == minVal);
    
    sdSingleUnits(i) = mean(allSdSnips(i,:));
    snrMinMaxSingleUnits(i) = (maxVal-minVal)./sqrt((allSdSnips(i,maxInd).^2 + allSdSnips(i,minInd).^2)/2);
    
    if abs(maxVal) > abs(minVal)
        snrMaxSingleUnits(i) = abs(maxVal)./allSdSnips(i,maxInd);
    else
        snrMaxSingleUnits(i) = abs(minVal)./allSdSnips(i,minInd);
    end
    
    snrMeanSingleUnits(i) = mean(allMeanSnips(i,:).^2)./mean(allSdSnips(i,:).^2);
    
    
end



% Single vs Multi unit histogram
figure(1);
[nvalsMulti, xvalsMulti] = hist(sdMultiUnits);
[nvalsSingle] = hist(sdSingleUnits, xvalsMulti);

bar(xvalsMulti, [nvalsMulti'./sum(nvalsMulti), nvalsSingle'./sum(nvalsSingle)]);
ylabel('Prob');
xlabel('Standard Deviation');
legend('Non Sorted', 'Sorted');



figure(2);
[nvalsMulti, xvalsMulti] = hist(snrMinMaxMultiUnits);
[nvalsSingle] = hist(snrMinMaxSingleUnits, xvalsMulti);


bar(xvalsMulti, [nvalsMulti'./sum(nvalsMulti), nvalsSingle'./sum(nvalsSingle)]);
ylabel('Prob');
xlabel('SNR Min Max');
legend('Non Sorted', 'Sorted');

figure(3);
[nvalsMulti, xvalsMulti] = hist(snrMaxMultiUnits);
[nvalsSingle] = hist(snrMaxSingleUnits, xvalsMulti);


bar(xvalsMulti, [nvalsMulti'./sum(nvalsMulti), nvalsSingle'./sum(nvalsSingle) ]);
ylabel('Prob');
xlabel('SNR Max');
legend('Non Sorted', 'Sorted');

figure(4);
[nvalsMulti, xvalsMulti] = hist(snrMeanMultiUnits);
[nvalsSingle] = hist(snrMeanSingleUnits, xvalsMulti);


bar(xvalsMulti, [nvalsMulti'./sum(nvalsMulti), nvalsSingle'./sum(nvalsSingle) ]);
ylabel('Prob');
xlabel('SNR Mean');
legend('Non Sorted', 'Sorted');

% Print out some summary results.
nG5MultiUnit = length(find(snrMinMaxMultiUnits > 5));
nG5SingleUnit = length(find(snrMinMaxSingleUnits > 5 ));

fprintf(1, 'Number of non sorted units with Min Max SNR > 5 %d/%d or %.1f%% \n', nG5MultiUnit, length(snrMinMaxMultiUnits), ...
    100.0*nG5MultiUnit/length(snrMinMaxMultiUnits) );
fprintf(1, 'Number of sorted units with Min Max SNR > 5 %d/%d or %.1f%% \n', nG5SingleUnit, length(snrMinMaxSingleUnits), ...
    100.0*nG5SingleUnit/length(snrMinMaxSingleUnits) );


