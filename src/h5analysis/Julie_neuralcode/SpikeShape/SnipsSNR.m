% Load all meanSnips and calculate the SNR...

% Load the data
load meanSnipsAll.mat


% Make a figure showing a couple of examples.
nsnips = size(allMeanSnips, 1);

indShow = randi(nsnips, [1, 9]);

figure(1)
for i=1:9
    subplot(3,3,i);
    
    plot(allMeanSnips(indShow(i),:),'k');
    hold on;
    plot(allMeanSnips(indShow(i),:)-allSdSnips(indShow(i),:),'k--');
    plot(allMeanSnips(indShow(i),:)+allSdSnips(indShow(i),:),'k--');
    title(allNames{indShow(i)});
    hold off;
end

snrSnips = zeros(1, nsnips);
snrSnips2 = zeros(1, nsnips);
snrSnips3 = zeros(1, nsnips);
sortColor = zeros(nsnips,1);
sortTypes = unique(allSortType);

for i=1:nsnips
    maxVal = max(allMeanSnips(i,:));
    maxInd = find(allMeanSnips(i,:) == maxVal);
    minVal = min(allMeanSnips(i,:));
    minInd = find(allMeanSnips(i,:) == minVal);
    
    snrSnips(i) = (maxVal-minVal)./sqrt((allSdSnips(i,maxInd).^2 + allSdSnips(i,minInd).^2)/2);
    
    if abs(maxVal) > abs(minVal)
        snrSnips2(i) = abs(maxVal)./allSdSnips(i,maxInd);
    else
        snrSnips2(i) = abs(minVal)./allSdSnips(i,minInd);
    end
    
    snrSnips3(i) = mean(allMeanSnips(i,:).^2)./mean(allSdSnips(i,:).^2);
    
    if strcmp(sortTypes{1}, allSortType{i})
        sortColor(i) = 1;
    elseif strcmp(sortTypes{2}, allSortType{i})
        sortColor(i) = 32;
    else
        sortColor(i) = 64;
    end
end

figure(2);
hist(snrSnips,20);
title('Max - Min SNR');

figure(3);
hist(snrSnips2,20);
title('Max SNR');

figure(4);
hist(snrSnips3,20);
title('Mean SNR');

figure(5);
subplot(1,3,1);
scatter(snrSnips, snrSnips2,6, sortColor, 'fill');
xlabel('Max - Min SNR');
ylabel('Max SNR');

subplot(1,3,2);
scatter(snrSnips, snrSnips3, 6, sortColor, 'fill');
xlabel('Max - Min SNR');
ylabel('Mean SNR');

subplot(1,3,3);
scatter(snrSnips2, snrSnips3, 6, sortColor, 'fill');
xlabel('Max SNR');
ylabel('Mean SNR');

colorVals = colormap();
text(30, 25, 'multi' , 'Color', colorVals(1,:));
text(30, 23, 'noise' , 'Color', colorVals(32,:));
text(30, 21, 'single' , 'Color', colorVals(64,:));


    