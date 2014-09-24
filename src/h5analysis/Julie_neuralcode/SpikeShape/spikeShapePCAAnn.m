% Calculate the PCA of all the mean spike snippets from Julie's data set

load meanSnipsAnn     % Loads the mat file that was made with getMeanSnips.m and then annotade by hand

nSnips = size(allSnips, 1);   % Number of mean snippets in our set
snipL = size(allSnips, 2);

% Lign up the snips at peak 
centeredSnips = zeros(nSnips, 2*snipL+1);
centeredSnipsNorm = zeros(nSnips, 2*snipL+1);
for i=1:nSnips
    [pkVal pkInd] = findpeaks(allSnips(i,4:end));    % Find peaks again
    centeredSnips(i, snipL+1-pkInd(1):2*snipL-pkInd(1)) = allSnips(i,:);
    centeredSnipsNorm(i,:) = centeredSnips(i,:)./pkVal(1);
end

% Calculate the PCA of the snipets first without normalization
[eigenVect, coords, eigenVals] = princomp(centeredSnips);


% Plot the coordinates of the spike shapes in the PCA space 
figure(1);
subplot(1,2,1);
plot(coords(:,1), coords(:,2), '+');
title('Mean Snippet PCs with Amplitude Information');

subplot(1,2,2);
ind_n = find(typeSnips =='n');
plot(coords(indPlot(ind_n),1), coords(indPlot(ind_n), 2), '+');
hold on;
ind_w = find(typeSnips == 'w');
plot(coords(indPlot(ind_w),1), coords(indPlot(ind_w), 2), '+r');
ind_l = find(typeSnips == 'l');
plot(coords(indPlot(ind_l),1), coords(indPlot(ind_l), 2), '+g');
legend('narrow', 'wide', 'large');



% Now repeat calculation with normalized snippets
allSnipsNorm = allSnips;

for i=1:nSnips
    [pkVal pkInd] = findpeaks(allSnipsNorm(i,4:end));    % Find peaks again
    allSnipsNorm(i,:) = allSnipsNorm(i,:)./pkVal(1);
end

% Calculate the PCA of the snipets first without normalization
[eigenVectNorm, coordsNorm, eigenValsNorm] = princomp(centeredSnipsNorm);


% Plot the coordinates of the spike shapes in the PCA space 
figure(2);
subplot(1,2,1);
plot(coordsNorm(:,1), coordsNorm(:,2), '+');
title('Mean Snippet PCs Amplitude Normalized');

subplot(1,2,2);
plot(coordsNorm(indPlot(ind_n),1), coordsNorm(indPlot(ind_n), 2), '+');
hold on;
plot(coordsNorm(indPlot(ind_w),1), coordsNorm(indPlot(ind_w), 2), '+r');
plot(coordsNorm(indPlot(ind_l),1), coordsNorm(indPlot(ind_l), 2), '+g');
legend('narrow', 'wide', 'large');






