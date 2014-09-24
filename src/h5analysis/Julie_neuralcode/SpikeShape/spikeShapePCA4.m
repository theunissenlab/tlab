% Calculate the PCA of all the mean spike snippets from Julie's data set

load meanSnipsAll     % Loads the mat file that was made with getMeanSnipsSS.m
tdt25k = 24414.0625;

nSnips = size(allMeanSnips, 1);   % Number of mean snippets in our set
snipL = size(allMeanSnips, 2);

snrMinMaxUnits = zeros(1, nSnips);


% Select all files that have min-max snr > 5
for i = 1:nSnips    
    maxVal = max(allMeanSnips(i,:));
    maxInd = find(allMeanSnips(i,:) == maxVal);
    minVal = min(allMeanSnips(i,:));
    minInd = find(allMeanSnips(i,:) == minVal);
    
    snrMinMaxUnits(i) = (maxVal-minVal)./sqrt((allSdSnips(i,maxInd).^2 + allSdSnips(i,minInd).^2)/2);  
end
indL5 = find(snrMinMaxUnits < 5 );
indG5 = find(snrMinMaxUnits >= 5);

G5Snips = allMeanSnips;
G5Snips(indL5,:) = [];

nSingle = size(G5Snips, 1);


powSnips = zeros(nSingle, snipL/2+1);

for i=1:nSingle
    fftOneSnip = fft(G5Snips(i,:)-mean(G5Snips(i,:)));
    powSnips(i,:) = abs(fftOneSnip(1:snipL/2+1));
    sumPow = sum(powSnips(i,:));
    powSnips(i,:) = powSnips(i,:)./sumPow;
end

% Plot 100 random snips
figure(1);
indPlot = randperm(nSingle, 100);
for i=1:100
    subplot(10,10,i);
    plot(powSnips(indPlot(i),:), 'k', 'LineWidth', 1);
    axis off
end
    

% Calculate the PCA of the snipets power first without normalization
[eigenVect, coords, eigenVals] = princomp(powSnips);

% Perform a k-means clustering to get the three groups.
[indKmeans,clusterMeans,sumD,distanceToMean] = kmeans(coords(:,1:3),4);

% Load the annotated snips - warning - this will overide allSnips
% load meanSnipsAnn
load meanSnipsAllAnn

% Plot the coordinates of the spike shapes in the PCA space 
figure(2);
subplot(1,3,1);
plot(coords(:,1), coords(:,2), '+');
xlabel('PCA 1 Coords');
ylabel('PCA 2 Coords');
title('Snippet PCs in Power Spectrum');
subplot(1,3,2);
ind_n = find(typeSnips =='n');
plot(coords(ind_n,1), coords(ind_n, 2), '+');
hold on;
ind_w = find(typeSnips == 'w');
plot(coords(ind_w,1), coords(ind_w, 2), '+r');
ind_l = find(typeSnips == 'l');
plot(coords(ind_l,1), coords(ind_l, 2), '+g');
legend('narrow', 'wide', 'large');

subplot(1,3,3);
plot(coords(indKmeans==1,1), coords(indKmeans==1, 2), '+');
hold on;
plot(coords(indKmeans==2,1), coords(indKmeans==2, 2), '+r');
plot(coords(indKmeans==3,1), coords(indKmeans==3, 2), '+g');
plot(coords(indKmeans==4,1), coords(indKmeans==4, 2), '+k');
plot(clusterMeans(1,1),clusterMeans(1,2),'x',...
     'MarkerSize',16,'LineWidth',4)
 plot(clusterMeans(2,1),clusterMeans(2,2),'rx',...
     'MarkerSize',16,'LineWidth',4)
 plot(clusterMeans(3,1),clusterMeans(3,2),'gx',...
     'MarkerSize',16,'LineWidth',4)
  plot(clusterMeans(4,1),clusterMeans(4,2),'kx',...
     'MarkerSize',16,'LineWidth',4)

legend('Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4');
hold off;

figure(3);
plot3(coords(indKmeans==1,1), coords(indKmeans==1,2), coords(indKmeans==1,3), '+');
hold on;
plot3(coords(indKmeans==2,1), coords(indKmeans==2,2), coords(indKmeans==2,3), '+r');
plot3(coords(indKmeans==3,1), coords(indKmeans==3,2), coords(indKmeans==3,3), '+g');
plot3(coords(indKmeans==4,1), coords(indKmeans==4,2), coords(indKmeans==4,3), '+k');
hold off;

figure(4);
subplot(1,3,1);
plot(mean(powSnips), 'k');
title('Mean');
subplot(1,3,2);
plot(eigenVect(:,1),'k');
title('PCA1');
subplot(1,3,3);
plot(eigenVect(:,2),'k');
title('PCA2');


% Figure 5 - plot 30 examples of each.
figure(5);
nPlot = 30;
titleName = {'unclassified'; 'large'; 'wide'; 'narrow'};
        

for nk=1:4  % for the four clusters
    subplot(1,4,nk);
    
    nCluster = sum(indKmeans == nk);
    plotInd = randperm(nCluster, nPlot);
    clusterInd = find(indKmeans == nk);
    for iPlot=1:nPlot
        snipToPlot = G5Snips(clusterInd(plotInd(iPlot)),:);
        
        % We are going to normalize by the first peak (usually in the first
        % half) and display going up...
        maxval = max(snipToPlot(1:10));   
        minval = min(snipToPlot(1:10));
        if abs(minval) > maxval
            snipToPlot = -snipToPlot;
        end
        maxval = max(snipToPlot(1:10)); 
        snipToPlot = snipToPlot./maxval;
        plot(snipToPlot, 'k-');
        hold on;
    end
    title(titleName{nk});
    hold off;
end

    
    


% The code below might have to be changed everytime the clustering is run
% because the cluster id (indKmeans) is random...

% Assigns the cluster type to the single spikes
% Find the mean width at have maximum for each group

spikeWidth = zeros(nSingle,1);
spikeType = cell(nSingle,1);
nSkip = 1;
for i=1:nSingle
    
    % Assign Cluster type
    if indKmeans(i) == 1
        spikeType{i} = 'unclassified';
    elseif indKmeans(i) == 2
        spikeType{i} = 'large';
    elseif indKmeans(i) == 3
        spikeType{i} = 'wide';
    else
        spikeType{i} = 'narrow';
    end
        
    meanSnip = allMeanSnips(indG5(i),:);
    [pkVal, pkInd] = findpeaks(meanSnip(1,nSkip:end));
    tEnd = snipL;
    for it=pkInd(1)+nSkip:snipL
        if meanSnip(it) < pkVal(1)./2
            tEnd = it;
            break;
        end
    end
    tBeg = 1;
    for it=pkInd(1)+nSkip-1:-1:1
        if meanSnip(it) < pkVal(1)./2
            tBeg = it;
            break
        end
    end
    spikeWidth(i) = 1000.0*(tEnd - tBeg)./tdt25k; % Spike width in ms
    if tEnd == tBeg
        fprintf(1,'Time end = %d Time beg = %d Amp = %f\n', tEnd, tBeg, pkVal(1));
    end
  
end

[meanSpikeWidth, semSpikeWidth, stdSpikeWidth] = grpstats(spikeWidth, indKmeans, {'mean'; 'sem'; 'std'});
[p, table, stats] = anova1(spikeWidth, indKmeans);
fprintf(1,'Mean Spike Widths at half maximum for each group:\n');
fprintf(1,'\tUnclassified %.2f +- %.3f sd: %.3f\n', meanSpikeWidth(1), semSpikeWidth(1), stdSpikeWidth(1));
fprintf(1,'\tLarge %.2f +- %.3f sd: %.3f\n', meanSpikeWidth(2), semSpikeWidth(2), stdSpikeWidth(2));
fprintf(1,'\tWide %.2f +- %.3f sd: %.3f\n', meanSpikeWidth(3), semSpikeWidth(3), stdSpikeWidth(3));
fprintf(1,'\tNarrow %.2f +- %.3f sd: %.3f\n', meanSpikeWidth(4), semSpikeWidth(4), stdSpikeWidth(4));


save spikeShapeResults allMeanSnips snrMinMaxUnits indL5 indG5 allNames allSortType spikeType spikeWidth





