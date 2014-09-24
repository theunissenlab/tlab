% Calculate the PCA of all the mean spike snippets from Julie's data set

load meanSnips     % Loads the mat file that was made with getMeanSnips.m

nSnips = size(allSnips, 1);   % Number of mean snippets in our set
snipL = size(allSnips, 2);

% % Modify sign of snippets so that the first peak, after 4 pts, is positive
% for i=1:nSnips
% 
%     [pkVal pkInd] = findpeaks(allSnips(i,4:end));    % Find peaks
%     [vaVal vaInd] = findpeaks(-allSnips(i,4:end));   % Find valleys
%     
%     if (isempty(pkVal))
%         allSnips(i,:) = -allSnips(i,:);  % Flip if there are no peaks
%     elseif (isempty(vaVal) )
%         ;
%     elseif vaInd(1) < pkInd(1)                 % Flip if valley is before peak
%         allSnips(i,:) = -allSnips(i,:); %#ok<SAGROW>  
%     end
%     
% end

% Lign up the snips at positive peak 
centeredSnips = zeros(nSnips, 2*snipL+1);
centeredSnipsNorm = zeros(nSnips, 2*snipL+1);
for i=1:nSnips
    [pkVal pkInd] = findpeaks(allSnips(i,1:end));    % Find peaks again
    if (isempty(pkVal))
        pkVal(1) = allSnips(i,snipL);
        pkInd(1) = snipL;
    end
    centeredSnips(i, snipL+1-pkInd(1):2*snipL-pkInd(1)) = allSnips(i,:);
    maxVal = max(abs(allSnips(i,:)));
    centeredSnipsNorm(i,:) = centeredSnips(i,:)./maxVal;
end

% Plot 100 random snips
figure(1);
indPlot = randperm(nSnips, 100);
for i=1:100
    subplot(10,10,i);
    plot(centeredSnips(indPlot(i),:), 'k', 'LineWidth', 1);
    axis off
end
    

% Calculate the PCA of the snipets first without normalization
[eigenVect, coords, eigenVals] = princomp(centeredSnips);

% Now repeat calculation with normalized snippets
allSnipsNorm = allSnips;


% Calculate the PCA of the snipets first without normalization
[eigenVectNorm, coordsNorm, eigenValsNorm] = princomp(centeredSnipsNorm);

% Load the annotated snips - warning - this will overide allSnips
load meanSnipsAnn

% Plot the coordinates of the spike shapes in the PCA space 
figure(2);
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


% Plot the coordinates of the spike shapes in the PCA space 
figure(3);
subplot(1,2,1);
plot(coordsNorm(:,1), coordsNorm(:,2), '+');
title('Mean Snippet PCs Amplitude Normalized');

subplot(1,2,2);
plot(coordsNorm(indPlot(ind_n),1), coordsNorm(indPlot(ind_n), 2), '+');
hold on;
plot(coordsNorm(indPlot(ind_w),1), coordsNorm(indPlot(ind_w), 2), '+r');
plot(coordsNorm(indPlot(ind_l),1), coordsNorm(indPlot(ind_l), 2), '+g');
legend('narrow', 'wide', 'large');





