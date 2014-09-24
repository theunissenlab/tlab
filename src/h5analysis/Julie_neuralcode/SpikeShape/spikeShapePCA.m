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

% Plot 100 random snips
figure(1);
indPlot = randperm(nSnips, 100);
for i=1:100
    subplot(10,10,i);
    plot(allSnips(indPlot(i),:), 'k', 'LineWidth', 1);
    axis off
end
    

% Calculate the PCA of the snipets first without normalization
[eigenVect, coords, eigenVals] = princomp(allSnips);


% Plot the coordinates of the spike shapes in the PCA space 
figure(2);
plot(coords(:,1), coords(:,2), '+');
title('Mean Snippet PCs with Amplitude Information');

% This will only work if cursor info was saved...
% We are plotting example snippets from each group
load cursorfig2
figure(2);
hold on;
for i=1:3
    figure(2);
    plot(cursor_info(i).Position(1), cursor_info(i).Position(2), 'kx', 'MarkerSize', 12);
    figure(3);
    subplot(1,3,i);
    plot(allSnips(cursor_info(i).DataIndex, :), 'k');
    axis off;
end
figure(2);
hold off;


% Now repeat calculation with normalized snippets
allSnipsNorm = allSnips;

for i=1:nSnips
    [pkVal pkInd] = findpeaks(allSnipsNorm(i,4:end));    % Find peaks again
    allSnipsNorm(i,:) = allSnipsNorm(i,:)./pkVal(1);
end

% Calculate the PCA of the snipets first without normalization
[eigenVectNorm, coordsNorm, eigenValsNorm] = princomp(allSnipsNorm);


% Plot the coordinates of the spike shapes in the PCA space 
figure(4);
plot(coordsNorm(:,1), coordsNorm(:,2), '+');
title('Mean Snippet PCs Amplitude Normalized');

load cursorfig4
figure(4);
hold on;
for i=1:6
    figure(4);
    plot(cursor_info(i).Position(1), cursor_info(i).Position(2), 'kx', 'MarkerSize', 14, 'LineWidth', 3);
    figure(5);
    subplot(1,6,i);
    plot(allSnips(cursor_info(i).DataIndex, :), 'k');
    axis off;
end
figure(4);
hold off;

% Look at slope of decay
slopeSnips = zeros(nSnips, 1);
for i=1:nSnips
    [pkVal pkInd] = findpeaks(allSnipsNorm(i,4:end));    % Find peaks again
    if pkInd(1) <= snipL-5
        dVal = (pkVal(1)-allSnipsNorm(i,pkInd(1)+4))./5;
    else
        dVal = (pkVal(1)-allSnipsNorm(i,snipL))./(snipL-pkInd(1));
    end
    slopeSnips(i) = dVal;
     
end

figure(6);
hist(slopeSnips,100);
axis([0 0.1 0 70]);
title('Histogram of Spike Slopes');



