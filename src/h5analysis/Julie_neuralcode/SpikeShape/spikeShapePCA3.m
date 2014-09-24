% Calculate the PCA of all the mean spike snippets from Julie's data set

load meanSnips     % Loads the mat file that was made with getMeanSnips.m

nSnips = size(allSnips, 1);   % Number of mean snippets in our set
snipL = size(allSnips, 2);


powSnips = zeros(nSnips, snipL/2+1);

for i=1:nSnips
    fftOneSnip = fft(allSnips(i,:)-mean(allSnips(i,:)));
    powSnips(i,:) = abs(fftOneSnip(1:snipL/2+1));
    sumPow = sum(powSnips(i,:));
    powSnips(i,:) = powSnips(i,:)./sumPow;
end

% Plot 100 random snips
figure(1);
indPlot = randperm(nSnips, 100);
for i=1:100
    subplot(10,10,i);
    plot(powSnips(indPlot(i),:), 'k', 'LineWidth', 1);
    axis off
end
    

% Calculate the PCA of the snipets power first without normalization
[eigenVect, coords, eigenVals] = princomp(powSnips);


% Load the annotated snips - warning - this will overide allSnips
load meanSnipsAnn

% Plot the coordinates of the spike shapes in the PCA space 
figure(2);
subplot(1,2,1);
plot(coords(:,1), coords(:,2), '+');
title('Snippet PCs in Power Spectrum');
subplot(1,2,2);
ind_n = find(typeSnips =='n');
plot(coords(indPlot(ind_n),1), coords(indPlot(ind_n), 2), '+');
hold on;
ind_w = find(typeSnips == 'w');
plot(coords(indPlot(ind_w),1), coords(indPlot(ind_w), 2), '+r');
ind_l = find(typeSnips == 'l');
plot(coords(indPlot(ind_l),1), coords(indPlot(ind_l), 2), '+g');
legend('narrow', 'wide', 'large');

figure(3);
plot3(coords(:,1), coords(:,2), coords(:,3), '+');

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

% Show some examples in wide box 1 and wide box 2 that corresponds
% approximately to negative and positive values of PC1 for values of PC2
% between -0.016 and -0.08

indBox1 = find(coords(:,2)<-0.016 & coords(:,2) > -0.08 & coords(:,1) < -0.07);
indBox2 = find(coords(:,2)<-0.016 & coords(:,2) > -0.08 & coords(:,1) > 0);


% plot 100 snips from each
figure(5);
indPlotBox1 = randperm(length(indBox1), 100);
for i=1:100
    subplot(10,10,i);
    plot(allSnips(indBox1(indPlotBox1(i)),:), 'k', 'LineWidth', 1);
    axis off
end

figure(6);
indPlotBox2 = randperm(length(indBox2), 100);
for i=1:100
    subplot(10,10,i);
    plot(allSnips(indBox2(indPlotBox2(i)),:), 'k', 'LineWidth', 1);
    axis off
end





