% Calculate the PCA of all the mean spike snippets from Julie's data set

load meanSnipsAll     % Loads the mat file that was made with getMeanSnipsSS.m

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

G5Snips = allMeanSnips;
G5Snips(indL5,:) = [];

nSingle = size(G5Snips, 1);

% Modify sign of snippets so that the first peak, after 4 pts, is positive
for i=1:nSnips

    [pkVal pkInd] = findpeaks(allMeanSnips(i,4:end));    % Find peaks
    [vaVal vaInd] = findpeaks(-allMeanSnips(i,4:end));   % Find valleys
    
    if (isempty(pkVal))
        allMeanSnips(i,:) = -allMeanSnips(i,:);  % Flip if there are no peaks
    elseif (isempty(vaVal) )
        ;
    elseif vaInd(1) < pkInd(1)                 % Flip if valley is before peak
        allMeanSnips(i,:) = -allMeanSnips(i,:); %#ok<SAGROW>  
    end
    
end

% Annotate 300 random G5 snips
figure(1);
nAnn = 300;
indPlot = randperm(nSingle, nAnn);
typeSnips = zeros(nSingle, 1);
for i=1:nAnn
    plot(G5Snips(indPlot(i), :), 'k');

    snipType = input('Type in the snippet Type n (narrow), w (wide), l (large)', 's');
    while ~(strcmp(snipType, 'n') || strcmp(snipType, 'w') || strcmp(snipType, 'l') )
        fprintf(1,'Wrong input. Please use n, w or l\n');
        snipType = input('Type in the snippet Type n (narrow), w (wide), l (large)', 's');
    end
    typeSnips(indPlot(i)) = snipType;
    
end

save meanSnipsAllAnn allMeanSnips G5Snips typeSnips indPlot nAnn
    
