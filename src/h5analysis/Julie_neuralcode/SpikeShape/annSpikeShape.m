% Calculate the PCA of all the mean spike snippets from Julie's data set

load meanSnips     % Loads the mat file that was made with getMeanSnips.m

nSnips = size(allSnips, 1);   % Number of mean snippets in our set
snipL = size(allSnips, 2);

% Modify sign of snippets so that the first peak, after 4 pts, is positive
for i=1:nSnips

    [pkVal pkInd] = findpeaks(allSnips(i,4:end));    % Find peaks
    [vaVal vaInd] = findpeaks(-allSnips(i,4:end));   % Find valleys
    
    if (isempty(pkVal))
        allSnips(i,:) = -allSnips(i,:);  % Flip if there are no peaks
    elseif (isempty(vaVal) )
        ;
    elseif vaInd(1) < pkInd(1)                 % Flip if valley is before peak
        allSnips(i,:) = -allSnips(i,:); %#ok<SAGROW>  
    end
    
end

% Annotate 300 random snips
figure(1);
nAnn = 300;
indPlot = randperm(nSnips, nAnn);
typeSnips = zeros(nSnips, 1);
for i=1:nAnn
    plot(allSnips(indPlot(i), :), 'k');

    snipType = input('Type in the snippet Type n (narrow), w (wide), l (large)', 's');
    while ~(strcmp(snipType, 'n') || strcmp(snipType, 'w') || strcmp(snipType, 'l') )
        fprintf(1,'Wrong input. Please use n, w or l\n');
        snipType = input('Type in the snippet Type n (narrow), w (wide), l (large)', 's');
    end
    typeSnips(i) = snipType;
    
end

save meanSnipsAnn allSnips typeSnips indPlot nAnn
    
