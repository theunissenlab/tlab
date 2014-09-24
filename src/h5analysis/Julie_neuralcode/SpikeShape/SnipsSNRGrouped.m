% Load all meanSnips and calculate the SNR...

% Get spike type information
spikeInfo;

nMultiUnits = length(multiUnitFiles);
nSingleUnits = length(singleUnitFiles);
nLargeUnits = length(largeUnitFiles);
nIsolatedUnits = length(isolatedUnitFiles);
nNoiseUnits = length(noiseUnitFiles);

% First store the information for the multi unit data
load meanSnipsNS.mat
nUnits = length(allNames);

sdMultiUnits = zeros(1, nMultiUnits);
snrMinMaxMultiUnits = zeros(1, nMultiUnits);
snrMaxMultiUnits = zeros(1, nMultiUnits);
snrMeanMultiUnits = zeros(1, nMultiUnits);

for iu = 1:nMultiUnits
    notFoundFlg = 1;
    for i = 1:nUnits
        if strcmp(allNames{i}, multiUnitFiles{iu})
            notFoundFlg = 0;
            maxVal = max(allMeanSnips(i,:));
            maxInd = find(allMeanSnips(i,:) == maxVal);
            minVal = min(allMeanSnips(i,:));
            minInd = find(allMeanSnips(i,:) == minVal);
            
            sdMultiUnits(iu) = mean(allSdSnips(i,:));
            snrMinMaxMultiUnits(iu) = (maxVal-minVal)./sqrt((allSdSnips(i,maxInd).^2 + allSdSnips(i,minInd).^2)/2);
            
            if abs(maxVal) > abs(minVal)
                snrMaxMultiUnits(iu) = abs(maxVal)./allSdSnips(i,maxInd);
            else
                snrMaxMultiUnits(iu) = abs(minVal)./allSdSnips(i,minInd);
            end
            
            snrMeanMultiUnits(iu) = mean(allMeanSnips(i,:).^2)./mean(allSdSnips(i,:).^2);
            
            % Show a few examples on figure 1
            if (iu <= 5)
                figure(1);
                subplot(5,5,iu);
                plot(allMeanSnips(i,:),'k');
                hold on;
                plot(allMeanSnips(i,:)-allSdSnips(i,:),'k--');
                plot(allMeanSnips(i,:)+allSdSnips(i,:),'k--');
                axis([0 20 -1e-4 1e-4]);
                if(iu == 1)
                    title('Multi units');
                end
                hold off;
            end
            
            break;
        end
    end
    if notFoundFlg
        fprintf(1, 'Warning: no file match of unit %s\n', multiUnitFiles{iu});
    end
end


% Load the spike sorted data and find other variables data
load meanSnipsAll.mat
nUnits = length(allNames);

sdSingleUnits = zeros(1, nSingleUnits);
snrMinMaxSingleUnits = zeros(1, nSingleUnits);
snrMaxSingleUnits = zeros(1, nSingleUnits);
snrMeanSingleUnits = zeros(1, nSingleUnits);

for iu = 1:nSingleUnits
    notFoundFlg = 1;
    for i = 1:nUnits
        if strcmp(allNames{i}, singleUnitFiles{iu})
            notFoundFlg = 0;
            maxVal = max(allMeanSnips(i,:));
            maxInd = find(allMeanSnips(i,:) == maxVal);
            minVal = min(allMeanSnips(i,:));
            minInd = find(allMeanSnips(i,:) == minVal);
            
            sdSingleUnits(iu) = mean(allSdSnips(i,:));
            snrMinMaxSingleUnits(iu) = (maxVal-minVal)./sqrt((allSdSnips(i,maxInd).^2 + allSdSnips(i,minInd).^2)/2);
            
            if abs(maxVal) > abs(minVal)
                snrMaxSingleUnits(iu) = abs(maxVal)./allSdSnips(i,maxInd);
            else
                snrMaxSingleUnits(iu) = abs(minVal)./allSdSnips(i,minInd);
            end
            
            snrMeanSingleUnits(iu) = mean(allMeanSnips(i,:).^2)./mean(allSdSnips(i,:).^2);
            
                        % Show a few examples on figure 1
             if (iu <= 5)
                figure(1);
                subplot(5,5,5+iu);
                plot(allMeanSnips(i,:),'k');
                hold on;
                plot(allMeanSnips(i,:)-allSdSnips(i,:),'k--');
                plot(allMeanSnips(i,:)+allSdSnips(i,:),'k--');
                axis([0 20 -1e-4 1e-4]);
                if(iu == 1)
                    title('Single units');
                end
                hold off;
            end
            break;
        end
    end
    if notFoundFlg
        fprintf(1, 'Warning: no file match of unit %s\n', singleUnitFiles{iu});
    end
end


sdLargeUnits = zeros(1, nLargeUnits);
snrMinMaxLargeUnits = zeros(1, nLargeUnits);
snrMaxLargeUnits = zeros(1, nLargeUnits);
snrMeanLargeUnits = zeros(1, nLargeUnits);

for iu = 1:nLargeUnits
    notFoundFlg = 1;
    for i = 1:nUnits
        if strcmp(allNames{i}, largeUnitFiles{iu})
            notFoundFlg = 0;
            maxVal = max(allMeanSnips(i,:));
            maxInd = find(allMeanSnips(i,:) == maxVal);
            minVal = min(allMeanSnips(i,:));
            minInd = find(allMeanSnips(i,:) == minVal);
            
            sdLargeUnits(iu) = mean(allSdSnips(i,:));
            snrMinMaxLargeUnits(iu) = (maxVal-minVal)./sqrt((allSdSnips(i,maxInd).^2 + allSdSnips(i,minInd).^2)/2);
            
            if abs(maxVal) > abs(minVal)
                snrMaxLargeUnits(iu) = abs(maxVal)./allSdSnips(i,maxInd);
            else
                snrMaxLargeUnits(iu) = abs(minVal)./allSdSnips(i,minInd);
            end
            
            snrMeanLargeUnits(iu) = mean(allMeanSnips(i,:).^2)./mean(allSdSnips(i,:).^2);
            
            % Show a few examples on figure 1
             if (iu <= 5)
                figure(1);
                subplot(5,5,10+iu);
                plot(allMeanSnips(i,:),'k');
                hold on;
                plot(allMeanSnips(i,:)-allSdSnips(i,:),'k--');
                plot(allMeanSnips(i,:)+allSdSnips(i,:),'k--');
                axis([0 20 -1e-4 1e-4]);
                if(iu == 1)
                    title('Large units');
                end
                hold off;
            end
            break;
        end
    end
    if notFoundFlg
        fprintf(1, 'Warning: no file match of unit %s\n', largeUnitFiles{iu});
    end
end

sdIsolatedUnits = zeros(1, nIsolatedUnits);
snrMinMaxIsolatedUnits = zeros(1, nIsolatedUnits);
snrMaxIsolatedUnits = zeros(1, nIsolatedUnits);
snrMeanIsolatedUnits = zeros(1, nIsolatedUnits);

for iu = 1:nIsolatedUnits
    notFoundFlg = 1;
    for i = 1:nUnits
        if strcmp(allNames{i}, isolatedUnitFiles{iu})
            notFoundFlg = 0;
            maxVal = max(allMeanSnips(i,:));
            maxInd = find(allMeanSnips(i,:) == maxVal);
            minVal = min(allMeanSnips(i,:));
            minInd = find(allMeanSnips(i,:) == minVal);
            
            sdIsolatedUnits(iu) = mean(allSdSnips(i,:));
            snrMinMaxIsolatedUnits(iu) = (maxVal-minVal)./sqrt((allSdSnips(i,maxInd).^2 + allSdSnips(i,minInd).^2)/2);
            
            if abs(maxVal) > abs(minVal)
                snrMaxIsolatedUnits(iu) = abs(maxVal)./allSdSnips(i,maxInd);
            else
                snrMaxIsolatedUnits(iu) = abs(minVal)./allSdSnips(i,minInd);
            end
            
            snrMeanIsolatedUnits(iu) = mean(allMeanSnips(i,:).^2)./mean(allSdSnips(i,:).^2);
            
            % Show a few examples on figure 1
             if (iu <= 5)
                figure(1);
                subplot(5,5,15+iu);
                plot(allMeanSnips(i,:),'k');
                hold on;
                plot(allMeanSnips(i,:)-allSdSnips(i,:),'k--');
                plot(allMeanSnips(i,:)+allSdSnips(i,:),'k--');
                axis([0 20 -1e-4 1e-4]);
                if(iu == 1)
                    title('Isolated units');
                end
                hold off;
             end
            
            break;
        end
    end
    if notFoundFlg
        fprintf(1, 'Warning: no file match of unit %s\n', isolatedUnitFiles{iu});
    end
end

sdNoiseUnits = zeros(1, nNoiseUnits);
snrMinMaxNoiseUnits = zeros(1, nNoiseUnits);
snrMaxNoiseUnits = zeros(1, nNoiseUnits);
snrMeanNoiseUnits = zeros(1, nNoiseUnits);

for iu = 1:nNoiseUnits
    notFoundFlg = 1;
    for i = 1:nUnits
        if strcmp(allNames{i}, noiseUnitFiles{iu})
            notFoundFlg = 0;
            maxVal = max(allMeanSnips(i,:));
            maxInd = find(allMeanSnips(i,:) == maxVal);
            minVal = min(allMeanSnips(i,:));
            minInd = find(allMeanSnips(i,:) == minVal);
            
            sdNoiseUnits(iu) = mean(allSdSnips(i,:));
            snrMinMaxNoiseUnits(iu) = (maxVal-minVal)./sqrt((allSdSnips(i,maxInd).^2 + allSdSnips(i,minInd).^2)/2);
            
            if abs(maxVal) > abs(minVal)
                snrMaxNoiseUnits(iu) = abs(maxVal)./allSdSnips(i,maxInd);
            else
                snrMaxNoiseUnits(iu) = abs(minVal)./allSdSnips(i,minInd);
            end
            
            snrMeanNoiseUnits(iu) = mean(allMeanSnips(i,:).^2)./mean(allSdSnips(i,:).^2);
            
                        % Show a few examples on figure 1
             if (iu <= 5)
                figure(1);
                subplot(5,5,20+iu);
                plot(allMeanSnips(i,:),'k');
                hold on;
                plot(allMeanSnips(i,:)-allSdSnips(i,:),'k--');
                plot(allMeanSnips(i,:)+allSdSnips(i,:),'k--');
                axis([0 20 -1e-4 1e-4]);
                if(iu == 1)
                    title('Noise units');
                end
                hold off;
             end
             
            break;
        end
    end
    if notFoundFlg
        fprintf(1, 'Warning: no file match of unit %s\n', noiseUnitFiles{iu});
    end
end

% Single vs Multi unit histogram
figure(2);
[nvalsMulti, xvalsMulti] = hist(sdMultiUnits);
[nvalsSingle] = hist(sdSingleUnits, xvalsMulti);
[nvalsLarge] = hist(sdLargeUnits, xvalsMulti);
[nvalsIsolated] = hist(sdIsolatedUnits, xvalsMulti);
[nvalsNoise] = hist(sdNoiseUnits, xvalsMulti);

bar(xvalsMulti, [nvalsMulti'./sum(nvalsMulti), nvalsSingle'./sum(nvalsSingle), ...
                 nvalsLarge'./sum(nvalsLarge), nvalsIsolated'./sum(nvalsIsolated), ...
                 nvalsNoise'./sum(nvalsNoise) ]);
ylabel('Prob');
xlabel('Standard Deviation');
legend('Multi', 'Single', 'Large', 'Isolated', 'Noise');



figure(3);
[nvalsMulti, xvalsMulti] = hist(snrMinMaxMultiUnits);
[nvalsSingle] = hist(snrMinMaxSingleUnits, xvalsMulti);
[nvalsLarge] = hist(snrMinMaxLargeUnits, xvalsMulti);
[nvalsIsolated] = hist(snrMinMaxIsolatedUnits, xvalsMulti);
[nvalsNoise] = hist(snrMinMaxNoiseUnits, xvalsMulti);

bar(xvalsMulti, [nvalsMulti'./sum(nvalsMulti), nvalsSingle'./sum(nvalsSingle), ...
                 nvalsLarge'./sum(nvalsLarge), nvalsIsolated'./sum(nvalsIsolated), ...
                 nvalsNoise'./sum(nvalsNoise) ]);
ylabel('Prob');
xlabel('SNR Min Max');
legend('Multi', 'Single', 'Large', 'Isolated', 'Noise');

figure(4);
[nvalsMulti, xvalsMulti] = hist(snrMaxMultiUnits);
[nvalsSingle] = hist(snrMaxSingleUnits, xvalsMulti);
[nvalsLarge] = hist(snrMaxLargeUnits, xvalsMulti);
[nvalsIsolated] = hist(snrMaxIsolatedUnits, xvalsMulti);
[nvalsNoise] = hist(snrMaxNoiseUnits, xvalsMulti);

bar(xvalsMulti, [nvalsMulti'./sum(nvalsMulti), nvalsSingle'./sum(nvalsSingle), ...
                 nvalsLarge'./sum(nvalsLarge), nvalsIsolated'./sum(nvalsIsolated), ...
                 nvalsNoise'./sum(nvalsNoise) ]);
ylabel('Prob');
xlabel('SNR Max');
legend('Multi', 'Single', 'Large', 'Isolated', 'Noise');

figure(5);
[nvalsMulti, xvalsMulti] = hist(snrMeanMultiUnits);
[nvalsSingle] = hist(snrMeanSingleUnits, xvalsMulti);
[nvalsLarge] = hist(snrMeanLargeUnits, xvalsMulti);
[nvalsIsolated] = hist(snrMeanIsolatedUnits, xvalsMulti);
[nvalsNoise] = hist(snrMeanNoiseUnits, xvalsMulti);

bar(xvalsMulti, [nvalsMulti'./sum(nvalsMulti), nvalsSingle'./sum(nvalsSingle), ...
                 nvalsLarge'./sum(nvalsLarge), nvalsIsolated'./sum(nvalsIsolated), ...
                 nvalsNoise'./sum(nvalsNoise) ]);
ylabel('Prob');
xlabel('SNR Mean');
legend('Multi', 'Single', 'Large', 'Isolated', 'Noise');

