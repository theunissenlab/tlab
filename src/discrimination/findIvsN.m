function[I, cH, iH] = findIvsN(spiketrain,winSize,N)
T = length(spiketrain);
I = zeros(N,N);
cH = zeros(N,N);
iH = zeros(N,N);
[I(1,1) cH(1,1) iH(1,1)] = findH1(spiketrain,winSize);

for i = 2:N,
    lenRemove = floor(length(spiketrain)/i);
    for j = 0:i-1,
       rt = spiketrain(:,j*lenRemove+1:j*lenRemove+lenRemove);
       [I(i,j+1) cH(i,j+1) iH(i,j+1)] = findH1(rt,winSize); 
    end
    % I(i) = mean(It);
    % cH(i) = mean(cHt);
    % iH(i) = mean(iHt);
end

I = I/(winSize*10^-3);

function [I cH iH] = findH1(rt,winSize)
y = 50; % Maximum No. of Spikes -1 in Window
[noTrials,t] = size(rt);
noWin = floor(t/winSize);

p = zeros(noWin,y);
for j = 1:noWin,
    q = zeros(1,y);
    for k = 1:noTrials,            
        x = sum(rt(k,j:j+winSize-1));
        if (x < y)
            q(x+1) = q(x+1) + 1;
        else
            q(y) = q(y) + 1;
            fprintf(1,'Warning in findH1: number of spikes %d exceeded limit %d\n', x+1, y);
        end
    end
    p(j,:) = q/noTrials;
end
[cH iH] = findH2(p,y,noWin);
I = iH - cH;

function [cH iH] = findH2(p,y,noWin)
x = 0;
for i = 1:noWin,
    for j = 1:y,
        if (p(i,j) > 0)
            x = x - p(i,j)*log2(p(i,j));
        end
    end
end
cH = x/noWin;

z = 0;
if (noWin > 1)
    q = sum(p)./noWin;
else
    q = p;
end
for i = 1:y,

    if(q(i) > 0)
        z = z - q(i)*log2(q(i));
    end
end
iH = z;

% lenChunk = floor(T/N)
% spiketrain = spiketrain(:,1:N*lenChunk);
% [I(1) cH(1) iH(1)] = findH1(spiketrain,winSize);
% for i = 1:N-1,
%     lenRemove = i*lenChunk
%     noRemoveChunks = floor(length(spiketrain)/lenRemove);
%     It = zeros(1,noRemoveChunks);
%     cHt = zeros(1,noRemoveChunks);
%     iHt = zeros(1,noRemoveChunks);
%     for j = 0:noRemoveChunks-1,
%         rt = spiketrain;
%         rt(:,j*lenRemove+1:j*lenRemove+lenRemove) = [];
%         [It(j+1) cHt(j+1) iHt(j+1)] = findH1(rt,winSize);
%     end
%     I(i+1) = mean(It);
%     cH(i+1) = mean(cHt);
%     iH(i+1) = mean(iHt);
% end

% x = 1/data size
% x = zeros(1,N);
% for i = 1:N,
% x(i) = 1/(length(spiketrain)-(i-1)*lenChunk);
% end