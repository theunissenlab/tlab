function [info infoInterval] = findInfoBiasCorr(spiketrain,winSize)
% this function finds I(S;R) for a given neuron
% two options: 1) without bias correction, 2) with bias correction
% 1) without bias correction

N = 10;
[I, cH, iH] = findIvsN(spiketrain,winSize,N);
infoValues = I(1:1);
x = 1;
for i = 2:N,
    x = cat(2,x,i*ones(1,i));
    infoValues = cat(2,infoValues,I(i,1:i));
end

[B,BINT,R,RINT,STATS] = regress(infoValues',[x' ones(length(infoValues),1)]);
info = B(2);
infoInterval = BINT(2,:);
plot(x,infoValues,'+');
title('Bias Correction Plot');
xlabel('1/Fraction Data Used');
ylabel('Information (Bits/sec)');
hold on;
plot(x,x*B(1) + B(2),'k');
hold off;
    
function [cH iH] = findH1(rt,winSize)
y = 50; % Maximum No. of Spikes -1 in Window
[noTrials,t] = size(rt);
noWin = floor(t/winSize);
p = zeros(noWin,y);
for j = 1:noWin,
    q = zeros(1,y);
    for k = 1:noTrials,            
        x = sum(rt(k,j:j+winSize-1));
        if ( x+1 > y)
            fprintf(1, 'Warning: number of spikes in bin exceeds limit set at %d', y);
            x = y-1;
        end
        q(x+1) = q(x+1) + 1;
    end
    p(j,:) = q/noTrials;
end
[cH iH] = findH2(p,y,noWin);

function [I,Q] = findMI(p,y,noWin)
Q = p./noWin;
[r c] = size(Q);
Qr1 = sum(Q,2);
Qr2 = sum(Q,1);

t = 0;
for i = 1:r,
    for j = 1:c,
      if (Q(i,j) > 0)
          t = t+ Q(i,j)*log2(Q(i,j)/(Qr1(i)*Qr2(j)));
      end
    end
end
I = t;

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
q = sum(p)./noWin;
for i = 1:y,
    if(q(i) > 0)
        z = z - q(i)*log2(q(i));
    end
end
iH = z;