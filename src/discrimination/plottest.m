clear
load zfnormal_fieldL.mat

ndat = length(all_outputs);

mizS = zeros (ndat, 4);
mizT = zeros (ndat, 4);

pzS = zeros (ndat, 4);
pzT = zeros (ndat, 4);


mizT2 = zeros (ndat, 4);

pzT2 = zeros (ndat, 4);


mizSA = zeros (ndat, 4);
mizTA = zeros (ndat, 4);


pzSA = zeros (ndat, 4);
pzTA = zeros (ndat, 4);

mizSB = zeros (ndat, 4);
mizTB = zeros (ndat, 4);
mizncSB = zeros (ndat, 4);
mizncTB = zeros (ndat, 4);

pzSB = zeros (ndat, 4);
pzTB = zeros (ndat, 4);



miConf = zeros (ndat, 4);
miConfA = zeros (ndat, 4);
miConfB = zeros (ndat, 4);

perCor = zeros(ndat, 4);
perCorA = zeros(ndat, 4);
perCorB = zeros(ndat, 4);

migamma = zeros(ndat, 1);

for idat = 1:ndat
    mizS(idat,:) = all_outputs{idat,2}.mizdistS;
    mizT(idat,:) = all_outputs{idat,2}.mizdistT;
%    mizT2(idat,:) = all_outputs{idat,2}.mizdistT2;    
  
    mizSB(idat,:) = all_outputs{idat,2}.mizdistSB;
    mizTB(idat,:) = all_outputs{idat,2}.mizdistTB;
    mizncSB(idat,:) = all_outputs{idat,2}.mizdistncSB;
    mizncTB(idat,:) = all_outputs{idat,2}.mizdistncTB;
    
    pzS(idat,:) = all_outputs{idat,2}.pzdistS;
    pzT(idat,:) = all_outputs{idat,2}.pzdistT;
%    pzT2(idat,:) = all_outputs{idat,2}.pzdistT2;    
 
    pzSB(idat,:) = all_outputs{idat,2}.pzdistSB;
    pzTB(idat,:) = all_outputs{idat,2}.pzdistTB;
    
    miConf(idat,:) = all_outputs{idat,2}.mi_confusion;
%    miConfB(idat,:) = all_outputs{idat,2}.mi_confusionB;
    
    perCor(idat,:) = all_outputs{idat,2}.percorrect;
%    perCorB(idat,:) = all_outputs{idat,2}.percorrectB;
    
    migamma(idat) = all_outputs{idat,2}.gamma_mutual_info;
    
    
end
  
figure(1);
subplot(2,2,1);
plot(mizSB(:,1), migamma);
plot(mizSB(:,1), migamma, '+');
xlabel('MI Z-distance Song (bit/s');
ylabel('MI Gamma (bits/s)');
title('10 ms');

subplot(2,2,2);
plot(mizSB(:,2), migamma);
plot(mizSB(:,2), migamma, '+');
xlabel('MI Z-distance Song (bit/s');
ylabel('MI Gamma (bits/s)');
title('30 ms');

subplot(2,2,3);
plot(mizSB(:,3), migamma);
plot(mizSB(:,3), migamma, '+');
xlabel('MI Z-distance Song (bit/s');
ylabel('MI Gamma (bits/s)');
title('50 ms');

subplot(2,2,4);
plot(mizSB(:,4), migamma);
plot(mizSB(:,4), migamma, '+');
xlabel('MI Z-distance Song (bit/s');
ylabel('MI Gamma (bits/s)');
title('100 ms');

figure(2);
plot(mizTB(:,2), migamma, '+');
xlabel('MI Z-distance Trial(bit/s)');
ylabel('MI Gamma (bits/s)');
title('10 ms');

figure(3);
plot(mizTB(:,2), mizSB(:,2), '+');
ylabel('MI Z-distance Song (bits/s)');
xlabel('MI Z-distance Trial(bit/s)');

figure(4);
plot(mizSB(:,2), mizncSB(:,2), '+');
xlabel('MI Z-distance Song (bit/s)');
ylabel('MI Z-distance Song No Correction (bits/s)');

figure(5);
plot(mizSB(:,2), mizncSB(:,2), '+');
xlabel('MI Z-distance Song (bits/s)');
ylabel('MI Z-distance Song No Correction (bit/s)');

figure(6);
plot(mizSB(:,2), miConf(:,2), '+');
xlabel('MI Z-distance Song (bits/s)');
ylabel('MI Confusion Matrix (bit/s)');

