% this m files finds the names of the stimulus used for each neuron

cd /auto/fdata/pgill/dataCD_withXF
dirout = dir('/auto/fdata/pgill/dataCD_withXF/maybe_L');
birdnames = {dirout(:).name};
len = length(birdnames);
birdnames = birdnames(:,3:len);
len = length(birdnames);
stimNamesAll = cell(1,len);

for i = 1:len,
    direct = fullfile('/auto/fdata/pgill/dataCD_withXF/maybe_L',birdnames{i},'conspecific');
    cd(direct)
    dirout1 = dir('.');
    stims = {dirout1(:).name};
    len1 = length(stims);
    stims = stims(:,2+len1/2:len1);
    lenS = length(stims);
    stimNamesBird = cell(1,lenS);
    for j = 1:lenS,
        fid = fopen(stims{j},'r');
        stim =char(fread(fid,'char')');
        stimNamesBird{1,j} = stim;
        fclose(fid);
    end
    stimNamesAll{1,i} = stimNamesBird;
end

I= zeros(len,len);
for i = 1:len,
s1 = stimNamesAll{1,i};
s1 = s1(:,1:min(length(s1),20));
    for j = i+1:len,
        s2 = stimNamesAll{1,j};
        s2 = s2(:,1:min(length(s2),20));
        I(i,j) = length(intersect(s1,s2));
    end
end

tI = (I >= 20);
k = 1;
while(sum(tI(k,:)) <= 0)
k = k + 1;
end

sIndex= tI(k,:);
sIndex = find(sIndex);
birds = cell(1,length(sIndex)+1);
birds{1,1} = birdnames{1,k};
for i = 2:length(sIndex)+1
    birds{1,i} = birdnames{1,sIndex(i-1)};
end
stimulus = stimNamesAll{1,k}

