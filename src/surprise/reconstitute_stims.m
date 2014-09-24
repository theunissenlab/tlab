tempDir = '/auto/fdata/pgill/Zador/Surprise/in_progress';
outDir = '/auto/fdata/pgill/Zador/Surprise/Stim_corp_only';
fiatdir(outDir);
names = dir(fullfile(tempDir,'*.mat'));
clear allStims

for iFile = 1:length(names)
    theName = [num2str(iFile) '.mat'];
    loaded = load(fullfile(tempDir,theName));
    for iStim = 1:length(loaded.out)
        if iFile == 1 
             allStims(iStim).stim_name = loaded.out(iStim).stim_name;
             allStims(iStim).louder = loaded.out(iStim).louder;
             allStims(iStim).quieter = loaded.out(iStim).quieter;
        else
             allStims(iStim).louder = [allStims(iStim).louder; loaded.out(iStim).louder];
             allStims(iStim).quieter = [allStims(iStim).quieter; loaded.out(iStim).quieter];
        end
    end

end

for iStim = 1:length(allStims)
    stimName = strrep(allStims(iStim).stim_name,'_',filesep);
    toWrite = fullfile(outDir,stimName);
    neededPath = fileparts(toWrite);
    fiatdir(neededPath);
    outstim = [allStims(iStim).quieter; allStims(iStim).louder];
    save(toWrite,'outstim');
end
