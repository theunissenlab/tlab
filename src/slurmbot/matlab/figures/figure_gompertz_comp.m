function allData = figure_gompertz_comp(cellNames)
    
    cinfo = get_cell_info();
    
    regions = {'mld', 'OV', 'L'};
    stimTypes = {'conspecific', 'bengalese', 'flatrip', 'songrip'};
    
    allData = struct;
    for m = 1:length(regions)
        for k = 1:length(stimTypes)
            allData.(regions{m}).(stimTypes{k}) = struct;
            allData.(regions{m}).(stimTypes{k}).stft.a = [];
            allData.(regions{m}).(stimTypes{k}).stft.b = [];
            allData.(regions{m}).(stimTypes{k}).stft.c = [];
            
            allData.(regions{m}).(stimTypes{k}).lyons.a = [];
            allData.(regions{m}).(stimTypes{k}).lyons.b = [];
            allData.(regions{m}).(stimTypes{k}).lyons.c = [];
        end
    end
    
    for k = 1:length(cellNames)
        
        cname = cellNames{k}; 
        reg = cinfo.(cname).region;
        if reg(1) == 'L'
            reg = 'L';
        end
        stypes = cinfo.(cname).stimTypes;        
        nstypes = length(stypes);
        
        for m = 1:nstypes
        
            stimType = stypes{m};
            outputDir = fullfile(cinfo.(cname).rootDir, stimType, 'output');
            
            stftFile = fullfile(outputDir, 'outputNL.gompertz.tfType_stft.txt');
            cvals = csvread(stftFile);
            
            allData.(reg).(stimType).stft.a(end+1) = cvals(1);
            allData.(reg).(stimType).stft.b(end+1) = cvals(2);
            allData.(reg).(stimType).stft.c(end+1) = cvals(3);            
            
            lyonsFile = fullfile(outputDir, 'outputNL.gompertz.tfType_lyons.txt');
            cvals = csvread(lyonsFile);
            
            allData.(reg).(stimType).lyons.a(end+1) = cvals(1);
            allData.(reg).(stimType).lyons.b(end+1) = cvals(2);
            allData.(reg).(stimType).lyons.c(end+1) = cvals(3);                        
        end
        
    end
    
    %% plot by region
    for m = 1:length(regions)
        reg = regions{m};
        aStft.(reg) = [];
        bStft.(reg) = [];
        cStft.(reg) = [];
        aLyons.(reg) = [];
        bLyons.(reg) = [];
        cLyons.(reg) = [];
    end    
    
    for k = 1:length(stimTypes)    
        for m = 1:length(regions)           
            reg = regions{m};
            aStft.(reg) = [aStft.(reg) allData.(reg).(stimTypes{k}).stft.a];
            bStft.(reg) = [bStft.(reg) allData.(reg).(stimTypes{k}).stft.b];
            cStft.(reg) = [cStft.(reg) allData.(reg).(stimTypes{k}).stft.c];
            
            aLyons.(reg) = [aLyons.(reg) allData.(reg).(stimTypes{k}).lyons.a];
            bLyons.(reg) = [bLyons.(reg) allData.(reg).(stimTypes{k}).lyons.b];
            cLyons.(reg) = [cLyons.(reg) allData.(reg).(stimTypes{k}).lyons.c];
        end
    end
    
    for m = 1:length(regions)
        reg = regions{m};
    
        figure('Name', reg); hold on;
        
        subplot(3, 2, 1); hold on;
        hist(aStft.(reg), 20);
        title(sprintf('a | stft | %0.3f +/- %0.4f', mean(aStft.(reg)), std(aStft.(reg))));
        
        subplot(3, 2, 3); hold on;
        hist(bStft.(reg), 20);
        title(sprintf('b | stft | %0.3f +/- %0.4f', mean(bStft.(reg)), std(bStft.(reg))));
        
        subplot(3, 2, 5); hold on;
        hist(cStft.(reg), 20);
        title(sprintf('c | stft | %0.3f +/- %0.4f', mean(cStft.(reg)), std(cStft.(reg))));
        
        
        subplot(3, 2, 2); hold on;
        hist(aLyons.(reg), 20);
        title(sprintf('a | lyons | %0.3f +/- %0.4f', mean(aLyons.(reg)), std(aLyons.(reg))));
        
        subplot(3, 2, 4); hold on;
        hist(bLyons.(reg), 20);
        title(sprintf('b | lyons | %0.3f +/- %0.4f', mean(bLyons.(reg)), std(bLyons.(reg))));
        
        subplot(3, 2, 6); hold on;
        hist(cLyons.(reg), 20);
        title(sprintf('c | lyons | %0.3f +/- %0.4f', mean(cLyons.(reg)), std(cLyons.(reg))));
        
    end
    
    
    