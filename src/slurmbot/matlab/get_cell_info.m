function cellinfoByName = get_cell_info()

    cellinfoByName = struct;

    dataDir = '/auto/fdata/mschachter/data';
    
    sbDir = '~/berkeley/tlab/trunk/src/slurmbot';
    sbDataDir = fullfile(sbDir, 'data');
    
    allCellsData = fullfile(sbDataDir, 'celldata.all.csv');
    fid = fopen(allCellsData, 'r');    
    cellData = textscan(fid, '%s %s %s', 'delimiter', ',');
    fclose(fid);
    
    for k = 1:length(cellData{1})
        
        bname = cellData{1}{k};
        cname = cellData{2}{k};
        reg = cellData{3}{k};
        cellName = sprintf('%s_%s', bname, cname);        
        cellinfoByName.(cellName) = struct;
        cellinfoByName.(cellName).region = reg;
        cellinfoByName.(cellName).rootDir = fullfile(dataDir, cellName);
    end
    clear cellData;
    
    stimTypesData = fullfile(sbDataDir, 'stimtypes.csv');
    fid = fopen(stimTypesData, 'r');
    stimData = textscan(fid, '%s %s %s %s', 'delimiter', ',');
    fclose(fid);
    
    for k = 1:length(stimData{1})
       
        cellName = stimData{1}{k};
        stimTypes = {};
        for m = 2:4           
            stype = stimData{m}{k};
            if ~isempty(stype)
                stimTypes{end+1} = stype;
            end
        end
        cellinfoByName.(cellName).stimTypes = stimTypes;        
    end
    clear stimData;
    
    
    
    
    