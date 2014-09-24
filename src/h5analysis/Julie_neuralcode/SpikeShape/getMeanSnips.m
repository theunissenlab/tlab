% addpath(genpath('/auto/fhome/julie/matlab/tlab'));
% addpath(genpath('/auto/k5/matlab714/toolbox/stats/stats'));
%cd /auto/fdata/julie/h5

% Get the birds from Julie h5 directory
input_dir = '/auto/k6/julie/h5';
Subjects = dir(input_dir);
ncells = 0;
allMeanSnips = [];
allSdSnips = [];


for ss=1:length(Subjects)   % Loop over birds
    Indiv=Subjects(ss).name;
    if length(Indiv)==11    % Bird names have 11 characters
        Idir=fullfile(input_dir, Indiv);
        h5files=dir(fullfile(Idir,'Site*ss*.h5'));
        
        
        LH=length(h5files);
        
        fprintf(1, 'Bird %s\n',Indiv);
        % figure();
        % title(sprintf('Bird %s\n', Indiv));
        for hh=1:LH
            h5file=fullfile(Idir, h5files(hh).name);
            
            fprintf(1,'Reading data from h5 file %s\n', h5file);
            unit = read_unit_h5file(h5file, 'r');
            fprintf(1,'...Done\n');
            % if strcmp(unit.sortType, 'single')
            [meanSnip, sdSnip] = meanSnip_h5(unit);
            allMeanSnips = [allMeanSnips; meanSnip];
            allSdSnips = [allSdSnips; sdSnip];
            ncells = ncells + 1;
            allNames{ncells} = h5file;
            allSortType{ncells} = unit.sortType;
            
        end
    else
        fprintf(1,'Warning: skipping %s - not a bird directory\n', Indiv);
    end
    
end

save meanSnipsAll.mat allMeanSnips allSdSnips allNames allSortType
