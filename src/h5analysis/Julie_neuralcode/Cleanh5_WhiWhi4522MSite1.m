%% This code is rewriting h5 files from site 1 WhiWhi4522M that come out...
... of tdt2hdf5 with a aborted recording session. It write new h5 files...
    ... without the data from the aborted recording session and add...
    ... "Old_Units" at the begining of the old h5.
addpath(genpath('/auto/fhome/julie/matlab/tlab'));
addpath(genpath('/auto/fhome/julie/matlab/STRFLab'));
addpath(genpath('/auto/k5/matlab714/toolbox/stats/stats'));
%cd /auto/fdata/julie/h5
cd /auto/k6/julie/h5/WhiWhi4522M
%cd /Users/elie/Documents/MATLAB/data/h5/WhiWhi4522M
Idir=pwd;
outdir=strcat(Idir, '/Old_Units');

h5filesall=dir(fullfile(Idir,'Site1*s*.h5'));
lh=length(h5filesall);
SS_Ind=zeros(lh,1);
for ff=1:lh
    if ~isempty(strfind(h5filesall(ff).name, 'ss'))
        SS_Ind(ff)=1;
    end
end
Indices=find(~SS_Ind);
LH=length(Indices);
for hh=1:LH
    %first read the unit
    nameh5=h5filesall(Indices(hh)).name;
    h5Path=fullfile(Idir, nameh5);
    unit = read_unit_h5file(h5Path, 'r');
    
    %then rename that old file
    movefile (h5Path, strcat(outdir, nameh5));
    fprintf(1,'the old %s is read and moved to Old_Units folder\n',h5Path);
                
    
    %then create a new one
    h5 = h5utils();
    
    % Create file
    unitPath=h5Path;
    fid = h5.create(unitPath);
    
    % Set top level attributes 
    % First the mandatory ones
    if (isfield(unit, 'electrode'))
        h5.set_attr(fid, '/', 'electrode', unit.electrode);
    else
        fprintf(1, 'Warning: unit structure in write_unit_h5file does not include site field\n');
    end
    if (isfield(unit, 'site'))
        h5.set_attr(fid, '/', 'site', unit.site);
    else
        fprintf(1, 'Warning: unit structure in write_unit_h5file does not include site field\n');
    end
    if (isfield(unit, 'source_directory'))
        h5.set_attr(fid, '/', 'source_directory', unit.source_directory);
    else
        fprintf(1, 'Warning: unit structure in write_unit_h5file does not include source_directory field\n');
    end
    
    % Then the more optional ones
    if (isfield(unit, 'ldepth'))
        if ~isempty(unit.ldepth)    
           h5.set_attr(fid, '/', 'ldepth', unit.ldepth);
        end
    end
    if (isfield(unit, 'rdepth'))
        if ~isempty(unit.rdepth)
           h5.set_attr(fid, '/', 'rdepth', unit.rdepth);
        end
    end
    if (isfield(unit, 'sortid'))
        h5.set_attr(fid, '/', 'sortid', unit.sortid);
    end
    if (isfield(unit, 'sortType'))
        h5.set_attr(fid, '/', 'sortType', unit.sortType);
    end

    % Create the class subgroups 
    iclass=2;
    h5.create_group(fid, ['/' unit.classes{iclass}]);
    responses = unit.class_responses.(unit.classes{iclass});
    nfiles = length(responses); 
    % Any attributes here?? */
        
    for ifile=1:nfiles
        fieldGroupName = ['/' unit.classes{iclass} '/' responses{ifile}.number];
        h5.create_group(fid, fieldGroupName);

        % Set the attributes of the file (stimulus) level
        responsesFieldNames = fieldnames(responses{ifile});
        nfields = length(responsesFieldNames);            
        for ifield=1:nfields
            if strcmp(responsesFieldNames{ifield}, 'number')
                continue;
            elseif strcmp(responsesFieldNames{ifield}, 'trials')
                continue;
            else
                h5.set_attr(fid, fieldGroupName, responsesFieldNames{ifield}, responses{ifile}.(responsesFieldNames{ifield}));
            end
        end

        % Now create the trials and write the attributes and the data
        ntrials = length(responses{ifile}.trials);                

        for j=1:ntrials

           trial = responses{ifile}.trials{j};

           % First create group
           respPath = sprintf('%s/%d', fieldGroupName, j);
           h5.create_group(fid, respPath);

           % Then set attributes...
           trialFieldNames = fieldnames(trial);
           nfields = length(trialFieldNames);            
           for ifield=1:nfields
               if strcmp(trialFieldNames{ifield}, 'spikeIds')
                   continue;
               elseif strcmp(trialFieldNames{ifield}, 'spikeTimes')
                   continue;
               elseif strcmp(trialFieldNames{ifield}, 'spikeWav')
                   continue;
               else
                   h5.set_attr(fid, respPath, trialFieldNames{ifield}, trial.(trialFieldNames{ifield}));
               end
           end


           % Finally the data

           % Our code for an empty trial
           if (isempty(trial.spikeIds))
               h5.set_ds(fid, respPath, 'spike_ids', -999);
               h5.set_ds(fid, respPath, 'spike_times', -999);
           else                  
               h5.set_ds(fid, respPath, 'spike_ids', trial.spikeIds');
               h5.set_ds(fid, respPath, 'spike_times', trial.spikeTimes');
           end

           % We could also save the snippet but are not doing that...             
        end    
    end
    
    % add the supplementary info
    pathPrefix = sprintf('/class_info/%s', unit.classes{iclass});
    h5.set_ds(fid, pathPrefix, 'mean_rate', unit.class_info.Call1c.meanSpikeRate);
    h5.set_ds(fid, [pathPrefix '/coherence'], 'cutoff_frequency', unit.class_info.Call1c.cutoffFreq);
    h5.set_ds(fid, [pathPrefix '/coherence'], 'frequencies', unit.class_info.Call1c.f);
    h5.set_ds(fid, [pathPrefix '/coherence'], 'upper', unit.class_info.Call1c.cUpper);
    h5.set_ds(fid, [pathPrefix '/coherence'], 'lower', unit.class_info.Call1c.cLower);
    h5.set_ds(fid, [pathPrefix '/coherence'], 'mean', unit.class_info.Call1c.c);

    h5.set_ds(fid, [pathPrefix '/coherence/info'], 'upper', unit.class_info.Call1c.infoUpper);
    h5.set_ds(fid, [pathPrefix '/coherence/info'], 'lower', unit.class_info.Call1c.infoLower);
    h5.set_ds(fid, [pathPrefix '/coherence/info'], 'mean', unit.class_info.Call1c.info); 

    pathPrefix = [sprintf('/extra_info/%s', unit.classes{iclass}) '/sections_data'];
    h5.set_ds(fid, pathPrefix, 'num_signif_sections_Bonferroni', unit.extra_info.Call1c.sections_data.nb_signif_sections_Bonf);
    h5.set_ds(fid, pathPrefix, 'num_signif_sections_BonferroniDeg', unit.extra_info.Call1c.sections_data.nb_signif_sections_BonfDeg);
    h5.set_ds(fid, pathPrefix, 'sections_pvalues', unit.extra_info.Call1c.sections_data.sections_pval);
    
    h5.close(fid);
end
    
    
    