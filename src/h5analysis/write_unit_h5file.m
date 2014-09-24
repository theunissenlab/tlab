function unit = write_unit_h5file(unitPath, unit)
% Write the h5 file for Theunissen Neurophysiology data using all the
% information in the unit structure.  For consistency write_unit_h5file and
% read_unit_h5_file should be consistant...

    h5 = h5utils();
    
    % Create file
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
    nclasses = length(unit.classes);
    for iclass=1:nclasses
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
                   h5.set_ds(fid, respPath, 'spike_ids', trial.spikeIds);
                   h5.set_ds(fid, respPath, 'spike_times', trial.spikeTimes);
               end
               
               % We could also save the snippet but are not doing that...             
            end    
        end        
    end
    
    h5.close(fid);
    
    