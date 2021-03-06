% Runs the neural discrimination for normal L
brainregion = 'sure_OV';

% All the data files pasted from the Excel file
ov_singleunits = { ...
    {'pupi0414', '1' }, ...
    {'pupi0414', '2'}, ...
    {'pupi0414', '3' }, ...
    {'pupi0414', '8' }, ...
    {'pupi0414', '9' }, ...
    {'pupi0414', '10' }, ...
    {'glb5202', '5' }, ...
    {'glb5202', '6' }, ...
    {'glb5202', '7' }, ...
    {'glb5656', '2' }, ... %note that this bird glb5656 has file names that are glb5202 from last bird
    {'glb5656', '3' }, ... %ditto for this cell and all cells for this bird
    {'glb5656', '4' }, ... %this cell had both S1 and M, and I think only con was played for M, so use S1...
    {'glb5656', '5' }, ... 
    {'blahp4108', '4' }, ...
    {'blahp4108', '5' }, ...
    {'obla1305', '8' }, ...
    {'obla17102', '2' }, ...
    {'obla17102', '7' }, ... %{'gr0869', '2' }, ... %this cell ran 2 blocks of same stimuli to capture more spikes, so didn't end up STRFing, so not in STRF database...
    {'gr0869', '5' }, ...
    {'gr0965', '4' }, ...
    {'gr0965', '6' }, ...
    {'gr1070', '1' }, ...
    {'gr1070', '2' }, ...
    {'gr1070', '3' }, ...
    {'gr1070', '4' }, ...
    {'blugr1612', '4' }, ... %both cells 4 and 5 seem crappy for STRFs, so wondering whether should leave 'em in the database...
    {'blugr1612', '5' }}; %both cells 4 and 5 seem crappy for STRFs, so wondering whether should leave 'em in the database...

ov_twounits = { ...
    {'pupi0414', '6' }, ...
    {'pupi0414', '7' }, ...
    {'glb5202', '9' }, ...
    {'obla1305', '4' }, ...
    {'obla1305', '5' }, ... 
    {'obla1305', '6' }, ...
    {'obla1305', '7' }, ...
    {'blahp9904', '6' }, ...
    {'gr0869', '8' }, ...
    {'gr0869', '9' }, ...
    {'gr0965', '1' }, ...
    {'gr0965', '2' }, ...
    {'gr0965', '5' }, ...
    {'gr1070', '5' }, ...
    {'blugr1612', '3' }};

ov_multiunits = { ...
    {'pupi0414', '5' }, ... %kicking in trial 14, but should be okay to run...
    {'blahp9904', '7' }, ...
    {'blahp9904', '8' }, ...
    {'blahp9904', '9' }, ...
    {'obla17102', '4' }, ...
    {'obla17102', '5' }, ...
    {'obla17102', '6' }, ...
    {'obla17102', '8' }, ...
    {'obla17102', '9' }, ...
    {'obla17102', '10' }, ... 
    ... %{'gr0869', '7' }, ... %this one only has 11 trials, so didn't STRF, so probably not in database...
    {'gr0869', '10' }, ...
    {'gr0869', '11' }, ...
    {'gr0965', '3' }, ...
    {'gr0965', '7' }, ...
    {'blugr1612', '2' }}; 

ov_allunits = { ...
    {'pupi0414', '1' }, ...
    {'pupi0414', '2'}, ...
    {'pupi0414', '3' }, ...
    {'pupi0414', '8' }, ...
    {'pupi0414', '9' }, ...
    {'pupi0414', '10' }, ...
    {'glb5202', '5' }, ...
    {'glb5202', '6' }, ...
    {'glb5202', '7' }, ...
    {'glb5656', '2' }, ... %note that this bird glb5656 has file names that are glb5202 from last bird
    {'glb5656', '3' }, ... %ditto for this cell and all cells for this bird
    {'glb5656', '4' }, ... %this cell had both S1 and M, and I think only con was played for M, so use S1...
    {'glb5656', '5' }, ... 
    {'blahp4108', '4' }, ...
    {'blahp4108', '5' }, ...
    {'obla1305', '8' }, ...
    {'obla17102', '2' }, ...
    {'obla17102', '7' }, ...
    ...%{'gr0869', '2' }, ... %this cell ran 2 blocks of same stimuli to capture more spikes, so didn't end up STRFing, so not in STRF database...
    {'gr0869', '5' }, ...
    {'gr0965', '4' }, ...
    {'gr0965', '6' }, ...
    {'gr1070', '1' }, ...
    {'gr1070', '2' }, ...
    {'gr1070', '3' }, ...
    {'gr1070', '4' }, ...
    {'blugr1612', '4' }, ... %both cells 4 and 5 seem crappy for STRFs, so wondering whether should leave 'em in the database...
    {'blugr1612', '5' }, ... %both cells 4 and 5 seem crappy for STRFs, so wondering whether should leave 'em in the database...
    {'pupi0414', '6' }, ...
    {'pupi0414', '7' }, ...
    {'glb5202', '9' }, ...
    {'obla1305', '4' }, ...
    {'obla1305', '5' }, ... 
    {'obla1305', '6' }, ...
    {'obla1305', '7' }, ...
    {'blahp9904', '6' }, ...
    {'gr0869', '8' }, ...
    {'gr0869', '9' }, ...
    {'gr0965', '1' }, ...
    {'gr0965', '2' }, ...
    {'gr0965', '5' }, ...
    {'gr1070', '5' }, ...
    {'blugr1612', '3' }, ...
    {'pupi0414', '5' }, ... %kicking in trial 14, but should be okay to run...
    {'blahp9904', '7' }, ...
    {'blahp9904', '8' }, ...
    {'blahp9904', '9' }, ...
    {'obla17102', '4' }, ...
    {'obla17102', '5' }, ...
    {'obla17102', '6' }, ...
    {'obla17102', '8' }, ...
    {'obla17102', '9' }, ...
    {'obla17102', '10' }, ...
    ...%{'gr0869', '7' }, ... %this one only has 11 trials, so didn't STRF, so probably not in database...
    {'gr0869', '10' }, ...
    {'gr0869', '11' }, ...
    {'gr0965', '3' }, ...
    {'gr0965', '7' }, ...
    {'blugr1612', '2' }}; 


% Allocate space for output values
number_cells = length(ov_allunits);

nfiles_con = zeros(1, number_cells);
zscore_con = zeros(1, number_cells);
stdz_con = zeros(1, number_cells);
dprime_con = zeros(1, number_cells);
stdd_con = zeros(1, number_cells);
p_con = zeros(1, number_cells);
rank_con = zeros(1, number_cells);
info_con =  zeros(1, number_cells);
stdinfo_con = zeros(1, number_cells);

% The next three loops could be merged into one but are left separate for
% clarity.
% First conspecific

% Then bengalese song
stimtype = 'conspecific';
for nc=1:number_cells
    [nfiles_con(nc) zscore_con(nc) stdz_con(nc) dprime_con(nc) stdd_con(nc) p_con(nc) rank_con(nc) info_con(nc) stdinfo_con(nc)] = ...
        neural_discrimination(ov_allunits{nc}{1}, brainregion, ov_allunits{nc}{2}, stimtype);
end

% Print the results to a text file
fname = 'Ov_conspecific.txt';
fid = fopen(fname, 'w');
fprintf(fid,'Bird\tCell\tNstim zf\tZscore zf\tZstd zf\tDprime Zf\tDstd zf\tDstd bf\tP zf\tRank zf\tInfo zf\tIstd zf\n');
for nc=1:number_cells
    fprintf(fid,'%s\t%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', ...
        ov_allunits{nc}{1}, ov_allunits{nc}{2}, ...
         nfiles_con(nc), ...
        zscore_con(nc),  stdz_con(nc), ...
         dprime_con(nc),  stdz_con(nc), ...
        p_con(nc),  rank_con(nc), ...
         info_con(nc),  stdinfo_con(nc));
end

fclose(fid);
