not_done = 1;
n_loop = 0;
njobs = 0;
% Add jobs to queue until done
    the_command = sprintf('cd(''/auto/fhome/smunro/matlab/Neural Discrimination/''); Normal_L_std');
    the_comment = sprintf('');
    njobs = njobs +1;
    qid_submitted(njobs) = dbaddqueuemaster(the_command, the_comment);
% 
% % Clean the queue
% for ij=1:njobs
%     silent_dbdeletequeue(qid_submitted(ij));
% end