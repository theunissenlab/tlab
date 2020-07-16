%This .m file is to calculate pairwise similarity of each STRF in the wn
%collection.
clear all;

addpath('/auto/fhome/pgill/clustering_STRFs');

load used_names;
load strfs;

differences = ones(length(used_names));
Y = [];
for jj = 1:length(used_names)
     %disp(['Doing ' num2str(length((jj+1):length(used_names))) 'comparisons with cell ' used_names{jj} '.']);
     for ii = (jj+1):length(used_names)
         strf1 = strfs{jj};
         strf2 = strfs{ii};
         the_diff = get_sim_score2(strf1,strf2);
         differences(ii,jj) = the_diff;
         differences(jj,ii) = the_diff;
         Y = [Y the_diff];
     end
end

save Y;