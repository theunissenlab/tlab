function out = faster_surprise_cond_p_7D_to_output(cond_p,which_cubes,what_position,use_stim,lengths,stim_names)
%  Takes the conditional probability matrix cond_p and the 7 X N arrays
%  what_cubes and what_position and returns the negative
%  log of the conditional probability of the target being what it is in the
%  following way.  First the conditional probability of that particular
%  response is determined, then it is determined if the target stimulus is
%  higher or lower than the mode stimulus expected.  If it's higher, than
%  the corresponding bin of higher_than_pred is set to -log(cond_p) +
%  log(p_max) and the bin of lower_than_pred is set to 0, where p_max is
%  the maximum probability of any particular target value given these
%  stimulus conditions.

higher_than_pred = zeros(1,size(which_cubes,2));
lower_than_pred = higher_than_pred;
use_max = 1;  %Whether to normalize my S_ml, the most likely stim

total = which_cubes + what_position;
%VI = interpn(1:15,1:4,1:6,1:6,1:6,1:6,1:15,cond_p,total(1,:),total(2,:),total(3,:),total(4,:),total(5,:),total(6,:),total(7,:));
VI = interpn(cond_p,total(1,:),total(2,:),total(3,:),total(4,:),total(5,:),total(6,:),total(7,:));
maxProb = zeros(size(cond_p,1),size(total,2));
for jj = 1:size(maxProb,1); 
    maxProb(jj,:) = interpn(cond_p,jj*ones(1,size(maxProb,2)),total(2,:),total(3,:),total(4,:),total(5,:),total(6,:),total(7,:)); 
end
[aOut,bOut] = max(maxProb);
isHigher = total(1,:)>=bOut;
surprise = -log(VI) + use_max*log(aOut);

lower_than_pred = surprise.*(1-isHigher);
higher_than_pred = surprise.*(isHigher);

%{

for t = find(use_stim)
    cube = which_cubes(:,t);
    domain = cond_p(:, cube(2) + (0:1),cube(3) + (0:1),cube(4) + (0:1),cube(5) + (0:1),cube(6) + (0:1),cube(7) + (0:1));
    the_p = recursive_interpolate(what_position(:,t),domain(cube(1) + (0:1),:,:,:,:,:,:,:,:,:,:));
    [the_max, index] = max(max(max(max(max(max(max(domain,[],2),[],3),[],4),[],5),[],6),[],7),[],1);

    if (cube(1)+what_position(1,t) )< index
        lower_than_pred(t) = -log(the_p) + use_max*log(the_max);
    else
        higher_than_pred(t) = -log(the_p) + use_max*log(the_max);
    end

end
%}



out = struct('stim_name',[],'louder',[],'quieter',[]);
start = 1;
for istim = 1:length(stim_names)
    active_index = start:(start -1 + lengths(istim));
    start = start + lengths(istim);
    if use_stim(active_index(1))
        out(end+1).stim_name = stim_names{istim};
        out(end).louder = higher_than_pred(active_index);
        out(end).quieter = lower_than_pred(active_index);
    end
end
out = out(2:end);  %Get rid of the empty template on the out struct.