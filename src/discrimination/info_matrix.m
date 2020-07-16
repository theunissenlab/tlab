function mi=info_matrix(prob_matrix)
% Calculates the mutual information for a prob of matrices
% This mutual information is not corrected for low probabilities.
% Correction is Panzeri et al.

% Check for sum = 1
sump = sum(sum(prob_matrix));
if ( sump < 0.9999 || sump > 1.0001)
    fprintf(1,'Error in info_matrix: input matrix sums to %f\n', sump);
end

zero_ind = prob_matrix == 0;
prob_matrix_for_entropy = prob_matrix;
prob_matrix_for_entropy(zero_ind) = 1;                   % Set 0 to 1 to insure that xlogx goes to zero as x goes to zero
tot_ent = sum(sum(-prob_matrix_for_entropy.*log2(prob_matrix_for_entropy)));

row_prob = sum(prob_matrix, 1);
zero_ind = row_prob == 0;
row_prob(zero_ind)= 1;  % Set 0 to 1 to insure that xlogx goes to zero as x goes to zero
row_ent = sum(-row_prob.*log2(row_prob));

col_prob = sum(prob_matrix, 2);
zero_ind = col_prob == 0;
col_prob(zero_ind) = 1;  % Set 0 to 1 to insure that xlogx goes to zero as x goes to zero
col_ent = sum(-col_prob.*log2(col_prob));

mi = row_ent + col_ent - tot_ent;

return