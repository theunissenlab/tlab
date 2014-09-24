function [pcc_diag, pcc_cat, mi_tot, mi_error, mi_error2, mi_tot2]=info_matrix_perarea(prob_matrix, cat)
% Calculates the mutual information for a matrix of probabilities within the
% diagonal, the category, and the error. cat is a cell of vectors, each
% vector contain the indices of the element belonging to the same category
% in the matrix.

% Check for sum = 1
sump = sum(sum(prob_matrix));
if ( sump < 0.9999 || sump > 1.0001)
    fprintf(1,'Error in info_matrix: input matrix sums to %f\n', sump);
end

% calculate mi for the whole matrix

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

mi_tot = row_ent + col_ent - tot_ent;

% calculate mi for the diagonal
Diag_mat = diag(prob_matrix);
zero_ind = Diag_mat == 0;
Diag_mat_for_entropy = Diag_mat;
Diag_mat_for_entropy(zero_ind) = 1; % Set 0 to 1 to insure that xlogx goes to zero as x goes to zero
mi_diag = sum(Diag_mat_for_entropy.*log2(Diag_mat_for_entropy./(row_prob'.*col_prob)));

% calculate mi for the categories
% loop through the categories
Ncat = length(cat);
mi_cat = 0;
for nc = 1:Ncat
    icat = cat{nc};
    ni = length(icat);
    for a = 1:ni-1
        aa = icat(a);
        for p = (a+1):ni
            pp = icat(p);
            prob_aapp = prob_matrix_for_entropy(aa,pp);
            prob_ppaa = prob_matrix_for_entropy(pp,aa);
            mi_ap = prob_aapp * (log2(prob_aapp/(row_prob(aa)*col_prob(pp))));
            mi_pa = prob_ppaa * (log2(prob_ppaa/(row_prob(pp)*col_prob(aa))));
            mi_cat = mi_cat + mi_pa + mi_ap;
        end
    end
end

% calculate MI when i different from j and both are not in the same
% category
nstim = size(prob_matrix,1);
Outside_cat=1;
mi_error = 0;
for ii = 1:(nstim-1)
    for jj = (ii+1) : nstim
        for nc = 1:Ncat
            icat = cat{nc};
            if ~isempty(intersect(ii, icat)) && ~isempty(intersect(jj, icat))
                Outside_cat=0;
                break
            end
        end
        if Outside_cat==1;
            prob_iijj = prob_matrix_for_entropy(ii,jj);
            prob_jjii = prob_matrix_for_entropy(jj,ii);
            mi_ij = prob_iijj * (log2(prob_iijj/(row_prob(ii)*col_prob(jj))));
            mi_ji = prob_jjii * (log2(prob_jjii/(row_prob(jj)*col_prob(ii))));
            mi_error = mi_error + mi_ij + mi_ji;
        end
    end
end

% deduct mi for the error (when actual and predicted don't belong to the
% same category
mi_error2 = mi_tot - mi_cat - mi_diag;

% alternative way to calculate mi_tot
mi_tot2 = 0;
for aa = 1:nstim
    for pp = 1:nstim
        prob_aapp = prob_matrix_for_entropy(aa,pp);
        mi_ap = prob_aapp * (log2(prob_aapp/(row_prob(aa)*col_prob(pp))));
        mi_tot2 = mi_tot2 + mi_ap;
    end
end
       

return