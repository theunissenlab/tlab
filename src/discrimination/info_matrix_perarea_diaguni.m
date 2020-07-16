function [ mi_tot, mi_diag_uni, mi_all_error_uni]=info_matrix_perarea_diaguni(prob_matrix,FigFlag)
% Calculates the mutual information for a matrix of probabilities within the
% diagonal, the category, and the error. cat is a cell of vectors, each
% vector contain the indices of the element belonging to the same category
% in the matrix.
if nargin==1
    FigFlag=0;
end

% Check for sum = 1
sump = sum(sum(prob_matrix));
if ( sump < 0.9999 || sump > 1.0001)
    fprintf(1,'Error in info_matrix: input matrix sums to %f\n', sump);
end

%% calculate mi for the whole matrix

zero_ind = prob_matrix == 0;
prob_matrix_for_entropy = prob_matrix;
prob_matrix_for_entropy(zero_ind) = 1;                   % Set 0 to 1 to insure that xlogx goes to zero as x goes to zero
tot_ent = sum(sum(-prob_matrix_for_entropy.*log2(prob_matrix_for_entropy)));

col_prob = sum(prob_matrix, 1);
zero_ind = col_prob == 0;
col_prob(zero_ind)= 1;  % Set 0 to 1 to insure that xlogx goes to zero as x goes to zero
col_ent = sum(-col_prob.*log2(col_prob));

row_prob = sum(prob_matrix, 2);
zero_ind = row_prob == 0;
row_prob(zero_ind) = 1;  % Set 0 to 1 to insure that xlogx goes to zero as x goes to zero
row_ent = sum(-row_prob.*log2(row_prob));

mi_tot = row_ent + col_ent - tot_ent;
EntropyValues.mi_tot.col_ent=col_ent;
EntropyValues.mi_tot.row_ent=row_ent;
EntropyValues.mi_tot.tot_ent=tot_ent;

%Plot the matrix
if FigFlag==1
    Clim=[0 max(max(prob_matrix))];
    createconfusionmatrix(prob_matrix, Clim)
    figure(1)
    title(sprintf('MI whole Matrix=%f',mi_tot));
end

%% calculate mi for the uniform matrix with diagonal intact
% Generate a matrix with uniform error probabilities within each row.
nstim = size(prob_matrix, 1);
prob_matrix_uni = zeros(size(prob_matrix));
for aa = 1:nstim
    pdiag = prob_matrix(aa,aa);
    prow = sum(prob_matrix(aa,:));
    puni = (prow - pdiag)./(nstim-1);
    prob_matrix_uni(aa,:) = repmat(puni, 1, nstim);
    prob_matrix_uni(aa,aa) = pdiag;
end

% Calculate mi for uniform matrix
zero_ind = prob_matrix_uni == 0;
prob_matrix_uni_for_entropy = prob_matrix_uni;
prob_matrix_uni_for_entropy(zero_ind) = 1;                   % Set 0 to 1 to insure that xlogx goes to zero as x goes to zero
tot_ent_uni = sum(sum(-prob_matrix_uni_for_entropy.*log2(prob_matrix_uni_for_entropy)));

col_prob_uni = sum(prob_matrix_uni, 1);
zero_ind = col_prob_uni == 0;
col_prob_uni(zero_ind)= 1;  % Set 0 to 1 to insure that xlogx goes to zero as x goes to zero
col_ent_uni = sum(-col_prob_uni.*log2(col_prob_uni));

row_prob_uni = sum(prob_matrix_uni, 2);
zero_ind = row_prob_uni == 0;
row_prob_uni(zero_ind) = 1;  % Set 0 to 1 to insure that xlogx goes to zero as x goes to zero
row_ent_uni = sum(-row_prob_uni.*log2(row_prob_uni));

mi_diag_uni = row_ent_uni + col_ent_uni - tot_ent_uni;

mi_all_error_uni = mi_tot - mi_diag_uni;

%plot the UniDiag Matrix
if FigFlag==1
    Clim=[0 max(max(prob_matrix))];
    createconfusionmatrix(prob_matrix_uni, Clim)
    figure(2)
    title(sprintf('MI UniDiag Matrix=%f',mi_diag_uni));
end


return