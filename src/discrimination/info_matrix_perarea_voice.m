function [ mi_tot, mi_diag_uni_cat, mi_real_error_uni]=info_matrix_perarea_voice(prob_matrix, CAT,FigFlag)
% Calculates the mutual information for a matrix of probabilities within the
% diagonal, the category, and the error. cat is a cell of cells of vectors, each
% vector contain the indices of the element belonging to the same category
% in the matrix.
if nargin==2
    FigFlag=0;
end

% Check for sum = 1
sump = sum(sum(prob_matrix));
if ( sump < 0.9999 || sump > 1.0001)
    fprintf(1,'Error in info_matrix: input matrix sums to %f\n', sump);
end

%% Get ready some output variables
mi_diag_uni_cat=nan(length(CAT),1);
mi_real_error_uni=nan(length(CAT),1);

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

%Plot the matrix
if FigFlag==1
    Clim=[0 max(max(prob_matrix))];
    createconfusionmatrix(prob_matrix, Clim)
    figure(1)
    title(sprintf('MI whole Matrix=%f',mi_tot));
end

%% calculate mi for the uniform matrix with diagonal + categories intact
% Generate a matrix with uniform error probabilities within each row for stim that are not in the same category only.
for cl=1:length(CAT)
    cat=CAT{cl};
    nstim = size(prob_matrix, 1);
    prob_matrix_uni_cat = prob_matrix;
    Ncat = length(cat);
    for aa = 1:nstim
        for nc = 1:Ncat
                icat = cat{nc};
                if ~isempty(intersect(aa, icat))
                    break
                end
        end
        outcat_pp = setdiff(1:nstim,icat);
        prow_outcat = sum(prob_matrix(aa,outcat_pp));
        puni = prow_outcat./length(outcat_pp);
        prob_matrix_uni_cat(aa,outcat_pp) = repmat(puni, 1, length(outcat_pp));
    end

    % Calculate MI for outside category uniform matrix
    zero_ind = prob_matrix_uni_cat == 0;
    prob_matrix_uni_cat_for_entropy = prob_matrix_uni_cat;
    prob_matrix_uni_cat_for_entropy(zero_ind) = 1;                   % Set 0 to 1 to insure that xlogx goes to zero as x goes to zero
    tot_ent_uni_cat = sum(sum(-prob_matrix_uni_cat_for_entropy.*log2(prob_matrix_uni_cat_for_entropy)));

    col_prob_uni_cat = sum(prob_matrix_uni_cat, 1);
    zero_ind = col_prob_uni_cat == 0;
    col_prob_uni_cat(zero_ind)= 1;  % Set 0 to 1 to insure that xlogx goes to zero as x goes to zero
    col_ent_uni_cat = sum(-col_prob_uni_cat.*log2(col_prob_uni_cat));

    row_prob_uni_cat = sum(prob_matrix_uni_cat, 2);
    zero_ind = row_prob_uni_cat == 0;
    row_prob_uni_cat(zero_ind) = 1;  % Set 0 to 1 to insure that xlogx goes to zero as x goes to zero
    row_ent_uni_cat = sum(-row_prob_uni_cat.*log2(row_prob_uni_cat));

    mi_diag_uni_cat(cl) = row_ent_uni_cat + col_ent_uni_cat - tot_ent_uni_cat;

    mi_real_error_uni(cl) = mi_tot - mi_diag_uni_cat(cl);

    %Plot the MI_UniDiagCat matrix and the position of the cell in the new MI
    %space
    if FigFlag==1
        Clim=[0 max(max(prob_matrix))];
        createconfusionmatrix(prob_matrix_uni_cat, Clim)
        figure(3)
        title(sprintf('MI UniDiagCat Matrix=%f/nProportion MI UniDiagCat=%f', mi_diag_uni_cat(cl), mi_diag_uni_cat(cl)/mi_tot));
        figure(4)
        plot(mi_tot, mi_diag_uni_cat(cl)/mi_tot,'ko', 'MarkerFaceColor','k');
        xlabel('Mutual Information of the call matrix mi_tot (bits)')
        ylabel('Proportion of MI in diagonal and category')
        axis([0 6 0 1])
        pause
    end
end
       

return