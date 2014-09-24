function [ mi_tot, mi_tot2, mi_diag, mi_error, mi_diag_uni, mi_all_error_uni, mi_diag_uni_cat, mi_real_error_uni]=info_matrix_perarea(prob_matrix, cat,FigFlag)
% Calculates the mutual information for a matrix of probabilities within the
% diagonal, the category, and the error. cat is a cell of vectors, each
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

% calculate mi for the whole matrix

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

%% calculate mi for the uniform matrix with diagonal + categories intact
% Generate a matrix with uniform error probabilities within each row for stim that are not in the same category only.
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

mi_diag_uni_cat = row_ent_uni_cat + col_ent_uni_cat - tot_ent_uni_cat;

mi_real_error_uni = mi_tot - mi_diag_uni_cat;

%Plot the MI_UniDiagCat matrix and the position of the cell in the new MI
%space
if FigFlag==1
    Clim=[0 max(max(prob_matrix))];
    createconfusionmatrix(prob_matrix_uni_cat, Clim)
    figure(3)
    title(sprintf('MI UniDiagCat Matrix=%f/nProportion MI UniDiagCat=%f', mi_diag_uni_cat, mi_diag_uni_cat/mi_tot));
    figure(4)
    plot(mi_tot, mi_diag_uni_cat/mi_tot,'ko', 'MarkerFaceColor','k');
    xlabel('Mutual Information of the call matrix mi_tot (bits)')
    ylabel('Proportion of MI in diagonal and category')
    axis([0 6 0 1])
end

%% alternative way to calculate mi_tot
mi_tot2 = 0;
mi_diag = 0;
for aa = 1:nstim
    for pp = 1:nstim
        prob_aapp = prob_matrix(aa,pp);
        if prob_aapp == 0   % this is a flag for zero probabilities
            mi_ap = 0;
        else
            mi_ap = prob_aapp * (log2(prob_aapp/(row_prob(aa)*col_prob(pp))));   % Watch out for row and col here!
        end
        mi_tot2 = mi_tot2 + mi_ap;
        if aa == pp
            mi_diag = mi_diag + mi_ap;
        end
    end
end

if (abs(mi_tot2-mi_tot) > 10^-5 )
    fprintf(1, 'There is something really wrong in the land of MI...');
end
    
mi_error = mi_tot - mi_diag;


%% OLD CODE

% % calculate mi for the categories
% % loop through the categories
% Ncat = length(cat);
% mi_cat = 0;
% for nc = 1:Ncat
%     icat = cat{nc};
%     ni = length(icat);
%     for a = 1:ni-1
%         aa = icat(a);
%         for p = (a+1):ni
%             pp = icat(p);
%             prob_aapp = prob_matrix_for_entropy(aa,pp);
%             prob_ppaa = prob_matrix_for_entropy(pp,aa);
%             mi_ap = prob_aapp * (log2(prob_aapp/(row_prob(aa)*col_prob(pp))));
%             mi_pa = prob_ppaa * (log2(prob_ppaa/(row_prob(pp)*col_prob(aa))));
%             mi_cat = mi_cat + mi_pa + mi_ap;
%         end
%     end
% end
% 
% % calculate MI when i different from j and both are not in the same
% % category
% nstim = size(prob_matrix,1);
% Outside_cat=1;
% mi_error = 0;
% for ii = 1:(nstim-1)
%     for jj = (ii+1) : nstim
%         for nc = 1:Ncat
%             icat = cat{nc};
%             if ~isempty(intersect(ii, icat)) && ~isempty(intersect(jj, icat))
%                 Outside_cat=0;
%                 break
%             end
%         end
%         if Outside_cat==1;
%             prob_iijj = prob_matrix_for_entropy(ii,jj);
%             prob_jjii = prob_matrix_for_entropy(jj,ii);
%             mi_ij = prob_iijj * (log2(prob_iijj/(row_prob(ii)*col_prob(jj))));
%             mi_ji = prob_jjii * (log2(prob_jjii/(row_prob(jj)*col_prob(ii))));
%             mi_error = mi_error + mi_ij + mi_ji;
%         end
%     end
% end
% 
% % deduct mi for the error (when actual and predicted don't belong to the
% % same category
% mi_error2 = mi_tot - mi_cat - mi_diag;


       

return