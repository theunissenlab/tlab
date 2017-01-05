function [ mi_tot, mi_tot_max, mi_diag_uni_cat, mi_real_error_uni, mi_diag_uni_cat_maxInv,mi_maxdiag_uni_cat_maxInv, EntropyValues, mi_diag_uni_cat_minInv]=info_matrix_perarea_ext(prob_matrix, cat,FigFlag, ForFig)
% if you want to restore the calculation of mi_diag_uni function [ mi_tot, mi_tot_max,mi_diag_uni, mi_all_error_uni, mi_diag_uni_cat, mi_real_error_uni, mi_diag_uni_cat_maxInv, mi_diag_uni_cat_minInv]=info_matrix_perarea_ext(prob_matrix, cat,FigFlag)
% Calculates the mutual information for a matrix of probabilities within the
% diagonal, the category, and the error. cat is a cell of vectors, each
% vector contain the indices of the element belonging to the same category
% in the matrix.
if nargin<3
    FigFlag=0;
    ForFig={};
end
if nargin==3
    fprintf('WARNING: please specify input for figure!!!\n')
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
    createconfusionmatrix(prob_matrix, Clim, ForFig.MIN, ForFig.MAX, ForFig.UVT)
    figure(1)
    title(sprintf('MI whole Matrix=%f',mi_tot));
end

% %% calculate mi for the uniform matrix with diagonal intact
% % Generate a matrix with uniform error probabilities within each row.
% nstim = size(prob_matrix, 1);
% prob_matrix_uni = zeros(size(prob_matrix));
% for aa = 1:nstim
%     pdiag = prob_matrix(aa,aa);
%     prow = sum(prob_matrix(aa,:));
%     puni = (prow - pdiag)./(nstim-1);
%     prob_matrix_uni(aa,:) = repmat(puni, 1, nstim);
%     prob_matrix_uni(aa,aa) = pdiag;
% end
% 
% % Calculate mi for uniform matrix
% zero_ind = prob_matrix_uni == 0;
% prob_matrix_uni_for_entropy = prob_matrix_uni;
% prob_matrix_uni_for_entropy(zero_ind) = 1;                   % Set 0 to 1 to insure that xlogx goes to zero as x goes to zero
% tot_ent_uni = sum(sum(-prob_matrix_uni_for_entropy.*log2(prob_matrix_uni_for_entropy)));
% 
% col_prob_uni = sum(prob_matrix_uni, 1);
% zero_ind = col_prob_uni == 0;
% col_prob_uni(zero_ind)= 1;  % Set 0 to 1 to insure that xlogx goes to zero as x goes to zero
% col_ent_uni = sum(-col_prob_uni.*log2(col_prob_uni));
% 
% row_prob_uni = sum(prob_matrix_uni, 2);
% zero_ind = row_prob_uni == 0;
% row_prob_uni(zero_ind) = 1;  % Set 0 to 1 to insure that xlogx goes to zero as x goes to zero
% row_ent_uni = sum(-row_prob_uni.*log2(row_prob_uni));
% 
% mi_diag_uni = row_ent_uni + col_ent_uni - tot_ent_uni;
% 
% mi_all_error_uni = mi_tot - mi_diag_uni;
% 
% %plot the UniDiag Matrix
% if FigFlag==1
%     Clim=[0 max(max(prob_matrix))];
%     createconfusionmatrix(prob_matrix_uni, Clim)
%     figure(2)
%     title(sprintf('MI UniDiag Matrix=%f',mi_diag_uni));
% end

%% calculate Max mi for the whole matrix given the sample size
% Generate a matrix with uniform 0 error probabilities within each row.
nstim = size(prob_matrix, 1);
prob_matrix_max = zeros(size(prob_matrix));
for aa = 1:nstim
    prob_matrix_max(aa,aa) = sum(prob_matrix(aa,:));
end

% Calculate mi for uniform matrix
zero_ind = prob_matrix_max == 0;
prob_matrix_max_for_entropy = prob_matrix_max;
prob_matrix_max_for_entropy(zero_ind) = 1;                   % Set 0 to 1 to insure that xlogx goes to zero as x goes to zero
tot_ent_max = sum(sum(-prob_matrix_max_for_entropy.*log2(prob_matrix_max_for_entropy)));

col_prob_max = sum(prob_matrix_max, 1);
zero_ind = col_prob_max == 0;
col_prob_max(zero_ind)= 1;  % Set 0 to 1 to insure that xlogx goes to zero as x goes to zero
col_ent_max = sum(-col_prob_max.*log2(col_prob_max));

row_prob_max = sum(prob_matrix_max, 2);
zero_ind = row_prob_max == 0;
row_prob_max(zero_ind) = 1;  % Set 0 to 1 to insure that xlogx goes to zero as x goes to zero
row_ent_max = sum(-row_prob_max.*log2(row_prob_max));

mi_tot_max = row_ent_max + col_ent_max - tot_ent_max;

%plot the ProbMax Matrix
if FigFlag==1
    Clim=[0 max(max(prob_matrix))];
    createconfusionmatrix(prob_matrix_max, Clim, ForFig.MIN, ForFig.MAX, ForFig.UVT)
    figure(2)
    title(sprintf('MI MaxProba Matrix=%f',mi_tot_max));
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
EntropyValues.mi_diag_uni_cat.col_ent=col_ent_uni_cat;
EntropyValues.mi_diag_uni_cat.row_ent=row_ent_uni_cat;
EntropyValues.mi_diag_uni_cat.tot_ent=tot_ent_uni_cat;

mi_real_error_uni = mi_tot - mi_diag_uni_cat;

%Plot the MI_UniDiagCat matrix and the position of the cell in the new MI
%space
if FigFlag==1
    Clim=[0 max(max(prob_matrix))];
    createconfusionmatrix(prob_matrix_uni_cat, Clim, ForFig.MIN, ForFig.MAX, ForFig.UVT)
    figure(3)
    title(sprintf('Inclusive Categorical Info Matrix=%f\nNormalized ICI=%f', mi_diag_uni_cat, mi_diag_uni_cat/mi_tot_max));
    figure(4)
    plot(mi_tot, mi_diag_uni_cat/mi_tot,'ko', 'MarkerFaceColor','k');
    xlabel('Mutual Information of the call matrix mi_tot (bits)')
    ylabel('Proportion of MI in diagonal and category')
    axis([0 6 0 1])
end

%% calculate mi for the matrix with a uniform distribution of proba outside categories...
... and a uniform distribution of proba inside categories to obtain the max invariant matrix...
    ......or the so called Exclusive Categorical Information
% Generate a matrix with uniform error probabilities within each row for
% stim that are not in the same category and a second value of proba for
% stim in the same category
nstim = size(prob_matrix, 1);
prob_matrix_uni_cat_maxInv = prob_matrix;
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
    prob_matrix_uni_cat_maxInv(aa,outcat_pp) = repmat(puni, 1, length(outcat_pp));
    prow_incat  = sum(prob_matrix(aa, icat));
    puni_cat = prow_incat./length(icat);
    prob_matrix_uni_cat_maxInv(aa,icat) = repmat(puni_cat, 1, length(icat));
end

% Calculate MI for outside category uniform matrix
zero_ind = prob_matrix_uni_cat_maxInv == 0;
prob_matrix_uni_cat_maxInv_for_entropy = prob_matrix_uni_cat_maxInv;
prob_matrix_uni_cat_maxInv_for_entropy(zero_ind) = 1;                   % Set 0 to 1 to insure that xlogx goes to zero as x goes to zero
tot_ent_uni_cat_maxInv = sum(sum(-prob_matrix_uni_cat_maxInv_for_entropy.*log2(prob_matrix_uni_cat_maxInv_for_entropy)));

col_prob_uni_cat_maxInv = sum(prob_matrix_uni_cat_maxInv, 1);
zero_ind = col_prob_uni_cat_maxInv == 0;
col_prob_uni_cat_maxInv(zero_ind)= 1;  % Set 0 to 1 to insure that xlogx goes to zero as x goes to zero
col_ent_uni_cat_maxInv = sum(-col_prob_uni_cat_maxInv.*log2(col_prob_uni_cat_maxInv));

row_prob_uni_cat_maxInv = sum(prob_matrix_uni_cat_maxInv, 2);
zero_ind = row_prob_uni_cat_maxInv == 0;
row_prob_uni_cat_maxInv(zero_ind) = 1;  % Set 0 to 1 to insure that xlogx goes to zero as x goes to zero
row_ent_uni_cat_maxInv = sum(-row_prob_uni_cat_maxInv.*log2(row_prob_uni_cat_maxInv));

mi_diag_uni_cat_maxInv = row_ent_uni_cat_maxInv + col_ent_uni_cat_maxInv - tot_ent_uni_cat_maxInv;

%Plot the MI_UniDiagCat_maxInv matrix
if FigFlag==1
    Clim=[0 max(max(prob_matrix))];
    createconfusionmatrix(prob_matrix_uni_cat_maxInv, Clim, ForFig.MIN, ForFig.MAX, ForFig.UVT)
    figure(5)
    title(sprintf('Exclusive Categorical Information Matrix=%f\n', mi_diag_uni_cat_maxInv));
end

%% calculate mi for the Max proba matrix with a uniform distribution of proba outside categories...
... and a uniform distribution of proba inside categories to obtain the max invariant matrix...
    ...or the so called maximum value of Exclusive Categorical Information given the structure fo the dataset
% Generate a matrix with uniform error probabilities within each row for
% stim that are not in the same category and a second value of proba for
% stim in the same category
nstim = size(prob_matrix_max, 1);
prob_maxmatrix_uni_cat_maxInv = prob_matrix_max;
Ncat = length(cat);
for aa = 1:nstim
    for nc = 1:Ncat
            icat = cat{nc};
            if ~isempty(intersect(aa, icat))
                break
            end
    end
    outcat_pp = setdiff(1:nstim,icat);
    prow_outcat = sum(prob_matrix_max(aa,outcat_pp));
    puni = prow_outcat./length(outcat_pp);
    prob_maxmatrix_uni_cat_maxInv(aa,outcat_pp) = repmat(puni, 1, length(outcat_pp));
    prow_incat  = sum(prob_matrix_max(aa, icat));
    puni_cat = prow_incat./length(icat);
    prob_maxmatrix_uni_cat_maxInv(aa,icat) = repmat(puni_cat, 1, length(icat));
end

% Calculate MI for outside category uniform matrix
zero_ind = prob_maxmatrix_uni_cat_maxInv == 0;
prob_maxmatrix_uni_cat_maxInv_for_entropy = prob_maxmatrix_uni_cat_maxInv;
prob_maxmatrix_uni_cat_maxInv_for_entropy(zero_ind) = 1;                   % Set 0 to 1 to insure that xlogx goes to zero as x goes to zero
tot_maxent_uni_cat_maxInv = sum(sum(-prob_maxmatrix_uni_cat_maxInv_for_entropy.*log2(prob_maxmatrix_uni_cat_maxInv_for_entropy)));

col_maxprob_uni_cat_maxInv = sum(prob_maxmatrix_uni_cat_maxInv, 1);
zero_ind = col_prob_uni_cat_maxInv == 0;
col_maxprob_uni_cat_maxInv(zero_ind)= 1;  % Set 0 to 1 to insure that xlogx goes to zero as x goes to zero
col_maxent_uni_cat_maxInv = sum(-col_maxprob_uni_cat_maxInv.*log2(col_maxprob_uni_cat_maxInv));

row_maxprob_uni_cat_maxInv = sum(prob_maxmatrix_uni_cat_maxInv, 2);
zero_ind = row_maxprob_uni_cat_maxInv == 0;
row_maxprob_uni_cat_maxInv(zero_ind) = 1;  % Set 0 to 1 to insure that xlogx goes to zero as x goes to zero
row_maxent_uni_cat_maxInv = sum(-row_maxprob_uni_cat_maxInv.*log2(row_maxprob_uni_cat_maxInv));

mi_maxdiag_uni_cat_maxInv = row_maxent_uni_cat_maxInv + col_maxent_uni_cat_maxInv - tot_maxent_uni_cat_maxInv;

%Plot the MI_UniDiagCat_maxInv matrix
if FigFlag==1
    Clim=[0 max(max(prob_matrix_max))];
    createconfusionmatrix(prob_maxmatrix_uni_cat_maxInv, Clim, ForFig.MIN, ForFig.MAX, ForFig.UVT)
    figure(6)
    title(sprintf('Max value of Exclusive Categorical Information=%f\nNormalized ECI=%f\n', mi_maxdiag_uni_cat_maxInv,mi_diag_uni_cat_maxInv./mi_maxdiag_uni_cat_maxInv));
end

%% calculate mi for the matrix with a uniform distribution of proba outside categories...
... and a diagonal distribution of proba inside categories to obtain the min invariant matrix
% Generate a matrix with uniform error probabilities within each row for
% stim that are not in the same category and sum the value of proba for
% stim in the same category to attribute that value in the diagonal
nstim = size(prob_matrix, 1);
prob_matrix_uni_cat_minInv = prob_matrix;
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
    prob_matrix_uni_cat_minInv(aa,outcat_pp) = repmat(puni, 1, length(outcat_pp));
    prow_incat  = sum(prob_matrix(aa, icat));
    prob_matrix_uni_cat_minInv(aa,icat) = zeros(length(icat),1);
    prob_matrix_uni_cat_minInv(aa,aa)= prow_incat;
end

% Calculate MI for outside category uniform matrix
zero_ind = prob_matrix_uni_cat_minInv == 0;
prob_matrix_uni_cat_minInv_for_entropy = prob_matrix_uni_cat_minInv;
prob_matrix_uni_cat_minInv_for_entropy(zero_ind) = 1;                   % Set 0 to 1 to insure that xlogx goes to zero as x goes to zero
tot_ent_uni_cat_minInv = sum(sum(-prob_matrix_uni_cat_minInv_for_entropy.*log2(prob_matrix_uni_cat_minInv_for_entropy)));

col_prob_uni_cat_minInv = sum(prob_matrix_uni_cat_minInv, 1);
zero_ind = col_prob_uni_cat_minInv == 0;
col_prob_uni_cat_minInv(zero_ind)= 1;  % Set 0 to 1 to insure that xlogx goes to zero as x goes to zero
col_ent_uni_cat_minInv = sum(-col_prob_uni_cat_minInv.*log2(col_prob_uni_cat_minInv));

row_prob_uni_cat_minInv = sum(prob_matrix_uni_cat_minInv, 2);
zero_ind = row_prob_uni_cat_minInv == 0;
row_prob_uni_cat_minInv(zero_ind) = 1;  % Set 0 to 1 to insure that xlogx goes to zero as x goes to zero
row_ent_uni_cat_minInv = sum(-row_prob_uni_cat_minInv.*log2(row_prob_uni_cat_minInv));

mi_diag_uni_cat_minInv = row_ent_uni_cat_minInv + col_ent_uni_cat_minInv - tot_ent_uni_cat_minInv;

%Plot the MI_UniDiagCat_maxInv matrix
if FigFlag==1
    Clim=[0 max(max(prob_matrix))];
    createconfusionmatrix(prob_matrix_uni_cat_minInv, Clim, ForFig.MIN, ForFig.MAX, ForFig.UVT)
    figure(7)
    title(sprintf('MI Min invariant matrix within category=%f\n', mi_diag_uni_cat_minInv));
end


% %% alternative way to calculate mi_tot
% mi_tot2 = 0;
% mi_diag = 0;
% for aa = 1:nstim
%     for pp = 1:nstim
%         prob_aapp = prob_matrix(aa,pp);
%         if prob_aapp == 0   % this is a flag for zero probabilities
%             mi_ap = 0;
%         else
%             mi_ap = prob_aapp * (log2(prob_aapp/(row_prob(aa)*col_prob(pp))));   % Watch out for row and col here!
%         end
%         mi_tot2 = mi_tot2 + mi_ap;
%         if aa == pp
%             mi_diag = mi_diag + mi_ap;
%         end
%     end
% end
% 
% if (abs(mi_tot2-mi_tot) > 10^-5 )
%     fprintf(1, 'There is something really wrong in the land of MI...');
% end
%     
% mi_error = mi_tot - mi_diag;


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