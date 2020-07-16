function [PCA_Tmat,E,D] = surprise_preprocess(spect,jj,w,d,b)
%  Function which goes from the spectrogram to a PCA + target format.
%  w is the width of the predictor domain, d is the lag between domain and
%  target, and b is the number of frequency bands away from T the
%  PCA-prediction looks.
%  PCA_Tmat has the form of a 3-D array with dimensions
%  2D+1,sum(stim_durs), n_bands - 2b.  The first dimension is of the form
%  T, eig1, eig2 ... where eig1 is the projection onto the weakest
%  eigenvalue of the array, eig2 is the second weakest, etc.
%  Version 3 does this calculation for only one stimulus band, saving
%  memory so that more stims can be used.

if ~exist('w','var')
    w = 3;
end

if ~exist('d','var')
    d = 3;
end

if~exist('b','var')
    b = 4;
end

first_column = spect(:,1);

to_prepend = first_column*ones(1,d+w);
orig_sz = size(spect);
spect = [to_prepend spect];  %  To make the spectrogram a bit bigget so the lengths will work out.

n_bands = size(spect,1);

PCA_Tmat = zeros(w*(2*b+1) + 1, orig_sz(2),1);

jj = (jj+b);
domain = [];
for kk = 1:w
    domain = [domain ; spect((jj - b):(jj + b), kk:(orig_sz(2)+kk-1))];
end
cov_mat = cov(domain');
[E,D] = eig(cov_mat);
%  Eigenvectors are the columns of E.
PCA_Tmat(:,:,1) = [spect(jj,(1+d+w):end) ; E'*domain];

