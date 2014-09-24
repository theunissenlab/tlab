function [sep_index, forward_sep, tx, fy, k] = calc_svd(forward)
% Calculates the separability index and returns the time frequency separable filter

[u, s, v] = svd(forward);
s_diag = diag(s);
sep_index = s_diag(1).^2/sum(s_diag.^2);
count = 1;
sfilt = s;
nsend = length(s_diag);
for j=count+1:nsend
       sfilt(j,j) = 0;
end
forward_sep = u*sfilt*v';
tx = v(1:end,1);
fy = u(1:end,1);
k = sfilt(1,1);

% Code to check...
% forward_sep_test = zeros(size(forward_sep));
% nt = length(tx);
% nf = length(fy);
% 
% for i=1:nt
%     for j=1:nf
%         forward_sep_test(j,i) = fy(j)*tx(i);
%     end
% end
