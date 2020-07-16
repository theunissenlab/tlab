function y=computeCOf(x,f)
% compute the correlation of x and \chi(kf) kf<length x
n=length(x);a=x-mean(x);
y=sum(a(f:f:n));
