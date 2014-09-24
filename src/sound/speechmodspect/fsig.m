function y = fsigplus(x, m)
% Returns a compress function of x. m is the maximun value

k = 10/m;
y = 1./(1+exp(k.*(-x+0.75*m)));

end