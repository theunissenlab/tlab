function P_est = recursive_interpolate(alpha,domain);
%  Uses a recursive algorithm to find the linear interpolation of a
%  hypercube with values at verticies given by domain, position inside the
%  hypercube given by alpha where all values of alpha are between 0 and 1. 
%  (Extrapolates linearly otherwise; no biggie.)

if any(size(domain) ~=2)
    error(strcat(['Error in function recursive_interpolate: must be given a hypercube; dimensions of input are ' num2str(size(domain)) '.']));
end
if length(size(domain)) ~= length(alpha)
    error(strcat(['Error in function recursive_interpolate: domain is ' num2str(length(size(domain))) '-D while alpha is ' num2str(length(alpha)) '-D.']));
end

lin_dom = domain(:);  %Puts domain in a vector for easier handling.

P_est = internal_recursive_interpolate(alpha,lin_dom);

function out = internal_recursive_interpolate(alpha,lin_dom)
if length(alpha) == 1  %base case
    out = (1-alpha)*lin_dom(1) + alpha*lin_dom(2);
else
    half_length = length(lin_dom)/2;
    out = (1-alpha(end))*internal_recursive_interpolate(alpha(1:(end-1)),lin_dom(1:half_length)) + (alpha(end))*internal_recursive_interpolate(alpha(1:(end-1)),lin_dom((half_length+1):(2*half_length)));
end

