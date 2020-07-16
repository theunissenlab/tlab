function y = take_out_repeats(x)
%takes in a vector x and takes out any repeated numbers;
if size(x,1)<size(x,2);
    %if x is a row vector, change it to a column vector;
    x=x';
end

x=sort(x);
index=find(diff(x))+1;
y = [min(x); x(index)];