function z=reformPC(x,time);
%x=x+mean(big_data,1)';
if size(x)~=2*time;
    disp('error');
else
    z=[(1:time)',x(1:time),x(time+1:end)];
end;