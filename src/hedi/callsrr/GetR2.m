function r2=GetR2(y,pred);

% Compute the R2 between data (y) and prediction 
% Very simple calculation

%r2=1-mean((y-pred).^2)/var(y);
r=regstats(y,pred,'linear','rsquare');
r2=r.rsquare;