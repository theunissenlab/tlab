function r_allign_comp = compress_raster(r, cmode)
% function  r_allign_comp = compress_raster(r, cmode);
%
% r is a [Raster trial count] x [SampleSize] matrix
% cmode determines whether continuous trials are in a row or
%   just compress for nonnan values

if ~exist('cmode', 'var')
	cmode = 1;
end

switch cmode
case 0
	r_allign = NaN + zeros(size(r));
	
	rc = 1;
	for ii=1:size(r,1)
		nonnan = find(~isnan(r(ii,:)));
		r_allign(ii,1:length(nonnan)) = r(ii,nonnan);
	end

	for ii=1:size(r_allign,2)
		thisnan = length(find(isnan(r_allign(:,ii))));
		if thisnan > size(r,1)*0.15
			break;
		end
	end
	repcount = ii-1;
	r_allign_comp = r_allign(:,1:repcount);

case 1

	r_allign = NaN + zeros(size(r));
	
	rc = 1;
	for ii=1:size(r,2)
		nonnan = find(~isnan(r(:,ii)));
		r_allign(nonnan,rc) = nanmean([r_allign(nonnan,rc) r(nonnan,ii)], 2);
		if length(find(~isnan(r_allign(:,rc))))/size(r,1) > 0.85 & max(nonnan)>size(r,1)*0.98
			rc = rc + 1;
		end
	end
	
	if length(find(~isnan(r_allign(:,rc))))/size(r,1) < 0.85
		rc = rc - 1;
	end
	repcount = rc;
	r_allign_comp = r_allign(:,1:repcount);

end
