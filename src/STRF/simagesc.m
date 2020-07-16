function out = simagesc(in1,in2,in3);
if nargin == 1
    to_plot = in1;
    maxval = max(max(abs(to_plot)));
    out = imagesc(to_plot);
    axis xy
    caxis([-maxval maxval]);
else
    to_plot = in3;
    maxval = max(max(abs(to_plot)));
    out = imagesc(in1,in2,to_plot);
    axis xy
    caxis([-maxval maxval]);
end