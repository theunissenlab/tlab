function wav=ComputeWaveFromFit(spl);

n=length(spl);

for i=1:n;
    waf=spl.wavf{i};
    orig=spl.splits{i};
end;