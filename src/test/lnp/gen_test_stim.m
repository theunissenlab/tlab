function stim = gen_test_stim(nDataPoints, numStimChannels, sampleRate)

sint = 1/sampleRate;
t = 0:sint:(nDataPoints-1)*sint;

freq = 2;

stimMask = repmat(cos(t*2*pi*freq), numStimChannels, 1)';

stimBase = randn(nDataPoints, numStimChannels);

stim = abs(stimMask .* stimBase);
