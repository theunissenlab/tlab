function test_gauss_vs_complex()

    increment = 24;
    winLength = 186;
    
    [input, samprate, depth] = wavread('/auto/fdata/mschachter/data/channing/clean/stim5.wav');
    
    [s1, to1, fo1, pg1] = GaussianSpectrum(input, increment, winLength, samprate);
    
    doFFTShift = 0;
    [s2, fo2, pg2] = ComplexSpectrum(input, 25, 190, 25000, doFFTShift);
    
    s1Size = size(s1)
    s2Size = size(s2)
    
    figure; hold on;
    subplot(2, 1, 1); hold on;
    imagesc(abs(s1)); axis tight; colorbar;    
    title('Gaussian');
    
    subplot(2, 1, 2); hold on;
    imagesc(abs(s2)); axis tight; colorbar;
    title('Complex');
