function [filter_fit,params] = fit_filter(calc_dir,do_plot,sparsest,verbose)

%% Input handling
if nargin < 2
    do_plot=0;
end

if nargin < 3
    sparsest = 0;
end

if nargin < 4
    verbose = 0;
end

%% Load all data
g = load(fullfile(calc_dir,'strflab_result.mat'));

% Get STRF
if sparsest
    % Get sparsest STRF from the family
    forward = sparsest_strf(g);
else
    if isfield(g,'modelParamsDF')
        forward = g.modelParamsDF.w1;
    else
        forward = g.modelParams.w1;
    end
end

[nBands,nt] = size(forward);
t = 0:nt-1;

% Get band normalization constants and re-normalize STRF before computing
band_stds = g.stimStds;
forward = bsxfun(@times,forward,band_stds');

% Get globals

p = g.preprocParams.tfParams;
samprate = p.samplerate;
f_range = p.high_freq - p.low_freq;

f_band = f_range/(nBands - 1);
f = (0:nBands-1) * f_band + p.low_freq;

%% Initialize output parameters
params.bw = f_band;
params.samprate = samprate;



%% Fit bandwidths

[cf bw forward_fpeak] = fit_cf_bw(f, forward);
filter_fit.cf = cf./1000;  % In kHz
filter_fit.bw = bw./1000;  % In kHz
filter_fit.lat = fit_lat(t, forward);
[bmf fpow pow] = calc_tms(forward, samprate);
filter_fit.bmf = bmf; % In Hz

%% Separate into excitatory and inhibitory
forward_exc = forward;
forward_exc(find(forward<0)) = 0.0;
[cf bw forward_fpeak] = fit_cf_bw(f, forward_exc);
filter_fit.ecf = cf./1000.0;
filter_fit.ebw = bw./1000.0;
filter_fit.elat = fit_lat(t, forward_exc);

forward_inh = forward;
forward_inh(find(forward>0)) = 0.0;
[cf bw forward_fpeak] = fit_cf_bw(f, forward_inh);
filter_fit.icf = cf./1000.0;
filter_fit.ibw = bw./1000.0;
filter_fit.ilat = fit_lat(t, forward_inh);

%% Fourier Analysis of Filters
[p_forward dwt dwf best_dwt best_dwf cm_dwt cm_dwf quad_assym sep_q1 sep_q2]= filter_transfer(forward, samprate, f_band);
filter_fit.mscmt = cm_dwt;
filter_fit.mscms = cm_dwf*1000.0;
filter_fit.msbmf = best_dwt;
filter_fit.msbms = best_dwf*1000;
filter_fit.quad_assym = quad_assym;
filter_fit.sep_q1 = sep_q1;
filter_fit.sep_q2 = sep_q2;

%% Separate filter into temporal and spectral functions and fit with gabor

[sep_index, forward_sep, tx, fy, k] = calc_svd(forward);

[A, X0, Sig, Bw, P, fit_tx, ssegabor] = fit_gabor(tx);
filter_fit.t_A = A;
filter_fit.t_X0 = X0;
filter_fit.t_Sig = Sig;
filter_fit.t_Bw = Bw;
filter_fit.t_P = rem(P,2*pi);
[A, X0, Bw, fit_tx_guassian, ssegaussian] = fit_gaussian(tx);
if (ssegaussian < ssegabor )
    filter_fit.t_A = A;
    filter_fit.t_X0 = X0;
    filter_fit.t_Sig = 0;
    filter_fit.t_Bw = Bw;
    filter_fit.t_P = 0;
    fit_tx = fit_tx_guassian;
end

[A, X0, Sig, Bw, P, fit_fy, ssegabor] = fit_gabor(fy);
filter_fit.f_A = A;
filter_fit.f_X0 = X0;
filter_fit.f_Sig = Sig./(f_band/1000.0);
filter_fit.f_Bw = Bw.*(f_band/1000.0);   % Bandwidth in kHz
filter_fit.f_P = rem(P,2*pi);
[A, X0, Bw, fit_fy_guassian, ssegaussian] = fit_gaussian(fy);

if (ssegaussian < ssegabor )
    filter_fit.f_A = A;
    filter_fit.f_X0 = X0;
    filter_fit.f_Sig = 0;
    filter_fit.f_Bw = Bw.*(f_band/1000.0);
    filter_fit.f_P = 0;
    fit_fy = fit_fy_guassian;
end

filter_fit.f_X0 = (p.low_freq + filter_fit.f_X0 * f_band)/1000.0;

nt = length(tx);
nf = length(fy);

forward_sep_fit = zeros(size(forward));
for i=1:nt
    for j=1:nf
        forward_sep_fit(j,i) = fit_fy(j)*fit_tx(i);
    end
end

%% Do plots if requested
if do_plot

    % Original STRF
    figure(1);
    colormap('Jet');
    imagesc(t, f./1000.0, forward);
    axis xy;
    maxabsforward = max(max(abs(forward)));
    caxis([-maxabsforward maxabsforward]);
    
    % MTF
    figure(2);
    pcolor(dwt, 1000.*dwf, p_forward);
    axis xy;
    axis([-100 100 0 2]);
    shading('interp');
    cmap = colormap;
    for i=1:10
        cmap(i,:) = [1 1 1];
    end
    colormap(cmap);
    
    % Separable fit with marginals
    figure(3);
    subplot(2,2,1)
    imagesc(t, f./1000.0, forward_sep);
    caxis([-maxabsforward maxabsforward]);
    axis xy;
    title('Separable');
    
    subplot(2,2,3);
    plot(t, tx, 'k');
    hold on;
    plot(t, fit_tx, 'k--');
    hold off;

    subplot(2,2,2);
    plot(fy, f./1000, 'k');
    hold on;
    plot(fit_fy, f./1000, 'k');
    hold off;

    % Gabor fit
    figure(4);
    imagesc(t, f./1000.0, forward_sep_fit);
    maxabssepfit = max(abs(max(forward_sep_fit(:))),abs(min(forward_sep_fit(:))))
    caxis([-maxabssepfit,maxabssepfit]);
    axis xy;
    title('Gabor fit');
end

%% Verbose fit results output
if verbose

    % Print output crap
    if ( filter_fit.t_Sig ~= 0 )
        fprintf(1, 'SigT=%f(ms) BWT=%f(ms)', 1/filter_fit.t_Sig, filter_fit.t_Bw);
    else
        fprintf(1, 'SigT=%f(ms) BWT=%f(ms)\n',0, filter_fit.t_Bw);
    end

    if ( filter_fit.f_Sig ~= 0 )
            fprintf(1, '\t SigF=%f(kHZ) BWF=%f(kHz)\n', ...
               1/filter_fit.f_Sig, filter_fit.f_Bw );
    else
            fprintf(1, '\t SigF=%f(kHZ) BWF=%f(kHz)\n', ...
               0, filter_fit.f_Bw );
    end
    fprintf(1, 'Cf = %f (kHz) bw = %f (kHz) lat = %f (ms) BMF = %f (Hz)\n', filter_fit.cf,  filter_fit.bw, ...
                                                      filter_fit.lat, filter_fit.bmf);
    fprintf(1, 'ECf = %f (kHz) Ebw = %f (kHz) Elat = %f (ms)\n',   filter_fit.ecf, filter_fit.ebw, ...
                                                  filter_fit.elat);
    fprintf(1, 'ICf = %f (kHz) Ibw = %f (kHz) Ilat = %f (ms)\n',   filter_fit.icf, filter_fit.ibw, ...
                                                  filter_fit.ilat);

    fprintf(1, 'BMF = %f (Hz) BMS = %f (cycles/kHz) Qassym = %f\n',   filter_fit.msbmf, filter_fit.msbms, ...
                                                    filter_fit.quad_assym);
end