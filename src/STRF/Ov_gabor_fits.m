% Extract all the Gabor filter parameters

clear all;
load 'allOvStrfs.mat'

nstrfs = length(allOvStrfs);
Ov_strf_type = cell(1, nstrfs);

Ov_filter_fit = repmat(struct( 'name', [], 'cf', 0.0, 'bw', 0.0, 'lat', 0, 'bmf', 0, ...
                              'ecf', 0.0, 'ebw', 0.0, 'elat', 0.0, ...
                              'icf', 0.0, 'ibw', 0.0, 'ilat', 0.0, ...
                              'msbmf', 0.0, 'msbms', 0.0, 'quad_assym', 0.0, 'sep_q1', 0.0, 'sep_q2', 0.0, ...
                              't_A', 0.0, 't_X0', 0.0, 't_sig', 0.0, 't_Bw', 0.0, 't_P', 0.0, ...
                              'f_A', 0.0, 'f_X0', 0.0, 'f_sig', 0.0, 'f_Bw', 0.0, 'f_P', 0.0), 1, nstrfs);                      
                                                    
% Get time and Frequency parameters
forward = allOvStrfs{1};
forward = forward(:,101:201);
nt = size(forward, 2);
nf = size(forward, 1);

f=1:nf;
f = f.*125; % frequency in Hz
t=0:nt-1;  % time in ms
samprate = 1000;   % sampling rate is 1 ms

           
for istrf = 1:nstrfs

    forward = allOvStrfs{istrf};
    forward = forward(:,101:201);
    % Extract cf, bw, lat 
    [cf bw forward_fpeak] = fit_cf_bw(f, forward);
    Ov_filter_fit(istrf).cf = cf./1000;  % In kHz
    Ov_filter_fit(istrf).bw = bw./1000;  % In kHz
    Ov_filter_fit(istrf).lat = fit_lat(t, forward);
    [bmf fpow pow] = calc_tms(forward, samprate);
    Ov_filter_fit(istrf).bmf = bmf; % In Hz
    
    % Separate into excitatory and inhibitory
    forward_exc = forward;
    forward_exc(find(forward<0)) = 0.0;
    [cf bw forward_fpeak] = fit_cf_bw(f, forward_exc);
    Ov_filter_fit(istrf).ecf = cf./1000.0;
    Ov_filter_fit(istrf).ebw = bw./1000.0;
    Ov_filter_fit(istrf).elat = fit_lat(t, forward_exc);

    forward_inh = forward;
    forward_inh(find(forward>0)) = 0.0;
    [cf bw forward_fpeak] = fit_cf_bw(f, forward_inh);
    Ov_filter_fit(istrf).icf = cf./1000.0;
    Ov_filter_fit(istrf).ibw = bw./1000.0;
    Ov_filter_fit(istrf).ilat = fit_lat(t, forward_inh);

    % Fourier Analysis of Filters
    [p_forward dwt dwf best_dwt best_dwf cm_dwt cm_dwf quad_assym sep_q1 sep_q2]= filter_transfer(forward, 1000, 125);
    Ov_filter_fit(istrf).msbmf = cm_dwt;
    Ov_filter_fit(istrf).msbms = cm_dwf.*1000.0;
    Ov_filter_fit(istrf).quad_assym = quad_assym;
    Ov_filter_fit(istrf).sep_q1 = sep_q1;
    Ov_filter_fit(istrf).sep_q2 = sep_q2;


% Separate filter into temporal and spectral functions and fit with
    % Gabor

    [sep_index, forward_sep, tx, fy, k] = calc_svd(forward);

    [A, X0, Sig, Bw, P, fit_tx, ssegabor] = fit_gabor(tx);
    Ov_filter_fit(istrf).name = sure_OV{istrf};
    Ov_filter_fit(istrf).t_A = A;
    Ov_filter_fit(istrf).t_X0 = X0;
    Ov_filter_fit(istrf).t_Sig = Sig;
    Ov_filter_fit(istrf).t_Bw = Bw;
    Ov_filter_fit(istrf).t_P = rem(P,2*pi);
    [A, X0, Bw, fit_tx_guassian, ssegaussian] = fit_gaussian(tx);
    if (ssegaussian < ssegabor )
        Ov_filter_fit(istrf).t_A = A;
        Ov_filter_fit(istrf).t_X0 = X0;
        Ov_filter_fit(istrf).t_Sig = 0;
        Ov_filter_fit(istrf).t_Bw = Bw;
        Ov_filter_fit(istrf).t_P = 0;
        fit_tx = fit_tx_guassian;
    end
    
    [A, X0, Sig, Bw, P, fit_fy, ssegabor] = fit_gabor(fy);
    Ov_filter_fit(istrf).f_A = A;
    Ov_filter_fit(istrf).f_X0 = X0;
    Ov_filter_fit(istrf).f_Sig = Sig./(125/1000.0);
    Ov_filter_fit(istrf).f_Bw = Bw.*(125/1000.0);   % Bandwidth in kHz
    Ov_filter_fit(istrf).f_P = rem(P,2*pi);
    [A, X0, Bw, fit_fy_guassian, ssegaussian] = fit_gaussian(fy);

    if (ssegaussian < ssegabor )
        Ov_filter_fit(istrf).f_A = A;
        Ov_filter_fit(istrf).f_X0 = X0;
        Ov_filter_fit(istrf).f_Sig = 0;
        Ov_filter_fit(istrf).f_Bw = Bw.*(125/1000.0);
        Ov_filter_fit(istrf).f_P = 0;
        fit_fy = fit_fy_guassian;
    end

    nt = length(tx);
    nf = length(fy);

    forward_sep_fit = zeros(size(forward));
    for i=1:nt
        for j=1:nf
            forward_sep_fit(j,i) = fit_fy(j)*fit_tx(i);
        end
    end

    % Uncoment following lines for visual check
    figure(1);
    
    subplot(1,3,1);
    imagesc(t, f./1000.0, forward); 
    title(sprintf('Cell %d - %s', istrf, sure_OV{istrf}));    
    axis xy;
    
    subplot(1,3,2);
    imagesc(t, f./1000.0, forward_sep);
    axis xy;
    title('Separable');
    subplot(1,3,3);
    imagesc(t, f./1000.0, forward_sep_fit);
    axis xy;
    title('Gabor fit');

    if ( Ov_filter_fit(istrf).t_Sig ~= 0 )
        fprintf(1, 'Cell number %d: SigT=%f(ms) BWT=%f(ms)\n', ...
            istrf, 1/Ov_filter_fit(istrf).t_Sig, Ov_filter_fit(istrf).t_Bw);
    else
        fprintf(1, 'Cell number %d: SigT=%f(ms) BWT=%f(ms)\n', ...
            istrf, 0, Ov_filter_fit(istrf).t_Bw);
    end
    
    if ( Ov_filter_fit(istrf).f_Sig ~= 0 )
            fprintf(1, '\t SigF=%f(kHZ) BWF=%f(kHz)\n', ...
               1/Ov_filter_fit(istrf).f_Sig, Ov_filter_fit(istrf).f_Bw );
    else
            fprintf(1, '\t SigF=%f(kHZ) BWF=%f(kHz)\n', ...
               0, Ov_filter_fit(istrf).f_Bw );
    end
    fprintf(1, 'Cf = %f (kHz) bw = %f (kHz) lat = %f (ms) BMF = %f (Hz)\n', Ov_filter_fit(istrf).cf,  Ov_filter_fit(istrf).bw, ...
                                                      Ov_filter_fit(istrf).lat, Ov_filter_fit(istrf).bmf);
    fprintf(1, 'ECf = %f (kHz) Ebw = %f (kHz) Elat = %f (ms)\n',   Ov_filter_fit(istrf).ecf, Ov_filter_fit(istrf).ebw, ...
                                                  Ov_filter_fit(istrf).elat);
    fprintf(1, 'ICf = %f (kHz) Ibw = %f (kHz) Ilat = %f (ms)\n',   Ov_filter_fit(istrf).icf, Ov_filter_fit(istrf).ibw, ...
                                                  Ov_filter_fit(istrf).ilat);
                                              
     figure(2);
     imagesc(dwt, dwf, p_forward);
     axis xy;
     title('Modulation Transfer function');
                                                  
    fprintf(1, 'BMF = %f (Hz) BMS = %f (cycles/kHz) Qassym = %f\n',   Ov_filter_fit(istrf).msbmf, Ov_filter_fit(istrf).msbms, ...
                                                    Ov_filter_fit(istrf).quad_assym);

    pause();

end

save('C:\Documents and Settings\Frederic\My Documents\Matlab\STRF\analysis\OV Paper\Ov_filter_fit.mat', 'Ov_filter_fit', 'nt', 'nf');

