function sound_syn = make_complex_sound(t_fund, fund, f_env, s_env, t_env, a_env, fs)
% Make a synthetic complex time varying sound from a harmonic stack.  The
% stack has a fundamental given by fund frequency values given in tfund.  The
% spectral enveloppe is given by senv sampled at frequencies fenv.  The
% amplitude enveloppe is given by aenv sampled at tenv.  The sampling
% rateof the output sound is fs.  The sound lenght is given by the minimum
% of t_fund and t_env.  Both t_fund and t_env are in ms.

% Set to 1 for debugging figure
debug_fig = 1;

% Find lenght of new sound
len_sound = min(max(t_fund), max(t_env));
fprintf(1, 'Lenght env %f Length fund %f Length sound %f', max(t_env), max(t_fund), len_sound);

npts_sound = fix(len_sound*fs/1000.0);                  % Length of sound in number of points


% Find upper range for frequency spectrum
fmax = max(f_env);
if (fmax > fs/2)
    fmax = fs/2;
end

sound_syn = zeros(1, npts_sound);
t_fund_current_id = 1;
t_env_current_id = 1;

fund_test = zeros(1, npts_sound);
amp_test = zeros(1, npts_sound);
senv_test = zeros(3, npts_sound);

for i=1:npts_sound
    
    % Current time in ms
    ti = (i-1)*1000.0/fs;
    
    % Current fundamental
    while (ti > t_fund(t_fund_current_id))
        
        if (t_fund_current_id == length(t_fund) )
            break;
        end
        t_fund_current_id = t_fund_current_id + 1;
    end
    if t_fund_current_id == 1
        if isnan(fund(t_fund_current_id))
            fti = 0;
        else
            fti = fund(t_fund_current_id);
        end
    else
        if isnan(fund(t_fund_current_id)) && isnan(fund(t_fund_current_id-1))
            fti = 0;
        elseif isnan(fund(t_fund_current_id))
            fti = fund(t_fund_current_id-1);
        elseif isnan(fund(t_fund_current_id-1))
            fti = fund(t_fund_current_id);
        else
            fti = fund(t_fund_current_id-1) + (ti-t_fund(t_fund_current_id-1)).*(fund(t_fund_current_id)-fund(t_fund_current_id-1))./(t_fund(t_fund_current_id)-t_fund(t_fund_current_id-1));
        end
    end
    if fti < 50.0
        countf = 1;
    else
        countf = floor(fmax/fti);
    end
    fund_test(i) = fti;
    
    % Current temporal amplitude enveloppe
    while (ti > t_env(t_env_current_id))      
        if (t_env_current_id == length(a_env) )
            break;
        end
        t_env_current_id = t_env_current_id + 1;
    end
    if t_env_current_id == 1
        t_amp = a_env(t_env_current_id);
    else
        t_amp = a_env(t_env_current_id-1) + (ti-t_env(t_env_current_id-1)).*(a_env(t_env_current_id)-a_env(t_env_current_id-1))./(t_env(t_env_current_id)-t_env(t_env_current_id-1));
    end
    amp_test(i) = t_amp;
    
    % Use the spectral enveloppe to get relative amplitude
    s_amp = zeros(1, countf);
    
    % countf = 1;   % for debugging purposes...
    for iif=1:countf
        fval = iif*fti;
        if fval <= 0
            s_amp(iif) = 0;
        else
            flow_ind = find(f_env <= fval, 1, 'last');
            if flow_ind == length(f_env)
                s_amp(iif) = s_env(flow_ind);
            else
                s_amp(iif) = s_env(flow_ind) + (fval - f_env(flow_ind)).*(s_env(flow_ind+1) - s_env(flow_ind))./(f_env(flow_ind+1)-f_env(flow_ind));
          
            end
        end
    end
    sum_s_amp = sum(s_amp);
    if (sum_s_amp ~= 0)
        s_amp = s_amp./sum_s_amp;
    end
    
    % Now generate harmonic sound
    sound_syn(i) = 0;
    for iif=1:countf
        fval = iif*fti;
        if fval ~= 0
            sound_syn(i)= sound_syn(i) + s_amp(iif)*cos(2*pi*fval*ti/1000.0);
        end
        if iif <=3
            senv_test(iif, i) = s_amp(iif);
        end
    end
    
    % ... and multiply by amplitude envelope
    sound_syn(i) = sound_syn(i)*t_amp;
end

if debug_fig
    figure(10);
    plot(amp_test);
    title('Amplitude');
    figure(11);
    plot(fund_test);
    title('Fundamental');
    figure(12);
    plot(senv_test(1,:),'b');
    hold on;
    plot(senv_test(2,:),'r');
    plot(senv_test(3,:),'k');
    hold off;
    figure(13);
    clf;
    spec(sound_syn, 125, fs, 50);
    pause;
end
