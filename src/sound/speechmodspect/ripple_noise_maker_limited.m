% Global variables

nsongs = 2;% Number of noise ripples to be made
ripple_type = 2; % 1 = song ripple noise 2 = flat ripple noise % 3 = specify on command line values for w_x and w_t
%Defining the cutoff for w_x and w_t for flat ripple noise
max_wx = 2.5; %in cycles/KHz
min_wx = 2;
max_wt = 70; %in Hz
min_wt = 30;
%save_dir = '/home/fet/matlab/stim/flatrip3_100/';
save_dir = 'C:/Documents and Settings/Administrator/My Documents/Matlab/make_stim/flatrip_limited';
load_dir = 'C:/Documents and Settings/Administrator/My Documents/Matlab/Amp Mod Power spectra/Data/Birdsong/dc/lin';

% The global parameters for calculating the ripple 
global f_low f_high f_step f_val ntones; 
global stim_dur samp_rate samp_ms tone_ramp amp_samp_rate amp_samp_ms; 
global no_of_ripples;
f_low=250;   
f_high=8000; 
f_step=50;    % 10 cycles per KHz Nyquist limit 
no_of_ripples=100; 
tone_ramp=25; % tone ramp in ms 
stim_dur=2000; % stimulus duration in ms 
samp_rate=32000; 
amp_samp_rate=1000;   % 500 Hz Nyquist limit
amp_samp_ms=amp_samp_rate/1000; % amplitude env sampling rate in ms
samp_ms=samp_rate/1000; 
f_val=0.0:f_step:f_high;
ntones = size(f_val,2);

 
%parameters required for selecting the power 
if (ripple_type == 1)
    power_tol=0.005; 
    cd /home/nandini/matlab/make_stim/new_analysis/125;
    stim_stat_pcolor_with_contour;
    new_stim_stat_fplot=zeros(size(stim_stat_fplot)); 
    %selecting the power from songs power spectrum to calculate ripple 
    max_power = max(max(abs(stim_stat_fplot))); 
    power_cutoff=power_tol*(max_power); 
    
    % this find command give the indices of (x,y) for power values > power_cutoff 
    [w_x,w_t]=find(abs(stim_stat_fplot)> power_cutoff); 
    sum_of_amp=0; 
    for k=1:length(w_x) 
        sum_of_amp=sum_of_amp + (abs(stim_stat_fplot(w_x(k),w_t(k)))); 
    end 
    
    %set up probability density distribution matrix 
    for k=1:length(w_x) 
        new_stim_stat_fplot(w_x(k),w_t(k))= (abs(stim_stat_fplot(w_x(k),w_t(k))))/sum_of_amp; 
    end 
    
    % set up cumulative probability density function 
    cum_power_dist=zeros(size(w_x)); 
    cum_power_dist(1)=new_stim_stat_fplot(w_x(1),w_t(1)); 
    for k=2:length(w_x) 
        cum_power_dist(k)=new_stim_stat_fplot(w_x(k),w_t(k)) + cum_power_dist(k-1); 
    end 
    
end 

%memory allocation for calculating the ensemble of song ripples 
amp_env_sum = zeros(nsongs, ntones,stim_dur*amp_samp_ms); 
mean_amp=zeros(1,ntones);
var_amp=zeros(1,ntones);
maxrip=zeros(1,nsongs);

%memory allocation for omega_x's and omega_t's selected
omx=zeros(nsongs,no_of_ripples);
omt=zeros(nsongs,no_of_ripples);
omxrip=zeros(nsongs,no_of_ripples);
omtrip=zeros(nsongs,no_of_ripples);
 

% Generate initial amplitude enveloppe for ripple songs
for ns=1:nsongs 
    i=1; 
    if (ripple_type == 1)
        while (i <= no_of_ripples) 
            ran_power_val=rand;
            for j=1:size(w_x) 
                if(ran_power_val < cum_power_dist(j)); 
                    omega_x(i)=(w_x(j)-nb+rand-0.5)*(1/(2*fstep*(nb-1)))*1000;% in cycles per KHz 
                    omega_t(i)=(w_t(j)-ntt+rand-0.5)*(1000/(2*(ntt-1))) ;
                    % Check for correct bounds: must be below Nyquist sampling rate for temp and spect modulation
                    % and within the omega_t*omega_x < i region.  omega_x is in cycles /kHz thus 1000
                    if (  abs(omega_x(i)) > 1000./(2*f_step)  | ...
                            abs(omega_t(i)) > amp_samp_rate./2 | ...
                            abs(omega_x(i)*omega_t(i)) >= 1000      ) 
                        i = i +0;
                    else                    
                        i=i+1;  
                    end
                    break;     
                end 
            end 
        end 
    elseif(ripple_type == 2)
        while (i <= no_of_ripples) 
            omega_x(i)=rand*2.0*max_wx - max_wx;
            omega_t(i)=rand*2.0*max_wt - max_wt;
            % Check for correct bounds: must be below Nyquist sampling rate for temp and spect modulation
            % and within the omega_t*omega_x < i region.  omega_x is in cycles /kHz thus 1000
            if (  abs(omega_x(i)) > 1000./(2*f_step)  | ...
                  abs(omega_t(i)) > amp_samp_rate./2 | ...
                  abs(omega_x(i)*omega_t(i)) >= 1000 | ...
                  abs(omega_x(i)) < min_wx | ...
                  abs(omega_t(i)) < min_wt ) 
                continue; 
            end
            i=i+1; 
        end 
    end
    
    
    for i=1:no_of_ripples      
        omx(ns,i)=omega_x(i);
        omt(ns,i)=omega_t(i);
        wrip_x = omega_x(i);
        wrip_t = omega_t(i);
        amp_env = make_amp_env(wrip_x,wrip_t,1);
        amp_env_sum(ns,:,:) = amp_env_sum(ns,:,:) + amp_env(1,:,:);
    end    
    mean_amp = mean_amp + mean(amp_env_sum(ns,:,:),3);
    
end


% Subtract the mean of envelopes and calculate the variance

mean_amp = mean_amp./nsongs;
% Make the songs have zero mean
for ns=1:nsongs
   for i=1:ntones
		amp_env_sum(ns,i,:)=amp_env_sum(ns,i,:) - mean_amp(i);
		var_amp(i) = var_amp(i) + sum(amp_env_sum(ns,i,:).*amp_env_sum(ns,i,:));
	end
end

var_amp = var_amp./(nsongs*stim_dur*amp_samp_ms);

    
% Load the song mean enveloppe and variance

load_file = sprintf('%s/stim_init_avg.avg', load_dir);
load(load_file);
load_file = sprintf('%s/stim_var.avg', load_dir);
load(load_file);
load_file = sprintf('%s/stim_init_count.avg', load_dir);
load(load_file);

flow = stim_init_count(4);
fhigh = stim_init_count(5);
fstep = stim_init_count(6);
f1=flow+fstep/2:fstep:fhigh;
amp_song=interp1(f1,stim_init_avg,f_val,'cubic');
var_song=interp1(f1,stim_var,f_val,'cubic');

    
% Match the variance and mean
if ripple_type == 1
    var_fix = sqrt(var_song./var_amp);
    for ns=1:nsongs
        for i=1:ntones
            amp_env_sum(ns,i,:)=amp_env_sum(ns,i,:)*var_fix(i) + amp_song(i);
        end
    end
else
    var_fix = sqrt(mean(var_song)/mean(var_amp));
    amp_fix = mean(amp_song);
    amp_env_sum = amp_env_sum*var_fix + amp_fix;
end

% Transform to exp after zero-ing out all negative values.
for ns=1:nsongs
    [indx,indy]= find(squeeze(amp_env_sum(ns,:,:))< 0.0);
    fprintf(1,'Found %d points below zero for song %d\n', length(indx), ns);
    if(length(indx) > 0)
        i=length(indx);
        amp_env_sum(ns,indx,indy) = 0.0;
    end
end
% amp_env_sum = exp(amp_env_sum)-1;

save_file = sprintf('%s/flatrip_limited.mat', save_dir);
save(save_file,'amp_env_sum'); 

figure;
subplot(2,1,1);
imagesc(0:stim_dur/1000, f_val/1000, squeeze(amp_env_sum(1,:,:)));
subplot(2,1,2);
plot(squeeze(amp_env_sum(1,20,:)));

% Plot the spectrograms for the stims
%figure(2);
%figure;
%if (ripple_type == 1)
 %   for ns=1:nsongs;
  %      nrows=nsongs/2;
   %     ncoll=2;
    %    subplot(nrows,ncoll,ns);
        %subplot(2,1,1)
     %   imagesc(0:stim_dur/1000,f_val/1000,squeeze(amp_env_sum(ns,:,:)));
     %   axis xy;
      %  axis([0 stim_dur/1000 f_low/1000 f_high/1000]);
      %  axis square;
        %xlabel('Time (secs)');
        %ylabel('Frequency (Khz)');
        %titlestring=sprintf('Amplitude envelope with wx=%f and wt=%f',wrip_x,wrip_t);
        %title(titlestring);
        %hold;
        % end
        %else
    %for ns=1:nsongs;
     %   nrows=nsongs/4;
      %  if (nrows < 1 )
       %     nrows = 1;
       % end
       % ncoll=4;
       % subplot(nrows,ncoll,ns);
        %subplot(2,1,1);
        %imagesc(0:stim_dur/1000,f_val/1000,squeeze(amp_env_sum(ns,:,:)));
        %axis xy;
        %axis square;
        %xlabel('Time (secs)');
        %ylabel('Frequency (Khz)');
        %titlestring=sprintf('Amplitude envelope for a song ripple with wx=%f and wt=%f',wrip_x,wrip_t);
        %title(titlestring);
        %hold;
        %end
        %end


%figure
%imagesc(0:stim_dur/1000,f_val/1000,squeeze(amp_env_sum(1,:,:)));

countrip=0;

for ns=1:nsongs   
    ns
    ripple = make_ripple(squeeze(amp_env_sum(ns,:,:)),1);
    maxr = max(ripple);
    minr = min(ripple);
    maxabs = max(maxr,abs(minr))
    if((ripple_type == 1) & (maxabs < 2.5e+10))
        countrip=countrip+1;
        ripple_file =sprintf('/home/nandini/matlab/make_stim/new_analysis/ripples/examples/songrip%d.mat',countrip);
        save(ripple_file,'ripple');
        maxrip(countrip) = max(maxr,abs(minr));
        for j=1:no_of_ripples
            omxrip(countrip,j)=omx(ns,j);
            omtrip(countrip,j)=omt(ns,j);
        end
    elseif(ripple_type == 2)
        countrip=countrip+1;
        ripple_file = sprintf('%s/flatrip%d.mat',save_dir, countrip);
        save(ripple_file,'ripple');
        maxrip(countrip) = max(maxr,abs(minr));
        for j=1:no_of_ripples
            omxrip(countrip,j)=omx(ns,j);
            omtrip(countrip,j)=omt(ns,j);
        end
    end
    
    %figure(3);
    %specgram(ripple,512,samp_rate); 
end
%figure
%specgram(ripple,512,samp_rate);

%axis xy;
ripmax=max(maxrip);


for ns=1:countrip
    if(ripple_type == 1)
        ripple_file= sprintf('/home/nandini/matlab/make_stim/new_analysis/ripples/examples/songrip%d.mat',ns);
        load(ripple_file);
        dbatt=0;
        wav_name = sprintf('/home/nandini/matlab/make_stim/new_analysis/ripples/examples/songrip%d.wav',ns);
        dcp_name = sprintf('/home/nandini/matlab/make_stim/new_analysis/ripples/examples/rip%d.dcp',ns);
        %plot(squeeze(omtrip(ns,:)),squeeze(omxrip(ns,:)),'b+');
        %hold on;
    elseif(ripple_type == 2)
        ripple_file= sprintf('%s/flatrip%d.mat',save_dir, ns);
        load(ripple_file);
        dbatt=10;
        wav_name = sprintf('%s/flatrip%d.wav',save_dir,ns);
        dcp_name = sprintf('%s/flatrip%d.dcp',save_dir,ns);
    end
    
    ripple = ripple./(ripmax*power(10.0,(dbatt/20))); 
    wavwrite(ripple, samp_rate, 16, wav_name);
    write_msong(dcp_name,samp_rate,samp_ms*stim_dur,ripple);
end

%if(ripple_type == 1 )
 %   print -dmfile /home/nandini/matlab/make_stim/new_analysis/ripples/song_ripples/figures/wx_wt_plot_song;
 %else
 %   print -dmfile /home/nandini/matlab/make_stim/new_analysis/ripples/flat_ripples/figures/wx_wt_plot_flat
%figure(3)
%subplot('position',[0.3 0.075 0.35 0.35])
%specgram(ripple,512,samp_rate);
%axis xy;
%axis([0 stim_dur/1000 0 f_high]);
%axis square;
%xlabel('Time (s)');
 %ylabel('Frequency (Hz)');
 %title('Spectrogram for flat ripple')
