
function [ripple] = make_ripple_sum(amp_env_sum, f_val, t_val, stim_dur, samp_rate, amp_samp_rate);
% Inverts spectrogram by simple sum of random sine waves. 

rip_tones = 1000;    % Number of random sine waves
random_flag = 1;     % Flag to set phase to random value

f_n = length(f_val);

ripple = zeros(1,samp_rate*stim_dur);
ratio=samp_rate/amp_samp_rate;

t = 1:stim_dur*samp_rate; 

for tones = 1:rip_tones

    % Get the cf from a particular distribution - Using constant distribution here
    fval = f_val(1) + (f_val(end) - f_val(1))*rand;
    for i=2:f_n
        if fval < f_val(i)
            f_frac = (fval-f_val(i-1))/(f_val(i)-f_val(i-1));
            x = amp_env_sum(i-1,:) + f_frac*(amp_env_sum(i,:)-amp_env_sum(i-1,:));
            break;
        end
    end 
   y_temp = interp(x,ratio);    
   y = zeros(1,length(t));
   nz = floor((length(y)-length(y_temp))/2);
   y(1+nz:nz+length(y_temp)) = y_temp;
   y(1:1+nz) = y_temp(1);
   y(nz+length(y_temp):end) = y_temp(end);
   
   if (random_flag) 
      phasewav = 2*pi*rand; 
   else 
      phasewav =0; 
   end 
   
    sinwav = sin((2*pi*fval.*t/samp_rate) + phasewav); 
    ripple = ripple + y.*sinwav;
   
end 



 