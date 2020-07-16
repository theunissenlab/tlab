function [me, sd, v] = cv2(spike, window)
%[me, sd, v] = cv2(spike, window)
%this is like cv1, but takes a spike train instead of the nth data output of dcp_show
wind1 = ones(1,window);

ntrials = size(spike,1);

spikeavg = mean(spike);
spkavg = conv(spikeavg, wind1)/window;
convspike = zeros(ntrials, length(spkavg));
varspike = zeros(ntrials, length(spkavg));
for i = 1:ntrials
   convspike(i, :) = conv(spike(i,:), wind1);
   %/size(wind1,2);
   %varspike(i,:) = abs(convspike(i,:) - spkavg);
end
%avg = sum(convspike,1)/ntrials;
%variance = sum(varspike,1)/ntrials;
me = mean(convspike);
sd=std(convspike);
v = var(convspike);

%figure
%plot(me, sd,'x');
%grid on;