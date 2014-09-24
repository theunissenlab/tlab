% Global variables
stimlen= 2.0;
fsamp=32000;
fsms = fsamp/1000;
fmax = 8000.0;
meanfund=700;
stdfund=100
meanripple=4000;
stdripple=3000;
loudness=[0.125 0.25 0.5 1.0];
nloud = size(loudness,2);
ntonefiles = 40;
rand('state',sum(100*clock));
randn('state',sum(100*clock));

tonelen = 95.0;
varlen = 66.0;
toneint = 37.0;
varint = 21.0;
toneramp = 25.0;

for nf=1:ntonefiles
	stim=zeros(1,fsamp*stimlen);
   fname = sprintf('randripple%i.dcp',nf);
	i = 1;
	start_sine=i;
	end_sine=start_sine+floor((tonelen+varlen*randn)*fsms);
	while end_sine <= start_sine
		end_sine=start_sine+floor((tonelen+varlen*randn)*fsms);
	end

	while end_sine <= stimlen*fsamp
		randfund = randn(1)*stdfund+ meanfund;  % This is the fundemental of the harmonic stack
		countf = floor(fmax/randfund);
		randtr = randn(1)*stdripple + meanripple;
		randph = rand(1);
		randld = floor(rand(1)*nloud + 1);
		for iif=1:countf
			amp(iif) = 0.5 + 0.5*cos(2*pi*(iif*randfund/randtr+randph));
		end
		for i=start_sine:end_sine                 % Loop trough time in sampling rate
			if (i - start_sine) < toneramp*fsms
				mult1 = (i-start_sine)/(fsms*toneramp);
				mult1 = 0.5*(1.0 - cos(pi*mult1));
			else
				mult1 = 1.0;
			end
			if (end_sine - i) < fsms*toneramp
				mult2 = (end_sine-i)/(fsms*toneramp);
				mult2 = 0.5*(1.0 - cos(pi*mult2));
			else
				mult2 = 1;
			end
			mult = min(mult1, mult2);                   % mult1 and mult2 are onset and offset ramps
			mult = mult * loudness(randld);
			for iif=1:countf                           % Loop over frequency number in harmonic stack iif = 1 is the fundamental, iff=2 the second harmonic
			   stim(i)= stim(i) + mult*amp(iif)*sin(2*pi*iif*randfund*i/fsamp);
            end
            % Add second loop for stim_inh
		end
		start_sine=end_sine+floor((toneint+varint*randn)*fsms);
		while start_sine <= end_sine
			start_sine=end_sine+floor((toneint+varint*randn)*fsms);
		end
		end_sine=start_sine+floor((tonelen+varlen*randn)*fsms);
		while end_sine <= start_sine
			end_sine=start_sine+floor((tonelen+varlen*randn)*fsms);
		end
   end

	write_song(fname,fsamp, fsamp*stimlen, stim);
	figure;
	specgram(stim,512,fsamp);
end


