function [spiketrain, rescaled, debugvec]=make_spikes3(meanrate, gammaconst, ntrials);
%[spiketrain, debugvec]=make_spikes2(meanrate, gammaconst, ntrials);
%same as make_spikes2, but faster.
%meanrate must not have any zeros for spline to work.
%in gamrnd, a is the gammaconst, (1/ba) is the meanrate.
%In comments, new times refer to the non-homogeneous time.
len=length(meanrate);
spiketrain=zeros(ntrials, len);
debugvec=[];
%only make a certain length for rescaled (saves time.)
rescaled = gamrnd(gammaconst, 1/mean(meanrate)/gammaconst, ntrials, floor(len*mean(meanrate)*1.5));
crescaled=cumsum(rescaled,2);
%if it doesn't reach len, add more
while ~isempty(find(crescaled(:,end)<len))
    rescaled = [rescaled gamrnd(gammaconst, 1/mean(meanrate)/gammaconst, ntrials, 100)];
crescaled=cumsum(rescaled,2);
end
%extra are the positions that need an extra spike.

%the following has it start at a random point in the first isi.
start=rand(ntrials,1).*rescaled(:,1);
rescaled=[rescaled(:,1)-start rescaled(:,2:end)];
crescaled=cumsum(rescaled,2);
%check again to make sure rescaled still adds up to len
while ~isempty(find(crescaled(:,end)<len))
    rescaled = [rescaled gamrnd(gammaconst, 1/mean(meanrate)/gammaconst, ntrials, 10)];
crescaled=cumsum(rescaled,2);
end
maxlength=0;
for n= 1:ntrials
    temp=min(find(crescaled(n,:)>len));
    if and(temp>maxlength, temp<=len)
        maxlength=temp;
    end
end
rescaled=rescaled(:,1:maxlength);
%start will be the vector of all spike times.
oldtimes=cumsum(rescaled,2);
rnew=cumsum(meanrate)/mean(meanrate);
%The following sorts oldtimes into rnew.  Thus it's placement will be next to
%its corresponding between its possible corresonding new time
%(non-homogeneous time).  
newtimes=interp1([0 rnew], 0:length(rnew), oldtimes);
newtimes=round(newtimes);
newtimes(find(newtimes<1))=1;
%the following accounts for unnatural peaks due to interp failure;
[outindexx, outindexy]=find(diff(newtimes)<0);
if ~isempty(outindexy)
    for ay=1:length(outindexy)
        if outindexx(ay)<length(meanrate)
            newtimes(outindexx(ay), outindexy(ay))=newtimes(outindexx(ay), outindexy(ay)+1);
        end
    end
end
newtimes(find(newtimes>length(meanrate)))=length(meanrate);



edges=1:len;
hvec=histc(newtimes', edges);

spiketrain=hvec';
