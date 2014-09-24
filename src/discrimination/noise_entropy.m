function [y, nentmat, h]=noise_entropy(spiketrain, wordlength, base);
%[y, nentmat, h]=noise_entropy(spiketrain, wordlength, base);
%base should be the maximum number in a bin of spiketrain.  This determines how many words possible.
ntrials=size(spiketrain,1);

%for k = 1:4
  %  if k ~=1
   % spiketrain=spiketrain(:,2:end);
   %end
 usablelength=size(spiketrain,2)-rem(size(spiketrain,2), wordlength);   
spiketrain2=spiketrain(:,1:usablelength);
clear spiketrain;
words=reshape(spiketrain2', wordlength, usablelength*ntrials/wordlength);
%now words are in columns.
maxbase=min([base 3]);
maxwords=maxbase^wordlength;

basevec=zeros(1, wordlength);
for bs=1:wordlength
   basevec(bs)=maxbase^(bs-1);
end
uncommonbasevec=zeros(1, wordlength);
for bs=1:wordlength
   uncommonbasevec(bs)=base^(bs-1);
end
wordsnew=words'*basevec';
%the following finds words with values >= maxbase, these won't be 
%represented uniquely with maxbase.

[ix iy]=find(max(words, [], 1)>=maxbase);
if ~isempty(iy)
wordsuncommon=words(:,iy);
%the following gets the values of unique uncommon words and adds maxwords
%to it to put it out of the numerical range of common words.

wordsnewuncommon=wordsuncommon'*uncommonbasevec'+maxwords;


wordsnew(iy)=wordsnewuncommon;
end

wordsnew=reshape(wordsnew, usablelength/wordlength, ntrials)';

clear words 
%edges=0:max(wordstotal);
if ~isempty(iy)
%totalpossiblewords=maxbase^wordlength+length(iy);
edges=[0:maxbase^wordlength-1 sort(take_out_repeats(wordsnewuncommon))'];
else
   %totalpossiblewords=maxbase^wordlength;
   edges=[0:maxbase^wordlength-1];
end
clear wordsnewuncommon

%the following histograms a little at a time to preserve memory.
h=[];
startvec=1:1000:size(wordsnew,2);
startvec=[startvec size(wordsnew,2)+1];
for sv=1:length(startvec)-1
h=histc(wordsnew(:,startvec(sv):startvec(sv+1)-1), edges);
%h=[h htemp];
p=zeros(size(h));
for hind=1:size(h,2)
  p(:,hind)=h(:,hind)/sum(h(:,hind));
end

%p=(h)./repmat(sum(h), size(h,1),1);
%pold=p;
p(find(p==0))=1;
nentmat(sv)=sum(-sum(p.*log2(p)));
end
y=sum(nentmat)/size(wordsnew,2);

