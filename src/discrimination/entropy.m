function [y, pold, h]=entropy(spiketrain, wordlength);
%[y, p]=entropy(spiketrain, wordlength);
ntrials=size(spiketrain,1);
wordstotal=[];
%for k = 1:4
  %  if k ~=1
   % spiketrain=spiketrain(:,2:end);
   %end
 usablelength=size(spiketrain,2)-rem(size(spiketrain,2), wordlength);   
spiketrain2=spiketrain(:,1:usablelength);
words=reshape(spiketrain2', wordlength, usablelength*ntrials/wordlength);
%now words are in columns.
base=max(max(spiketrain))+1;
basevec=zeros(1, wordlength);
maxbase=min([max(max(spiketrain))+1 3]);
maxwords=maxbase^wordlength;
for bs=1:wordlength
   basevec(bs)=maxbase^(bs-1);
end
wordsnew=words'*basevec';


uncommonbasevec=zeros(1, wordlength);
for bs=1:wordlength
   uncommonbasevec(bs)=base^(bs-1);
end
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
%end


%edges=0:max(wordstotal);
if ~isempty(iy)
totalpossiblewords=maxbase^wordlength+length(iy);
edges=[0:maxbase^wordlength-1 sort(take_out_repeats(wordsnewuncommon))'];
else
   totalpossiblewords=maxbase^wordlength;
   edges=[0:maxbase^wordlength-1];
end
h=histc(wordsnew, edges);
p=h/(sum(h));
pold=p;
p(find(p==0))=1;
y=-sum(p.*log2(p));

