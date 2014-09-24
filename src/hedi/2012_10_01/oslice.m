function s=oslice(wavfile,slice);
%% optimizing wav data

% load target wave

fname=wavfile;
[wav_target,fs,nbits]=wavread(fname);
spec_target=timefreq(wav_target,44100,'stft');
tmax=spec_target.t(end) % max in s
slice_length=0.125; % slices of 20 ms
nslices= 1; %floor(tmax/slice_length)
slice_l=floor(min(length(spec_target.t)/nslices,length(spec_target.t)))

s=[]
fil=[0 0 0 1000 10000 12000 5000 500 500 0];
X0=[3000,3000,2,3,1,0,fil];
 ns=0;  
%for ns=1:(nslices-1);
%ns=slice;
    offset=1+ns*slice_l;
    slice=spec_target.spec(:,offset:(offset+slice_l-1));
    target=slice;
    g=@(X)LMSslice(X,target,slice_length,nslices,ns);
    options=optimset('MaxIter',200);
    X=fminsearch(g,X0,options);
    X0=X
    %X=simulannealbnd(g,X0);
    s=[s;X];
    
%end;

