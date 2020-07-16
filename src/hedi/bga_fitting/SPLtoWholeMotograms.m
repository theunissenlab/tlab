function mspl=SPLtoWholeMotograms(spl)
mspl=spl;
mspl.wspec=htimefreq(spl.w,44100,'stft');
nf=length(mspl.wspec.t);
mspl.alpha=zeros(nf,1);
mspl.mu=zeros(nf,1);
mspl.sigma1=zeros(nf,1);
mspl.sigma2=zeros(nf,1);
nl=length(spl.wavF);
for n=1:nl;
   moto=mspl.wavF{n}.b;
   sp_start=mspl.splits{n}.start;
   sp_end=sp_start+length(moto(:,1))-1;
   mspl.alpha(sp_start:sp_end)=moto(:,2);
   mspl.mu(sp_start:sp_end)=moto(:,3);
   mspl.sigma1(sp_start:sp_end)=moto(:,4);
   mspl.sigma2(sp_start:sp_end)=moto(:,5);
end;
% 
%     mspl.s_alpha=zeros(nf,1);
%     mspl.s_mu=zeros(nf,1);
%     mspl.s_sigma1=zeros(nf,1);
%     mspl.s_sigma2=zeros(nf,1);
%     nl=length(spl.smooth);
%     for n=1:nl;
%         s_moto=mspl.smooth{n}.moto;
%         sp_start=mspl.splits{n}.start;
%         sp_end=sp_start+length(s_moto(:,1))-1;
%         mspl.s_alpha(sp_start:sp_end)=s_moto(:,2);
%         mspl.s_mu(sp_start:sp_end)=s_moto(:,3);
%         mspl.s_sigma1(sp_start:sp_end)=s_moto(:,4);
%         mspl.s_sigma2(sp_start:sp_end)=s_moto(:,5);
%     end;
%     
    
nsp=htimefreq(mspl.nw,44100,'stft');
subplot(4,1,1);
imagesc(mspl.wspec.t,mspl.wspec.f,mspl.wspec.spec);axis xy;

subplot(4,1,2);
imagesc(nsp.t,nsp.f,nsp.spec);axis xy;

subplot(4,1,3);

plot(mspl.wspec.t,mspl.mu,'k',mspl.wspec.t,mspl.mu+mspl.sigma1,'-r',mspl.wspec.t,mspl.mu-mspl.sigma2,'-r')
axis([0 mspl.wspec.t(end) 0 20000]);

subplot(4,1,4);
plot(mspl.wspec.t,-min(mspl.alpha,0.4)-0.4)

axis([0 mspl.wspec.t(end) 0.0 0.2]);

% sp=htimefreq(mspl.snw,44100,'stft');
% subplot(4,2,2);
% imagesc(mspl.wspec.t,mspl.wspec.f,mspl.wspec.spec);axis xy;
% 
% subplot(4,2,4);
% imagesc(nsp.t,sp.f,sp.spec);axis xy;
% 
% subplot(4,2,6);
% 
% plot(mspl.wspec.t,mspl.s_mu,'k',mspl.wspec.t,mspl.s_mu+mspl.s_sigma1,'-r',mspl.wspec.t,mspl.s_mu-mspl.s_sigma2,'-r')
% axis([0 mspl.wspec.t(end) 0 20000]);
% 
% subplot(4,2,8);
% plot(mspl.wspec.t,-min(mspl.s_alpha,0.4)-0.4)
% 
% axis([0 mspl.wspec.t(end) 0.0 0.2]);


