% try to fit a whole call

x=wavread(fname);
tx=timefreq(x,44100,'stft');
tmax=tx.t(end);
n=length(tx.t);
%xpa=[];
X0=[-0.2553, -0.1187, 3.724e3,1.383e3];


% for t=(n/2-1):-1:1;
%     t
%     spec=tx.spec(:,t);
%     [pmax,imax]=max(spec);
%     spec=spec/pmax;   
%     mz=findAlphaBetaX0(tx.f,spec,X0);
%     xpa=[xpa;t mz];
%     X0=mz(1:4)
% end;
X0=[-0.2553, -0.1187, 3.724e3,1.383e3];
for t=n/2:n;
    t
    spec=tx.spec(:,t);
    [pmax,imax]=max(spec);
    spec=spec/pmax;   
    mz=findAlphaBetaX0(tx.f,spec,X0);
    xpa=[xpa;t mz];
    X0=mz(1:4)
end;

