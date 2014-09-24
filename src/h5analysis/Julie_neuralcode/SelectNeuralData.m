Res=load(MatfilePath);
StimType=unique(Res.VocType);
NstimType=length(StimType);
%% Create the data set for the type of comparison you want:
% here with Selector = type of stim (AggC, DisC Tet...)
SpikeTrains=cell(NstimType,1);
%Sections_len=zeros(Nsections*length(Res.Trials{1},1));
Sections_len=zeros(NstimType,1);

for vt=1:NstimType
    st=StimType(vt);
    selector=strcmp(Res.VocType, st);
    selectorInd=find(selector);
    NselectorInd=length(selectorInd);
    %Sections_len{vt}=zeros(NselectorInd*length(Res.Trials{selectorInd(1)}));
    ISec = 1;
    for NSI = 1:NselectorInd
        SI=selectorInd(NSI);
        SpikeTrains{vt} = [SpikeTrains{vt} ; Res.Trials{SI}];
        NSec = length(Res.Trials{SI});
        Sections_len(vt)=min(Res.SectionLength(SI));
        ISec=ISec+NSec;
    end
end

% then you can use the output s such in info_distanceB for instance:
%[pc, mi_conf, zdT, pzdT, mi_zdT, mi_zdT_nc, zdS, pzdS, mi_zdS, mi_zdS_nc, confusion_matrix] = info_distanceB(NstimType, SpikeTrains , Sections_len, @VR_distanceB, winSize(is));
       