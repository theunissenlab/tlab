function morphs=MgramMorphTimber(mgStart,mgEnd,nstep);

[startDur,npa]=size(mgStart);
[endDur,npa]=size(mgEnd);
morphs={};
%morphs{1}=mgStart;
%morphs{nstep}=mgEnd;
Dur=startDur;
if startDur<endDur;
    Dur=endDur;
end;

tStart=interp1(mgStart(:,1),mgStart(:,1),linspace(1,startDur,Dur));
tEnd=interp1(mgEnd(:,1),mgEnd(:,1),linspace(1,endDur,Dur));

astart=interp1(mgStart(:,1),mgStart(:,2),tStart);
aend=interp1(mgEnd(:,1),mgEnd(:,2),tEnd);
mstart=interp1(mgStart(:,1),mgStart(:,3),tStart);
mend=interp1(mgEnd(:,1),mgEnd(:,3),tEnd);
sstart=interp1(mgStart(:,1),mgStart(:,4),tStart);
send=interp1(mgEnd(:,1),mgEnd(:,4),tEnd);
pstart=interp1(mgStart(:,1),mgStart(:,5),tStart);
pend=interp1(mgEnd(:,1),mgEnd(:,5),tEnd);
i=1;
for t=linspace(0,1,nstep/2);    
    tin=(1-t/2)*tStart+(t/2)*tEnd;
    ain=astart;
    min=(1-t)*mstart+t*mend;
    sin=(1-t)*sstart+t*send;
    pin=(1-t/2)*pstart+t/2*pend;
    mg=[tin',ain',min',sin',pin'];
    morphs{i}=mg;
    i=i+1;    
end;

for t=linspace(0,1,nstep/2);    
    tin=(1/2-t/2)*tStart+(1/2+t/2)*tEnd;
    ain=(1-t)*astart+t*aend;
    min=mend;
    sin=send;
    pin=(1/2-t/2)*pstart+(1/2+t/2)*pend;
    mg=[tin',ain',min',sin',pin'];
    morphs{i}=mg;
    i=i+1;    
end;
