function vspl=SmoothSPLtoWave(spl)
vspl=spl;
vspl.snw=0*spl.w;
nl=length(spl.splits);
for n=1:nl;   
   sp_start=floor(vspl.splits{n}.start/1000*44100)
   sp_end=sp_start+length(vspl.smooth{n}.nw)-1
   vspl.snw(sp_start:sp_end)=vspl.smooth{n}.nw; 
end;
