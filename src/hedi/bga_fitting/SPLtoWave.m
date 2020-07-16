function vspl=SPLtoWave(spl)
vspl=spl;
vspl.nw=0*spl.w
nl=length(spl.splits);
for n=1:nl;
   vspl.r_wav{n}=ComputeWaveFromMoto(spl.r_moto{n}.b,vspl.splits{n});
   sp_start=floor(vspl.splits{n}.start/1000*44100);
   sp_end=sp_start+length(vspl.r_wav{n})-1;
   vspl.nw(sp_start:sp_end)=vspl.r_wav{n}; 
end;
