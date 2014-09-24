function nl=GetLongestNF(sex,wave_dir)
% Compute the longest call in ms in the list sent
% INPUT: 
%        stims: array of nstims wavnames
%        stim_dir: directory of all calls 
% OUTPUT: nl longest call in ms

nl=-1;
wild_card=sprintf('%s/NF2.%cbird*.txt',wave_dir,sex);
dwaves=dir(wild_card);
nstims=length(dwaves);
for ns=1:nstims;
    wave_name=dwaves(ns).name;
    a=load([wave_dir,'/',wave_name]);
    [nu,nf]=size(a);   
    if nu > nl;
        nl=nu;
    end;
end;


