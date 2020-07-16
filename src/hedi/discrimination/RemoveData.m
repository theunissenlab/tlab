function mgram=RemoveData(motodata,type);

mgram=motodata;
[ncalls,ndims]=size(mgram);
nsize=ndims/4;
if strcmp(type,'syrinx');
    start=nsize+1;
    mgram=mgram(:,start:end);
else
    inx=[1:nsize,3*nsize+1:ndims];
    mgram=mgram(:,inx);
end;
