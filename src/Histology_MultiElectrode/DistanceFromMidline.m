fid = fopen('/Users/elie/Documents/ManipBerkeley/Brain_Histology/LblBlu2028M/RH_distanceTomidlines.txt');
header = textscan(fid, '%s\t%s', 1);
data = textscan(fid, '%f\t%s');
Med=find(data{1}<0.4);
Lat=find(data{1}>0.4);
MedMean=mean(data{1}(Med))
LatMean=mean(data{1}(Lat))
ElecDist=LatMean-MedMean
RealMedMean=MedMean/ElecDist*0.5
RealLatMean=LatMean*0.5/ElecDist
