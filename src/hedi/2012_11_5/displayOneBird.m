function [syr,fil,o_spec,f_spec,res]=displayOneBird(sex,bird,call);
data_dir='../data';
output_dir='../output';
image_dir='../output/image';

fname=sprintf('%s/res.pAv4.%cbird%dcall%d.mat',output_dir,sex,bird,call)
ori_wave=sprintf('%s/%cbird%dcall%d.wav',data_dir,sex,bird,call)
fin_wave=sprintf('%s/pAv3.%cbird%dcall%d.wav',output_dir,sex,bird,call)

load(fname);
syr=lowP(res.data(:,2),res.data(:,1),5);
fi=getFilter(res.data(:,1),res.data(:,3:end));        
fil=lowP(fi(:,1),res.data(:,1),5);

o_spec=plotSpectrum(ori_wave);
f_spec=plotSpectrum(fin_wave);
        
