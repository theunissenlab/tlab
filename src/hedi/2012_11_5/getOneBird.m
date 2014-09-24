function bird_data=getOneBird(sex,birdn);
data_dir='../data';
output_dir='../output';
image_dir='../output/image';

bird_data={};
bird_data.sex=sex;
bird_data.name=birdn;
bird_data.calls=[]
for call=1:16;
    fname=sprintf('%s/res.pAv4.%cbird%dcall%d.mat',output_dir,sex,birdn,call);
    if exist(fname,'file');
        [syr,fil,o_spec,f_spec]=displayOneBird(sex,birdn,call);
        bird={};        
        bird.call=call;
        bird.syr=syr;
        bird.fil=fil;
        bird.o_spec=o_spec;
        bird.f_spec=f_spec;
        bird_data.calls=[bird_data.calls;bird];        
    end;
end;
