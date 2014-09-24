
addpath /auto/k2/share/fairslurm/;
in_dir='/auto/k6/hsoula/stress/wavs/';
out_dir='/auto/k6/hsoula/stress/new_output/';
dirs_calls=dir([in_dir,'*.wav']);
n=length(dirs_calls);

for i=1:n;
  name=dirs_calls(i).name;
  fname=[in_dir,name];
  out_name=[out_dir,'/','fit.',name(1:(end-4)),'.mat'];
  if ~exist(out_name,'file'); 
  num_com='';%%maxNumCompThreads(1);' ;
  fit_cm=['fit_calls(''',fname,''');'];
  
  params={};
  params.out=['/auto/k6/hsoula/stress/',name,'.txt'];
  params.memory=10000;
  params.err='/auto/k6/hsoula/stress/err.file';
  params.cpus=2;
  fprintf(['launching',' ',num_com,' ',fit_cm,'\n'])
  slurm_sbatch([num_com,fit_cm],params);
   end;
end;
