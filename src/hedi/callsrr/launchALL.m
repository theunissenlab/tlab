function launchALL(bird);

h5_path='/Users/hsoula/Research/AudioNeuro/SolveigData/h5/';
path_files=sprintf('%s/%s/*.h5',h5_path,bird);
d=dir(path_files);
n=length(d);
for i=1:n;
    file=d(i).name;
    fprintf('reading %s',sprintf('%s/%s/%s',h5_path,bird,file));
    unit=read_unit_h5file(sprintf('%s/%s/%s',h5_path,bird,file),'r');
    site=unit.site;
    site_val=str2num(site(end));
    stim_dir=sprintf('/Users/hsoula/Research/AudioNeuro/SolveigData/%s/',bird);
    motograms_dir='../output';
    ParseUnit(unit,stim_dir,motograms_dir);
end;