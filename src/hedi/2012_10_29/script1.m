addpath('../sound');
addpath('../data');
addpath('../output');
a = dir('../data/01M01*')
for x=1:length(a);   
    fname=a(x).name
    ProcessFileFitv3(fname,-0.51,-0.49);
end;