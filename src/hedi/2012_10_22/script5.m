a = dir('../data/01M05*')
for x=1:length(a);   
    fname=a(x).name
    ProcessFileFitv3(fname,-0.54,-0.49);
end;