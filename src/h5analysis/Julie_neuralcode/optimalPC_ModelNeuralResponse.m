cd /Users/elie/Documents/MATLAB/data/matfile/CalMatfile
input_dir=pwd;
Subjects = dir(input_dir);
Subj=[3 4 5 7 8];
NBPC=1:10:110;
PCnb=cell(5,1);
for SS=1:5
    ss=Subj(SS);
    Idir=fullfile(input_dir, Subjects(ss).name);
    matfiles=dir(fullfile(Idir,'Cal*.mat'));
    LM=length(matfiles);
    fprintf('Reading Calfile of %s\n',Subjects(ss).name);
    PCnb{SS}=zeros(LM,1);
    for hh=1:LM
        MatfilePath=fullfile(Idir, matfiles(hh).name);
        fprintf('Reading %s\n', MatfilePath);
        Cal=load(MatfilePath);
        R2A=Cal.R2A;
        figure(1)
        plot(NBPC, R2A(1,:), 'rs-', NBPC, R2A(2,:), 'co-', NBPC, R2A(3,:), 'g*-');
        axis([0 100 0 1])
        PCnb{SS}(hh)=input('Optimal number of PC is? (1 11 21 31 41 51 61 71 81 91 or 101)');
    end
end
        