function run_neural_discr(Birdname)

%% This function runs all spike sorted files from one bird subject (Birdname) using strfinator

if nargin==0
    Birdname = input('Which bird subject should be analyzed? ', 's');
end

cd /auto/k8/fdata/solveig/tdt_h5
cd(sprintf('%s', Birdname))
input_dir = pwd;

Units = dir(input_dir);
nUnits = length(Units);

if (Birdname == 'GreWhi2513F' | Birdname == 'BlaLbl1986M') % birds not tested with Syn stimuli
    for i=1:nUnits
        h5Path = Units(i).name;
        if strfind(h5Path, 'ss') & strfind(h5Path, '.h5')
            neural_discrimination_h5_slvgNoSyn(h5Path, Birdname) % run neural_discrimination code
        end
    end
else
    for i=1:nUnits
        h5Path = Units(i).name;
        if strfind(h5Path, 'ss') & strfind(h5Path, '.h5')
            neural_discrimination_h5_slvg(h5Path, Birdname) % run neural_discrimination code
        end
    end
end

return