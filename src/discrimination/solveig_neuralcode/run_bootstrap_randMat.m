function run_bootstrap_randMat()
% This script runs the script RandomPredictConfusionMatrices_SM for all units from the same subject

% if nargin==0
%     Birdname = input('Which bird subject should be analyzed? ', 's');
% end

Birds = {'BluBlu1221M', 'GraBla1602F', 'GraBla1703F', 'WhiWhi4826M', 'WhiWhi1009F', 'GreWhi2513F', 'BlaLbl1986M', 'GraGre1804M'};
nBirds = length(Birds);

if ismac()
        [status username] = system('who am i');
        if strcmp(strtok(username), 'Solveig')
            Path = fullfile('/Users','Solveig','PhD','ELECTROPHY','h5_Sorted_(m.u.)', Birds{jj}, [Birds{jj} '_mat']);
        end
else
    Path = fullfile('/auto','k8','fdata', 'solveig','tdt_h5', Birds{jj});
end

for jj = 1:nBirds
    fprintf('Dealing with bird %s...\n', Birds{jj});
    cd(Path)
    input_dir = pwd;
    
    Files = dir(input_dir);
    nFiles = length(Files);
    for i=1:nFiles
        FileName = Files(i).name;
        if strfind(FileName, '.mat')
            RandomPredictConfusionMatrices_SM(FileName, Birds{jj});
        end
    end
end

return