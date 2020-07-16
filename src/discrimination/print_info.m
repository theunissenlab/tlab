% Print the results to a text file

data_dir = '/auto/fhome/fet/matlab/Neural Discrimination/fieldL/';
data_file = 'zfnormal_fieldL.mat';

% Load data file
load(strcat(data_dir, data_file));


fid = fopen(fname, 'w');
fprintf(fid,'Bird\tCell\t');
fprintf(fid,'Zscore Ben\tZscore Zeb\tZscore Ml\t');
fprintf(fid,'Dprime Ben\tDprime Zeb\tDprime Ml\t');
fprintf(fid,'GammaInfo Ben\tGammaInfo Zeb\tGammaInfo Ml\t');
fprintf(fid,'GammaNoise Ben\tGammaNoise Zeb\tGammaNoise Ml\t');
fprintf(fid,'GammaRate Ben\tGammaRate Zeb\tGammaRate Ml\t');
fprintf(fid,'GammaCnst Ben\tGammaCnst Zeb\tGammaCnst Ml\t');
fprintf(fid,'Bandwidth Ben\tBandwdith Zeb\tBandwidth Ml\t');
fprintf(fid,'RateGamma Ben\tRateGamma Zeb\tRateGamma Ml\t');
fprintf(fid,'Pcorrect Ben\tPcorrect Zeb\tPcorrect ML\t');
fprintf(fid,'ConfInfo Ben\tConfInfo Zeb\tConfInfo ML\t');
fprintf(fid,'zdist Ben\tzdist Zeb\tzdist ML\t\n');



for nc=1:number_cells
    
    % This needs to be changed depending on the brain region
    fprintf(fid, '%s\t%s\t', fieldL_normal{nc}{1}, fieldL_normal{nc}{2});
    
    % Check the order of the stimulus in each loaded file to make sure it
    % is Bengalese, Conspecific and Flatrip
    fprintf(fid,'%f\t%f\t%f\t', alloutputs{nc,1}.avgzscore, alloutputs{nc,2}.avgzscore,alloutputs{nc,3}.avgzscore);
fprintf(fid,'Dprime Ben\tDprime Zeb\tDprime Ml\t');
fprintf(fid,'%f\t%f\t%f\t', alloutputs{nc,1}.avgdprime, alloutputs{nc,2}.avgdprime,alloutputs{nc,3}.avgdprime);
fprintf(fid,'GammaInfo Ben\tGammaInfo Zeb\tGammaInfo Ml\t');
fprintf(fid,'%f\t%f\t%f\t', alloutputs{nc,1}.gamma_mutual_info, alloutputs{nc,2}.gamma_mutual_info,loutputs{nc,3}.gamma_mutual_info);
fprintf(fid,'GammaNoise Ben\tGammaNoise Zeb\tGammaNoise Ml\t');
fprintf(fid,'GammaRate Ben\tGammaRate Zeb\tGammaRate Ml\t');
fprintf(fid,'GammaCnst Ben\tGammaCnst Zeb\tGammaCnst Ml\t');
fprintf(fid,'Bandwidth Ben\tBandwdith Zeb\tBandwidth Ml\t');
fprintf(fid,'RateGamma Ben\tRateGamma Zeb\tRateGamma Ml\t');
fprintf(fid,'Pcorrect Ben\tPcorrect Zeb\tPcorrect ML\t');
fprintf(fid,'ConfInfo Ben\tConfInfo Zeb\tConfInfo ML\t');
fprintf(fid,'zdist Ben\tzdist Zeb\tzdist ML\t\n');

    
end