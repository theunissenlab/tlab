srcRootDir='$SRC_ROOT_DIR'

tlabDir=fullfile(srcRootDir, 'tlab', 'src')

pystrfsDir=fullfile(tlabDir, 'pystrfs', 'matlab')
addpath(pystrfsDir)

strflabDir=fullfile(srcRootDir, 'strflab')
addpath(genpath(strflabDir))

output_dir = '$OUTPUT_DIR'
file_list = '$FILE_LIST'
output_files = textread(file_list, '%s');

for k = 1:length(output_files) compute_output_nl_for_figure(fullfile(output_dir, output_files{k})); end

exit

