srcRootDir='$SRC_ROOT_DIR'

tlabDir=fullfile(srcRootDir, 'tlab', 'src')

pystrfsDir=fullfile(tlabDir, 'pystrfs', 'matlab')
addpath(pystrfsDir)

strflabDir=fullfile(srcRootDir, 'strflab')
addpath(genpath(strflabDir))

unitFile = '$UNIT_FILE'

compute_coherence_and_info_bound_individual(unitFile)
exit

