srcRootDir='$SRC_ROOT_DIR'

tlabDir=fullfile(srcRootDir, 'tlab', 'src')

pystrfsDir=fullfile(tlabDir, 'pystrfs', 'matlab')
addpath(pystrfsDir)

strflabDir=fullfile(srcRootDir, 'strflab')
addpath(genpath(strflabDir))

unitFile = '$UNIT_FILE'
preprocFile = '$PREPROC_FILE'
stimClass = '$STIM_CLASS'

trainingGroups = $TRAINING_GROUPS
validationGroups = $VALIDATION_GROUPS

strfLength = 60
dfResponseFile = '$DF_RESPONSE_FILE'
dfDistsResponseFile = '$DF_DISTS_RESPONSE_FILE'
dfSplineResponseFile = '$DF_SPLINE_RESPONSE_FILE'

transforms = {$TRANSFORMS}

run_directfit(unitFile, preprocFile, stimClass, trainingGroups, validationGroups, strfLength, transforms, dfResponseFile)
compute_response_info_coherence(dfResponseFile)

run_outputnl(dfResponseFile, dfDistsResponseFile, 'dists')
compute_response_info_coherence(dfDistsResponseFile)

run_outputnl(dfResponseFile, dfSplineResponseFile, 'spline')
compute_response_info_coherence(dfSplineResponseFile)

