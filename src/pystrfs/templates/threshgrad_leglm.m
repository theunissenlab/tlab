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
outputFile = '$RESPONSE_FILE'

transforms = {$TRANSFORMS}

run_threshgrad_leglm(unitFile, preprocFile, stimClass, trainingGroups, validationGroups, strfLength, transforms, $THRESHOLD, outputFile, '$FAMILY_NAME')
compute_response_info_coherence(outputFile)

