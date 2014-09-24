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
earlyStoppingGroups = $EARLY_STOPPING_GROUPS

strfLength = 60
tgResponseFile = '$TG_RESPONSE_FILE'
tgDistsResponseFile = '$TG_DISTS_RESPONSE_FILE'
tgSplineResponseFile = '$TG_SPLINE_RESPONSE_FILE'

transforms = {$TRANSFORMS}

run_threshgrad_cv(unitFile, preprocFile, stimClass, trainingGroups, validationGroups, earlyStoppingGroups, strfLength, transforms, $THRESHOLD, tgResponseFile, '$MODEL_TYPE')
compute_response_info_coherence(tgResponseFile)

if $COMPUTE_OUTPUTNL
   run_outputnl(tgResponseFile, tgDistsResponseFile, 'dists')
   compute_response_info_coherence(tgDistsResponseFile)

   run_outputnl(tgResponseFile, tgSplineResponseFile, 'spline')
   compute_response_info_coherence(tgSplineResponseFile)
end
