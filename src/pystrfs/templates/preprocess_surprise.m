srcRootDir='$SRC_ROOT_DIR'

tlabDir=fullfile(srcRootDir, 'tlab', 'src')

pystrfsDir=fullfile(tlabDir, 'pystrfs', 'matlab')
addpath(pystrfsDir)

surpriseDir=fullfile(tlabDir, 'surprise')
addpath(surpriseDir)

strflabDir=fullfile(srcRootDir, 'strflab')
addpath(genpath(strflabDir))

params = struct
params.domainFrequencyWidth = $DOMAIN_FREQUENCY_WIDTH
params.domainTimeWidth = $DOMAIN_TIME_WIDTH
params.domainGap = $DOMAIN_GAP

transforms = {'log'}

[stimFiles, md5s] = read_stimdata_from_file('$STIMDATA_FILE')
preprocess_surprise(md5s, '$PREPROCESS_FILE', params, transforms, '$OUTPUT_FILE')
exit

