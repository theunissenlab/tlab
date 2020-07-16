srcRootDir='$SRC_ROOT_DIR'

tlabDir=fullfile(srcRootDir, 'tlab', 'src')

pystrfsDir=fullfile(tlabDir, 'pystrfs', 'matlab')
addpath(pystrfsDir)

atDir = fullfile(srcRootDir, 'tlab', 'AuditoryToolbox')
addpath(atDir)

strflabDir=fullfile(srcRootDir, 'strflab')
addpath(genpath(strflabDir))

tfParams = struct
tfParams.low_freq = 100
tfParams.high_freq = 15000
tfParams.log = 0
tfParams.earQ = $LYONS_EARQ
tfParams.agc = $LYONS_AGC
tfParams.differ = $LYONS_AGC
tfParams.tau = 3
tfParams.step = $LYONS_STEP        

[stimFiles, md5s] = read_stimdata_from_file('$STIMDATA_FILE')
preprocess_lyons(stimFiles, md5s, '$OUTPUT_FILE', tfParams)
exit

