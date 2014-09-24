srcRootDir='$SRC_ROOT_DIR'

tlabDir=fullfile(srcRootDir, 'tlab', 'src')

pystrfsDir=fullfile(tlabDir, 'pystrfs', 'matlab')
addpath(pystrfsDir)

strflabDir=fullfile(srcRootDir, 'strflab')
addpath(genpath(strflabDir))

tfParams = struct
tfParams.nstd = 6
tfParams.fband = 125
tfParams.high_freq = 8000
tfParams.low_freq = 375         
tfParams.log = 0       
tfParams.dbnoise = 80           
tfParams.refpow = 6

[stimFiles, md5s] = read_stimdata_from_file('$STIMDATA_FILE')
preprocess_spectrogram(stimFiles, md5s, '$OUTPUT_FILE', tfParams)
exit

