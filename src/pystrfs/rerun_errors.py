import os
import subprocess


cprts = os.path.split(__file__)
curScript = cprts[1]
print 'Current script: %s' % curScript

slurmDir = '/auto/k6/mschachter/pystrfs/slurm/*'
cmdStr = 'grep -l \'Error in\' %s' % slurmDir
proc = subprocess.Popen(cmdStr, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
proc.wait()
lns = proc.stdout.readlines()
for ln in lns:
    fpath = ln.strip()
    fprts = os.path.split(fpath)
    rootDir = fprts[0]
    fileName = fprts[1]

    if fileName != curScript:

        sprts = fileName.split('.')
        rprts = sprts[0:len(sprts)-3]
        if len(rprts) > 3:
            rootName = '.'.join(rprts)
            sbatchName = '%s.sbatch' % rootName
            sbatchPath = os.path.join(rootDir, sbatchName)

            print 'Running script %s...' % sbatchName
            subprocess.call(['sbatch', '-v', sbatchPath])
            os.remove(fpath)
