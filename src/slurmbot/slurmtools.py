import subprocess
import string
import time
import tempfile
import os
import re


def parse_time(timeStr):
   """Parses an squeue-like time, returns time in seconds"""

   sleft = timeStr.split('-')

   ndays = 0
   nhours = 0
   nmins = 0
   nsecs = 0

   if len(sleft) > 1:
      ndays = int(sleft[0])
      sleft = sleft[1]
   else:
      sleft = sleft[0]
   
   sleft = sleft.split(':')
   if len(sleft) == 3:
      nhours = int(sleft[0])
      nmins = int(sleft[1])
      nsecs = int(sleft[2])
   elif len(sleft) == 2:
      nmins = int(sleft[0])
      nsecs = int(sleft[1])
   
   totalTime = ndays*86400 + nhours*3600 + nmins*60 + nsecs
   return totalTime


def slurm_squeue(**keywords):

    sqcmd = ['squeue', '-h', '-o', '%N %P %Q %u %M %T %i %C']
    if 'username' in keywords:
        sqcmd.append('-u')
        sqcmd.append(keywords['username'])  

    proc = subprocess.Popen(sqcmd, stdout=subprocess.PIPE)
    outStr = proc.communicate()[0]
   
    jobsInfo = []
    for l in string.split(outStr, '\n'):
        if len(l) > 0:
            jstr = string.split(l)
            if len(jstr) < 8:
               jstr.insert(0, 'None')
            jinfo = {}
            jinfo['NODELIST'] = jstr[0]
            jinfo['PARTITION'] = jstr[1]
            jinfo['PRIORITY'] = int(jstr[2])
            jinfo['USER'] = jstr[3]
            jinfo['TIME'] = jstr[4]
            jinfo['STATE'] = jstr[5]
            jinfo['ID'] = int(jstr[6])
            jinfo['CPUS'] = int(jstr[7])
            jobsInfo.append(jinfo)

    return jobsInfo
    

def slurm_sbatch(cmdList, **sbatchParams):
    """Run a command using sbatch, returning the job id.
    
    Keyword arguments:
    partition -- the name of the partition to run on
    depends -- a comma-separated list of job id dependencies
    out -- file path to write stdout to
    err -- file path to write stderr to
    nodes -- the # of nodes this job will require
    qos -- the quality-of-service for this job
    cpus -- the # of CPUs required by this job
    script -- the name of the script file that will be written, defaults to a temporary file

    Returns the file path to the script that was executed.
    """

    sbatchCmds = []
    if 'partition' in sbatchParams:
        sbatchCmds.append('-p')
        sbatchCmds.append(sbatchParams['partition'])
    if 'depends' in sbatchParams:
        sbatchCmds.append('-d')
        sbatchCmds.append(sbatchParams['depends'])
    if 'out' in sbatchParams:
        sbatchCmds.append('-o')
        sbatchCmds.append(sbatchParams['out'])
    if 'err' in sbatchParams:
        sbatchCmds.append('-e')
        sbatchCmds.append(sbatchParams['err'])
    if 'nodes' in sbatchParams:
        sbatchCmds.append('-N')
        sbatchCmds.append('%d' % sbatchParams['nodes'])
    if 'qos' in sbatchParams:
        sbatchCmds.append('--qos=%d' % sbatchParams['qos'])
    if 'cpus' in sbatchParams:
        sbatchCmds.append('-c')
        sbatchCmds.append('%d' % sbatchParams['cpus'])

    #Create a temp batch script
    if 'script' in sbatchParams:
       ofname = sbatchParams['script']
       ofd = open(ofname, 'w')
    else:
       (ofid, ofname) = tempfile.mkstemp(prefix='slurmtools_')
       ofd = open(ofname, 'w')

    #write the temp batch script
    ofd.write('#!/bin/sh\n')
    sbHeader = '#SBATCH %s\n' % ' '.join(sbatchCmds)
    ofd.write(sbHeader)
    ofd.write(' '.join(cmdList))
    ofd.write('\n')
    ofd.close()

    jobId = slurm_sbatch_from_file(ofname)
    return jobId


def slurm_sbatch_from_file(fileName):
   
    finalCmds = ['sbatch', '-v', fileName]

    #run the sbatch process
    proc = subprocess.Popen(finalCmds, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    jobId = -1
    attempts = 0
    maxAttempts = 100
    while attempts < maxAttempts:

        time.sleep(0.1)

        odata = proc.stdout.read()
        m = re.search('batch\\sjob\\s\\d*', odata)
        if m is not None:
            mstr = m.group()
            jobId = int(mstr.split()[2])
            break

        attempts += 1

    return jobId
   
   

def slurm_srun(cmdList, **srunParams):
    """Run a command using srun, returning the job id.
    
    Keyword arguments:
    partition -- the name of the partition to run on
    depends -- a comma-separated list of job id dependencies
    out -- file path to write stdout to
    err -- file path to write stderr to
    nodes -- the # of nodes this job will require
    qos -- the quality-of-service for this job
    """

    srunCmds = ['srun', '-v']
    if 'partition' in srunParams:
        srunCmds.append('-p')
        srunCmds.append(srunParams['partition'])
    if 'depends' in srunParams:
        srunCmds.append('-d')
        srunCmds.append(srunParams['depends'])
    if 'out' in srunParams:
        srunCmds.append('-o')
        srunCmds.append(srunParams['out'])
    if 'err' in srunParams:
        srunCmds.append('-e')
        srunCmds.append(srunParams['err'])
    if 'nodes' in srunParams:
        srunCmds.append('-N')
        srunCmds.append(srunParams['nodes'])
    if 'qos' in srunParams:
        srunCmds.append('--qos=%d' % srunParams['qos'])
    
    for c in cmdList:
        srunCmds.append(c)

    (ofid, ofname) = tempfile.mkstemp(prefix='slurmtools_')
    proc = subprocess.Popen(srunCmds, stdout=ofid, stderr=ofid)
    
    jobId = -1
    attempts = 0
    maxAttempts = 30
    while attempts < maxAttempts:

        time.sleep(0.5)

        f = open(ofname)
        odata = f.read()
        f.close()    

        retCode = proc.poll()
        if retCode is not None:
            print 'srun may have failed, output below:'
            print odata
            return -1     

        m = re.search('jobid\\s\\d*', odata)
        if m is not None:
            mstr = m.group()
            jobId = int(mstr.split()[1])
            break

        attempts += 1

    return jobId
                    
   
def get_user_level_stats():

    jobInfos = slurm_squeue()
    userStats = {}
    for jinfo in jobInfos:
      
       uname = jinfo['USER']
       if uname not in userStats:
          userStats[uname] = {'PENDING':0,'RUNNING':0,'SUSPENDED':0,'TOTAL':0,
                              'PENDING_CPUS':0,'RUNNING_CPUS':0,'SUSPENDED_CPUS':0,'TOTAL_CPUS':0}
       ustat = userStats[uname]
       js = jinfo['STATE']
      
       if js in ustat:
          ncpus = jinfo['CPUS']
          jsc = '%s_CPUS' % js
          ustat[js] += 1
          ustat[jsc] += 1*ncpus
          ustat['TOTAL'] += 1
          ustat['TOTAL_CPUS'] += 1*ncpus

    return userStats




if __name__ == '__main__':
    
    """
    cmdList = ['sleep', '20']
    outFile = '/auto/fhome/mschachter/slurm_test_out.txt'
    errFile = '/auto/fhome/mschachter/slurm_test_err.txt'

    jobId = slurm_sbatch(cmdList, out=outFile, err=errFile, partition='all')
    print jobId

    #jobsInfo = slurm_squeue()
    """

    #jinfo = slurm_squeue(username='mschachter')
    #print jinfo

    ustats = get_user_level_stats()


    for (uname, ustat) in ustats.iteritems():

        print '%s:' % uname
        print '\tRunning CPUs: %d' % ustat['RUNNING_CPUS']
        print '\tSuspended CPUs: %d' % ustat['SUSPENDED_CPUS']
        print '\tPending CPUs: %d' % ustat['PENDING_CPUS']


