import sys
import string
import os
import time

from pystrfs.slurmtools import *


class SlurmJob:

    def __init__(self):
        self.cmds = None
        self.params = None
        self.completion_file = None
        self.running = False
        self.start_times = 0
        self.completed = False
        self.job_id = -1


class SlurmBot:

    def __init__(self):
        self.jobs = []
        self.poll_interval = 37.0
        self.max_jobs = 10


    def add(self, cmdStrs, sbatchParams, completion_fileName = None):
        sj = SlurmJob()
        sj.cmds = cmdStrs
        sj.params = sbatchParams
        sj.completion_file = completion_fileName
        self.jobs.append(sj)


    def clear(self):
        del self.jobs[0:len(self.jobs)]


    def get_queued_jobs(self):
        return filter(lambda sj: not sj.completed and not(sj.running), self.jobs)


    def get_running_jobs(self):
        return filter(lambda sj: sj.running and not(sj.completed), self.jobs)


    def get_uncompleted_jobs(self):
        ujList = []
        for sj in self.jobs:
            if not sj.completed or not os.path.exists(sj.completion_file):
                ujList.append(sj)
        return ujList

    def mark_completed_jobs(self):
        for sj in self.jobs:
            if sj.completion_file != None and os.path.exists(sj.completion_file):
                sj.completed = True

    def update_running_jobs(self):

        jobInfo = slurm_squeue()
        for rj in self.get_running_jobs():
            jinfo = filter(lambda j: j['ID'] == rj.job_id, jobInfo)
            if len(jinfo) == 0:
                rj.running = False
                if rj.completion_file is None:
                    print 'Marking job %d as completed (no completion file specified)' % rj.job_id
                    rj.completed = True
                elif rj.completion_file != None and os.path.exists(rj.completion_file):
                    print 'Marking job %d as completed' % rj.job_id
                    rj.completed = True


    def run_and_wait(self, ignoreCompleted = False):

        if not ignoreCompleted:
            self.mark_completed_jobs()
            uj = self.get_uncompleted_jobs()
            nComp = len(self.jobs) - len(uj)
            print 'Marked %d jobs as already completed' % nComp

        self.update_running_jobs()
        queuedJobs = self.get_queued_jobs()
        runningJobs = self.get_running_jobs()

        while len(queuedJobs) > 0 or len(runningJobs) > 0:

            self.update_running_jobs()

            queuedJobs = self.get_queued_jobs()
            runningJobs = self.get_running_jobs()

            print '# of queued jobs: %d' % len(queuedJobs)
            print '# of running jobs: %d' % len(runningJobs)

            nJobsAvailable = min(len(queuedJobs), self.max_jobs - len(runningJobs))
            print '# of available job space to run: %d' % nJobsAvailable
            if nJobsAvailable > 0 and len(queuedJobs) > 0:

                for k in range(nJobsAvailable):
                    #Run jobs
                    sj = queuedJobs[k]
                    job_id = slurm_sbatch(sj.cmds, **sj.params)
                    sj.running = True
                    sj.job_id = job_id
                    if sj.start_times == 0:
                        print 'Started job with id %d' % job_id
                    else:
                        print 'Restarted job with id %d' % job_id
                    sj.start_times += 1

            time.sleep(self.poll_interval)


    def run(self):

        for sj in self.get_queued_jobs():
            cmdStrs = sj.cmds
            slurmParams = sj.params
            job_id = slurm_sbatch(cmdStrs, **slurmParams)
            sj.job_id = job_id
            print 'Started job with id %d' % job_id


if __name__ == '__main__':

    sb = SlurmBot()
    sb.max_jobs = 5

    cmdStr = ['sleep', '10s']
    sparams = {'out': '/auto/fhome/mschachter/slurm_%j_out.txt',
               'err': '/auto/fhome/mschachter/slurm_%j_err.txt',
               'partition':'all'}

    print 'adding processes'
    for k in range(10):
        sb.add(cmdStr, sparams)

    print 'running and waiting'
    sb.run_and_wait()

    print 'done'
