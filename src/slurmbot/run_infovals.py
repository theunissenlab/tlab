import os
import sys
import string
import numpy as np

import cellinfo
import matlabtemplate
import slurmbot
from modelcomps import *


def add_infoval(cellName, slurmBot, stimType):
   
   mt = matlabtemplate.MatlabTemplate()
   mt.template_from_file('scripts/groupinfo.template')
   
   fileDesc = 'infovals.%s' % stimType

   sout = os.path.join(SLURM_DIR, '%s.%s.%%j.out.txt' % (cellName, fileDesc))
   serr = sout
   scriptName = os.path.join(SLURM_DIR, '%s.%s.sbatch' % (cellName, fileDesc))

   sparams = {'out': sout,
              'err': serr,
              'partition':'all',
              'cpus': 1,
              'script': scriptName}

   params = {'ROOT_DIR':ROOT_DIR, 'CELL_NAME':cellName, 'STIM_TYPE':stimType}

   add_matlab_run(slurmBot, mt, params, sparams)      


if __name__ == '__main__':

   slurmBot = slurmbot.SlurmBot()
   slurmBot.pollInterval = 20
   slurmBot.maxJobs = 10

   cellInfos = cellinfo.get_all_cellinfo()
   
   for cinfo in cellInfos:
      for stype in cinfo.stimTypes:
         add_infoval(cinfo.cellName, slurmBot, stype)
   
   print '%d jobs to be run...' % len(slurmBot.jobs)
   slurmBot.run_and_wait()



