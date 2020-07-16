import os
import sys
import string
import shutil
import numpy as np

import cellinfo
import matlabtemplate
import slurmbot
from modelcomps import *


def add_directfit_stft_single(cellName, slurmBot, nlType, preprocType, stimType):
   
   fileDesc = '%s.%s.%s' % (preprocType, stimType, nlType)

   #choose the MATLAB template and set up slurm parameters
   mt = matlabtemplate.MatlabTemplate()
   mt.template_from_file('scripts/outputnl.sim.template')
   
   sout = os.path.join(SLURM_DIR, '%s.outputnl.sim.%s.%%j.out.txt' % (cellName, fileDesc))
   serr = sout
   scriptName = os.path.join(SLURM_DIR, '%s.outputnl.sim.%s.sbatch' % (cellName,fileDesc))

   sparams = {'out': sout,
              'err': serr,
              'partition':'all',
              'cpus': 2,
              'script': scriptName}

   outputFilePath = os.path.join(DATA_DIR, 'outputNL', 'output')

   params = {'ROOT_DIR':ROOT_DIR, 'DATA_DIR':DATA_DIR, 'CELL_NAME':cellName, 'STIM_TYPE':stimType, 
             'PREPROC_TYPE':preprocType, 'NL_TYPE':nlType, 'OUTPUT_FILE_PATH':outputFilePath}

   add_matlab_run(slurmBot, mt, params, sparams)
   


if __name__ == '__main__':
   
   if len(sys.argv) < 3:
      print 'Usage:'
      print '  python run_outputnl.py run <cell1,cell2,...>'
      exit()

   cmd = sys.argv[1]

   if cmd == 'run':      
      cellNamesStr = sys.argv[2]
      cellNames = cellNamesStr.split(',')
      
      slurmBot = slurmbot.SlurmBot()
      slurmBot.maxJobs = 1

      cellInfos = cellinfo.get_all_cellinfo()

      # add optimizations      
      for cellName in cellNames:
         print 'Adding optimizations for cell %s...' % cellName
         cinfo = filter(lambda ci: ci.cellName == cellName, cellInfos)[0]
         preprocTypes = ['stft', 'lyons', 'surprise']
         stimTypes = ['conspecific', 'flatrip']
         nlTypes = ['exponential']
         for stimType in stimTypes:
            for preprocType in preprocTypes:
               for nlType in nlTypes:
                  add_directfit_stft_single(cinfo.cellName, slurmBot, nlType, preprocType, stimType)

      print 'Running %d optimizations...' % len(slurmBot.jobs)
      slurmBot.run()
      #slurmBot.run_and_wait()
      print 'Done!'
      
      exit()

