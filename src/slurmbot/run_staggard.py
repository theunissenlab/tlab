import os
import sys
import string
import numpy as np

import cellinfo
import matlabtemplate
import slurmbot
from modelcomps import *

  

def add_staggard_stft_single(cellName, slurmBot, stimType, surprise=False):
   #split into training and validation groups (validation are top 2 w/ regards to info)
   (ts, vs) = get_cv_groups(cellName, 1, stimType)
   
   if surprise:
      fileDesc = 'tfType_stft.surprise.%s' % stimType
   else:
      fileDesc = 'tfType_stft.%s' % stimType

   #choose the MATLAB template and set up slurm parameters
   mt = matlabtemplate.MatlabTemplate()
   mt.template_from_file('scripts/staggard.stft.template')
   
   sout = os.path.join(SLURM_DIR, '%s.staggard.%s.%%j.out.txt' % (cellName, fileDesc))
   serr = sout
   scriptName = os.path.join(SLURM_DIR, '%s.staggard.%s.sbatch' % (cellName,fileDesc))

   sparams = {'out': sout,
              'err': serr,
              'partition':'all',
              'cpus': 2,
              'script': scriptName}

   vgStr = ('%s' % vs.ravel()).strip().replace('\n', '')
   tgStr = ('%s' % ts.ravel()).strip().replace('\n', '')

   outputFileName = 'strflab.staggard.%s.mat' % fileDesc

   preprocDesc = 'stft.default'

   doSurprise = 0
   if surprise:
      doSurprise = 1

   params = {'ROOT_DIR':ROOT_DIR, 'CELL_NAME':cellName, 'STRF_LENGTH':STRF_LENGTH,
             'TRAINING_GROUPS':tgStr, 'PREPROC_DESC':preprocDesc, 'STIM_TYPE':stimType, 
             'VALIDATION_GROUPS':vgStr, 'OUTPUT_FILE_NAME':outputFileName, 'SURPRISE':doSurprise}

   fullOutputFileName = os.path.join(DATA_DIR, cellName, stimType, 'output', outputFileName)

   add_matlab_run(slurmBot, mt, params, sparams, fullOutputFileName)

   
def add_staggard_lyons_single(cellName, slurmBot, stimType):
   
   #get best lyons params
   lyonsBestParamFile = os.path.join(DATA_DIR, cellName, stimType, 'output', 'best.lyons.params.txt')
   f = open(lyonsBestParamFile, 'r')
   lstr = f.read()
   f.close()
   larr = lstr.split(',')
   earQ = int(larr[0])
   step = float(larr[1])
   agc = int(larr[2])

   lyonsDesc = 'lyons.earq_%d.step_%0.2f.agc_%d' % (earQ, step, agc)

   #split into training and validation groups (validation are top 2 w/ regards to info)
   (ts, vs) = get_cv_groups(cellName, 1, stimType)
   
   fileDesc = 'tfType_lyons.%s' % stimType

   #choose the MATLAB template and set up slurm parameters
   mt = matlabtemplate.MatlabTemplate()
   mt.template_from_file('scripts/staggard.lyons.template')
   
   sout = os.path.join(SLURM_DIR, '%s.staggard.%s.%%j.out.txt' % (cellName, fileDesc))
   serr = sout
   scriptName = os.path.join(SLURM_DIR, '%s.staggard.%s.sbatch' % (cellName,fileDesc))

   sparams = {'out': sout,
              'err': serr,
              'partition':'all',
              'cpus': 2,
              'script': scriptName}

   vgStr = ('%s' % vs.ravel()).strip().replace('\n', '')
   tgStr = ('%s' % ts.ravel()).strip().replace('\n', '')

   outputFileName = 'strflab.staggard.%s.mat' % fileDesc

   preprocDesc = lyonsDesc

   params = {'ROOT_DIR':ROOT_DIR, 'CELL_NAME':cellName, 'STRF_LENGTH':STRF_LENGTH,
             'TRAINING_GROUPS':tgStr, 'PREPROC_DESC':preprocDesc, 'STIM_TYPE':stimType, 
             'VALIDATION_GROUPS':vgStr, 'OUTPUT_FILE_NAME':outputFileName,
             'LYONS_EARQ':earQ, 'LYONS_STEP':step, 'LYONS_AGC':agc}

   fullOutputFileName = os.path.join(DATA_DIR, cellName, stimType, 'output', outputFileName)

   add_matlab_run(slurmBot, mt, params, sparams, fullOutputFileName)
   



if __name__ == '__main__':
   
   if len(sys.argv) < 3:
      print 'Usage:'
      print '  python run_staggard.py run <cell1,cell2,...>'
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
         for stype in cinfo.stimTypes:
            add_staggard_lyons_single(cellName, slurmBot, stype)
            add_staggard_stft_single(cellName, slurmBot, stype, False)
            add_staggard_stft_single(cellName, slurmBot, stype, True)
            

      print 'Running %d optimizations...' % len(slurmBot.jobs)
      slurmBot.run()
      print 'Done!'
      
      exit()


