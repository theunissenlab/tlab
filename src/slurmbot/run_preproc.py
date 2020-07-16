import os
import sys
import string
import numpy as np

import cellinfo
import matlabtemplate
import slurmbot
from modelcomps import *


def add_preproc_stft(cellName, slurmBot, stimType):
   
   mt = matlabtemplate.MatlabTemplate()
   mt.template_from_file('scripts/preproc.stft.template')
   
   fileDesc = 'preproc.stft.%s' % stimType

   sout = os.path.join(SLURM_DIR, '%s.%s.%%j.out.txt' % (cellName, fileDesc))
   serr = sout
   scriptName = os.path.join(SLURM_DIR, '%s.%s.sbatch' % (cellName, fileDesc))

   sparams = {'out': sout,
              'err': serr,
              'partition':'all',
              'cpus': 1,
              'script': scriptName}

   preprocDesc = 'stft.default'
   params = {'ROOT_DIR':ROOT_DIR, 'CELL_NAME':cellName, 'PREPROC_DESC':preprocDesc, 'STIM_TYPE':stimType}

   add_matlab_run(slurmBot, mt, params, sparams)      


def add_preproc_surprise(cellName, slurmBot, stimType):
   
   mt = matlabtemplate.MatlabTemplate()
   mt.template_from_file('scripts/preproc.surprise.template')
   
   fileDesc = 'preproc.surprise.%s' % stimType

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


def add_preproc_lyons(cellName, slurmBot, stimType):

   steps = [0.5, 0.25]
   earQs = [4, 6, 8]

   for step in steps:
      for earQ in earQs:
         tfParams = {}
         tfParams['TF_TYPE'] = 'lyons'
         tfParams['LYONS_STEP'] = step
         tfParams['LYONS_EARQ'] = earQ
         tfParams['LYONS_AGC'] = 0         
         add_preproc_lyons_single(cellName, slurmBot, stimType, tfParams)
         tfParams['LYONS_AGC'] = 1         
         add_preproc_lyons_single(cellName, slurmBot, stimType, tfParams)


def add_preproc_lyons_single(cellName, slurmBot, stimType, lyonsParams):
   
   mt = matlabtemplate.MatlabTemplate()
   mt.template_from_file('scripts/preproc.lyons.template')
      
   lyonsDesc = 'lyons.earq_%d.step_%0.2f.agc_%d' % (lyonsParams['LYONS_EARQ'],
                                                    lyonsParams['LYONS_STEP'],
                                                    lyonsParams['LYONS_AGC'])
   fileDesc = 'preproc.%s.%s' % (lyonsDesc, stimType)

   sout = os.path.join(SLURM_DIR, '%s.%s.%%j.out.txt' % (cellName, fileDesc))
   serr = sout
   scriptName = os.path.join(SLURM_DIR, '%s.%s.sbatch' % (cellName, fileDesc))

   sparams = {'out': sout,
              'err': serr,
              'partition':'all',
              'cpus': 1,
              'script': scriptName}

   params = {'ROOT_DIR':ROOT_DIR, 'CELL_NAME':cellName, 'PREPROC_DESC':lyonsDesc, 'STIM_TYPE':stimType}
   params.update(lyonsParams)

   add_matlab_run(slurmBot, mt, params, sparams)



if __name__ == '__main__':

   slurmBot = slurmbot.SlurmBot()
   slurmBot.pollInterval = 20
   slurmBot.maxJobs = 1

   cellInfos = cellinfo.get_all_cellinfo()

   if len(sys.argv) > 1:
      cellStr = sys.argv[1]
      cellNames = cellStr.split(',')
   else:
      cellNames = [cinfo.cellName for cinfo in cellInfos]
           
   for cname in cellNames:

      cinfo = filter(lambda x: x.cellName == cname, cellInfos)
      cinfo = cinfo[0]

      for stype in cinfo.stimTypes:
         #add_preproc_stft(cinfo.cellName, slurmBot, stype)
         #add_preproc_lyons(cinfo.cellName, slurmBot, stype)
         add_preproc_surprise(cinfo.cellName, slurmBot, stype)
   
   print '%d jobs to be run...' % len(slurmBot.jobs)
   slurmBot.run()



