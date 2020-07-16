import os
import sys
import string
import numpy as np

import matlabtemplate
import slurmbot
from modelcomps import *


threshVals = [0.0, 0.25, 0.5, 0.75, 1]

def add_threshgrad_optimizations(cellName, nlType, slurmBot):
   
   #break training groups into 3 partitions after holding out the
   #top 2 sets in terms of info values
   nFolds = 3
   (trainingGroups, validationGroups) = get_cv_groups(cellName, nFolds)

   #choose the MATLAB template and set up slurm parameters
   mt = matlabtemplate.MatlabTemplate()
   mt.template_from_file('scripts/threshgrad.template')
   
   #do K-fold CV on each lambda2 value
   for threshVal in threshVals:
   
      for k in range(nFolds):
         
         ts = trainingGroups[k]
         vs = validationGroups[k]

         vgStr = ('%s' % vs.ravel()).strip().replace('\n', '')
         tgStr = ('%s' % ts.ravel()).strip().replace('\n', '')

         fileDesc = 'nlType_%s.thresh_%0.2f.JN_%d' % (nlType, threshVal, k)
         scriptName = os.path.join(SLURM_DIR, ('%s.threshgrad.%s.sbatch' % (cellName,fileDesc)))

         sout = os.path.join(SLURM_DIR, '%s.threshgrad.%s.%%j.txt' % (cellName, fileDesc))
         serr = sout

         sparams = {'out': sout,
                    'err': serr,
                    'partition':'highprio',
                    'cpus': 2,
                    'script':scriptName}

         params = {'ROOT_DIR':ROOT_DIR, 'CELL_NAME':cellName, 'STRF_LENGTH':STRF_LENGTH,
                   'THRESHOLD':threshVal, 'NL_TYPE':nlType,
                   'TRAINING_GROUPS':tgStr,
                   'VALIDATION_GROUPS':vgStr,
                   'FILE_DESC':fileDesc,
                   'EARLY_STOP': 1,
                   'MAX_ITERS': 8500}

         add_matlab_run(slurmBot, mt, params, sparams)

def add_threshgrad_final_optimization(cellName, nlType, threshVal, slurmBot, maxIters):
   
   #break training groups into 3 partitions after holding out the
   #top 2 sets in terms of info values
   nFolds = 3
   (ts, vs) = get_cv_groups(cellName, 1)

   #choose the MATLAB template and set up slurm parameters
   mt = matlabtemplate.MatlabTemplate()
   mt.template_from_file('scripts/threshgrad.template')
   
   vgStr = ('%s' % vs.ravel()).strip().replace('\n', '')
   tgStr = ('%s' % ts.ravel()).strip().replace('\n', '')

   fileDesc = 'best.nlType_%s' % (nlType)
   scriptName = os.path.join(SLURM_DIR, ('%s.threshgrad.%s.sbatch' % (cellName,fileDesc)))

   sout = os.path.join(SLURM_DIR, '%s.threshgrad.%s.%%j.txt' % (cellName, fileDesc))
   serr = sout

   sparams = {'out': sout,
              'err': serr,
              'partition':'all',
              'cpus': 2,
              'script':scriptName}

   params = {'ROOT_DIR':ROOT_DIR, 'CELL_NAME':cellName, 'STRF_LENGTH':STRF_LENGTH,
             'THRESHOLD':threshVal, 'NL_TYPE':nlType,
             'TRAINING_GROUPS':tgStr,
             'VALIDATION_GROUPS':vgStr,
             'FILE_DESC':fileDesc,
             'EARLY_STOP': 0,
             'MAX_ITERS': maxIters}

   
   add_matlab_run(slurmBot, mt, params, sparams)



if __name__ == '__main__':
   
   if len(sys.argv) < 3:
      print 'Usage:'
      print '  python run_threshgrad.py run <cell1,cell2,...>'
      print '  python run_threshgrad.py findbest <cell1,cell2,...>'
      print '  python run_threshgrad.py trainbest <cell1,cell2,...>'
      
      exit()

   cmd = sys.argv[1]


   if cmd == 'run':      
      cellNamesStr = sys.argv[2]
      cellNames = cellNamesStr.split(',')
      
      slurmBot = slurmbot.SlurmBot()

      # add optimizations      
      for cellName in cellNames:
         print 'Adding optimizations for cell %s...' % cellName
         add_threshgrad_optimizations(cellName, 'linear', slurmBot)
         add_threshgrad_optimizations(cellName, 'exponential', slurmBot)
         
      print 'Running optimizations...'
      slurmBot.run()
      print 'Done!'
      
      exit()

   if cmd == 'findbest':

      cellNamesStr = sys.argv[2]
      cellNames = cellNamesStr.split(',')
      
      nFolds = 3
      hyperParamVals = threshVals
      fileDescPrefix = 'strflab.stft.threshgrad'

      for cellName in cellNames:         
         print '-------------------\nHyperparams for cell %s' % cellName
         fileDescFmt = 'nlType_linear.thresh_%0.2f.JN_%d'
         outputFileName = 'bestparams.threshgrad.nlType_linear.txt'
         pick_best_hyperparams(cellName, nFolds, hyperParamVals, fileDescPrefix, fileDescFmt, outputFileName)
         fileDescFmt = 'nlType_exponential.thresh_%0.2f.JN_%d'
         outputFileName = 'bestparams.threshgrad.nlType_exponential.txt'
         pick_best_hyperparams(cellName, nFolds, hyperParamVals, fileDescPrefix, fileDescFmt, outputFileName)
      

      exit()

   if cmd == 'trainbest':

      slurmBot = slurmbot.SlurmBot()

      cellNamesStr = sys.argv[2]
      cellNames = cellNamesStr.split(',')

      for cellName in cellNames:
         cellDir = os.path.join(DATA_DIR, cellName)

         nlTypes = ['linear', 'exponential']
         threshVals = []
         for nlType in nlTypes:
            fname = 'bestparams.threshgrad.nlType_%s.txt' % nlType
            bestFile = os.path.join(cellDir, 'conspecific', 'output', fname)
            bfd = open(bestFile)
            bstr = bfd.read()
            bfd.close()
            bvals = bstr.split(',')
            maxIters = int(bvals[0])
            tval = float(bvals[1])
            print tval
            add_threshgrad_final_optimization(cellName, nlType, tval, slurmBot, maxIters)

      print 'About to run %d jobs...' % len(slurmBot.queue)
      slurmBot.run()
      exit()
