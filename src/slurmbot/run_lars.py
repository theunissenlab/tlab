import os
import sys
import string
import numpy as np

import matlabtemplate
import slurmbot
from modelcomps import *

lambda2Vals = [0.001, 0.01, 0.1, 0.5]

def add_lars_final_optimization(cellName, slurmBot, nlType, lambda2, maxIters):

   #split into training and validation groups (validation are top 2 w/ regards to info)
   (ts, vs) = get_cv_groups(cellName, 1)

   #choose the MATLAB template and set up slurm parameters
   mt = matlabtemplate.MatlabTemplate()
   mt.template_from_file('scripts/lars.template')

   vgStr = ('%s' % vs.ravel()).strip().replace('\n', '')
   tgStr = ('%s' % ts.ravel()).strip().replace('\n', '')

   fileDesc = 'best.nlType_%s' % nlType

   scriptName = os.path.join(SLURM_DIR, ('%s.lars.%s.sbatch' % (cellName,fileDesc)))

   sout = os.path.join(SLURM_DIR, '%s.lars.%s.%%j.txt' % (cellName,fileDesc))
   serr = sout

   sparams = {'out': sout,
              'err': serr,
              'partition':'all',
              'cpus': 2,
              'script':scriptName}

   params = {'ROOT_DIR':ROOT_DIR, 'CELL_NAME':cellName, 'STRF_LENGTH':STRF_LENGTH,
             'LAMBDA2':lambda2, 
             'MAX_ITERS': maxIters,
             'NL_TYPE': nlType,
             'TRAINING_GROUPS':tgStr,
             'VALIDATION_GROUPS':vgStr,
             'FILE_DESC':fileDesc,
             'FIND_BEST': 0}

   add_matlab_run(slurmBot, mt, params, sparams)



def add_lars_optimizations(cellName, nlType, slurmBot, outputDesc=None):
   
   #break training groups into 3 partitions after holding out the
   #top 2 sets in terms of info values
   nFolds = 3
   (trainingGroups, validationGroups) = get_cv_groups(cellName, nFolds)

   #choose the MATLAB template and set up slurm parameters
   mt = matlabtemplate.MatlabTemplate()
   mt.template_from_file('scripts/lars.template')
   
   #do K-fold CV on each lambda2 value
   for lambda2 in lambda2Vals:
   
      for k in range(nFolds):
         
         ts = trainingGroups[k]
         vs = validationGroups[k]

         skipRun = False
         if outputDesc is not None:
            #output file checking
            outputName = outputDesc % (lambda2, k)
            outputPath = os.path.join(DATA_DIR, cellName, 'conspecific', 'output', outputName)
            skipRun = os.path.exists(outputPath)
            if not skipRun:
               print 'Output does not exist: %s' % outputPath


         if not skipRun:
            vgStr = ('%s' % vs.ravel()).strip().replace('\n', '')
            tgStr = ('%s' % ts.ravel()).strip().replace('\n', '')

            fileDesc = 'nlType_%s.lambda_%0.3f.JN_%d' % (nlType, lambda2, k)         
            scriptName = os.path.join(SLURM_DIR, ('%s.lars.%s.sbatch' % (cellName,fileDesc)))

            sout = os.path.join(SLURM_DIR, '%s.lars.%s.%%j.txt' % (cellName,fileDesc))
            serr = sout

            sparams = {'out': sout,
                       'err': serr,
                       'partition':'highprio',
                       'cpus': 2,
                       'script':scriptName}

            params = {'ROOT_DIR':ROOT_DIR, 'CELL_NAME':cellName, 'STRF_LENGTH':STRF_LENGTH,
                      'LAMBDA2':lambda2, 
                      'MAX_ITERS': 400,
                      'NL_TYPE': nlType,
                      'TRAINING_GROUPS':tgStr,
                      'VALIDATION_GROUPS':vgStr,
                      'FILE_DESC':fileDesc,
                      'FIND_BEST': 1}

            add_matlab_run(slurmBot, mt, params, sparams)





if __name__ == '__main__':
   
   if len(sys.argv) < 3:
      print 'Usage:'
      print '  python run_lars.py run <cell1,cell2,...>'
      print '  python run_lars.py findbest <cell1,cell2,...>'
      print '  python run_lars.py trainbest <cell1,cell2,...>'
      exit()

   cmd = sys.argv[1]

   if cmd == 'run':      
      cellNamesStr = sys.argv[2]
      cellNames = cellNamesStr.split(',')
      
      slurmBot = slurmbot.SlurmBot()

      # add optimizations      
      for cellName in cellNames:
         print 'Adding optimizations for cell %s...' % cellName
         add_lars_optimizations(cellName, 'linear', slurmBot, 'strflab.stft.LARS.nlType_linear.lambda_%0.3f.JN_%d.txt')
         add_lars_optimizations(cellName, 'exponential', slurmBot, 'strflab.stft.LARS.nlType_exponential.lambda_%0.3f.JN_%d.txt')

      print 'Queuing optimizations (%d jobs)...' % len(slurmBot.queue)
      slurmBot.run()
      print 'Done!'
      
      exit()

   if cmd == 'findbest':

      cellNamesStr = sys.argv[2]
      cellNames = cellNamesStr.split(',')

      nFolds = 3
      fileDescPrefix = 'strflab.stft.LARS'

      for cellName in cellNames:
         print 'Best Hyperparams for Cell %s, linear:' % cellName
         pick_best_hyperparams(cellName, nFolds, lambda2Vals, fileDescPrefix, 'nlType_linear.lambda_%0.3f.JN_%d', 'bestparams.LARS.nlType_linear.txt')
         print ''
         print 'Best Hyperparams for Cell %s, exponential:' % cellName
         pick_best_hyperparams(cellName, nFolds, lambda2Vals, fileDescPrefix, 'nlType_exponential.lambda_%0.3f.JN_%d', 'bestparams.LARS.nlType_exponential.txt')
         print ''   
   
      exit()
   
   if cmd == 'trainbest':
      slurmBot = slurmbot.SlurmBot()
      slurmBot.maxJobs = 1

      cellNamesStr = sys.argv[2]
      cellNames = cellNamesStr.split(',')
      
      for cellName in cellNames:
         cellDir = os.path.join(DATA_DIR, cellName)

         nlTypes = ['linear', 'exponential']
         lambda2Vals = []
         for nlType in nlTypes:
            fname = 'bestparams.LARS.nlType_%s.txt' % nlType
            bestFile = os.path.join(cellDir, 'conspecific', 'output', fname)
            bfd = open(bestFile)
            bstr = bfd.read()
            bfd.close()
            bvals = bstr.split(',')
            maxIters = int(bvals[0])
            lambda2 = float(bvals[1])
            add_lars_final_optimization(cellName, slurmBot, nlType, lambda2, maxIters)

      slurmBot.run()

