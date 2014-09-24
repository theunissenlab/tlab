import os
import sys
import string
import shutil
import numpy as np

import cellinfo
import matlabtemplate
import slurmbot
from modelcomps import *


def copy_lyons_best(cellName, stimType):

   outputDir = os.path.join(DATA_DIR, cellName, stimType, 'output')
   lyonsBestFile = os.path.join(outputDir, 'best.lyons.txt')
   
   f = open(lyonsBestFile, 'r')
   lyonsPathTxt = f.read().strip()
   f.close()

   (odir, lyonsNameTxt) = os.path.split(lyonsPathTxt)
   lyonsNameMat = string.replace(lyonsNameTxt, '.txt', '.mat')

   outputPrefix = 'strflab.tfType_lyons.%s' % stimType
   
   outputMat = '%s.mat' % outputPrefix
   f1 = os.path.join(outputDir, lyonsNameMat)
   f2 = os.path.join(outputDir, outputMat)
   print 'Copying %s to %s...' % (f1, f2)
   shutil.copy2(f1, f2)

   outputTxt = '%s.txt' % outputPrefix
   f1 = os.path.join(outputDir, lyonsNameTxt)
   f2 = os.path.join(outputDir, outputTxt)
   print 'Copying %s to %s...' % (f1, f2)
   shutil.copy2(f1, f2)      


def score_stft_and_lyons(cellName, stimType):

   outputDir = os.path.join(DATA_DIR, cellName, stimType, 'output')
   
   stftName = 'strflab.tfType_stft.%s.txt' % stimType
   
   lyonsBestFile = os.path.join(outputDir, 'best.lyons.txt')
   
   f = open(lyonsBestFile, 'r')
   lyonsName = f.read().strip()
   f.close()
   
   (stftTrainPerf, stftValidPerf) = read_perf_file(os.path.join(outputDir, stftName))
   (lyonsTrainPerf, lyonsValidPerf) = read_perf_file(lyonsName)

   stftScore = 0.3*stftTrainPerf + 0.7*stftValidPerf
   lyonsScore = 0.3*lyonsTrainPerf + 0.7*lyonsValidPerf

   return (stftScore, lyonsScore)

def compare_lyons_vs_stft(cellNames, region, stimType, allCellInfos):

   regCellInfos = filter(lambda x: x.region == region, allCellInfos)

   stftScores = []
   lyonsScores = []

   for cname in cellNames:
      clst = filter(lambda x: x.cellName == cname, regCellInfos)
      if len(clst) > 0:
         cinfo = clst[0]
         if stimType in cinfo.stimTypes:         
            (stftScore, lyonsScore) = score_stft_and_lyons(cinfo.cellName, stimType)
            stftScores.append(stftScore)
            lyonsScores.append(lyonsScore)
   
   stftScores = np.array(stftScores)
   lyonsScores = np.array(lyonsScores)

   diffs = lyonsScores - stftScores

   return diffs
   

def find_best_lyons(cellName, stimType):

   steps = [0.5, 0.25]
   earQs = [4, 6, 8]

   bestScore = -1
   bestParams = {'LYONS_STEP': -1, 'LYONS_EARQ': -1, 'LYONS_AGC': 0,'OUTPUT_FILE': None}

   print 'Lyons Data | %s (%s)' % (cellName, stimType)

   agcPerfSum = 0.0

   for step in steps:
      for earQ in earQs:
         agcInc = 0.0
         for agc in range(2):
            lyonsDesc = 'lyons.earq_%d.step_%0.2f.agc_%d' % (earQ, step, agc)
            fileDesc = 'tfType_%s.%s' % (lyonsDesc, stimType)   
            outputFileName = 'strflab.%s.txt' % fileDesc

            fullOutputFileName = os.path.join(DATA_DIR, cellName, stimType, 'output', outputFileName)

            if not os.path.exists(fullOutputFileName):
               print 'Missing file: %s' % fullOutputFileName

            else:

               f = open(fullOutputFileName, 'r')
               fstr = f.read().strip()
               f.close()
               
               darr = fstr.split(',')
               if len(darr) < 3:
                  print 'Corrupt file: %s' % fullOutputFileName 
               else:
                  trainingPerf = float(darr[1])
                  validationPerf = float(darr[2])
                  score = 0.3*trainingPerf + 0.7*validationPerf
                   
                  agcInc += (2*agc - 1)*score
  
                  print '\tEarQ=%d, Step=%0.2f, AGC=%d' % (earQ, step, agc)
                  print '\t%s' % fullOutputFileName
                  print '\t\ttraining=%0.4f | validation=%0.4f' % (trainingPerf, validationPerf)
                  print '\t\tScore = %0.4f' % score

                  if score > bestScore:
                     bestScore = score
                     bestParams['LYONS_EARQ'] = earQ
                     bestParams['LYONS_STEP'] = step
                     bestParams['LYONS_AGC'] = agc
                     bestParams['OUTPUT_FILE'] = fullOutputFileName

         agcPerfSum += agcInc

   agcGain = agcPerfSum / 6
   print 'Best: EarQ = %d  | step = %0.2f  | agc = %d  | Score = %0.4f' % (bestParams['LYONS_EARQ'], bestParams['LYONS_STEP'], bestParams['LYONS_AGC'], bestScore)
   print 'AGC Gain: %0.3f' % agcGain
   print '--------------------------'

   bestOutputFileName = os.path.join(DATA_DIR, cellName, stimType, 'output', 'best.lyons.txt')
   f = open(bestOutputFileName, 'w')
   f.write('%s\n' % bestParams['OUTPUT_FILE'])
   f.close()
   
   bestOutputParamsFileName = os.path.join(DATA_DIR, cellName, stimType, 'output', 'best.lyons.params.txt')
   f = open(bestOutputParamsFileName, 'w')
   f.write('%d,%0.2f,%d\n' % (bestParams['LYONS_EARQ'], bestParams['LYONS_STEP'], bestParams['LYONS_AGC']))
   f.close()



def add_directfit_lyons_all(cellName, slurmBot, stimType):

   steps = [0.5, 0.25]
   earQs = [4, 6, 8]

   for step in steps:
      for earQ in earQs:
         tfParams = {}
         tfParams['TF_TYPE'] = 'lyons'
         tfParams['LYONS_STEP'] = step
         tfParams['LYONS_EARQ'] = earQ
         tfParams['LYONS_AGC'] = 0         
         add_directfit_lyons_single(cellName, slurmBot, stimType, tfParams)
         tfParams['LYONS_AGC'] = 1         
         add_directfit_lyons_single(cellName, slurmBot, stimType, tfParams)


def add_directfit_lyons_single(cellName, slurmBot, stimType, lyonsParams):
   #split into training and validation groups (validation are top 2 w/ regards to info)
   (ts, vs) = get_cv_groups(cellName, 1, stimType)
   
   lyonsDesc = 'lyons.earq_%d.step_%0.2f.agc_%d' % (lyonsParams['LYONS_EARQ'],
                                                    lyonsParams['LYONS_STEP'],
                                                    lyonsParams['LYONS_AGC'])
   fileDesc = 'tfType_%s.%s' % (lyonsDesc, stimType)

   #choose the MATLAB template and set up slurm parameters
   mt = matlabtemplate.MatlabTemplate()
   mt.template_from_file('scripts/directfit.lyons.template')
   
   sout = os.path.join(SLURM_DIR, '%s.directfit.%s.%%j.out.txt' % (cellName, fileDesc))
   serr = sout
   scriptName = os.path.join(SLURM_DIR, '%s.directfit.%s.sbatch' % (cellName,fileDesc))

   sparams = {'out': sout,
              'err': serr,
              'partition':'all',
              'cpus': 2,
              'script': scriptName}

   vgStr = ('%s' % vs.ravel()).strip().replace('\n', '')
   tgStr = ('%s' % ts.ravel()).strip().replace('\n', '')

   outputFileName = 'strflab.%s.mat' % fileDesc

   preprocDesc = lyonsDesc

   params = {'ROOT_DIR':ROOT_DIR, 'CELL_NAME':cellName, 'STRF_LENGTH':STRF_LENGTH,
             'TRAINING_GROUPS':tgStr, 'PREPROC_DESC':preprocDesc, 'STIM_TYPE':stimType, 
             'VALIDATION_GROUPS':vgStr, 'OUTPUT_FILE_NAME':outputFileName}
   params.update(lyonsParams)
   fullOutputFileName = os.path.join(DATA_DIR, cellName, stimType, 'output', outputFileName)

   add_matlab_run(slurmBot, mt, params, sparams, fullOutputFileName)
   

def add_directfit_stft_single(cellName, slurmBot, stimType):
   #split into training and validation groups (validation are top 2 w/ regards to info)
   (ts, vs) = get_cv_groups(cellName, 1, stimType)
   
   fileDesc = 'tfType_stft.%s' % stimType

   #choose the MATLAB template and set up slurm parameters
   mt = matlabtemplate.MatlabTemplate()
   mt.template_from_file('scripts/directfit.stft.template')
   
   sout = os.path.join(SLURM_DIR, '%s.directfit.%s.%%j.out.txt' % (cellName, fileDesc))
   serr = sout
   scriptName = os.path.join(SLURM_DIR, '%s.directfit.%s.sbatch' % (cellName,fileDesc))

   sparams = {'out': sout,
              'err': serr,
              'partition':'all',
              'cpus': 2,
              'script': scriptName}

   vgStr = ('%s' % vs.ravel()).strip().replace('\n', '')
   tgStr = ('%s' % ts.ravel()).strip().replace('\n', '')

   outputFileName = 'strflab.%s.mat' % fileDesc

   preprocDesc = 'stft.default'

   params = {'ROOT_DIR':ROOT_DIR, 'CELL_NAME':cellName, 'STRF_LENGTH':STRF_LENGTH,
             'TRAINING_GROUPS':tgStr, 'PREPROC_DESC':preprocDesc, 'STIM_TYPE':stimType, 
             'VALIDATION_GROUPS':vgStr, 'OUTPUT_FILE_NAME':outputFileName}

   fullOutputFileName = os.path.join(DATA_DIR, cellName, stimType, 'output', outputFileName)

   add_matlab_run(slurmBot, mt, params, sparams, fullOutputFileName)
   

def add_directfit_surprise_single(cellName, slurmBot, stimType):
   #split into training and validation groups (validation are top 2 w/ regards to info)
   (ts, vs) = get_cv_groups(cellName, 1, stimType)
   
   fileDesc = 'tfType_surprise.%s' % stimType

   #choose the MATLAB template and set up slurm parameters
   mt = matlabtemplate.MatlabTemplate()
   mt.template_from_file('scripts/directfit.surprise.template')
   
   sout = os.path.join(SLURM_DIR, '%s.directfit.%s.%%j.out.txt' % (cellName, fileDesc))
   serr = sout
   scriptName = os.path.join(SLURM_DIR, '%s.directfit.%s.sbatch' % (cellName,fileDesc))

   sparams = {'out': sout,
              'err': serr,
              'partition':'all',
              'cpus': 2,
              'script': scriptName}

   vgStr = ('%s' % vs.ravel()).strip().replace('\n', '')
   tgStr = ('%s' % ts.ravel()).strip().replace('\n', '')

   outputFileName = 'strflab.%s.mat' % fileDesc

   preprocDesc = 'stft.default'

   params = {'ROOT_DIR':ROOT_DIR, 'CELL_NAME':cellName, 'STRF_LENGTH':STRF_LENGTH,
             'TRAINING_GROUPS':tgStr, 'PREPROC_DESC':preprocDesc, 'STIM_TYPE':stimType, 
             'VALIDATION_GROUPS':vgStr, 'OUTPUT_FILE_NAME':outputFileName}

   fullOutputFileName = os.path.join(DATA_DIR, cellName, stimType, 'output', outputFileName)

   add_matlab_run(slurmBot, mt, params, sparams, fullOutputFileName)



if __name__ == '__main__':
   
   if len(sys.argv) < 3:
      print 'Usage:'
      print '  python run_directfit.py run <cell1,cell2,...>'
      print '  python run_directfit.py lyonsbest <cell1,cell2,...>'
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
            #add_directfit_stft_single(cellName, slurmBot, stype)
            #add_directfit_lyons_all(cellName, slurmBot, stype)
            add_directfit_surprise_single(cellName, slurmBot, stype)

      print 'Running %d optimizations...' % len(slurmBot.jobs)
      slurmBot.run()
      print 'Done!'
      
      exit()

   elif cmd == 'lyonsbest':

      cellNamesStr = sys.argv[2]
      cellNames = cellNamesStr.split(',')
      
      cellInfos = cellinfo.get_all_cellinfo()

      for cellName in cellNames:
         cinfo = filter(lambda ci: ci.cellName == cellName, cellInfos)[0]
         for stype in cinfo.stimTypes:
            find_best_lyons(cinfo.cellName, stype)

      exit()

   elif cmd == 'lyonscompare':
      cellNamesStr = sys.argv[2]
      cellNames = cellNamesStr.split(',')
      
      regions = ['mld', 'OV', 'L']
      stimTypes = ['conspecific', 'flatrip', 'songrip']

      allCellInfos = cellinfo.get_all_cellinfo()

      for reg in regions:
         for stype in stimTypes:
            print '---------------------'
            print 'Region=%s | Stim=%s' % (reg, stype)
            diffs = compare_lyons_vs_stft(cellNames, reg, stype, allCellInfos)
            print diffs
            print 'Avg Diff: %0.4f +/- %0.4f' % (diffs.mean(), diffs.std())

   elif cmd == 'lyonscopybest':
      cellNamesStr = sys.argv[2]
      cellNames = cellNamesStr.split(',')

      allCellInfos = cellinfo.get_all_cellinfo()

      for cname in cellNames:
         cinfo = filter(lambda x: x.cellName == cname, allCellInfos)[0]
         for stype in cinfo.stimTypes:
            copy_lyons_best(cname, stype)
      
      








      

