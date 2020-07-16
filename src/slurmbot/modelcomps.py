import os
import glob
import string
import numpy as np

import slurmbot
import slurmtools
import matlabtemplate

ROOT_DIR = '/auto/fdata/mschachter'
DATA_DIR = os.path.join(ROOT_DIR, 'data')
SLURM_DIR = os.path.join(ROOT_DIR, 'slurm')
STRF_LENGTH = 75


def read_perf_file(filePath):

   tperf = 0.0
   vperf = 0.0
   if os.path.exists(filePath):
      f = open(filePath, 'r')
      fstr = f.read().strip()
      f.close()
      darr = fstr.split(',')
      if len(darr) == 3:
         tperf = float(darr[1])
         vperf = float(darr[2])
   return (tperf, vperf)
   


def add_matlab_run(slurmBot, matlabTemplate, matlabParams, slurmParams, completionFileName = None):
   
   mstr = '\"%s\"' % matlabTemplate.fill_template(matlabParams)      

   cmds = ['matlabbg', mstr]
   slurmBot.add(cmds, slurmParams, completionFileName)


def read_info_vals(dataDir, cellName, stimType):

   odir = os.path.join(dataDir, cellName, stimType, 'output')
   ifile = os.path.join(odir, 'info_vals.txt')
   iinfo = {}
   if os.path.exists(ifile):      
      ifd = open(ifile)
      iflines = ifd.readlines()
      ifd.close()      
      for il in iflines:
         ivals = il.split('=')
         iinfo[ivals[0]] = ivals[1].strip()            

   else:
      print 'No such directory %s' % ifile

   return iinfo         
   

def get_cv_groups(cellName, nFolds, stimType):
   # get info values for each pair in the dataset
   infoVals = read_info_vals(DATA_DIR, cellName, stimType)
   groupInfos = []
   for k, val in enumerate(infoVals['meanGroupInfoVals'].split(',')):
      groupInfos.append( [k, float(val)] )
   
   # sort the info values from greatest to least
   groupInfos.sort(key=lambda v: v[1], reverse=True)

   #pick hold out groups as top two pairs
   gi = np.array(groupInfos)
   gi += 1
   
   holdOutGroups = gi[:2, 0]
   remainingGroups = gi[2:, 0]

   #use K-fold cross-validation on the remaining groups
   #to determine training and validation sets for LARS
   pSize = int(len(remainingGroups) / nFolds)

   if nFolds == 1:
      trainingGroups = remainingGroups
      validationGroups = holdOutGroups
   
   else:

      trainingGroups = []
      validationGroups = []
      for k in range(nFolds):
         a = k*pSize
         b = a + pSize
         vs = remainingGroups[a:b]
         ts = np.hstack((remainingGroups[:a], remainingGroups[b:]))
         trainingGroups.append(ts)
         validationGroups.append(vs)

   return (trainingGroups, validationGroups)


def pick_best_hyperparams(cellName, nFolds, hyperParamVals, fileDescPrefix, fileDescFmt, outputFileName):

   outputDir = os.path.join(DATA_DIR, cellName, 'conspecific', 'output')
   basefname = '%s.%s.txt' % (fileDescPrefix, fileDescFmt)

   trainPerfMeans = []
   trainPerfDevs = []
   validPerfMeans =[]
   validPerfDevs = []
   bestIterMeans = []
   scores = []

   for m,hval in enumerate(hyperParamVals):

      trainPerfRatios = []
      validPerfRatios = []
      bestIters = []
      for k in range(nFolds):
         fname = basefname % (hval, k)
         fpath = os.path.join(outputDir, fname)

         fid = open(fpath)
         fstr = fid.read()
         fid.close()
         vals = fstr.strip().split(',')
         
         bestIter = float(vals[0])
         tpr = float(vals[1])
         vpr = float(vals[2])
         
         bestIters.append(bestIter)
         trainPerfRatios.append(tpr)
         validPerfRatios.append(vpr)
         #print '%s: bestIter=%d, perf=%0.2f, valid=%0.2f' % (fname, bestIter, vpr, tpr)
         
         
      bestIters = np.array(bestIters)
      trainPerfRatios = np.array(trainPerfRatios)
      validPerfRatios = np.array(validPerfRatios)
      
      iterMean = int(bestIters.mean())
      iterDev = int(np.sqrt(bestIters.var()))
      trainMean = trainPerfRatios.mean()
      trainDev = np.sqrt(trainPerfRatios.var())
      validMean = validPerfRatios.mean()
      validDev = np.sqrt(validPerfRatios.var())

      score = 0.30*trainMean + 0.70*validMean

      numTrainingIters = iterMean + int(0.5*iterDev)

      bestIterMeans.append(numTrainingIters)
      trainPerfMeans.append(trainMean)
      trainPerfDevs.append(trainDev)
      validPerfMeans.append(validMean)
      validPerfDevs.append(validDev)
      scores.append( [m, score] )

      print 'hyperParam = %0.5f' % hval
      print '\tTraining Info: %0.2f +/- %0.4f' % (trainMean, trainDev)
      print '\tValidation Info: %0.2f +/- %0.4f' % (validMean, validDev)
      print '\t# of Iterations: %d +/- %0.0f' % (iterMean, iterDev)
      print '\tScore: %0.5f' % score
      print ''

   scores.sort(key=lambda v: v[1], reverse=True)
   bestIndex = scores[0][0]
   bestScore = scores[0][1]
   bestHyperParam = hyperParamVals[bestIndex]
   bestIteration = int(bestIterMeans[bestIndex])
   print 'Best Score: %f with hyperparam=%f, iters=%d' % (bestScore, bestHyperParam, bestIteration)
   
   bestFilePath = os.path.join(outputDir, outputFileName)
   ofd = open(bestFilePath, 'w')
   ofd.write('%d,%f\n' % (bestIteration, bestHyperParam))
   ofd.close()

if __name__ == '__main__':

   cellName = 'yy1617_4_A'
   (tg, vg) = get_cv_groups(cellName, 3)
      
   print tg
   print vg



   
