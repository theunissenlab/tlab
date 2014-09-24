import sys
import string

from os import system
import os.path

import tempfile

import shutil

import numpy as np

from pyevolve import G1DList
from pyevolve import GSimpleGA
from pyevolve import Selectors
from pyevolve import Statistics
from pyevolve import DBAdapters
from pyevolve import Initializators
import pyevolve

import slurmbot


class LyonsObjectiveFunction:

   slurmBot = slurmbot.SlurmBot()
   
   paramNames = ['EARQ', 'EAR_STEP']

   def __init__(self, birdName, cellName):
      self.slurmBot.setTemplateFile('lyons_objfunc.template')
      self.birdName = birdName
      self.cellName = cellName

   def getNumParams(self):
      return len(self.paramNames)
   
   def score(self, chromosome):

      if len(chromosome) != len(self.paramNames):
         print 'Cannot score, length of input params does not match up to expected value: %d != %d' % (len(chromosome), len(self.paramNames))
         return -1.0

      paramSet = {}
      for k,val in enumerate(chromosome):
         print val
         pname = self.paramNames[k]
         paramSet[pname] = val

      paramSet['FILE_DESC'] = '%s_%s' % (self.birdName, self.cellName)

      scoreFileName = 'score.txt'
      paramSet['SCORE_FILE_NAME'] = scoreFileName

      tempDir = tempfile.mkdtemp()
      paramSet['TEMP_DIR'] = tempDir

      paramSet['BIRD_NAME'] = self.birdName
      paramSet['CELL_NAME'] = self.cellName

      print paramSet
      
      self.slurmBot.setParameterSets([paramSet,])
      self.slurmBot.run()

      scorePath = os.path.join(tempDir, scoreFileName)

      if not os.path.isfile(scorePath):
         print 'Can not locate score path file: %s' % scorePath
         return -1.0
      f = open(scorePath, 'r')
      fstr = f.read()
      f.close()
      scoreVal = float(fstr.strip())
      print 'scoreVal=%f' % scoreVal

      shutil.rmtree(tempDir)
      
      return scoreVal


def FixedStepInitializator(genome, **args):
   
   genome.clearList()
   
   valRanges = genome.getParam("valueRanges")
   stepSizes = genome.getParam("stepSizes")
   for k,valRange in enumerate(valRanges):
      step = stepSizes[k]
      rngLen = valRange[1] - valRange[0]
      numPoss = int(round(rngLen / step))
      rnum = np.random.randint(0, numPoss)
      rval = rnum*step + valRange[0]
      genome.append(rval)
     

def FixedStepMutator(genome, **args):
   
   numMutations = 0
   valRanges = genome.getParam("valueRanges")
   mutationRates = genome.getParam("mutationRates")
   stepSizes = genome.getParam("stepSizes")

   for k,paramVal in enumerate(genome):
      mrate = mutationRates[k]
      rng = valRanges[k]
      rnum = np.random.rand()
      if rnum < mrate:
         stepSize = stepSizes[k]
         if np.random.rand() < 0.5:
            stepSize = -stepSize
         oldVal = genome[k]
         genome[k] = genome[k] + stepSize
         if genome[k] > rng[1]:
            genome[k] = rng[1]            
         if genome[k] < rng[0]:
            genome[k] = rng[0]
         numMutations += (abs(oldVal - genome[k]) > 0)

   return numMutations


if __name__ == '__main__':

   # Enable the pyevolve logging system
   #pyevolve.logEnable()

   lyonsObjFunc = LyonsObjectiveFunction('pipi1112', '1_A')
   numParams = lyonsObjFunc.getNumParams()

   #set mutation rates
   mRates=0.35*np.squeeze(np.ones([1, numParams]))

   #set parameter ranges
   valRanges = np.ones([numParams, 2])
   valRanges[0, :] = [4.0, 13.0]
   valRanges[1, :] = [0.20, 0.80]

   #set step sizes
   stepSzs = np.array([0.5, 0.05])

   # Genome instance
   genome = G1DList.G1DList(numParams)
   genome.setParams(stepSizes=stepSzs,
                    mutationRates=mRates,
                    valueRanges=valRanges)
   genome.initializator.set(FixedStepInitializator)
   genome.mutator.set(FixedStepMutator)
   genome.evaluator.set(lyonsObjFunc.score)
      
   # Genetic Algorithm Instance
   ga = GSimpleGA.GSimpleGA(genome)
   ga.setPopulationSize(5)
   ga.selector.set(Selectors.GRouletteWheel)   
   ga.setGenerations(60)
   ga.terminationCriteria.set(GSimpleGA.ConvergenceCriteria)

   sqlite_adapter = DBAdapters.DBSQLite(identify="ex1", resetDB=True)
   ga.setDBAdapter(sqlite_adapter)

   ga.evolve(freq_stats=1)

   print ga.bestIndividual()

