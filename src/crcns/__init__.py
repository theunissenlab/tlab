import os
import re
import glob


CRCNS_DIR = '/auto/k6/tlab/crcns'   
CELLS_DIR = os.path.join(CRCNS_DIR, 'all_cells')
STIMS_DIR = os.path.join(CRCNS_DIR, 'all_stims')


class CRCNSBird:
   def __init__(self):
      self.name = None
      self.blocks = {}

   def __repr__(self):
      bstr = '\n'.join([str(x) for x in self.blocks.values()])
      return '%s:\n%s' % (self.name, bstr)


class CRCNSBlock:
   def __init__(self):
      self.name = None
      self.protocol = None
      self.units = {}
   def __repr__(self):
      return '%s: [%s]' % (self.name, ','.join(self.units.keys()))
       

class CRCNSWavFile:
   def __init__(self):
      self.fullPath = None
      self.md5 = None
      self.nSamples = None
      self.depth = None
      self.sampleRate = None
      self.stimType = None

   def __repr__(self):
      return self.md5


class CRCNSProtocol:
   def __init__(self):
      self.name = None
      self.stims = []
      self.stimTypes = []
      self.stimNumbers = []

   def __repr__(self):      
      return ','.join(self.stims)

   def __eq__(self, other):
      if isinstance(other, self.__class__):
         stStr1 = ','.join(self.stimTypes)
         stStr2 = ','.join(other.stimTypes)
         return str(self) == str(other) and stStr1 == stStr2
      else:
         return False

   def __ne__(self, other):
      return not self.__eq__(other)
   

class CRCNSCell:
   def __init__(self):
      self.path = None
      self.bird = None
      self.unit = None
      self.stimTypes = []
      self.region = None
      self.protocol = None
      self.spikeFiles = []

   def __repr__(self):
      pname = 'None'
      if self.protocol is not None:
         pname = self.protocol.name
      return '%s: region=%s, stimTypes=%s, protocol=%s' % (self.name(), self.region, ','.join(self.stimTypes), pname)
   def name(self):
      return '%s_%s' % (self.bird, self.unit)

class CRCNSSpikeFile:
   def __init__(self):
      self.name = None      
      self.trials = []

   def __repr__(self):
      return '%s: # trials=%d' % (self.name, len(self.trials))

   def __eq__(self, other):
      
      if not isinstance(other, self.__class__):
         return False
      
      if len(self.trials) != len(other.trials):
         return False

      for k,trials1 in enumerate(self.trials):
         stimes1 = trials1.spikes
         stimes2 = other.trials[k].spikes
         if len(stimes1) != len(stimes2):
            return False
         for n,ti1 in enumerate(stimes1):
            ti2 = stimes2[n]
            if ti1 != ti2:
               return False
      return True

   def __ne__(self, other):
      return not self.__eq__(other)


class CRCNSTrial:
   def __init__(self):
      self.number = None
      self.spikes = []


def get_wavdata():
   """ Get a list of CRCNSWavFile objects corresponding to each stimulus used """

   wdata = []
   dataFile = os.path.join(CRCNS_DIR, 'stim_data.csv')
   f = open(dataFile, 'r')
   for ln in f.readlines():
      ln = ln.strip()
      if len(ln) > 0:
         vals = ln.split(',')
         fname = vals[0]
         sampleRate = float(vals[1])
         depth = int(vals[2])
         nSamples = int(vals[3])
         stimType = vals[4]

         fvals = fname.split('.')
         md5Name = fvals[0]

         wfile = CRCNSWavFile()
         wfile.fullPath = os.path.join(STIMS_DIR, fname)
         wfile.md5 = md5Name
         wfile.nSamples = nSamples
         wfile.depth = depth
         wfile.sampleRate = sampleRate
         wfile.stimType = stimType

         wdata.append(wfile)

   return wdata


def get_all_cells():
   """ Get a list of all CRCNSCell objects on the filesystem """

   cellList = {}

   cellStimFile = os.path.join(CRCNS_DIR, 'cell_stim_classes.csv')

   f = open(cellStimFile, 'r')
   for ln in f.readlines():
      ln = ln.strip()
      if len(ln) > 0:
         svals = ln.split(',', 1)
         cellName = svals[0]
         stimTypes = svals[1].split(',')

         svals = cellName.split('_', 1)
         birdName = svals[0]
         unitName = svals[1]

         cell = CRCNSCell()
         cell.bird = birdName
         cell.unit = unitName
         cell.stimTypes = stimTypes
         cell.path = os.path.join(CRCNS_DIR, 'all_cells', cellName)

         cellList[cellName] = cell
   f.close()
   

   cellRegFile = os.path.join(CRCNS_DIR, 'cell_regions.csv')
   f = open(cellRegFile, 'r')
   for ln in f.readlines():
      ln = ln.strip()
      if len(ln) > 0:
         svals = ln.split(',')
         cellName = svals[0]
         region = svals[1]
         cellList[cellName].region = region      
   f.close()

   return cellList


def get_stim_md5s_from_dir(dirName, stimType):
   """ Get a list of MD5s of stimulus files for a given directory """

   gstr = os.path.join(dirName, 'stim*')
   slist = {}
   for fullName in glob.glob(gstr):
      (rootDir, stimName) = os.path.split(fullName)
      sindx = int(stimName[4:])
      spikeName = 'spike%d' % sindx
      if not os.path.exists(os.path.join(rootDir, spikeName)):
         print 'No spike file found for stim file %s' % fullName
      slist[sindx] = (stimName, spikeName)
   skeys = sorted(slist)

   md5Names = []
   stimTypes = []
   spikeFiles = []

   for k in skeys:
      stimName = slist[k][0]
      spikeName = slist[k][1]
      fname = os.path.join(dirName, stimName)
      
      f = open(fname, 'r')
      wavFile = f.read().strip()
      svals = wavFile.split('.')
      md5Name = svals[0]
      md5Names.append(md5Name)
      stimTypes.append(stimType)
      spikeFiles.append(spikeName)

   return (md5Names, stimTypes, spikeFiles)


def get_spike_times_from_file(fname):
   stimes = []
   f = open(fname, 'r')
   for ln in f.readlines():
      ln = ln.strip()
      if len(ln) > 0:
         tstrs = ln.split(' ')
         ti = [float(x)*1e-3 for x in tstrs]
         stimes.append(ti)
            
   f.close()
   return stimes


def identify_unique_protocols():
   """ Get a list of CRCNSProtocol objects that correspond to each cell """

   cellList = get_all_cells()      
   
   protocols = {}

   for (cellName, cell) in cellList.iteritems():

      cell.protocol = CRCNSProtocol()

      for stype in cell.stimTypes:

         stimDir = os.path.join(cell.path, stype)
         (stimMd5s, stimTypes, spikeFiles) = get_stim_md5s_from_dir(stimDir, stype)
         for k,stimMd5 in enumerate(stimMd5s):
            stype = stimTypes[k]
            cell.protocol.stims.append(stimMd5)
            cell.protocol.stimTypes.append(stype)
            cell.protocol.stimNumbers.append(k+1)

      pkey = str(cell.protocol)
      if pkey not in protocols:
         protocols[pkey] = cell.protocol

   print 'Identified %d unique protocols out of %d cells' % (len(protocols), len(cellList))

   #set names of protocols
   for k, proto in enumerate(protocols.values()):
      proto.name = 'CRCNS Protocol %d' % k

   #reset the protocols in each cell
   for cell in cellList.values():
      pkey = str(cell.protocol)
      cell.protocol = protocols[pkey]

   return (cellList, protocols.values())


def parse_unit_id(unitName):
   """ Given a unit name like 4_A, return ('4', 'A') """ 

   m = re.search('[\d]*', unitName)
   if m is not None:
      blockId = m.group(0)
   else:
      print 'Error: no block id for unit %s_%s!' % (birdName, unitName)
      
   m = re.search('[A-Z]', unitName)
   if m is not None:
      pairId = m.group(0)
   else:
      pairId = None

   #print '%s | %s | %s' % (unitName, blockId, pairId)

   return (blockId, pairId)


def get_bird_hierarchy():
   """ Get a list of CRCNSBird objects with their associated protocols """

   (cellList, protocols) = identify_unique_protocols()

   birdList = {}
  
   
   for (cellName, cell) in cellList.iteritems():

      birdName = cell.bird
      unitName = cell.unit

      #set up blocks and protcols for each bird/cell
      (blockId, pairId) = parse_unit_id(unitName)

      if birdName not in birdList:
         b = CRCNSBird()
         b.name = birdName
         birdList[birdName] = b

      b = birdList[birdName]

      if blockId in b.blocks:
         #check for matching protocols
         for otherUnit in b.blocks[blockId].units.values():
            if otherUnit.protocol != cell.protocol:
               print 'ERROR: Protocol mismatch within same block:'
               print '\t%s' % str(otherUnit)
               print '\t%s' % str(cell)
               exit()
      
      #re-check for existence of blockId (could have been changed)
      if blockId not in b.blocks:         
         blk = CRCNSBlock()
         blk.name = blockId
         blk.protocol = cell.protocol
         b.blocks[blockId] = blk
      
      if pairId is None:
         pairId = 'A'

      if pairId not in b.blocks[blockId].units:
         b.blocks[blockId].units[pairId] = cell

      protocolSpikeFiles = [None]*len(cell.protocol.stims)

      #get responses for cell
      for stype in cell.stimTypes:
         stimDir = os.path.join(cell.path, stype)
         (sfiles, stimFiles, spikeFiles) = get_stim_md5s_from_dir(stimDir, stype)
         for k,sfile in enumerate(spikeFiles):
            stimMd5 = sfiles[k]
            spikeFile = CRCNSSpikeFile()
            spikeFile.name = sfile
            spikeFilePath = os.path.join(stimDir, sfile)
            stimes = get_spike_times_from_file(spikeFilePath)
            for k,ti in enumerate(stimes):
               trial = CRCNSTrial()
               trial.number = k+1
               trial.spikes = ti
               spikeFile.trials.append(trial)

            stmp = [n for n,md5val in enumerate(cell.protocol.stims) if md5val == stimMd5]
            if len(stmp) == 0:
               print 'Cannot find stim corresponding to spike file %s' % spikeFilePath
            else:
               if len(stmp) == 1:
                  spikeFileIndex = stmp[0]
               else:
                  #in this case the stimulus protocol has repeats in the stimulus used,
                  #so we need to handle that (which sucks!)
                  #print 'Stimulus protocol repeats for cell %s, stimType=%s' % (str(cell), stype)
                  ntmp = [n for n,sval in enumerate(protocolSpikeFiles) if sval is None]
                  firstNoneIndex = ntmp[0]
                  ttmp = [sindx for sindx in stmp if sindx >= firstNoneIndex]
                  spikeFileIndex = ttmp[0]
            
               protocolSpikeFiles[spikeFileIndex] = spikeFile

      cell.spikeFiles = protocolSpikeFiles
   print 'Got responses for all cells...'
            
   return birdList


def find_stimulus_repeats():

   birdList = get_bird_hierarchy()

   cellsToExclude = []

   for bird in birdList.values():   
      for (blkId, blk) in bird.blocks.iteritems():         
         for (pairId, cell) in blk.units.iteritems():
            #aggregate spike files from similar stimuli
            spikeFilesForStimuli = {}
            
            for k,stimMd5 in enumerate(cell.protocol.stims):
               if stimMd5 not in spikeFilesForStimuli:
                  spikeFilesForStimuli[stimMd5] = []
               spikeFilesForStimuli[stimMd5].append(cell.spikeFiles[k])

            #check for duplicates
            for stimMd5,spikeFiles in spikeFilesForStimuli.iteritems():
               if len(spikeFiles) > 1:
                  areEqual = True
                  for k,sfile in enumerate(spikeFiles):
                     for m in range(k,len(spikeFiles)):
                        areEqual = areEqual and (sfile == spikeFiles[m])
                  if areEqual:
                     print 'Duplicate stimuli (%d) for cell %s with equal spike files!' % (len(spikeFiles), cell)
                  else:
                     print 'Duplicate stimuli (%d) for cell %s' % (len(spikeFiles), cell)
                  if cell.name() not in cellsToExclude:
                     cellsToExclude.append(cell.name())

   return cellsToExclude

