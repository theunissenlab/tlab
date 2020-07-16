import sys
import string
import os
from math import isnan

class CellInfo:
   
   def __init__(self):      
      self.birdName = None
      self.cellName = None
      self.region = None
      self.infoMean = None
      self.infoUpper = None
      self.infoLower = None
      self.spikeRate = None


class CellOrganizer:

   def __init__(self):

      self.knownRegions = ['L', 'L1', 'L2a', 'L2b', 'L3', 'OV', 'mld', 'CM']
      self.cellsByRegion = {}
      for r in self.knownRegions:
         self.cellsByRegion[r] = []

   def mergeInfoData(self, dataDir, cellDataFile):
      
      stimTypes = ['conspecific', 'zfsongs']

      f = open(cellDataFile, 'r')
      for l in f.readlines():         
         [birdName, cellName, region] = string.split(l.strip(), ',')
      
         dirName = '%s_%s' % (birdName, cellName)   
         fullPath = os.path.join(dataDir, dirName)

         cellInfo = CellInfo()
         cellInfo.birdName = birdName
         cellInfo.cellName = cellName
         cellInfo.region = region

         found = False
         for stype in stimTypes:

            outputDir = os.path.join(fullPath, stype, 'output')
            if os.path.isdir(outputDir):
               ifile = os.path.join(outputDir, 'info_vals.csv')
               if os.path.isfile(ifile):
                  found = True
                  ifd = open(ifile)
                  istr = ifd.read()                  
                  ifd.close()

                  ivals = string.split(istr.strip(), ',')
                  if len(ivals) == 4:
                     try: 
                        cellInfo.spikeRate = float(ivals[0])
                     except ValueError as e:
                        cellInfo.spikeRate = -1.0
                     try:
                        cellInfo.infoLower = float(ivals[1])
                     except ValueError as e:
                        cellInfo.infoLower = -1.0
                     try:
                        cellInfo.infoMean = float(ivals[2])   
                     except ValueError as e:
                        cellInfo.infoMean = -1.0
                     try:
                        cellInfo.infoUpper = float(ivals[3])         
                     except ValueError as e:
                        cellInfo.infoUpper = -1.0

         if not found:
            print 'Could not find info file for %s_%s' % (cellInfo.birdName, cellInfo.cellName)
            cellInfo.spikeRate = -1.0
            cellInfo.infoMean = -1.0
            cellInfo.infoUpper = -1.0
            cellInfo.infoLower = -1.0

         if isnan(cellInfo.spikeRate):
            cellInfo.spikeRate = -1.0
         if isnan(cellInfo.infoMean):
            cellInfo.infoMean = -1.0
         if isnan(cellInfo.infoUpper):
            cellInfo.infoUpper = -1.0
         if isnan(cellInfo.infoLower):
            cellInfo.infoLower = -1.0
      
   
         #c = cellInfo
         #print '%s,%s,%s,%f,%f\n' % (c.birdName, c.cellName, c.region, c.spikeRate, c.infoMean)

         self.cellsByRegion[region].append(cellInfo)

      f.close()


   def createCellList(self, dataCDDir, allowedCertainties = ['sure']):
      
      for baseDir in os.listdir(dataCDDir):
         
         fullPath = os.path.join(dataCDDir, baseDir)
         if os.path.isdir(fullPath):
            vals = string.split(baseDir, '_')

            if len(vals) > 1:
               certainty = vals[0]
               region = vals[1]

               if certainty in allowedCertainties:               
                  for dname in os.listdir(fullPath):

                     [birdName, cellName] = string.split(dname, '_', 1)
                     cellInfo = CellInfo()
                     cellInfo.cellName = cellName
                     cellInfo.birdName = birdName  
                     cellInfo.region = region               
                     
                     self.cellsByRegion[region].append(cellInfo)


   def getSortedCellsFromRegions(self, regionList):
      
      cellInfos = []
      for reg in regionList:
         clist = self.cellsByRegion[reg]
         cellInfos.extend(clist)

      sortedCellInfos = sorted(cellInfos, key=lambda c: c.infoMean, reverse=True)

      return sortedCellInfos



def printUsage():
   print 'Usage: python identify_morphs.py find <data CD directory>'
   print 'Usage: python identify_morphs.py copy <data CD directory> <celldata csv file> <destination directory>'
   print 'Usage: python identify_morphs.py mergeinfo <data directory> <celldata csv file>'

if __name__ == '__main__':

   if len(sys.argv) < 2:
      printUsage()
      exit()

   cmd = sys.argv[1]

   if cmd == 'find':

      dataCDDir = sys.argv[2]

      cellOrg = CellOrganizer()
      cellOrg.createCellList(dataCDDir)
    
      mldFile = open('celldata.mld.csv', 'w')
      mldList = cellOrg.getSortedCellsFromRegions(['mld'])
      for c in mldList:
         mldFile.write('%s,%s,%s\n' % (c.birdName, c.cellName, c.region))      
      mldFile.close()
    
      ovFile = open('celldata.ov.csv', 'w')
      ovList = cellOrg.getSortedCellsFromRegions(['OV'])
      for c in ovList:
         ovFile.write('%s,%s,%s\n' % (c.birdName, c.cellName, c.region))
      ovFile.close()
      
      lFile = open('celldata.l.csv', 'w')
      lList = cellOrg.getSortedCellsFromRegions(['L', 'L1', 'L2a', 'L2b', 'L3'])
      for c in lList:
         lFile.write('%s,%s,%s\n' % (c.birdName, c.cellName, c.region))      
      lFile.close()

      cmFile = open('celldata.cm.csv', 'w')
      cmList = cellOrg.getSortedCellsFromRegions(['CM'])
      for c in cmList:
         cmFile.write('%s,%s,%s\n' % (c.birdName, c.cellName, c.region))      
      cmFile.close()

   elif cmd == 'copy':

      dataCDDir = sys.argv[2]
      cellDataFile = sys.argv[3]
      destDir = sys.argv[4]

      f = open(cellDataFile, 'r')
      for l in f.readlines():         
         [birdName, cellName, region] = string.split(l.strip(), ',')
      
         dirName = 'sure_%s/%s_%s' % (region, birdName, cellName)
         fullPath = os.path.join(dataCDDir, dirName)

         sysCmd = 'cp -R %s %s' % (fullPath, destDir)
         print sysCmd
         os.system(sysCmd)

      f.close()

   elif cmd == 'mergeinfo':

      dataDir = sys.argv[2]
      cellDataFile = sys.argv[3]
      regListCsv = sys.argv[4]

      regList = string.split(regListCsv, ',')

      cellOrg = CellOrganizer()
      cellOrg.mergeInfoData(dataDir, cellDataFile)

      fileParts = string.split(cellDataFile, '.')
      fileName = '%s+info.%s.csv' % (fileParts[0], fileParts[1])

      outFile = open(fileName, 'w')
      cellList = cellOrg.getSortedCellsFromRegions(regList)
      for c in cellList:
         outFile.write('%s,%s,%s,%f,%f,%f,%f\n' % (c.birdName, c.cellName, c.region, c.spikeRate, c.infoLower, c.infoMean, c.infoUpper))      
      outFile.close()    

   else:
      print 'Unknown command %s, options are: find copy mergeinfo' % (cmd)
      exit()
    

