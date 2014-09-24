import string
import os


DATA_ROOT_DIR = '/auto/fdata/mschachter/data'

if __name__ == '__main__':

   cellNames = []
   for cname in os.listdir(DATA_ROOT_DIR):
      cellDir = os.path.join(DATA_ROOT_DIR, cname)
      if os.path.isdir(cellDir):
         if cname != 'all_stims':
            cellNames.append(cname)

   cellStimDirs = {}
   for cname in cellNames:
      cellDir = os.path.join(DATA_ROOT_DIR, cname)
      stimTypes = []
      for sdir in os.listdir(cellDir):
         stimDir = os.path.join(cellDir, sdir)
         if os.path.isdir(stimDir):
            stimTypes.append(sdir)
      cellStimDirs[cname] = stimTypes

   outputFile = 'data/stimtypes.csv'
   fout = open(outputFile, 'w')
   for cname,stypes in cellStimDirs.iteritems():
      fout.write('%s,%s\n' % (cname, ','.join(stypes)))
   fout.close()



