import string
import os
import sys
import glob
import shutil

if __name__ == '__main__':

   if len(sys.argv) < 3:
      print 'Usage :python clean_datadir.py all <data root dir>'
      print 'Usage :python clean_datadir.py output <data root dir> globPattern'
      exit()

   cmd = sys.argv[1]

   if cmd == 'output':

      DATA_ROOT_DIR = sys.argv[2]
      print 'Data Root Dir: %s' % DATA_ROOT_DIR

      globPattern = sys.argv[3]

      if not os.path.exists(DATA_ROOT_DIR):
         print 'Data root directory does not exist!'
         exit()
      
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

      for cname,stimTypes in cellStimDirs.iteritems():
         cellDir = os.path.join(DATA_ROOT_DIR, cname)
         for stype in stimTypes:
            dataRootDir = os.path.join(cellDir, stype)
            outputDir = os.path.join(dataRootDir, 'output')

            if os.path.exists(outputDir):
               outputGlob =  os.path.join(outputDir, globPattern)
            
               print 'Removing %s' % outputGlob
               for ofile in glob.glob(outputGlob):
                  if os.path.isdir(ofile):                     
                     pass
                  else:
                     os.remove(ofile)



   elif cmd == 'all':

      DATA_ROOT_DIR = sys.argv[2]
      print 'Data Root Dir: %s' % DATA_ROOT_DIR

      if not os.path.exists(DATA_ROOT_DIR):
         print 'Data root directory does not exist!'
         exit()

      print 'This will remove ALL output and preprocessing files from the data directory, are you sure you want to do this?'
      doThis1 = raw_input('[Y/n] ')
      if doThis1 != 'Y':
         exit()
      print 'Are you REALLY REALLY sure?!?'
      doThis2 = raw_input('[Y/n] ')
      if doThis2 != 'Y': 
         exit()

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

      for cname,stimTypes in cellStimDirs.iteritems():
         cellDir = os.path.join(DATA_ROOT_DIR, cname)
         for stype in stimTypes:
            dataRootDir = os.path.join(cellDir, stype)
            preprocDir = os.path.join(dataRootDir, 'preproc')
            outputDir = os.path.join(dataRootDir, 'output')
            if os.path.exists(preprocDir):
               preprocGlob = os.path.join(preprocDir, '*')
               print 'Removing %s' % preprocGlob
               for pfile in glob.glob(preprocGlob):
                  if os.path.isdir(pfile):
                     shutil.rmtree(pfile)
                  else:
                     os.remove(pfile)
            if os.path.exists(outputDir):
               outputGlob =  os.path.join(outputDir, '*')
            
               print 'Removing %s' % outputGlob
               for ofile in glob.glob(outputGlob):
                  if os.path.isdir(ofile):
                     shutil.rmtree(ofile)
                  else:
                     os.remove(ofile)
