import sys
import os
import slurmbot
import matlabtemplate


def merge_info(dataDir, cellDataFile):
   
   f = open(cellDataFile)
   flines = f.readlines()
   f.close()

   cellInfos = []
   for l in flines:

      cdata = l.split(',')
      cdir = '%s_%s' % (cdata[0], cdata[1])

      cinfo = {}
      cinfo['bird'] = cdata[0].strip()
      cinfo['cell'] = cdata[1].strip()
      cinfo['region'] = cdata[2].strip()
      cinfo['dir'] = os.path.join(dataDir, cdir)
      
      cellInfos.append(cinfo)

   for cinfo in cellInfos:

      cdir = cinfo['dir']
      stimdir = os.path.join(cdir, 'conspecific')
      if os.path.exists(stimdir):
         
         outdir = os.path.join(stimdir, 'output')
         ifile = os.path.join(outdir, 'info_vals.txt')
         ifd = open(ifile)
         iflines = ifd.readlines()
         ifd.close()
         iinfo = {}
         for il in iflines:
            ivals = il.split('=')
            iinfo[ivals[0]] = ivals[1].strip()
         
         cinfo['rate'] = float(iinfo['meanSpikeRate'])
         cinfo['info'] = float(iinfo['meanInfoVal'])

      else:
         print 'Dir does not exist: %s' % stimdir
         cinfo['rate'] = -1.0
         cinfo['info'] = -1.0

   ofile = 'celldata+info.all.csv'
   ofd = open(ofile, 'w')
   for cinfo in cellInfos:
      ofd.write('%s,%s,%s,%0.3f,%0.3f\n' % (cinfo['bird'], cinfo['cell'], cinfo['region'], cinfo['rate'], cinfo['info']))
   ofd.close()

   


if __name__ == '__main__':

   if len(sys.argv) < 2:
      print 'Usage:'
      print '  python run_groupinfo.py run'
      print '  python run_groupinfo.py merge'
      exit()

   cmd = sys.argv[1]

   if cmd == 'merge':
      dataDir = '/auto/fdata/mschachter/data'
      cellDataFile = 'data/celldata.all.csv'
      merge_info(dataDir, cellDataFile)


   if cmd == 'run':

      f = open('data/celldata.all.csv')
      flines = f.readlines()
      f.close()
      
      cellDirs = []
      for l in flines:
         cdata = l.split(',')
         cdir = '%s_%s' % (cdata[0], cdata[1])
         cellDirs.append(cdir)

      sb = slurmbot.SlurmBot()
      sb.maxJobs = 10

      rootDir = '/auto/fdata/mschachter'
      outputDir = os.path.join(rootDir, 'slurm')
      
      mt = matlabtemplate.MatlabTemplate()
      mt.template_from_file('scripts/groupinfo.template')

      sparams = {'out': os.path.join(outputDir, 'slurm_%j_out.txt'),
                 'err': os.path.join(outputDir, 'slurm_%j_err.txt'),
                 'partition':'all'}

      for cellName in cellDirs:
         params = {'CLEAN_OUTPUT_FOLDER':0, 'CELL_NAME':cellName, 'ROOT_DIR':rootDir}
         mstr = mt.fill_template(params);
         
         logFile = os.path.join(outputDir, ('groupinfo_%s.txt' % cellName))

         cmds = ['matlabbg2', mstr, logFile]
         sb.add(cmds, sparams)
      
      sb.run_and_wait()
      exit()
   

   
      
