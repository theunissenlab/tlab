import os
import string

class MatlabTemplate:

   template = None

   def __init__(self):
      pass

   def template_from_file(self, fname):
      if not os.path.isfile(fname):
         print 'Template file does not exist: %s' % fname
         return
   
      f = open(fname, 'r')
      tLines = f.readlines()
      f.close()

      self.template = ';'.join(l.strip() for l in tLines)


   def fill_template(self, params):
   
      scriptStr = self.template      
      for name,val in params.iteritems():
         nStr = '$%s' % name
         scriptStr = string.replace(scriptStr, nStr, str(val))

      return scriptStr


if __name__ == '__main__':

   mt = MatlabTemplate()
   mt.template_from_file('scripts/groupinfo.template')
   
   params = {'CLEAN_OUTPUT_FOLDER':1, 'CELL_NAME':'pupu1203_4_A', 'FDATA_DIR':'/auto/fdata/mschachter/data'}

   mstr = mt.fill_template(params)
   print mstr
