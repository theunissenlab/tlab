import string

class CellInfo:

   def __init__(self):
      self.cellName = None
      self.stimTypes = []
      self.region = None
      pass

   def __repr__(self):
      stypes = ','.join(self.stimTypes)
      cstr = 'Cell: %s | Region: %s | (%s)' % (self.cellName, self.region, stypes)
      return cstr


def get_all_cellinfo():

   cellInfos = []
   
   cellDataPath = 'data/celldata.all.csv'   
   f = open(cellDataPath, 'r')
   lns = f.readlines()
   for ln in lns:
      darr = ln.strip().split(',')
      if len(darr) > 0:
         bname = darr[0]
         cname = darr[1]
         reg = darr[2]
         cellName = '%s_%s' % (bname, cname)
         cinfo = CellInfo()
         cinfo.cellName = cellName
         cinfo.region = reg         
         cellInfos.append(cinfo)
      
   stimTypesPath = 'data/stimtypes.csv'
   f = open(stimTypesPath, 'r')
   lns = f.readlines()
   f.close()
   for ln in lns:
      darr = ln.strip().split(',')
      cname = darr[0]
      stypes = darr[1:]
      clst = filter(lambda ci: ci.cellName == cname, cellInfos)
      if len(clst) > 0:
         clst[0].stimTypes = stypes
      else:
         print 'No such name found: %s' % cname

   return cellInfos


if __name__ == '__main__':
   
   cinfos = get_all_cellinfo()

   for ci in cinfos:
      print ci



   
   
   



