from os import path
from numpy import array as narray, fromfile
import array

__all__ = ['BinFile','DlmFile']

class BinFile(object):
    ''' Binary data files '''

    def __init__(self,filename,dtype):
        self.dtype = dtype
        self.name = filename

    def __len__(self):
        return int(path.getsize(self.name)/narray([],self.dtype).itemsize)
        
    def __getattr__(self,name):
        if name == 'data':
            return self.getdata()
        elif name == 'name':
            return self.filename
        else:
            raise AttributeError, "\'%s\' object has no attribute \'%s\'" % (self.__class__, name)
        
    def getdata(self):
        #data = array.array(self.typecode)
        nbytes = path.getsize(self.name)
        #nvals = nbytes/data.itemsize
        nvals = nbytes/self.itemsize
        try:
            fo = file(self.name,'rb')
            #data.fromfile(fo,nvals)
            #return narray(data)
            return fromfile(fo,dtype=self.dtype,count=nvals)
        finally:
            fo.close()
    
    def get_chunk(self,start,stop):
        #data = array.array(self.typecode)
        nvals = stop-start+1
        try:
            fo = open(self.name,'rb')
            fo.seek(self.itemsize*start)
            #data.fromfile(fo,nvals)
            #return narray(data)
            return fromfile(fo,dtype=self.dtype,count=nvals)
        finally:
            fo.close()
            
    def _itemsize(self):
        return narray([],dtype=self.dtype).itemsize
    
    itemsize = property(_itemsize)
    
class DlmFile(object):
    ''' Files with a fixed delimiter; returns as ASCII text'''
    
    def __init__(self,filename,delimiter):
        self.name = filename
        self.delimiter = delimiter
        
    def __getattr__(self,name):
        if name == 'data':
            return self.getdata()
        elif name == 'name':
            return self.filename
        else:
            raise AttributeError, "\'%s\' object has no attribute \'%s\'" % (self.__class__, name)

    def _rawdata(self):
        fo = file(self.name,'r')
        data = fo.read()
        fo.close
        return data

    def getdata(self):
        return self._rawdata().split(self.delimiter)
