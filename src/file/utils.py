from os import makedirs, path, rmdir
from shutil import rmtree
from tempfile import mkdtemp

__all__ = ['fiatdir']

def fiatdir(dirname,mode = 0755):
    if not path.exists(dirname):
        makedirs(dirname,mode)
        
class TempDir(object):
    
    def __enter__(self,path=None):
        self.path = path
    
    def __enter__(self):
        self.path = mkdtemp()
        return self.path
    
    def __exit__(self,*args):
        self.delete()
        
    def delete(self):
        rmtree(self.path)
        self.path = None
