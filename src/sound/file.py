from __future__ import with_statement, absolute_import
from wave import open as wopen, Wave_read
from scipy.io.wavfile import read as wread
from .data import PCM, si as data_si

__all__ = ['Wave','DCP']

class WaveRead(Wave_read):
    '''
    Extension of Python Wave_read class to include
    __enter__ and __exit__ methods
    '''
    
    def __init__(self,filename):
        Wave_read.__init__(self,filename)
        
    def __enter__(self):
        return self
        
    def __exit__(self,type,value,traceback):
        self.close()
        return type is None
        
class Sound(object):
    '''
    Base class for soundfiles
    '''
    
    def __init__(self,filename):
        self.name = filename

class Wave(Sound):
    '''
    Wave file class that persists and wraps file io methods
    '''
    
    _reader = WaveRead
    
    def __init__(self,filename):
        # Removed because I want these to be usable
        # without calling __init__; SQLAlchemy query retrieval
        # does not call __init__ for some reason.
        #self.reader = None
        Sound.__init__(self,filename)

    def __len__(self):
        self.open()
        return self.reader.getnframes()
    
    def __enter__(self):
        self.open()
        return self()
    
    def __exit__(self,type,value,traceback):
        self.close()
        return type is None
    
    def open(self):
        # Need this to work without having __init__ set null reader
        #if self.reader is None:
        if (not hasattr(self,'reader')) or self.reader is None:
            self.reader = self._reader(self.name)
        return self.reader
        
    def close(self):
        if self.reader is not None:
            self.reader.close()
            self.reader = None
    
    def get_samplerate(self):
        with self.open() as r:
            return r.getframerate()
        
    def get_data(self):
        return wread(self.name)[1]
    
    def get_len(self):
        with self.open() as r:
            return r.getnframes()
    
    def get_nchannels(self):
        with self.open() as r:
            return r.getnchannels()
    
    def get_samplesize(self):
        with self.open() as r:
            return r.getsampwidth()
        
    def get_sound(self):
        return PCM(self.get_data(),self.get_samplerate(),file=self)
            
class DCP(Sound):
    pass

def open(soundfile):
    soundfile.open()
    
def close(soundfile):
    soundfile.close()
    
def read(soundfile):
    return soundfile.get_sound()

def si(*files,**kwargs):
    return data_si(*(read(f) for f in files),**kwargs)

def si_from_name(*paths):
    return si(*(Wave(p) for p in paths))