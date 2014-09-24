from __future__ import with_statement
from numpy import array
from os import getenv, path

class Waveforms(object):
    
    def __init__(self,data,nsamples):
        if data.ndim != 2:
            raise ValueError, "Spike waveform data input must be 2-D"
        
        if data.shape[0] == nsamples:
            self.data = data.transpose()
        elif data.shape[1] == nsamples:
            self.data = data
        else:
            raise ValueError, "Spike waveform data must have one\
                                dimension of length %d" % nsamples
        
    def _pca(self):
        pass
    
    def sort_with_klustakwik(self):
        pass
    
    def to_kk_file(self):
        filename = path.join(getenv('TMPDIR'),'spikes.fet.1')
        with KKFile(filename,'w') as kkfile:
            kkfile.write(self.data)
            
        return kkfile.name
    
class KKFile(file):
    
    def write(self,array_data):
        if array_data.ndim != 2:
            raise ValueError, "Input data must be 2-D"
        
        # Write data size as first line in file
        file.write(self,'%d\n' % array_data.shape[1])
        
        for row in array_data:
            row.tofile(self,' ','%f')
            
def test(channel):
    from tlab.lib.file.data import BinFile
    from numpy import fromfile, mean, matrix, uint8
    from numpy.linalg import eig
    datadir = '/Users/channing/Documents/GrayGray1615/'
    datablockname = 'loc 7 depth 2500'
    datafilename = path.join(datadir,datablockname + ' waves.f32')
    
    # Get nsamples
    spikefilename = path.join(datadir,datablockname + ' spikes.txt')
    spikedata = BinFile(datafilename,'float32').getdata()
    with open(spikefilename) as spikesfile:
        nspikes = len(spikesfile.readlines())
    timestamps = fromfile(spikefilename,sep=' ',)
    timestamps = timestamps.reshape(nspikes,len(timestamps)/nspikes)
    channels = uint8(timestamps[:,1])
    sortcodes = uint8(timestamps[:,3])
    timestamps = timestamps[:,0]
    nsamples = len(spikedata)/nspikes
    
    spikedata = Waveforms(spikedata.reshape(nspikes,nsamples),nsamples)
    clusterfilename = path.join(datadir,'spikes.clu.1')
    clusters = fromfile(clusterfilename,dtype='uint8',sep=' ')
    
    sd = Waveforms(spikedata.data[channels==10,:],nsamples)
    
    datameans = mean(sd.data,0)
    meansub = sd.data - datameans
    acorr = matrix(meansub).transpose() * matrix(meansub)
    w,v = eig(acorr)
    scores = matrix(meansub) * v
    return spikedata,sd,clusters,w,v,scores