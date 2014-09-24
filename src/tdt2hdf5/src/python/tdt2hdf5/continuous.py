import os
import struct
import gc
import h5py

import numpy as np

from convert import Tank, Block


class ContinuousBlock(Block):

    def __init__(self):
        Block.__init__(self)
        self.lfp_by_electrode = dict()
        self.lfp_sample_rate = 381.4697

    def initialize(self):
        self.read_epoch_file()
        self.read_spike_file()
        self.aggregate_spikes()
        self.read_lfp_files()

    def aggregate_spikes(self):
        electrodes = self.spikes_data[:, 1].astype('int')
        sorts = self.spikes_data[:, -1].astype('int')
        self.unique_electrodes = range(1, 33)
        self.unique_sorts = np.unique(sorts)

        self.aggregate_spike_data = dict()

        for electrode in self.unique_electrodes:
            for sort in self.unique_sorts:
                index = (electrodes == electrode) & (sorts == sort)
                if index.sum() > 0:
                    spike_times = self.spikes_data[index, 0]
                    spike_ids = self.spikes_data[index, 2].astype('int')
                    self.aggregate_spike_data[(electrode, sort)] = (spike_times, spike_ids)
                else:
                    #print 'No spikes found for electrode/sort combination %d/%d' % (electrode, sort)
                    pass
        if len(self.aggregate_spike_data) == 0:
            raise Exception('No spikes were found for block %s!' % self.name)
        print '# of multiunit recordings: %d' % len(self.aggregate_spike_data)
        del self.spikes_data

    def read_lfp_files(self):

        for electrode in self.unique_electrodes:
            file_name = os.path.join(self.tank_dir, '%s LFPs %d stream.f32' % (self.name, electrode))
            if not os.path.exists(file_name):
                print 'Cannot locate LFP file for electrode %d: %s' % (electrode, file_name)
                continue

            lfp_data = read_f32_file(file_name)
            self.lfp_by_electrode[electrode] = lfp_data

    def write_snippets(self, parent_group, snippet_length=18):
        file_name = os.path.join(self.tank_dir, '%s waves.f32' % self.name)
        snippets_raw = read_f32_file(file_name)
        slen = len(snippets_raw)
        num_snippets = slen / snippet_length
        if slen % snippet_length != 0:
            print 'Length of snippets file is not a multiple of %d for block %s! Last snippet will be truncated.' % (snippet_length, self.name)
        parent_group['snippets'] = snippets_raw.reshape(num_snippets, snippet_length)
        del snippets_raw

    def write(self, grp):

        grp.attrs['name'] = self.name
        grp.attrs['site'] = self.site
        grp.attrs['protocol'] = self.protocol
        grp.attrs['rdepth'] = self.rdepth
        grp.attrs['ldepth'] = self.ldepth
        grp.attrs['start_date'] = self.start_date
        grp.attrs['end_date'] = self.end_date

        #write LFP data
        lfp_grp = grp.create_group('lfp')
        lfp_grp.attrs['sample_rate'] = self.lfp_sample_rate
        for e,lfp in self.lfp_by_electrode.iteritems():
            lfp_grp['%d' % e] = self.lfp_by_electrode[e]

        #write spike data
        spike_grp = grp.create_group('spikes')
        for (electrode,sort),(spike_times,spike_ids) in self.aggregate_spike_data.iteritems():
            es_grp = spike_grp.create_group('%d,%d' % (electrode,sort))
            es_grp.attrs['electrode'] = electrode
            es_grp.attrs['sort'] = sort
            es_grp['spike_times'] = spike_times
            es_grp['spike_ids'] = spike_ids

        #write stim data
        grp['stim_data'] = self.epoch_data

        #write snippets data
        self.write_snippets(grp)

    def plot_lfp(self):

        """
        for k in range(8):
            subplot(8, 2, k+1)
            plot(b.lfp_by_electrode[k+1])
            subplot(8, 2, k+1+8)
            plot(b.lfp_by_electrode[k+1+8])
            plt.title('Electrode %d' %
        """


def read_f32_file(fname):

    if not os.path.exists(fname):
        raise Exception('.f32 file does not exist: %s' % fname)

    statinfo = os.stat(fname)
    fsize = statinfo.st_size

    if fsize % 4 != 0:
        print 'File %s has a size that is not a multiple of 4! Will truncate the last 4 bytes.'

    nfloats = fsize / 4  # truncates the last 4 bytes for an abnormal file size

    f = open(fname, 'rb')
    rawdata = f.read()
    data = np.array(struct.unpack('f'*nfloats, rawdata))
    f.close()

    return data


class ContinuousTank(Tank):

    def __init__(self):
        Tank.__init__(self)
        self.block_class = ContinuousBlock

    def write(self, hf):

        sgrp = hf.create_group('stims')
        for stimnum,stim in self.stims.iteritems():

            ssgrp = sgrp.create_group('%d' % stim.number)

            ssgrp.attrs['number'] = stim.number
            ssgrp.attrs['original_wavfile'] = str(stim.original_wavfile)
            ssgrp.attrs['tdt_wavfile'] = str(stim.tdt_wavfile)
            ssgrp.attrs['stim_class'] = str(stim.stim_class)
            ssgrp.attrs['type'] = str(stim.type)
            ssgrp.attrs['source'] = str(stim.source)
            ssgrp.attrs['sex'] = str(stim.source_sex)
            for pname,pval in stim.parameters.iteritems():
                ssgrp.attrs[pname] = pval

    def write_block(self, bgrp, b):
        bname = '%s_%s_L%dR%d' % (b.site, b.protocol, b.ldepth, b.rdepth)
        bbgrp = bgrp.create_group(bname)
        b.write(bbgrp)


def import_tank_continuous(tank_dir, stim_file, output_file, block_mapping_file=None):
    """ Create a list of blocks to be written out """

    #open blocks file
    block_file = os.path.join(tank_dir, 'blocks.txt')
    if not os.path.exists(block_file):
        print 'Cannot find the blocks file, quitting: %s' % block_file
        return

    tank = ContinuousTank()
    tank.read_stim_protocol_file(stim_file)
    tank.read_blocks_file(block_file, tank_dir, block_mapping_file=block_mapping_file)

    hf = h5py.File(output_file, 'w')

    bird_name = tank.name[:-1]
    bird_gender = tank.name[-1]
    hf.attrs['bird_name'] = bird_name
    hf.attrs['bird_gender'] = bird_gender

    bgrp = hf.create_group('sites')

    print 'Writing stims to file...'
    tank.write(hf)
    for k,block in enumerate(tank.blocks):
        print ''
        print 'Initializing block %s, %d to go...' % (block.name, len(tank.blocks)-k-1)
        try:
            block.initialize()
        except ValueError:
            print 'Block %s is bad!' % block.name
            block.is_bad = True
        if not block.is_bad:
            tank.write_block(bgrp, block)
            del block
            gc.collect()
    hf.close()
    tank.blocks = None
    del tank
    hf = None
    gc.collect()