import hashlib
import os
import csv
import re
import h5py
import wave

import numpy as np

from scipy.stats import t as tdist
import time


class Stim(object):
    def __init__(self):
        self.number = None
        self.original_wavfile = None
        self.tdt_wavfile = None
        self.stim_class = None
        self.type = None
        self.source = None
        self.source_sex = None
        self.parameters = None


class Cell(object):
    def __init__(self):
        self.responses = []
        self.electrode = None
        self.sort = None
        self.site = None


class Response(object):
    def __init__(self):
        self.stim_number = None
        self.num_trials = None
        self.trials = []
        self.zscore = None
        self.tstat = None
        self.pval = None

    def compute_zscore(self):
        #get background and peri rates
        bg_rates = np.array([t.bg_rate for t in self.trials])
        peri_rates = np.array([t.peri_rate for t in self.trials])
        bg_counts = np.array([t.bg_count for t in self.trials])
        peri_counts = np.array([t.peri_count for t in self.trials])
        rate_diff = peri_rates - bg_rates
        count_diff = peri_counts - bg_counts

        if peri_counts.sum() + bg_counts.sum() < len(self.trials):
            pval = 0.5 #kludge
            z = 0.0
            tstat = 0.0
        else:
            rate_diff_std = rate_diff.std(ddof=1)
            count_diff_std = count_diff.std(ddof=1)
            if rate_diff_std == 0.0:
                print 'Very strange that this happenend, rate_diff_std=%0.3f, count_diff_std=%0.3f, stim_num=%d' % (rate_diff_std, count_diff_std, self.stim_number)
                rate_diff_std = 1.0
            z = rate_diff.mean() / rate_diff_std
            tstat = z*np.sqrt(len(bg_rates))
            pval = (1.0 - tdist.cdf(np.abs(tstat), len(bg_rates)-1))*2 #two-tailed t-test pvalue

        self.zscore = z
        self.tstat = tstat
        self.pval = pval


class Trial(object):
    def __init__(self):
        self.spike_ids = []
        self.spike_times = []
        self.bg_rate = 0.0
        self.peri_rate = 0.0
        self.bg_count = 0
        self.peri_count = 0
        self.pre_time = 0.0
        self.post_time = 0.0
        self.tank_dir = None


class Block(object):
    def __init__(self):
        self.name = None
        self.site = None
        self.protocol = None
        self.rdepth = None
        self.ldepth = None
        self.tank = None
        self.start_date = None
        self.end_date = None
        self.tank_dir = None
        self.epoch_data = None
        self.spikes_data = None
        self.electrodes = []
        self.is_bad = False

    def initialize(self):
        t1 = time.time()
        try:
            self.read_epoch_file()
        except ValueError:
            self.is_bad = True
            return
        t2 = time.time() - t1
        print 'read_epoch_file: %d seconds' % int(t2)

        t1 = time.time()
        self.process_epoch_data()
        t2 = time.time() - t1
        print 'process_epoch_data: %d seconds' % int(t2)

        t1 = time.time()
        self.read_spike_file()
        t2 = time.time() - t1
        print 'read_spike_file: %d seconds' % int(t2)

        t1 = time.time()
        self.process_spike_data()
        t2 = time.time() - t1
        print 'process_spike_data: %d seconds' % int(t2)

        t1 = time.time()
        self.compute_epoch_rates()
        t2 = time.time() - t1
        print 'compute_epoch_rates2: %d seconds' % int(t2)

        t1 = time.time()
        self.create_electrodes()
        t2 = time.time() - t1
        print 'create_electrodes: %d seconds' % int(t2)

    def read_epoch_file(self):
        file_name = os.path.join(self.tank_dir, '%s epoc Stim.txt' % self.name)
        if not os.path.exists(file_name):
            new_file_name = os.path.join(self.tank_dir, '%s epoc Stm+.txt' % self.name)
            print 'Cannot find epoch file: %s, trying %s' % (file_name, new_file_name)
            file_name = new_file_name
            if not os.path.exists(file_name):
                raise ValueError('Cannot find epoch file: %s, giving up you fucking fuck!!' % file_name)

        epoch_data = []
        f = open(file_name, 'r')
        creader = csv.reader(f, delimiter='\t')
        for row in creader:
            stim_num = float(row[0])
            start_time = float(row[1]) #in seconds
            end_time = float(row[2]) #in seconds
            epoch_data.append( [stim_num, start_time, end_time] )
        self.epoch_data = np.array(epoch_data)

        f.close()

    def process_epoch_data(self):

        stim_nums = self.epoch_data[:, 0].astype('int')
        min_stim_num = stim_nums.min()
        max_stim_num = stim_nums.max()

        stim_num_of_trials = []

        for stim_num in range(min_stim_num, max_stim_num+1):
            #count the # of times the stim appears
            num_trials = (stim_nums == stim_num).sum()
            if num_trials > 0:
                #only add good stimuli
                stim_num_of_trials.append( [stim_num, num_trials] )

        self.stim_num_of_trials = np.array(stim_num_of_trials)

        #find the trial number for each epoch
        new_epoch_data = []
        trial_num_map = {} #contains the next index for a trial, stim_num is key
        for (stim_num,start_time,end_time) in self.epoch_data:
            if stim_num not in trial_num_map:
                trial_num_map[stim_num] = 0
            trial_num = trial_num_map[stim_num]
            new_epoch_data.append( [stim_num, start_time, end_time, trial_num] )
            trial_num_map[stim_num] += 1
        #replace epoch_data with new one
        self.epoch_data = np.array(new_epoch_data)

        estart = self.epoch_data[:, 1]
        eend = self.epoch_data[:, 2].tolist()
        eend.insert(0, 0.0)
        inter_epoch_times = estart - np.array(eend[:-1])
        print 'Block %s: inter epoch time=%0.2f +/- %0.2f' % (self.name, inter_epoch_times.mean(), inter_epoch_times.std(ddof=1))
        self.inter_epoch_times = inter_epoch_times
        self.inter_epoch_mean = inter_epoch_times.mean()
        self.inter_epoch_std = inter_epoch_times.std(ddof=1)

    def compute_epoch_rates(self):

        bg_counts = {} #key = (epoch_index, electrode) value = spike count
        peri_counts = {} #key = (epoch_index, electrode) value = spike count

        #record start of background, and start and end time for each epoch for easy lookup
        epoch_prepost_time = [] #the amount of time considered "background" and post-stimulus
        epoch_times = []
        for epoch_index in range(self.epoch_data.shape[0]):
            epoch_start = self.epoch_data[epoch_index, 1]
            epoch_end = self.epoch_data[epoch_index, 2]

            pre_epoch_interval = self.inter_epoch_times[epoch_index]
            rel_background_start = pre_epoch_interval / 2.0  #take half this time
            abs_background_start = epoch_start - rel_background_start
            epoch_times.append( (epoch_start, epoch_end, abs_background_start) )

            rel_posttime = 0.0
            if epoch_index < self.epoch_data.shape[0]-1:
                next_epoch_index = epoch_index + 1
                next_epoch_interval = self.inter_epoch_times[next_epoch_index]
                rel_posttime = next_epoch_interval / 2.0

            epoch_prepost_time.append( (rel_background_start, rel_posttime))

        for spike_time, electrode, spike_id, sort, stim_num, trial_num, rel_time, epoch_index in self.spikes_data:

            epoch_index = int(epoch_index)
            #find out if this spike belongs in the background of an epoch
            (epoch_start, epoch_end, epoch_bg_start) = epoch_times[epoch_index]
            if spike_time < epoch_start and spike_time >= epoch_bg_start:
                bg_key = (epoch_index, electrode)
                if bg_key not in bg_counts:
                    bg_counts[bg_key] = 0
                bg_counts[bg_key] += 1
            elif spike_time >= epoch_start and spike_time <= epoch_end:
                #add spike if it happens during epoch
                peri_key = (epoch_index, electrode)
                if peri_key not in peri_counts:
                    peri_counts[peri_key] = 0
                peri_counts[peri_key] += 1


        espike_counts = {}
        espike_rates = {}
        #initialize everything to zero
        for epoch_index in range(self.epoch_data.shape[0]):
            for enum in self.unique_electrodes.keys():
                key = (epoch_index, enum)
                espike_counts[key] = [0.0, 0.0]
                espike_rates[key] = [0.0, 0.0]

        #append background rates
        for (epoch_index,electrode),count in bg_counts.iteritems():
            key = (epoch_index,electrode)
            bg_dur = epoch_times[epoch_index][0] - epoch_times[epoch_index][2]
            espike_counts[key][0] = count
            espike_rates[key][0] = float(count) / bg_dur

        #append peri-stimulus rates
        for (epoch_index,electrode),count in peri_counts.iteritems():
            key = (epoch_index,electrode)
            peri_dur = epoch_times[epoch_index][1] - epoch_times[epoch_index][0]
            espike_counts[key][1] = count
            espike_rates[key][1] = float(count) / peri_dur

        self.epoch_prepost_times = np.array(epoch_prepost_time)
        self.epoch_rates = espike_rates
        self.epoch_spike_counts = espike_counts

    def read_spike_file(self):

        file_name = os.path.join(self.tank_dir, '%s spikes.txt' % self.name)
        if not os.path.exists(file_name):
            print 'Cannot find spikes file: %s' % file_name
            return

        spikes_data = []
        f = open(file_name, 'r')
        creader = csv.reader(f, delimiter='\t')
        for row in creader:
            spike_time = float(row[0])
            electrode = float(row[1])
            spike_id = float(row[2])
            sort = float(row[3])

            spikes_data.append( [spike_time, electrode, spike_id, sort] )
        self.spikes_data = np.array(spikes_data)

        f.close()

    def process_spike_data(self):

        new_spikes_data = []
        epoch_onset_index = 1
        epoch_offset_index = 2
        self.unique_electrodes = {}

        for (spike_time, electrode, spike_id, sort) in self.spikes_data:
            spike_onset = spike_time - self.epoch_data[:, epoch_onset_index]
            spike_offset = spike_time - self.epoch_data[:, epoch_offset_index]
            if electrode not in self.unique_electrodes:
                self.unique_electrodes[electrode] = True

            #find the index of the last epoch onset that occurs before spike
            epoch_index = None
            epochs_before_spike = (spike_onset >= 0).nonzero()[0]
            if len(epochs_before_spike) == 0:
                #in this case, the spike time occurs before the first epoch
                epoch_index = 0
            else:
                last_epoch_before_spike = epochs_before_spike[-1]

            #find the index of the first epoch offset that occurs after a spike
            epochs_after_spike = (spike_offset <= 0).nonzero()[0]
            if len(epochs_after_spike) == 0:
                #in this case, spike time occurs after last epoch
                epoch_index = self.epoch_data.shape[0] - 1
            else:
                first_epoch_after_spike = epochs_after_spike[0]

            if epoch_index is None:
                if last_epoch_before_spike == first_epoch_after_spike:
                    #in this case, the spike falls within an epoch
                    epoch_index = last_epoch_before_spike
                elif first_epoch_after_spike > last_epoch_before_spike and (first_epoch_after_spike - last_epoch_before_spike == 1):
                    #in this case the spike time occurs in between two epochs. first we
                    #compute the distance between the spike time and the epoch preceding it (left)
                    #and the epoch following it (right)
                    right_epoch_diff = self.epoch_data[first_epoch_after_spike, epoch_onset_index] - spike_time
                    left_epoch_diff = spike_time - self.epoch_data[last_epoch_before_spike, epoch_offset_index]
                    #choose the epoch that is closest to the spike time
                    if right_epoch_diff < left_epoch_diff:
                        #choose right
                        epoch_index = first_epoch_after_spike
                    else:
                        #choose left
                        epoch_index = last_epoch_before_spike
                else:
                    raise Exception('Epoch data corrupted, potentially not sorted properly! first_epoch_after_spike=%d, last_epoch_before_spike=%d' % \
                                    (first_epoch_after_spike, last_epoch_before_spike))

            stim_num = self.epoch_data[epoch_index, 0]
            trial_num = self.epoch_data[epoch_index, 3]
            rel_time = spike_time - self.epoch_data[epoch_index, epoch_onset_index]

            new_spikes_data.append( [spike_time, electrode, spike_id, sort, stim_num, trial_num, rel_time, epoch_index] )

        self.spikes_data = np.array(new_spikes_data)

    def create_electrodes(self):

        electrodes = self.spikes_data[:, 1].squeeze()
        unique_electrodes = self.unique_electrodes.keys()
        print 'Unique electrodes: %s' % ','.join(['%d' % ev for ev in unique_electrodes])

        for electrode_num in unique_electrodes:

            #get submatrix that corresponds to only this electrode
            if len(unique_electrodes) == 1:
                espike_data = self.spikes_data
            else:
                espike_data = self.spikes_data[electrodes == electrode_num, :]
            sorts = espike_data[:, 3].squeeze()
            unique_sorts = np.unique(sorts)
            print '  Electrode %d, unique sorts: %s' % (electrode_num, ','.join(['%d' % sn for sn in unique_sorts]))

            for sort in unique_sorts:

                #get submatrix that corresponds to only this electrode/sort combination
                if len(unique_sorts) == 1:
                    esort_data = espike_data
                else:
                    esort_data = espike_data[sorts == sort, :]
                #print 'esort_data.shape=',esort_data.shape

                electrode = Electrode()
                electrode.id = int(electrode_num)
                electrode.sort = int(sort)

                if esort_data.shape[0] > 1:
                    stim_nums = esort_data[:, 4].squeeze()
                else:
                    stim_nums = np.array(esort_data[:, 4].squeeze()).reshape([1])
                #print 'stim_nums.shape=',stim_nums.shape
                unique_stim_nums = np.unique(self.epoch_data[:, 0])

                print '    Stim Nums: %s' % ','.join(['%d' % sn for sn in unique_stim_nums])

                for stim_num in unique_stim_nums:

                    num_stim_trials = (self.epoch_data[:, 0] == stim_num).sum()

                    #get submatrix for this electrode/sort/stimulus
                    if len(unique_stim_nums) == 1:
                        estim_data = esort_data
                    else:
                        estim_data = esort_data[stim_nums == stim_num, :]

                    if len(estim_data) == 0:
                        #this means there were no spikes for this stimulus number across trials
                        print '      No spikes for stim number %d' % stim_num
                        r = Response()
                        r.num_trials = num_stim_trials
                        r.stim_number = int(stim_num)
                        for tnum in range(r.num_trials):
                            t = Trial()
                            t.spike_ids = np.array([])
                            t.spike_times = np.array([])
                            r.trials.append(t)
                            r.compute_zscore() #now that all the trials are added, compute z-score for stimulus
                        electrode.responses.append(r)
                        continue

                    #print 'estim_data.shape=',estim_data.shape
                    if len(estim_data.shape) < 2:
                        estim_data = estim_data.reshape([1, len(estim_data)])
                    trials = estim_data[:, 5].squeeze()
                    #print 'trials.shape=',trials.shape
                    unique_trials = np.arange(0, num_stim_trials)

                    r = Response()
                    r.num_trials = len(unique_trials)
                    r.stim_number = int(stim_num)

                    print '      Stim number %d has %d trials' % (stim_num, num_stim_trials)

                    for trial in unique_trials:
                        t = Trial()

                        #find the epoch index that corresponds to this stimulus/trial
                        epoch_indices = (self.epoch_data[:, 0] == stim_num) & (self.epoch_data[:, 3] == trial)
                        if epoch_indices.sum() == 0:
                            #no such trial
                            print 'This shouldn\'t happen! No trial associated with electrode %d, sort %d, stim num %d, trial %d' % \
                                  (electrode_num, sort, stim_num, trial)
                            t.spike_ids = np.array([])
                            t.spike_times = np.array([])
                            r.trials.append(t)
                            continue

                        if epoch_indices.sum() > 1:
                            print 'Multiple epochs found for stimulus number %d and trial number %d, something is wrong!' % (stim_num, trial)
                        epoch_index = epoch_indices.nonzero()[0][0]

                        #find background rates that correspond to electrode/stimulus/trial
                        key = (epoch_index, electrode_num)
                        if key not in self.epoch_rates:
                            print 'Cannot find background/peri spike rates for epoch index=%d, electrode=%d, something is wrong!' % (epoch_index, electrode_num)
                        (bg_rate, peri_rate) = self.epoch_rates[ (epoch_index, electrode_num) ]
                        (bg_count, peri_count) = self.epoch_spike_counts[ (epoch_index, electrode_num) ]
                        r.stim_duration = self.epoch_data[epoch_index, 2] - self.epoch_data[epoch_index, 1]
                        (pre_bg_time, post_bg_time) = self.epoch_prepost_times[epoch_index]
                        t.pre_time = pre_bg_time
                        t.post_time = post_bg_time
                        t.bg_rate = bg_rate
                        t.peri_rate = peri_rate
                        t.bg_count = bg_count
                        t.peri_count = peri_count

                        #get submatrix for this electrode/sort/stimulus/trial
                        if len(trials.shape) > 0:
                            et_index = trials == trial
                        else:
                            if trials.squeeze() == trial:
                                et_index = np.array([True], dtype='bool')
                            else:
                                et_index = np.array([], dtype='bool')
                        if et_index.sum() > 0:
                            etrial_data = estim_data[et_index, :]
                        else:
                            etrial_data = np.array([])
                        #print 'et_index.shape=',et_index.shape
                        #print 'etrial_data.shape=',etrial_data.shape

                        for (xspike_time, xelectrode, xspike_id, xsort_num, xstim_num, xtrial_num, xrel_time, xepoch_index) in etrial_data:
                            t.spike_ids.append(xspike_id)
                            t.spike_times.append(xrel_time)

                        r.trials.append(t)

                    r.compute_zscore() #now that all the trials are added, compute z-score for stimulus
                    electrode.responses.append(r)

                self.electrodes.append(electrode)


class Tank(object):
    def __init__(self):
        self.directory = None
        self.blocks = []
        self.stims = {}
        self.name = None
        self.block_class = Block

    def read_blocks_file(self, block_file, tank_dir, block_mapping_file=None):

        (root_dir, tname) = os.path.split(tank_dir)
        self.name = tname

        block_name_maps = {}
        if block_mapping_file is not None:
            print 'Reading block mapping file from %s' % block_mapping_file
            f = open(block_mapping_file, 'r')
            creader = csv.reader(f, delimiter=',')
            for row in creader:
                if len(row) > 0:
                    orig_block_name = row[0]
                    formatted_block_name = row[1]
                    block_name_maps[orig_block_name] = formatted_block_name
            f.close()

        self.directory = tank_dir
        f = open(block_file, 'r')
        creader = csv.reader(f, delimiter='\t')
        for row in creader:
            block = self.block_class()

            formatted_block_name = row[1]
            if formatted_block_name in block_name_maps:
                print 'Mapping block name from %s to %s' % (formatted_block_name, block_name_maps[formatted_block_name])
                formatted_block_name = block_name_maps[formatted_block_name]

            site,protocol,rdepth,ldepth = parse_block_name(formatted_block_name)
            block.site = site
            block.protocol = protocol
            block.rdepth = rdepth
            block.ldepth = ldepth

            block.tank_dir = tank_dir
            block.tank = row[0]

            block.name = row[1]
            block.start_date = row[2]
            block.end_date = row[3]
            self.blocks.append(block)
        f.close()

    def read_stim_protocol_file(self, protocol_file):

        f = open(protocol_file, 'r')
        creader = csv.reader(f, delimiter='\t')
        for row in creader:
            stim_obj = Stim()
            stim_obj.number = int(row[0])
            stim_obj.original_wavfile = parse_na_to_none(row[1])
            stim_obj.tdt_wavfile = parse_na_to_none(row[2])
            stim_obj.stim_class = parse_na_to_none(row[3])
            stim_obj.type = parse_na_to_none(row[4])
            stim_obj.source = parse_na_to_none(row[5])
            stim_obj.source_sex = parse_na_to_none(row[6])
            stim_obj.parameters = parse_stim_parameters(row[7])

            self.stims[stim_obj.number] = stim_obj
        f.close()


def parse_block_name(bname):

    blist = bname.split('_')
    if len(blist) < 3:
        print 'Invalid block name, should be Site_Protocol_RxxxLxxx: %s' % bname
        return bname, None, None, None
    site = blist[0]
    protocol = blist[1]
    (rdepth, ldepth) = parse_depth(blist[2])

    return site,protocol,rdepth,ldepth


def parse_depth(dstr):

    rdepth = None
    rm = re.search('R[0-9]*', dstr)
    if rm is not None:
        rdepth = int(rm.group(0)[1:])

    ldepth = None
    lm = re.search('L[0-9]*', dstr)
    if lm is not None:
        ldepth = int(lm.group(0)[1:])

    return rdepth, ldepth


def parse_na_to_none(str):
    if str.lower() in ['na', 'unused']:
        return None
    return str


def parse_stim_parameters(pstr):
    pstr = parse_na_to_none(pstr)
    if pstr is None:
        return {}
    #split on spaces
    params = {}
    plist = pstr.split(' ')
    for nvp in plist:
        (name, value) = nvp.split('=')
        params[name] = value
    return params


class Electrode(object):
    def __init__(self):
        self.id = None
        self.sort = None
        self.responses = []


def import_tank(tank_dir, stim_file, block_mapping_file=None):
    """ Create a list of blocks to be written out """

    #open blocks file
    block_file = os.path.join(tank_dir, 'blocks.txt')
    if not os.path.exists(block_file):
        print 'Cannot find the blocks file, quitting: %s' % block_file
        return

    tank = Tank()
    tank.read_stim_protocol_file(stim_file)
    tank.read_blocks_file(block_file, tank_dir, block_mapping_file=block_mapping_file)
    #print 'WARNING: IN DEBUG MODE'
    #tank.blocks = tank.blocks[:1] #DEBUG

    for k,block in enumerate(tank.blocks):
        print 'Initializing block %s, %d to go...' % (block.name, len(tank.blocks)-k+1)
        try:
            block.initialize()
        except ValueError:
            print '!!!!!!!! Block %s is bad!' % block.name
            block.is_bad = True

    return tank


def write_tank(tank, output_dir, sort_codes='*', check_zscores=True, pval_thresh=0.05):
    """ Write a list of blocks to an HDF5 file """

    if sort_codes == '*':
        sorts = range(10)
    else:
        sorts = [int(x) for x in sort_codes.split(',')]

    good_blocks = [b for b in tank.blocks if not b.is_bad]

    #aggregate electrodes across blocks where sites and depths are same
    cell_map = {} #map a unique cell to a dict where the key is a protocol name, and the value is an electrode
    for block in good_blocks:
        for electrode in block.electrodes:
            if electrode.sort not in sorts:
                continue

            #a cell is uniquely identified by a block site, depth, electrode id, and sort code
            ekey = block.site,block.rdepth,block.ldepth,electrode.id,electrode.sort
            if ekey not in cell_map:
                cell_map[ekey] = {}
            cell_map[ekey][block.protocol] = electrode

    bad_cells = {}
    if check_zscores:
        #do Bonferroni correction on pvalue and check to see if responses are stimulus conditioned
        for (site,rdepth,ldepth,electrode_id,sort_code),cell_protocols in cell_map.iteritems():
            ekey = (site,rdepth,ldepth,electrode_id,sort_code)
            all_pvals = []
            for protocol_name in cell_protocols:
                electrode = cell_map[ekey][protocol_name]
                for resp in electrode.responses:
                    all_pvals.append(resp.pval)
            all_pvals = np.array(all_pvals)
            bf_pval_thresh = pval_thresh / len(all_pvals)
            #if there are any pvalues below the corrected pvalue, we're good
            num_below = (all_pvals < bf_pval_thresh).sum()
            if num_below == 0:
                print 'Found non-stimulus conditioned cell: site=%s_R%sL%s electrode=%d' % (site, str(rdepth), str(ldepth), electrode_id)
                bad_cells[ekey] = True

    #write the cells to a file
    for (site,rdepth,ldepth,electrode_id,sort_code),cell_protocols in cell_map.iteritems():

        ekey = (site,rdepth,ldepth,electrode_id,sort_code)
        if ekey in bad_cells:
            continue

        rdepth_str = ''
        if rdepth is not None:
            rdepth_str = 'R%d' % rdepth
        ldepth_str = ''
        if ldepth is not None:
            ldepth_str = 'L%d' % ldepth

        cell_name = '%s_%s%s_e%d_s%d' % (site, ldepth_str, rdepth_str, electrode_id,sort_code)

        #initialize the HDF5 file
        file_name = os.path.join(output_dir, '%s.h5' % cell_name)
        f = h5py.File(file_name, 'w')
        f.attrs['tank_name'] = tank.name
        f.attrs['source_directory'] = tank.directory
        f.attrs['site'] = site
        if rdepth is not None:
            f.attrs['rdepth'] = rdepth
        if ldepth is not None:
            f.attrs['ldepth'] = ldepth
        f.attrs['electrode'] = electrode_id
        f.attrs['sortid'] = sort_code
        f.attrs['sortType'] = 'tdt'

        for protocol,electrode in cell_protocols.iteritems():
            #create a group for the protocol
            protocol_grp = f.create_group(protocol)
            for k,resp in enumerate(electrode.responses):

                stim = tank.stims[resp.stim_number]
                (md5, stim_dur) = read_stim_info(stim.tdt_wavfile)

                grp_name = '%d' % stim.number
                try:
                    resp_grp = protocol_grp.create_group(grp_name)
                except ValueError:
                    print 'ValueError creating group %s, cell=%s, protocol=%s, electrode=%d' % \
                            (grp_name, cell_name, protocol, electrode.id)
                    continue
                resp_grp.attrs['stim_md5'] = md5
                resp_grp.attrs['stim_duration'] = stim_dur
                resp_grp.attrs['original_wavfile'] = str(stim.original_wavfile)
                resp_grp.attrs['tdt_wavfile'] = str(stim.tdt_wavfile)
                resp_grp.attrs['stim_source'] = str(stim.source)
                resp_grp.attrs['stim_source_sex'] = str(stim.source_sex)
                resp_grp.attrs['stim_class'] = str(stim.stim_class)
                resp_grp.attrs['stim_type'] = str(stim.type)
                resp_grp.attrs['zscore'] = resp.zscore
                resp_grp.attrs['tstat'] = resp.tstat
                resp_grp.attrs['pvalue'] = resp.pval
                for pname,pval in stim.parameters.iteritems():
                    resp_grp.attrs[pname] = pval

                for m,trial in enumerate(resp.trials):
                    tgrp_name = '%d' % (m+1)
                    trial_grp = resp_grp.create_group(tgrp_name)
                    trial_grp.attrs['bg_rate'] = trial.bg_rate
                    trial_grp.attrs['peri_rate'] = trial.peri_rate
                    trial_grp.attrs['pre_time'] = trial.pre_time
                    trial_grp.attrs['post_time'] = trial.post_time
                    st = trial.spike_times
                    if len(st) == 0:
                        st = -999
                    si = trial.spike_ids
                    if len(si) == 0:
                        si = -999
                    trial_grp['spike_times'] = st
                    trial_grp['spike_ids'] = si

        f.close()


def read_stim_info(wav_file):
    wf = wave.open(wav_file, 'r')
    samp_rate = wf.getframerate()
    num_frames = wf.getnframes()
    duration = float(num_frames) / samp_rate
    wf.close()
    md5 = md5_for_file(wav_file)

    return md5,duration


def md5_for_file(file_name, block_size=2**20):
    f = open(file_name, 'r')
    md5 = hashlib.md5()
    while True:
        data = f.read(block_size)
        if not data:
            break
        md5.update(data)
    f.close()
    return md5.hexdigest()
