import glob
import os
import re
import gc
import h5py
import numpy as np
import time


class SpikeSorter(object):
    def __init__(self, h5file):
        pass

    def sort(self, output_file, output_group):
        pass


class PreExistingRFSpikeSorter(SpikeSorter):
    """
        Takes spikes that have already been sorted by the random forest
        implementation within individual spike files and applies it.
    """

    def __init__(self, h5file, sorted_file_directory):
        SpikeSorter.__init__(self, h5file)

        self.h5file = h5file
        self.sorted_file_directory = sorted_file_directory
        self.sort_type = dict()

        #build a table of files keyed by site/electrode/sort
        self.file_map = dict()
        for ssfile in glob.glob(os.path.join(self.sorted_file_directory, '*ss*.h5')):
            hf = h5py.File(ssfile, 'r')

            site = hf.attrs['site']
            ldepth = int(hf.attrs['ldepth'])
            rdepth = int(hf.attrs['rdepth'])

            electrode = int(hf.attrs['electrode'])
            sort = int(hf.attrs['sortid'])
            sort_type = hf.attrs['sortType']

            #parse for TDT sort code
            tdt_sort = 0
            (ssroot,ssname) = os.path.split(ssfile)
            lm = re.search('_s[0-9]_', ssname)
            if lm is not None:
                tdt_sort = int(lm.group(0)[2])

            key = (site, ldepth, rdepth, electrode, tdt_sort)
            self.sort_type[key] = sort_type
            if key not in self.file_map:
                self.file_map[key] = list()
            self.file_map[key].append(ssfile)
            hf.close()

    def sort(self, output_hf, output_group):

        #iterate through electrode combinations in continuous file and separate spikes
        hf = h5py.File(self.h5file, 'r')
        for site_key,site_group in hf['sites'].iteritems():

            unit_count_by_electrode = dict()  #unique identifier for each sorted unit on an electrode

            site = site_group.attrs['site']
            ldepth = site_group.attrs['ldepth']
            rdepth = site_group.attrs['rdepth']
            protocol = site_group.attrs['protocol']
            start_date = site_group.attrs['start_date']
            end_date = site_group.attrs['end_date']

            output_site_name = '%s_%s_L%dR%d' % (site, protocol, ldepth, rdepth)
            print '[sort] Working on site: %s' % output_site_name
            output_site_group = output_group.create_group(output_site_name)
            output_site_group.attrs['site'] = site
            output_site_group.attrs['ldepth'] = ldepth
            output_site_group.attrs['rdepth'] = rdepth
            output_site_group.attrs['protocol'] = protocol
            output_site_group.attrs['start_date'] = start_date
            output_site_group.attrs['end_date'] = end_date
            
            for ekey,egroup in site_group['spikes'].iteritems():

                #check for spike-sorted files for this given electrode/sort combo
                electrode = int(egroup.attrs['electrode'])
                sort = int(egroup.attrs['sort'])
                key = (site, ldepth, rdepth, electrode, sort)
                if key not in self.file_map:
                    print 'No ss sort files found for %s_L%dR%d_e%d_s%d' % (site, ldepth, rdepth, electrode, sort)
                    continue
                sort_type = self.sort_type[key]

                print '\tElectrode/Sort: %d/%d' % (electrode, sort)

                if electrode not in unit_count_by_electrode:
                    unit_count_by_electrode[electrode] = 1

                #count the units by electrode, ignore sort code
                unit_number = unit_count_by_electrode[electrode]

                #get the continuous spikes and times for this electrode
                spike_ids = np.array(egroup['spike_ids'], dtype='int')
                spike_times = np.array(egroup['spike_times'], dtype='float')

                #create a map of spike ids to spike times
                spikeid2time = dict()
                for k,spike_time in enumerate(spike_times):
                    spike_id = spike_ids[k]
                    spikeid2time[spike_id] = spike_time
                del spike_times

                #find each sorted unit associated with the electrode/sort code pair, merge the spikes across stimuli
                #and trials, and then write the merged data to the output file
                for sort_file in self.file_map[key]:

                    stime = time.time()
                    print '\t\tSort File: %s' % sort_file
                    #construct list to contain all the spikes for this unit
                    unit_spike_ids = list()

                    sf = h5py.File(sort_file, 'r')
                    pgroup = sf[protocol]
                    #iterate through each stimulus in the protocol
                    for stim_key,stim_group in pgroup.iteritems():
                        #iterate through each trial in the stimulus
                        for trial_key,trial_group in stim_group.iteritems():
                            trial_spike_ids = np.array(trial_group['spike_ids'], dtype='int')
                            if trial_spike_ids.shape[0] == 0:
                                continue

                            #update list of unsorted spikes
                            spike_ids = np.setdiff1d(spike_ids, trial_spike_ids)
                            #deal with single spike trials
                            if trial_spike_ids.shape[0] == 1:
                                if trial_spike_ids[0] == -999:
                                    continue
                                trial_spike_ids = [int(trial_spike_ids[0, 0])]
                            else:
                                trial_spike_ids = trial_spike_ids.squeeze()
                            #append to list of spike ids for this unit
                            unit_spike_ids.extend(trial_spike_ids)

                    #write an output group to the file that contains the spike ids and times for the sorted unit
                    ugrp = self.write_unit_spikes(electrode, ldepth, rdepth, output_site_group, site,sort, sort_type, spikeid2time, unit_number, unit_spike_ids, sort_file)
                    output_hf.flush()
                    ugrp = None
                    del unit_spike_ids
                    unit_number += 1

                    sf.close()
                    etime = time.time() - stime
                    print '\t\tTime to process sort file: %f seconds' % etime

                #write the spikes that are left over to their own unit
                if len(spike_ids) > 0:
                    ugrp = self.write_unit_spikes(electrode, ldepth, rdepth, output_site_group, site, -1, 'unsorted', spikeid2time, unit_number, spike_ids, 'None')
                    output_hf.flush()
                    ugrp = None
                    del spike_ids
                    unit_number += 1
                del spikeid2time
                gc.collect()

                unit_count_by_electrode[electrode] = unit_number
        hf.close()

    def write_unit_spikes(self, electrode, ldepth, rdepth, output_site_group, site,sort, sort_type, spikeid2time, unit_number, unit_spike_ids, ss_file_name):

        if len(unit_spike_ids) == 0:
            return None

        unit_group_name = '%s_L%dR%d_e%d_%d' % (site,ldepth,rdepth,electrode,unit_number)
        if unit_group_name in output_site_group.keys():
            raise Exception('Error, already existing group name for unit! %s' % unit_group_name)
        unit_group = output_site_group.create_group(unit_group_name)
        unit_group.attrs['sort_type'] = sort_type
        unit_group.attrs['sort_file_name'] = ss_file_name
        unit_group.attrs['sort'] = sort
        unit_group.attrs['electrode'] = electrode
        for sid in unit_spike_ids:
            if type(sid).__name__ not in ['int64', 'int']:
                print 'Strange spike id: ',sid
        unit_group['spike_ids'] = unit_spike_ids
        unit_group['spike_times'] = [spikeid2time[int(spike_id)] for spike_id in unit_spike_ids]
        return unit_group