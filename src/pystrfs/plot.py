import os
import sys
import csv
import operator

import numpy as np

import h5py

import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg

from browser.model import meta
from browser.model.analysis.model import Model
from browser.model.analysis.strflab.models import StrflabLinearModel, StrflabSeparableNLModel, \
                                                  StrflabGlmModel, StrflabLogExpGlmModel
                                                  
from pystrfs import map_dcp_unit_to_fs, PYSTRFS_DIR, get_unit_by_old_id
from pystrfs.analysis import SingleModelData, AggregateModel, OutputNL


class UnitResponse:
    def __init__(self, unit, stimClass):        
        self.stimClass = stimClass
        self.unitFile = map_dcp_unit_to_fs(unit)        
        self.create_response()
        
    def create_response(self):
        
        f = h5py.File(self.unitFile, 'r')
        
        scgrp = f['class_responses'][self.stimClass]
        
        group_index = np.array([], dtype='int32')
        nstims = int(scgrp.attrs['num_stims'])
        all_psths = np.array([], dtype='float')
        num_trials = []
        total_len = 0
        for k in range(nstims):
            
            respgrp = scgrp[str(k+1)]
            stimdur = int(respgrp.attrs['trial_duration']*1000)
            total_len += stimdur
             
            ntrials = respgrp.attrs['num_trials']
            num_trials.append(ntrials)            
            psth_bins = np.zeros([1, stimdur+1], dtype='int32')
            
            gindx=[k+1]*(stimdur+1)
            group_index = np.append(group_index, gindx)
                        
            for j in range(ntrials):
                stimes = np.array(respgrp[str(j+1)]).squeeze()
                if len(stimes.shape) > 0:
                    stimes *= 1000
                    stimes = stimes.astype('int32')
                    sindx = stimes[(stimes >=0) & (stimes <= stimdur)]
                    psth_bins[0, sindx] += 1
            
            psth = psth_bins.astype('float') / float(ntrials)                               
            all_psths = np.append(all_psths, psth)
            
        f.close()
        
        self.groupIndex = group_index
        self.numTrials = np.array(num_trials)
        self.psth = all_psths
        

class ModelPlot(AggregateModel):

    def __init__(self, unit, output_file_template, nfolds = 5):
        output_files = []
        for k in range(nfolds):
            output_files.append(output_file_template % k)
        AggregateModel.__init__(self, output_files)
        
    def plot_model_strf(self, fold, output_file_prefix):
        md = self.modelData[fold]
        self.plot_strf(md, output_file_prefix)
        
    def plot_agg_strf(self, output_file_prefix):
        self.plot_strf(self.aggregateModel, output_file_prefix)

    def plot_strf(self, md, output_file_prefix):
        
        if not md.is_surprise:
            absstrf = np.abs(md.strf)
            strfabsmax = absstrf.max()
            
            ax = plt.subplot(111)
            ext = [md.timeLags.min(), md.timeLags.max(), md.low_freq/1000.0, md.high_freq/1000.0]
            self.subplot_strf(ax, ext, md.strf, strfabsmax)
             
        else:
            halfChans = int(md.numChans / 2)
            strf_louder = md.strf[0:halfChans, :]
            strf_quieter = md.strf[halfChans:, :]
            ext = [md.timeLags.min(), md.timeLags.max(), md.low_freq/1000.0, md.high_freq/1000.0]
            
            ax1 = plt.subplot(211)
            absstrf = np.abs(strf_louder)
            strfabsmax = absstrf.max()
            self.subplot_strf(ax1, ext, strf_louder, strfabsmax, title='Louder', show_xlabel=False)
            
            ax1 = plt.subplot(212)
            absstrf = np.abs(strf_quieter)
            strfabsmax = absstrf.max()
            self.subplot_strf(ax1, ext, strf_quieter, strfabsmax, title='Quieter', show_xlabel=True)
        
        fig = plt.gcf()
        fig.set_facecolor("#FFFFFF")     
        save_plot(fig, output_file_prefix)                
        plt.clf()
    
    def subplot_strf(self, ax, extent, strf, strfabsmax, title=None, show_xlabel=True):
        plt.imshow(strf, cmap=cm.jet, extent=extent, interpolation=None, \
                   origin='lower', vmin=-strfabsmax, vmax=strfabsmax, aspect='auto')
        
        if title is None:
            ax.set_title('STRF')
        else:
            ax.set_title(title)
        if show_xlabel:
            ax.set_xlabel('Time Lags (ms)')
        ax.set_ylabel('Frequency (kHz)')
        plt.colorbar()
        
    def plot_model_nl(self, fold, output_file_prefix):
        md = self.modelData[fold]
        self.plot_nl(md, output_file_prefix)
        
    def plot_agg_nl(self, output_file_prefix):
        self.plot_nl(self.aggregateModel, output_file_prefix)
    
    def plot_nl(self, md, output_file_prefix):        
        ax = plt.subplot(111)
        
        x = md.output_nl.domain
        y = md.output_nl.range
        y[np.isnan(y)] = 0.0
        y[np.isinf(y)] = 0.0
        y[y > 1000.0] = 1000.0

        ax.plot(x, y, 'k-', linewidth=3.0)
        ax.set_title('Output NL')
        ax.set_xlabel('Linear Output')
        ax.set_ylabel('Prediction')
        ax.set_ylim(0, md.output_nl.range.max())
        
        fig = plt.gcf()
        fig.set_facecolor("#FFFFFF")     
        save_plot(fig, output_file_prefix)    
        plt.clf()
        
    def plot_agg_pvals(self, output_file_prefix):
        am = self.aggregateModel
        ax = plt.subplot(111)
        ext = [am.timeLags.min(), am.timeLags.max(), am.low_freq/1000.0, am.high_freq/1000.0]
        
        plt.imshow(am.smoothedStrfPvals, cmap=cm.hot, interpolation=None, \
                   extent=ext, origin='lower', aspect='auto')
        plt.colorbar()
        ax.set_title('STRF p-values')
        ax.set_xlabel('Time Lags (ms)')
        ax.set_ylabel('Frequency (kHz)')
        fig = plt.gcf()
        fig.set_facecolor("#FFFFFF")     
        save_plot(fig, output_file_prefix)    
        plt.clf()
    
    def plot_response(self, fold, grp, output_file_prefix, unit_response=None):
        md = self.modelData[fold]
        
        gindx = (md.groupIndex == grp).squeeze()        
        mresp = md.response.squeeze()
        mresp = mresp[gindx]
                
        ax = plt.subplot(111)
        if unit_response is not None:
            if md.is_count_response:
                mresp /= unit_response.numTrials.mean()            
            ax.plot(unit_response.psth[gindx], 'k-', linewidth=2.0, label='Actual')
            ax.plot(mresp, 'r-', linewidth=2.0, alpha=0.75, label='Model')            
            ax.legend( ('Actual', 'Model') )
        else:        
            ax.plot(mresp, 'r-', linewidth=2.0)
        ax.set_xlim(0, len(mresp))        
        ax.set_ylim(0.0, 1.0)
        ax.set_title('Model Response')
        ax.set_xlabel('Time (ms)')
        ax.set_ylabel('PSTH')
        
        fig = plt.gcf()
        fig.set_facecolor("#FFFFFF")     
        save_plot(fig, output_file_prefix)    
        plt.clf()
        
    
def save_plot(fig, output_file_prefix):
    canvas = FigureCanvasAgg(fig)        
    png_path = output_file_prefix + '.png'
    canvas.print_png(png_path, dpi=72)        
    svg_path = output_file_prefix + '.svg'
    canvas.print_svg(svg_path, dpi=150)

def plot_info_bound_hist(info_bound_file):
    
    info_bounds_by_region = {}
    
    max_info = 0
    min_info = 300
    f = open(info_bound_file, 'r')
    cr = csv.reader(f, delimiter=',')
    for row in cr:
        reg = row[1]
        info = float(row[2])
        
        if info > max_info:
            max_info = info
        if info < min_info:
            min_info = info
        
        if reg not in info_bounds_by_region:
            info_bounds_by_region[reg] = []
        info_bounds_by_region[reg].append(info)
            
    f.close()
    
    for k,(reg,infos) in enumerate(info_bounds_by_region.iteritems()):
        ax = plt.subplot(111)
        infos = np.array(infos)
        imean = infos.mean()
        istd = infos.std()
        ax.hist(infos, bins=20)
        
        ax.set_xlim(min_info, max_info)
        ax.set_title('%s: mean=%0.1f +/- %0.1f bits/s' % (reg, imean, istd))
        ax.set_xlabel('Info (bits/s)')
        ax.set_ylabel('Count')
        
        ofile_prefix = os.path.join(PYSTRFS_DIR, 'web', 'images', 'info_bound_%s' % reg)
        
        fig = plt.gcf()
        fig.set_facecolor("#FFFFFF")     
        save_plot(fig, ofile_prefix)    
        plt.clf()
        
def plot_stimulus(preprocFile, md5, outputFile):
    
    f = h5py.File(preprocFile, 'r')
    
    is_stft = preprocFile.find('stft') > -1    
    is_lyons = preprocFile.find('lyons') > -1
    is_surprise = preprocFile.find('surprise') > -1
        
    spec = np.transpose(np.array(f[md5]['spectrogram']))
    
    all_specs = []
    all_titles = []
        
    if is_stft:

        dbnoise = float(f[md5].attrs['dbnoise'])
        low_freq = float(f[md5].attrs['low_freq'])
        high_freq = float(f[md5].attrs['high_freq'])
        
        refpow = spec.max()            
        spec = 20*np.log10(spec/refpow) + dbnoise;        
        spec[spec < 0] = 0
        
        all_specs.append(spec)
        all_titles.append('Spectrogram: %s' % md5)
    elif is_lyons:
        low_freq = float(f[md5].attrs['low_freq'])
        high_freq = float(f[md5].attrs['high_freq'])
        
        all_specs.append(spec)
        all_titles.append('Cochleagram: %s' % md5)
    
    elif is_surprise:
        low_freq = float(f[md5].attrs['low_freq'])
        high_freq = float(f[md5].attrs['high_freq'])
        
        halfChans = int(spec.shape[0] / 2)
        strf_louder = spec[0:halfChans, :]
        strf_quieter = spec[halfChans:, :]
        
        all_specs = [strf_louder, strf_quieter]        
        all_titles = ['Louder', 'Quieter']        
    
    else:
        print 'Unknown preproc type!'
        return
    
    f.close()
            
    plt.clf()
    plnum = (100*len(all_specs)) + 11
    for k,spec in enumerate(all_specs):        
        ax = plt.subplot(plnum + k)
        plt.imshow(spec, origin='lower', aspect='auto', extent=[0, spec.shape[1]-1, low_freq/1e3, high_freq/1e3])
        ax.set_xlabel('Time (ms)')
        ax.set_ylabel('Frequency (kHz)')
        ax.set_title(all_titles[k])
        if not is_surprise:
            plt.colorbar()
    
    if is_surprise:
        plt.subplots_adjust(hspace=0.35)
    
    fig = plt.gcf()
    fig.set_facecolor("#FFFFFF")     
    save_plot(fig, outputFile)    
    plt.clf()
    
            
        
def plot_neurogram(preprocFile, md5):
    
    f = h5py.File(preprocFile, 'r')
    
    unitNames = f[md5].keys()
    numUnits = len(unitNames)
    
    utlen = -1     
    unitTraces = {}
    for un in unitNames:
        unit = get_unit_by_old_id(un)
        unitRegions = [a.name for a in unit.recsite.areas]
        if len(unitRegions) > 1:
            if 'L' in unitRegions:
                region = 'L'
            else:
                region = unitRegions[0]
        else:
            region = unitRegions[0]
        if region not in unitTraces:
            unitTraces[region] = []
        ut = np.array(f[md5][un]['best']).squeeze()
        ut[np.isnan(ut)] = 0
        ut /= ut.max()
        utlen = len(ut)
        unitTraces[region].append(ut)        
    
    f.close()
    
    regs = ['MLd', 'OV', 'L', 'CM']    
    ngram = []    
    for reg in regs:
        for ut in unitTraces[reg]:            
            ngram.append(ut)
    ngram = np.array(ngram)
    
    plt.clf()
    ax = plt.subplot(111)
    plt.imshow(ngram, cmap=cm.jet, interpolation='nearest', \
                   origin='lower', vmin=0, vmax=1, aspect='auto')
    
    x = np.arange(utlen)
    y = np.ones(utlen)
    yinc = int(numUnits / len(regs))
    for k in range(len(regs)-1):
        ax = plt.gca()
        yln = yinc*(k+1)*y
        ax.plot(x, yln, 'w-')        
    
    ax.set_xlim(0, utlen-1)
    ax.set_ylim(0, numUnits-1)
    ax.set_title('Neurogram')
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('Cell #')
    
    fig = plt.gcf()
    fig.set_facecolor("#FFFFFF")     
    save_plot(fig, '/home/cheese63/Desktop/neurogram')
    plt.clf()
    
