
import os
import re
import glob

from sqlalchemy import Table, select, and_
from sqlalchemy.orm import mapper, relation, backref, deferred
from browser.model import meta

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

from pystrfs.modelrun import *
from pystrfs.analysis import *

class StepSize(object):
    def __init__(self):
        pass

step_sizes = Table('step_sizes', meta.metadata,
                   autoload=True, schema='analysis')

smapper = mapper(StepSize,step_sizes,
                 properties={'id':step_sizes.c.id,
                             'preproc':step_sizes.c.preproc,
                             'model':step_sizes.c.model,
                             'threshold':step_sizes.c.threshold,
                             'step_size':step_sizes.c.step_size,
                             'unit':step_sizes.c.unit})


def parse_output_file_name(outputFile):
    (rootDir, fname) = os.path.split(outputFile)
    fparts = fname.split('.')
    model = fparts[1]
    unit = fparts[2]
    preproc = fparts[3]

    m = re.search('thresh_[0-9].[0-9]*', fname)
    if m is None:
        print 'Problem getting threshold from file name: %s' % outputFile
        return None
    tstr = m.group(0)
    ts = tstr.split('_')
    thresh = float(ts[1])

    return (model, preproc, thresh, unit)

def parse_step_sizes_from_file(outputFile):

    if os.path.exists(outputFile):
        print 'Getting step sizes for %s' % outputFile
        (model, preproc, thresh, unit) = parse_output_file_name(outputFile)

        if model == 'lin+outputnl':
            model = 'linear'

        f = open(outputFile)
        for ln in f.readlines():
            step = None
            inum = None
            m = re.search('step: [0-9].[0-9]*', ln)
            if m is not None:
                mstr = m.group(0)
                t = mstr.split(':')
                step = float(t[1])

            m = re.search('iter: [0-9]*', ln)
            if m is not None:
                mstr = m.group(0)
                t = mstr.split(':')
                inum = int(t[1])

            if inum is not None and step is not None:
                ss = StepSize()
                ss.model = model
                ss.preproc = preproc
                ss.threshold = thresh
                ss.unit = unit
                ss.step_size = step
                ss.iter = inum
                meta.Session.add(ss)
        meta.Session.commit()
    else:
        print 'Output file does not exist: %s' % outputFile

def parse_step_sizes():

    flist = glob.glob('/auto/k6/mschachter/pystrfs/slurm/threshgrad/*.txt')
    for ofile in flist:
        parse_step_sizes_from_file(ofile)

def plot_step_sizes(lowIter = 0):

    preproc = ['stft', 'lyons', 'surprise']
    models = ['linear', 'poisson', 'binomial', 'leglm']
    thresh = [0.00, 0.25, 0.5, 0.75, 1.00]

    #for p in preproc:
    for t in thresh:
        stepSizes = []
        #sq = meta.Session.query(StepSize).filter_by(threshold=t, preproc=p).filter(StepSize.iter <= 100)
        sq = meta.Session.query(StepSize).filter_by(threshold=t).filter(StepSize.iter >= lowIter)
        for s in sq:
            stepSizes.append(s.step_size)
        stepSizes = np.array(stepSizes)
        if len(stepSizes) > 0:
            qval = sp.stats.scoreatpercentile(stepSizes, 25)
            mval = sp.stats.scoreatpercentile(stepSizes, 50)
            print '%0.2f: numsamples=%d, mean=%f +/- %f, quartile=%f, median=%f, min=%f' % (t, len(stepSizes), stepSizes.mean(), stepSizes.std(), qval, mval, stepSizes.min())
            """
            if len(stepSizes) > 0:
               fig = plt.figure()
               ax = fig.add_subplot(111)
               indx = np.abs(stepSizes) > (stepSizes.mean() + 2*stepSizes.std())
               (n, bins, patches) = ax.hist(stepSizes[indx], 50, normed = False, facecolor='green', alpha=0.75)
               plt.title('%0.2f: val=%f +/- %f' % (t, stepSizes.mean(), stepSizes.std()))
            """

        #plt.show()


def show_v2_performance(unit):

    (maxScore1, scoreValues1, scoreCube1, infoTrainCube1, infoValidCube1, modelCube1) = create_aggregate_score_cube(unit, 1)
    (maxScore2, scoreValues2, scoreCube2, infoTrainCube2, infoValidCube2, modelCube2) = create_aggregate_score_cube(unit, 2)

    preprocTypes = ['rawspectrogram', 'lyons', 'surprise']
    modelTypes = ['linear', 'sepnl_spline', 'sepnl_dists', 'poisson', 'binomial', 'leglm']
    thresholds = [0.0, 0.25, 0.5, 0.75, 1.0]

    scoreDiffs = []

    for pt in preprocTypes:
        for mtype in modelTypes:
            for tval in thresholds:
                if tval not in scoreValues1[pt][mtype] or tval not in scoreValues2[pt][mtype]:
                    print 'Missing value for preprocType=%s, modelType=%s' % (pt, mtype)
                    nscore = -1.0
                else:
                    scoreDiff = scoreValues2[pt][mtype][tval] - scoreValues1[pt][mtype][tval]
                    scoreDiffs.append(scoreDiff)

    scoreDiffs = np.array(scoreDiffs)
    scoreDiffsNeg = scoreDiffs[scoreDiffs < 0.0]
    scoreDiffsPos = scoreDiffs[scoreDiffs > 0.0]
    print 'Perf increase: %0.3f +/- %0.4f' % (scoreDiffs.mean(), scoreDiffs.std())
    print '# that did worse: %d, mean=%0.4f +/- %0.4f' % (len(scoreDiffsNeg), scoreDiffsNeg.mean(), scoreDiffsNeg.std())
    print '# that did better: %d, mean=%0.4f +/- %0.4f' % (len(scoreDiffsPos), scoreDiffsPos.mean(), scoreDiffsPos.std())


def show_v4_performance(unit):

    (maxScore1, scoreValues1, scoreCube1, infoTrainCube1, infoValidCube1, modelCube1) = create_aggregate_score_cube(unit, 2)
    (maxScore2, scoreValues2, scoreCube2, infoTrainCube2, infoValidCube2, modelCube2) = create_aggregate_score_cube(unit, 4)

    preprocTypes = ['rawspectrogram', 'lyons', 'surprise']
    modelTypes = ['linear', 'sepnl_spline', 'sepnl_dists', 'poisson', 'binomial', 'leglm']
    thresholds = [0.0, 0.25, 0.5, 0.75, 1.0]

    scoreDiffs = []
    timeDiffs = []

    for p,pt in enumerate(preprocTypes):
        for m,mtype in enumerate(modelTypes):
            for t,tval in enumerate(thresholds):
                if tval not in scoreValues1[pt][mtype]:
                    print 'Missing v2 value for preprocType=%s, modelType=%s' % (pt, mtype)
                    nscore = -1.0
                if tval not in scoreValues2[pt][mtype]:
                    print 'Missing v4 value for preprocType=%s, modelType=%s' % (pt, mtype)
                    nscore = -1.0
                else:
                    scoreDiff = scoreValues2[pt][mtype][tval] - scoreValues1[pt][mtype][tval]
                    scoreDiffs.append(scoreDiff)
                    rt1 = modelCube1[p][m][t].optimization.runtime
                    rt2 = modelCube2[p][m][t].optimization.runtime
                    timeDiffs.append(rt2 - rt1)

    scoreDiffs = np.array(scoreDiffs)
    scoreDiffsNeg = scoreDiffs[scoreDiffs < 0.0]
    scoreDiffsPos = scoreDiffs[scoreDiffs > 0.0]
    print 'Perf increase: %0.3f +/- %0.4f' % (scoreDiffs.mean(), scoreDiffs.std())
    print '# that did worse: %d, mean=%0.4f +/- %0.4f' % (len(scoreDiffsNeg), scoreDiffsNeg.mean(), scoreDiffsNeg.std())
    print '# that did better: %d, mean=%0.4f +/- %0.4f' % (len(scoreDiffsPos), scoreDiffsPos.mean(), scoreDiffsPos.std())

    timeDiffs = np.array(timeDiffs) / 60
    timeDiffsNeg = timeDiffs[timeDiffs < 0.0]
    timeDiffsPos = timeDiffs[timeDiffs > 0.0]
    print 'Avg. Runtime decrease for v4 vs v2: %0.1f min +/- %0.1f' % (timeDiffs.mean(), timeDiffs.std())
    print '# that did worse: %d, mean=%0.1f +/- %0.1f' % (len(timeDiffsPos), timeDiffsPos.mean(), timeDiffsPos.std())
    print '# that did better: %d, mean=%0.1f +/- %0.1f' % (len(timeDiffsNeg), timeDiffsNeg.mean(), timeDiffsNeg.std())

    return (scoreDiffs, timeDiffs)
