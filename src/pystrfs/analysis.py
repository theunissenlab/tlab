import re
import sys
import csv
from operator import itemgetter,attrgetter

import numpy as np

from scipy.signal import convolve2d
from scipy.interpolate import interp1d
from scipy.stats.distributions import t as tdist

import tables

from sqlalchemy import and_

from browser.model.analysis.model import Model, ResponseSet, find_response_set, ModelPerformance
from browser.model.analysis.strflab.models import StrflabLinearModel,StrflabSeparableNLModel, StrflabGlmModel, StrflabLogExpGlmModel
from browser.model.analysis.optimization import DirectFitOptimization,ThreshgradLSOptimization

from pystrfs import *
from pystrfs.preprocess import get_preprocessor_from_filename
from pystrfs.modelrun import *

def create_goodcells_file():

    regions = ['L', 'MLd', 'OV', 'CM']
    regUnits = {'L':[], 'MLd':[], 'OV':[], 'CM':[]}

    infoThresh = 10.0

    for reg in regions:
        fpath = os.path.join(UNITS_DIR, 'performance.%s.Con.csv' % reg)
        f = open(fpath, 'r')
        for ln in f.readlines():
            ln = ln.strip()
            if len(ln) > 0:
                vals = ln.split(',')
                info = float(vals[1])
                if info >= infoThresh:
                    regUnits[reg].append(vals[0])
    maxLen = 0
    for reg in regions:
        maxLen = max([maxLen, len(regUnits[reg])])
    fout = os.path.join(UNITS_DIR, 'allgoodcells.csv')
    f = open(fout, 'w')
    for k in range(maxLen):
        for reg in regions:
            if len(regUnits[reg]) > k:
                f.write('%s\n' % regUnits[reg][k])

    f.close()



def create_response_performance_files():

    regionExprs = ['L*', 'MLd', 'OV', 'CM']
    regions = ['L', 'MLd', 'OV', 'CM']
    regUnits = {'L':[], 'MLd':[], 'OV':[], 'CM':[]}

    umap = get_all_dcp_units(True)

    for uId,u in umap.iteritems():
        if len(u.recsite.areas) > 0:
            reg = u.recsite.areas[0].name
            rkey = None
            for k,rexpr in enumerate(regionExprs):
                m = re.match(rexpr, reg)
                if m is not None and len(m.group(0)) > 0:
                    rkey = regions[k]
                    break
            if rkey is not None:
                #print '%s associated with %s' % (uId, rkey)
                regUnits[rkey].append(u)
            else:
                print 'Cannot find regexp associated with region %s' % reg

        else:
            print 'Unit has no region: %s' % u.recsite.old_id

    for regName,ulist in regUnits.iteritems():
        conInfoMap = {}
        for unit in ulist:
            for classResp in unit.class_resps:
                stimGroup = stimClassToGroup[classResp.presentation.class_name]
                if stimGroup == 'Con':
                    conInfoMap[unit.recsite.old_id] = classResp.performance.info_lower
        infoSorted = sorted(conInfoMap.iteritems(), key=itemgetter(1), reverse=True)

        fpath = os.path.join(UNITS_DIR, 'performance.%s.Con.csv' % regName)
        f = open(fpath, 'w')
        for isort in infoSorted:
            f.write('%s,%f\n' % (isort[0], isort[1]))
        f.close()

def compute_output_nl_fig_data(unit, slurmBot, neurogram=False, min_runs=False):

    mr = ModelRunner()
    if neurogram:
        mruns = mr.create_threshgrad_runs_all_neurogram(unit)
    else:
        if min_runs:
            mruns = mr.create_threshgrad_runs_min(unit)
        else:
            mruns = mr.create_threshgrad_runs_all(unit)

    output_files = []
    template_file = os.path.join(TEMPLATE_DIR, 'output_nl_fig.m')

    output_dir = ''
    for mrun in mruns:
        ofiles = mrun.get_output_file()

        for ofile in ofiles:
            if not os.path.exists(ofile):
                print 'File not found: %s' % ofile
            else:
                (output_dir, fname) = os.path.split(ofile)
                output_files.append(fname)

    file_list_name = os.path.join(output_dir, 'output_file_list.csv')
    f = open(file_list_name, 'w')
    for ofile in output_files:
        f.write('%s\n' % ofile)
    f.close()

    ofparam = ','.join(['\'%s\'' % s for s in output_files])
    tparams = {'SRC_ROOT_DIR': SRC_ROOT_DIR}
    tparams['OUTPUT_DIR'] = output_dir
    tparams['FILE_LIST'] = file_list_name
    mt = MatlabTemplate()
    mt.template_from_file(template_file)
    mp = mt.to_matlab_process(tparams)
    
    fileDesc = 'output_nl_fig_%s' % unit.recsite.old_id
    sout = os.path.join(SLURM_DIR, '%s.%%j.out.txt' % fileDesc)
    serr = sout
    scriptName = os.path.join(SLURM_DIR, '%s.sbatch' % fileDesc)

    sparams = {'out': sout,
               'err': serr,
               'partition':'all',
               'cpus': 1,
               'script': scriptName}

    slurmBot.add(mp.get_commands(), sparams)


def create_opt_df_from_h5file(f):
    opt = None
    optMethod = ''.join(f['opt'].attrs['method'])
    if optMethod == 'direct_fit':
        infoFreqCutoff = f['opt']['info_freq_cutoff'][0][0]
        infoWindowSize = f['opt']['info_window_size'][0][0]
        spStrs = ','.join([str(x[0]) for x in f['opt']['sparsenesses']])
        tolStrs = ','.join([str(x[0]) for x in f['opt']['tolerances']])

        opt = DirectFitOptimization()
        opt.info_freq_cutoff = infoFreqCutoff
        opt.info_window_size = infoWindowSize
        opt.sparsenesses = spStrs
        opt.tolerances = tolStrs
        opt.bias_high_tolerance = 0

    return opt

def create_opt_tg_from_h5file(f):
    opt = None
    optMethod = ''.join(f['opt'].attrs['method'])
    if optMethod == 'threshgrad':
        numItersCv = f['opt']['num_iters'][0]
        maxIters = int(f['opt']['max_iter'][0][0])
        runTimes = f['opt']['run_times'][0]
        threshold = float(f['opt']['threshold'][0])
        if 'version' in f['opt']:
            version = int(f['opt']['version'][0])
        else:
            version = 1

        if 'final_run_time' in f['opt']:
            finalRunTime = float(f['opt']['final_run_time'][0])
        else:
            #estimate the final run time
            rt = np.array(runTimes)
            ni = np.array(numItersCv)
            timePerIter = rt / ni
            print 'Avg. Time per Iteration: %0.2f +/- %0.2f' % (timePerIter.mean(), timePerIter.std())
            finalRunTime = maxIters * (timePerIter.mean() + 0.5*timePerIter.std())

        numItersCvStr = ','.join(['%d' % x for x in numItersCv])
        runtimesStr = ','.join(['%d' % x for x in runTimes])

        opt = ThreshgradLSOptimization()
        opt.num_iters_cv = numItersCvStr
        opt.runtimes_cv = runtimesStr
        opt.runtime = finalRunTime
        opt.max_iters = maxIters
        opt.num_cv_partitions = len(numItersCv)
        opt.threshold = threshold
        opt.version = version

    return opt



def get_or_create_response_set(resps):

    idList = [r.id for r in resps]
    rset = find_response_set(idList)
    if rset is None:
        rset = ResponseSet()
        for r in resps:
            rset.members.append(r)
        meta.Session.add(rset)
        meta.Session.commit()
    return rset


def map_output_to_model(unit, outputFile):

    if not os.path.exists(outputFile):
        print 'Output file not found: %s' % outputFile

    #check if model already exists
    mq = meta.Session.query(Model).filter_by(directory=outputFile)
    if mq.count() > 0:
        return mq.one()

    (rootDir, fileName) = os.path.split(outputFile)
    fileParts = fileName.split('.')

    #read file and create corresponding Preprocessor, ModelInstance,
    #OptimizationInstance, and ResponseSet objects
    try:
        f = h5py.File(outputFile, 'r')
    except IOError:
        print 'Error opening file %s' % outputFile
        return None

    #get preprocessor
    preprocFile = ''.join(f['data'].attrs['preproc_file'])
    preproc = get_preprocessor_from_filename(preprocFile)

    #create optimization instance
    opt = None
    optMethod = ''.join(f['opt'].attrs['method'])
    if 'type' not in f['model'].attrs:
        modelType = 'linear'
    else:
        modelType = ''.join(f['model'].attrs['type'])

    if optMethod == 'direct_fit':
        if modelType in ['nl_dists', 'nl_spline']:
            origFile = ''.join(f['model'].attrs['original_response_file'])
            forig = h5py.File(origFile, 'r')
            opt = create_opt_df_from_h5file(forig)
            forig.close()
        else:
            opt = create_opt_df_from_h5file(f)
    elif optMethod == 'threshgrad':
        if modelType in ['nl_dists', 'nl_spline']:
            origFile = ''.join(f['model'].attrs['original_response_file'])
            forig = h5py.File(origFile, 'r')
            opt = create_opt_tg_from_h5file(forig)
            forig.close()
        else:
            opt = create_opt_tg_from_h5file(f)

    else:
        print 'Unknown optimization method %s' % optMethod
        return None


    #create or get ResponseSet instances for training and validation
    stimGroup = ''.join(f['data'].attrs['stim_class'])
    stimClass = stimGroupToClass[stimGroup]
    trainingResps = []
    validationResps = []
    trainingGroups = [int(x[0]) for x in f['data']['training_groups']]
    validationGroups = [int(x[0]) for x in f['data']['validation_groups']]

    classResp = None
    for cp in unit.class_resps:
        if cp.class_name == stimClass:
            classResp = cp
            break

    resps = [r for r in cp.repeated]

    for tg in trainingGroups:
        trainingResps.append(resps[tg-1])
    for vg in validationGroups:
        validationResps.append(resps[vg-1])

    trainingSet = get_or_create_response_set(trainingResps)
    validationSet = get_or_create_response_set(validationResps)

    #create ModelInstance
    if modelType == 'linear':
        strfLen = int(f['model'].attrs['strf_length'])
        delaysExpr = '%d:%d' % (0, strfLen-1)
        numChans = f['model']['weights'].shape[1]

        model = StrflabLinearModel()
        model.delays_expr = delaysExpr
        model.num_chans = numChans
        model.output_nl = 'linear'
        model.freq_dom = 0

    elif modelType in ['binomial', 'poisson']:
        strfLen = int(f['model'].attrs['strf_length'])
        delaysExpr = '%d:%d' % (0, strfLen-1)
        numChans = f['model']['weights'].shape[1]

        model = StrflabGlmModel()
        model.delays_expr = delaysExpr
        model.num_chans = numChans
        if modelType == 'poisson':
            model.output_nl_expr = 'exp(x)'
        elif modelType == 'binomial':
            model.output_nl_expr = '1 ./ (1 + exp(-x))'
        model.dispersion = 1.0
        model.family = modelType

    elif modelType == 'leglm':
        strfLen = int(f['model'].attrs['strf_length'])
        delaysExpr = '%d:%d' % (0, strfLen-1)
        numChans = f['model']['weights'].shape[1]
        m = f['model']['m'][0][0]

        model = StrflabLogExpGlmModel()
        model.delays_expr = delaysExpr
        model.num_chans = numChans
        model.dispersion = 1.0
        model.family = modelType
        model.m = m

    elif modelType in ['nl_dists', 'nl_spline']:
        matFile = ''.join(f['model'].attrs['mat_file'])

        model = StrflabSeparableNLModel()
        model.mat_file = matFile
        model.fit_type = modelType

    else:
        print 'Cannot handle model type %s' % modelType
        return None

    model.preprocesser = preproc
    model.optimization = opt
    model.training_set = trainingSet
    model.validation_set = validationSet
    model.directory = outputFile


    #add model to db
    meta.Session.add(opt)
    meta.Session.add(model)
    meta.Session.commit()

    #map model performance
    trainInfoLower = f['model']['performance']['coherence']['training']['info_lower'][0][0]
    trainInfoUpper = f['model']['performance']['coherence']['training']['info_upper'][0][0]
    trainInfoMean = f['model']['performance']['coherence']['training']['info_mean'][0][0]

    trainBoundInfoLower = f['model']['performance']['coherence']['bound']['training']['info_lower'][0][0]
    trainBoundInfoUpper = f['model']['performance']['coherence']['bound']['training']['info_upper'][0][0]
    trainBoundInfoMean = f['model']['performance']['coherence']['bound']['training']['info_mean'][0][0]

    validInfoLower = f['model']['performance']['coherence']['validation']['info_lower'][0][0]
    validInfoUpper = f['model']['performance']['coherence']['validation']['info_upper'][0][0]
    validInfoMean = f['model']['performance']['coherence']['validation']['info_mean'][0][0]

    validBoundInfoLower = f['model']['performance']['coherence']['bound']['validation']['info_lower'][0][0]
    validBoundInfoUpper = f['model']['performance']['coherence']['bound']['validation']['info_upper'][0][0]
    validBoundInfoMean = f['model']['performance']['coherence']['bound']['validation']['info_mean'][0][0]

    f.close()

    mpt = ModelPerformance()
    mpt.response_set = model.training_set
    mpt.model = model
    mpt.info_lower = trainInfoLower
    mpt.info_upper = trainInfoUpper
    mpt.info_mean = trainInfoMean
    mpt.info_bound_lower = trainBoundInfoLower
    mpt.info_bound_upper = trainBoundInfoUpper
    mpt.info_bound_mean = trainBoundInfoMean

    meta.Session.add(mpt)

    mpv = ModelPerformance()
    mpv.response_set = model.validation_set
    mpv.model = model
    mpv.info_lower = validInfoLower
    mpv.info_upper = validInfoUpper
    mpv.info_mean = validInfoMean
    mpv.info_bound_lower = validBoundInfoLower
    mpv.info_bound_upper = validBoundInfoUpper
    mpv.info_bound_mean = validBoundInfoMean

    meta.Session.add(mpv)

    meta.Session.commit()

    return model


def load_unit_output_into_db_df(unit):

    mr = ModelRunner()
    mruns = mr.create_directfit_runs(unit)

    for mrun in mruns:
        ofiles = mrun.get_output_file()

        for ofile in ofiles:
            print 'Mapping %s' % ofile
            map_output_to_model(unit, ofile)

def load_unit_output_into_db_tg(unit, min_runs=False):

    mr = ModelRunner()
    if min_runs:
        mruns = mr.create_threshgrad_runs_min(unit)
    else:
        mruns = mr.create_threshgrad_runs_all(unit)

    for mrun in mruns:
        ofiles = mrun.get_output_file()

        for ofile in ofiles:
            if not os.path.exists(ofile):
                print 'File not found: %s' % ofile
            else:
                print 'Mapping %s' % ofile
                map_output_to_model(unit, ofile)
                
def load_unit_output_into_db_tg_neurogram(unit):

    mr = ModelRunner()
    mruns = mr.create_threshgrad_runs_all_neurogram(unit)

    for mrun in mruns:
        ofiles = mrun.get_output_file()

        for ofile in ofiles:
            if not os.path.exists(ofile):
                print 'File not found: %s' % ofile
            else:
                print 'Mapping %s' % ofile
                map_output_to_model(unit, ofile)


def recompute_response_info(unit, sb = None):

    templateFile = os.path.join(TEMPLATE_DIR, 'compute_response_info.m')

    if sb is None:
        sb = SlurmBot()
    mr = ModelRunner()
    mruns = mr.create_threshgrad_runs_all(unit)
    for mrun in mruns:
        ofiles = mrun.get_output_file()
        for of in ofiles:
            if os.path.exists(of):
                mt = MatlabTemplate()
                mt.template_from_file(templateFile)
                tparams = {'SRC_ROOT_DIR':SRC_ROOT_DIR, 'RESPONSE_FILE': of}
                mp = mt.to_matlab_process(tparams)
                (root, fname) = os.path.split(of)
                fileDesc = '%s.%s.txt' % (unit.recsite.old_id, fname)
                sout = os.path.join(SLURM_DIR, 'recompute_response_info.%s.%%j.out.txt' % fileDesc)
                serr = sout
                scriptName = os.path.join(SLURM_DIR, '%s.sbatch' % fileDesc)

                sparams = {'out': sout,
                           'err': serr,
                           'partition':'all',
                           'cpus': 1,
                           'script': scriptName}

                sb.add(mp.get_commands(), sparams)

    return sb


def get_model_list(unit, criteria):

    lstr = '%%%s%%' % unit.recsite.old_id
    allCriteria = [Model.directory.like(lstr)]
    if criteria is not None:
        allCriteria.extend(criteria)

    mq = meta.Session.query(Model).filter_by(unit_id=unit.id).filter(and_(*allCriteria))

    return mq.all()


def get_model_list_kw(unit, modelClass, criteria):
    lstr = '%%/%s/%%' % unit.recsite.old_id
    mq = meta.Session.query(modelClass).filter(modelClass.directory.like(lstr)).filter_by(**criteria)
    return mq.all()


def find_best_lyons(unit, debug = False):

    crit = [Model.preprocess_class == 'lyons', Model.model_class == 'strflab_linear']
    mlist = get_model_list(unit, crit)

    #aggregate by earq and step
    mmap = {}
    for m in mlist:
        k = (m.preprocesser.earq, m.preprocesser.step)
        if k not in mmap:
            mmap[k] = []
        mmap[k].append(m)

    mstats = {}

    #quantify average performance
    mperfs = []
    for k,models in mmap.iteritems():
        plst = [m.validation_performance.info_mean for m in models]
        perfs = np.array(plst)
        mstats[k] = perfs
        mperfs.append((k, perfs.mean()))
        if debug:
            print 'earq=%d, step=%0.2f: info=%0.1f +/- %0.2f' % (k[0], k[1], perfs.mean(), perfs.std())

    #sort
    mperfs.sort(key=itemgetter(1), reverse=True)

    bestEarq = mperfs[0][0][0]
    bestStep = mperfs[0][0][1]

    if debug:
        print 'Best Lyons params: earq=%d, step=%0.2f' % (bestEarq, bestStep)

    #write to a file
    unitDir = get_unit_dir(unit)
    outputFile = os.path.join(unitDir, 'output', 'bestlyons.txt')
    f = open(outputFile, 'w')
    f.write('earq=%d\nstep=%0.2f\n' % (bestEarq, bestStep))
    f.close()

    return mstats

def get_lyons_stats(units):

    allStats = {}
    for u in units:
        mstats = find_best_lyons(u, False)
        for k,perfs in mstats.iteritems():
            if k not in allStats:
                allStats[k] = np.array([])
            allStats[k] = np.append(allStats[k], perfs)

    for k,stats in allStats.iteritems():
        print 'earq=%d, step=%0.2f: info=%0.1f +/- %0.1f' % (k[0], k[1], stats.mean(), stats.std())

    return allStats


def get_all_models_tg(unit, preprocClass, optVersion = 5, modelTypes=None, thresholds=None):
    models = {'linear': get_model_list_kw(unit, StrflabLinearModel, {'optimization_class':'threshgradls', 'preprocess_class':preprocClass}),
              'sepnl_spline': get_model_list_kw(unit, StrflabSeparableNLModel, {'optimization_class':'threshgradls', 'fit_type':'nl_spline', 'preprocess_class':preprocClass}),
              'sepnl_dists': get_model_list_kw(unit, StrflabSeparableNLModel, {'optimization_class':'threshgradls', 'fit_type':'nl_dists', 'preprocess_class':preprocClass}),
              'poisson': get_model_list_kw(unit, StrflabGlmModel, {'optimization_class':'threshgradls', 'family':'poisson', 'preprocess_class':preprocClass}),
              'binomial': get_model_list_kw(unit, StrflabGlmModel, {'optimization_class':'threshgradls', 'family':'binomial', 'preprocess_class':preprocClass}),
              'leglm': get_model_list_kw(unit, StrflabLogExpGlmModel, {'optimization_class':'threshgradls', 'preprocess_class':preprocClass})}

    for mname,mlist in models.iteritems():        
        finalList = []
        for k,m in enumerate(mlist):
            if m.optimization.optimization_class == 'threshgradls':
                if m.optimization.version == optVersion:
                    if thresholds is None or m.optimization.threshold in thresholds:
                        finalList.append(m)
        del models[mname]
        models[mname] = finalList
        
    final_models = {}
    for mname,mlist in models.iteritems():
        if modelTypes is None or mname in modelTypes:
            final_models[mname] = mlist
    

    return final_models

def create_model_cube_tg(unit, optVersion = 5, extraPreprocs=[],
                         preprocTypes=None, modelTypes=None, thresholds=None):
    if preprocTypes is None:
        preprocTypes = ['rawspectrogram', 'lyons', 'surprise']
    if modelTypes is None:
        modelTypes = ['linear', 'sepnl_spline', 'sepnl_dists', 'poisson', 'binomial', 'leglm']
    if thresholds is None:
        thresholds = [0.0, 0.25, 0.5, 0.75, 1.0]
    
    for ep in extraPreprocs:
        preprocTypes.append(ep)

    #initialize empty model cube
    modelCube = [None] * len(preprocTypes)
    for k in range(len(preprocTypes)):
        modelCube[k] = [None] * len(modelTypes)
        for p in range(len(modelTypes)):
            modelCube[k][p] = [None] * len(thresholds)
            for q in range(len(thresholds)):
                modelCube[k][p][q] = []

    #fill in the model cube
    for k,pt in enumerate(preprocTypes):
        models = get_all_models_tg(unit, pt, optVersion, modelTypes=modelTypes, thresholds=thresholds)
        for p,mtype in enumerate(modelTypes):
            mlist = models[mtype]
            for m in mlist:
                tval = m.optimization.threshold
                q = thresholds.index(tval)
                modelCube[k][p][q].append(m)

    return modelCube


def create_scoring_cubes_tg(unit, optVersion = 5, extraPreprocs=[],
                            preprocTypes=None, modelTypes=None, thresholds=None):

    if preprocTypes is None:
        preprocTypes = ['rawspectrogram', 'lyons', 'surprise']
    for ep in extraPreprocs:
        preprocTypes.append(ep)
    if modelTypes is None:
        modelTypes = ['linear', 'sepnl_spline', 'sepnl_dists', 'poisson', 'binomial', 'leglm']
    if thresholds is None:
        thresholds = [0.0, 0.25, 0.5, 0.75, 1.0]

    modelCube = create_model_cube_tg(unit, optVersion=optVersion, extraPreprocs=extraPreprocs,
                                     preprocTypes=preprocTypes, modelTypes=modelTypes, thresholds=thresholds)

    scoreCube = np.zeros([len(preprocTypes), len(modelTypes), len(thresholds)])
    infoTrainCube = np.zeros([len(preprocTypes), len(modelTypes), len(thresholds)])
    infoValidCube = np.zeros([len(preprocTypes), len(modelTypes), len(thresholds)])

    maxScore = 0.0    
    
    scoreValues = {}
    
    #compute score across all types
    for k,pt in enumerate(preprocTypes):
        scoreValues[pt] = {}
        for p,mtype in enumerate(modelTypes):
            scoreValues[pt][mtype] = {}

            for q,tval in enumerate(thresholds):

                trats = []
                vrats = []
                tivals = []
                vivals = []
                for m in modelCube[k][p][q]:
                    try:
                        tivals.append(m.training_performance.info_mean)
                        vivals.append(m.validation_performance.info_mean)
                        trat = m.training_performance.info_mean / m.training_performance.info_bound_mean
                        vrat = m.validation_performance.info_mean / m.validation_performance.info_bound_mean
                        trats.append(trat)
                        vrats.append(vrat)
                    except:
                        print 'Problem obtaining perf info: unit=%s, preproc=%s, model=%s, thresh=%0.2f' %\
                              (unit.recsite.old_id, pt, mtype, tval)
                        print sys.exc_info()[0]
                        print sys.exc_info()[1]
                        print sys.exc_info()[2]

                
                meanTrainingInfo = 0.0
                meanValidationInfo = 0.0
                vperfs = []
                meanScore = 0.0
                
                if len(modelCube[k][p][q]) > 0:
                    trats = np.array(trats)
                    vrats = np.array(vrats)
                    tivals = np.array(tivals)
                    vivals = np.array(vivals)

                    mscore = vrats.mean()
                    if mscore > maxScore:
                        maxScore = mscore
                    
                    vperfs = vrats
                    meanScore = mscore
                    meanTrainingInfo = tivals.mean()
                    meanValidationInfo = vivals.mean()                    
                
                scoreValues[pt][mtype][tval] = vperfs
                #print (pt, mtype, tval, k, p, q) #good for debugging
                scoreCube[k][p][q] = meanScore
                infoTrainCube[k][p][q] = meanTrainingInfo
                infoValidCube[k][p][q] = meanValidationInfo

    return (maxScore, scoreValues, scoreCube, infoTrainCube, infoValidCube, modelCube)


def create_info_bound_file(outputFile):
    
    f = open(outputFile, 'w')
    units = get_good_units()
    for u in units:        
        stimGroups = [stimClassToGroup[cr.presentation.class_name] for cr in u.class_resps]         
        stimIndx = stimGroups.index('Con')
        classResp = u.class_resps[stimIndx]
        
        info_mean = classResp.performance.info_mean
        unitRegions = [a.name for a in u.recsite.areas]
        if len(unitRegions) > 1:
            if 'L' in unitRegions:
                region = 'L'
            else:
                region = unitRegions[0]
        else:
            region = unitRegions[0]
            
        f.write('%s,%s,%0.6f\n' % (u.recsite.old_id, region, info_mean))
        
    f.close()
        
    

def examine_tg_runtimes(unit):
    mcrit = [Model.model_class == 'strflab_linear', Model.optimization_class == 'threshgradls']
    mlist = get_model_list(unit, mcrit)

    for m in mlist:
        rtcvStr = m.optimization.runtimes_cv
        rtcv = np.array([float(x) for x in rtcvStr.split(',')])
        totalRuntimes = rtcv.sum()
        numIter = m.optimization.max_iters
        runTime = m.optimization.runtime
        (rootDir, fileName) = os.path.split(m.directory)
        print '%s: numiter=%d, runtime=%0.0fm, total=%0.1f h' % (fileName, numIter, runTime / 60, (totalRuntimes + runTime) / 3600)


def create_anova_data(units, outputFile, optVersion = 5, extraPreprocs=[]):

    preprocTypes = ['rawspectrogram', 'lyons', 'surprise']
    modelTypes = ['linear', 'sepnl_spline', 'sepnl_dists', 'poisson', 'binomial', 'leglm']
    thresholds = [0.0, 0.25, 0.5, 0.75, 1.0]
    
    for ep in extraPreprocs:
        preprocTypes.append(ep)

    f = open(outputFile, 'w')
    f.write('unit,region,preproc,model,thresh,score,info_bound,num_iters,runtime\n')

    for u in units:

        (maxScore, scoreValues, scoreCube, infoTrainCube, infoValidCube, modelCube) = create_scoring_cubes_tg(u, optVersion, extraPreprocs=extraPreprocs)

        region = get_unit_str(u)
            
        stimGroups = [stimClassToGroup[cr.presentation.class_name] for cr in u.class_resps]         
        stimIndx = stimGroups.index('Con')
        classResp = u.class_resps[stimIndx]
        
        info_bound_mean = classResp.performance.info_mean

        for k,pt in enumerate(preprocTypes):
            for p,mt in enumerate(modelTypes):
                for q,tval in enumerate(thresholds):
                    svals = scoreValues[pt][mt][tval]
                    if len(svals) > 0:
                        
                        avg_iters = []
                        runtimes = []
                        for m in modelCube[k][p][q]:
                            nis = [float(s) for s in m.optimization.num_iters_cv.split(',')]
                            nis = np.array(nis)
                            avg_iters.append(nis.mean())
                            runtimes.append(m.optimization.runtime)
                        avg_iters = np.array(avg_iters)            
                        runtimes = np.array(runtimes)
                        score = svals.mean()
                        runtime = runtimes.mean()
                        print '%s,%s,%s,%s,%0.2f,%0.6f,%0.6f,%0.2f,%0.2f' % (u.recsite.old_id, region, pt, mt, tval, score, info_bound_mean, avg_iters.mean(), runtimes.mean())
                        f.write('%s,%s,%s,%s,%0.2f,%0.6f,%0.6f,%0.2f,%0.2f\n' % (u.recsite.old_id, region, pt, mt, tval, score, info_bound_mean, avg_iters.mean(), runtimes.mean()))

    f.close()


def get_unit_str(u):
    unitRegions = [a.name for a in u.recsite.areas]
    if len(unitRegions) > 1:
        if 'L' in unitRegions:
            region = 'L'
        else:
            region = unitRegions[0]
    else:
        region = unitRegions[0]
    return region


def get_strf_and_nl(output_file):
        
    output_nl_range = None
    output_nl_domain = None
        
    of = h5py.File(output_file, 'r')
    if 'original_response_file' in of['model'].attrs:
        output_nl_range = np.array(of['model']['output_nl']['range'])
        output_nl_domain = np.array(of['model']['output_nl']['domain'])
        orig_file = ''.join(of['model'].attrs['original_response_file'])
        f = h5py.File(orig_file)
        of.close()
    else:
        output_nl_range = np.array(of['model']['output_nl']['range'])
        output_nl_domain = np.array(of['model']['output_nl']['domain'])
        f = of                    
    
    biases = np.array(f['model']['cv_bias'])
    bias = biases.mean()
    strfs = np.array(f['model']['cv_weights'])
    strf = np.transpose(strfs.mean(axis=2))
    f.close()
    
    return (strf, bias, output_nl_range, output_nl_domain)   


class Struct:
    def __init__(self):
        pass

class OutputNL:
    def __init__(self):
        self.domain = None
        self.range = None

class SingleModelData:

    def __init__(self):
        self.unit = None
        self.model = None
        self.output_nl = None
        self.preprocFile = None
        self.stimClass = None
        self.strf = None
        self.timeLags = None
        self.freqs = None
        self.aggregate = False
        self.strfStd = None
        self.trainingIndex = None
        self.validationIndex = None
        self.groupIndex = None
        self.response = None

class AggregateModel:

    def __init__(self, unit, output_files):

        self.failed = False
        self.modelData = []
        self.unit = unit
        self.unitFile = map_dcp_unit_to_fs(unit)

        for output_file in output_files:            
            self.read_output_file(output_file)

        if len(self.modelData) > 0:
            self.build_aggregate_model()            
        else:
            self.failed = True

    def read_output_file(self, output_file):

        if not os.path.exists(output_file):
            print 'File does not exist: %s' % output_file
            return
        f = h5py.File(output_file, 'r')
        smd = SingleModelData()
        smd.outputFile = output_file
        smd.unit = self.unit
        try:        
            smd.preprocFile = ''.join(f['data'].attrs['preproc_file'])
            smd.stimClass = ''.join(f['data'].attrs['stim_class'])            
            
            smd.is_surprise = smd.preprocFile.find('surprise') >= 0
            smd.is_count_response = (smd.outputFile.find('poisson') >= 0) or (smd.outputFile.find('leglm') >= 0)
            smd.is_np_nl = (smd.outputFile.find('nl_dist') >= 0) or (smd.outputFile.find('nl_spline') >= 0)
            
            smd.trainingIndex = f['data']['training_index']
            smd.validationIndex = f['data']['validation_index']
            smd.trainingGroups = np.array(f['data']['training_groups']).squeeze()
            smd.validationGroups = np.array(f['data']['validation_groups']).squeeze()
            smd.groupIndex = np.array(f['data']['group_index'])
            smd.response = np.array(f['model']['response'])
    
            pf = h5py.File(smd.preprocFile, 'r')
            k1 = pf.keys()[0] #any element will do
            try:
                high_freq = pf[k1].attrs['high_freq']
                low_freq = pf[k1].attrs['low_freq']
            except:
                print 'Problem getting high_freq and low_freq from %s' % smd.preprocFile
                high_freq = -1
                low_freq = -1
            pf.close()        
    
            # get model
            mq = meta.Session.query(Model).filter_by(directory=output_file)
            if mq.count() == 0:
                print 'Cannot locate model corresponding to output file: %s' % output_file
            else:
                smd.model = mq.one()
    
            # get STRF properties
            if smd.is_np_nl:
                #open up original response file
                orfname = ''.join(f['model'].attrs['original_response_file'])
                orf = h5py.File(orfname, 'r')
                strf = np.transpose(np.array(orf['model']['weights']))
                orf.close()            
            else:
                strf = np.transpose(np.array(f['model']['weights']))
            smd.strf = strf
    
            numTimePts = strf.shape[1]
            numChans = strf.shape[0]
            if high_freq == -1:
                high_freq = numChans
                low_freq = 1
            smd.timeLags = np.array(range(0, numTimePts))
            smd.numChans = numChans
            smd.high_freq = high_freq
            smd.low_freq = low_freq
            if not smd.is_surprise:
                smd.freqs = np.linspace(low_freq, high_freq, numChans)
            else:
                hc = int(numChans / 2)
                freq = np.linspace(low_freq, high_freq, numChans)
                smd.freqs = np.array([freq, freq])
    
            output_nl = OutputNL()
            output_nl.domain = np.array(f['model']['output_nl']['domain'])
            output_nl.range = np.array(f['model']['output_nl']['range'])
    
            smd.output_nl = output_nl
    
            self.modelData.append(smd)
        except:  
            print 'Problem with model file %s' % output_file
            print sys.exc_info()[0]
            print sys.exc_info()[1]
            print sys.exc_info()[2]
        try:
            f.close()
        except:
            pass


    def build_aggregate_model(self):        
        aggModel = SingleModelData()
        aggModel.nfolds = len(self.modelData)
        aggModel.aggregate = True
        aggModel.unit = self.unit
        aggModel.preprocFile = self.modelData[0].preprocFile
        aggModel.stimClass = self.modelData[0].stimClass
        aggModel.freqs = self.modelData[0].freqs
        aggModel.numChans = self.modelData[0].numChans
        aggModel.timeLags = self.modelData[0].timeLags
        aggModel.high_freq = self.modelData[0].high_freq
        aggModel.low_freq = self.modelData[0].low_freq
        aggModel.is_surprise = self.modelData[0].is_surprise
        aggModel.is_count_response = self.modelData[0].is_count_response
        aggModel.output_nls = [md.output_nl for md in self.modelData]

        #compute aggregate STRF params
        strfs = np.array([smd.strf for smd in self.modelData])
        aggModel.strf = strfs.mean(axis=0).squeeze()
        aggModel.strfStd = strfs.std(axis=0).squeeze()
        
        smoothedStrfs = []
        g1 = gaussian_2d_kernel(1)
        for strf in strfs:
            sstrf = convolve2d(strf, g1, mode='same')
            smoothedStrfs.append(sstrf)
        smoothedStrfs = np.array(smoothedStrfs)
        
        aggModel.smoothedStrf = smoothedStrfs.mean(axis=0).squeeze()
        aggModel.smoothedStrfStd = smoothedStrfs.std(axis=0).squeeze()
                
        strf_tstat = np.abs(aggModel.smoothedStrf / aggModel.smoothedStrfStd)
        df = len(self.modelData) - 1
        strf_pvals = (1 - tdist.cdf(strf_tstat, df))*2
        aggModel.smoothedStrfPvals = strf_pvals

        #compute aggregate output nl
        minx_vals = []
        maxx_vals = []
        for nl in aggModel.output_nls:
            minx_vals.append(np.min(nl.domain))
            maxx_vals.append(np.max(nl.domain))
        minx_vals = np.array(minx_vals)
        maxx_vals = np.array(maxx_vals)

        minx = minx_vals.max()
        maxx = maxx_vals.min()
        avg_x = np.linspace(minx, maxx, 200)
        avg_x = avg_x[1:-2]
        y = np.zeros([len(aggModel.output_nls), len(avg_x)])
        
        if minx < maxx:
            for k,nl in enumerate(aggModel.output_nls):            
                xnl = nl.domain.squeeze()
                ynl = nl.range.squeeze()
                if len(xnl.shape) > 0 and len(ynl.shape) > 0:
                    f = interp1d(xnl, ynl)
                    y[k, :] = f(avg_x)

        agg_nl = OutputNL()
        agg_nl.domain = avg_x
        agg_nl.range = y.mean(axis=0)
        agg_nl.range_std = y.std(axis=0)
        aggModel.output_nl = agg_nl

        self.aggregateModel = aggModel


def gaussian_2d_kernel(size):
    """ Adapted from: http://www.scipy.org/Cookbook/SignalSmooth """
    size = int(size)
    x, y = np.mgrid[-size:size+1, -size:size+1]
    g = np.exp(-(x**2/float(size)+y**2/float(size)))
    return g / g.sum()


class ModelRow(tables.IsDescription):
    
    cell_name = tables.StringCol(100)
    region = tables.StringCol(10)
    preproc_type = tables.StringCol(30)
    is_best = tables.BoolCol()
    
    info_bound_min = tables.Float64Col()
    info_bound_max = tables.Float64Col()
    info_bound_mean = tables.Float64Col()
    
    train_info_min = tables.Float64Col()
    train_info_max = tables.Float64Col()
    train_info_mean = tables.Float64Col()
    
    train_info_min_std = tables.Float64Col()
    train_info_max_std = tables.Float64Col()
    train_info_mean_std = tables.Float64Col()
    
    valid_info_min = tables.Float64Col()
    valid_info_max = tables.Float64Col()
    valid_info_mean = tables.Float64Col()
    
    valid_info_min_std = tables.Float64Col()
    valid_info_max_std = tables.Float64Col()
    valid_info_mean_std = tables.Float64Col()
    
    train_perf = tables.Float64Col()
    train_perf_std = tables.Float64Col()
    
    valid_perf = tables.Float64Col()
    valid_perf_std = tables.Float64Col()
    
    
def create_strf_library(units, output_file):
    
    preprocTypes = ['rawspectrogram', 'lyons', 'surprise']
    mt = 'sepnl_spline'
    tval = 0.75
    
    f = tables.openFile(output_file, mode='w', title='Model Library')
    model_grp = f.createGroup('/', 'models', 'All Models')
    perf_table = f.createTable(model_grp, 'performance', ModelRow, "Model Performance")
    
    for u in units:

        cell_name = u.recsite.old_id
        print 'Processing cell %s...' % cell_name
        perf_row = perf_table.row

        (maxScore, scoreValues, scoreCube, infoTrainCube, infoValidCube, modelCube) = \
                create_scoring_cubes_tg(u, preprocTypes=preprocTypes, modelTypes=[mt], thresholds=[tval])

        region = get_unit_str(u)            
        stimGroups = [stimClassToGroup[cr.presentation.class_name] for cr in u.class_resps]         
        stimIndx = stimGroups.index('Con')
        classResp = u.class_resps[stimIndx]
        
        all_data = {}

        best_preproc_type = None
        best_valid_perf = -np.Inf

        for k,pt in enumerate(preprocTypes):    
            svals = scoreValues[pt][mt][tval]
            if len(svals) > 0:                
                s = Struct()
                output_files = [m.directory for m in modelCube[k][0][0]]
                s.train_info_means = np.array([m.training_performance.info_mean for m in modelCube[k][0][0]])
                s.valid_info_means = np.array([m.validation_performance.info_mean for m in modelCube[k][0][0]])
                s.train_info_uppers = np.array([m.training_performance.info_upper for m in modelCube[k][0][0]])
                s.valid_info_uppers = np.array([m.validation_performance.info_upper for m in modelCube[k][0][0]])
                s.train_info_lowers = np.array([m.training_performance.info_lower for m in modelCube[k][0][0]])
                s.valid_info_lowers = np.array([m.validation_performance.info_lower for m in modelCube[k][0][0]])
                
                s.train_perf_ratios = np.array([m.training_performance.info_mean / m.training_performance.info_bound_mean for m in modelCube[k][0][0]])
                s.valid_perf_ratios = np.array([m.validation_performance.info_mean / m.validation_performance.info_bound_mean for m in modelCube[k][0][0]])

                am = AggregateModel(u, output_files)
                s.agg_model = am.aggregateModel
                all_data[pt] = s
        
                vperf = s.valid_perf_ratios.mean()
                if vperf > best_valid_perf:
                    best_valid_perf = vperf
                    best_preproc_type = pt        
        
        if len(all_data) == 0:
            continue
        
        param_grp = f.createGroup(model_grp, '%s' % cell_name)
        
        for k,pt in enumerate(preprocTypes):
            
            preproc_param_grp = f.createGroup(param_grp, '%s' % pt)
            
            if pt not in all_data:
                continue
            
            s = all_data[pt]
            
            perf_row['cell_name'] = cell_name
            perf_row['region'] = region
            perf_row['preproc_type'] = pt
            perf_row['is_best'] = pt == best_preproc_type
            
            perf_row['info_bound_min'] = classResp.performance.info_lower
            perf_row['info_bound_max'] = classResp.performance.info_upper
            perf_row['info_bound_mean'] = classResp.performance.info_mean
                        
            perf_row['train_info_min'] = s.train_info_lowers.mean()
            perf_row['train_info_max'] = s.train_info_uppers.mean()
            perf_row['train_info_mean'] = s.train_info_means.mean()
            
            perf_row['train_info_min_std'] = s.train_info_lowers.std()
            perf_row['train_info_max_std'] = s.train_info_uppers.std()
            perf_row['train_info_mean_std'] = s.train_info_means.std()
            
            perf_row['valid_info_min'] = s.valid_info_lowers.mean()
            perf_row['valid_info_max'] = s.valid_info_uppers.mean()
            perf_row['valid_info_mean'] = s.valid_info_means.mean()
            
            perf_row['valid_info_min_std'] = s.valid_info_lowers.std()
            perf_row['valid_info_max_std'] = s.valid_info_uppers.std()
            perf_row['valid_info_mean_std'] = s.valid_info_means.std()
                        
            perf_row['train_perf'] = s.train_perf_ratios.mean()
            perf_row['train_perf_std'] = s.train_perf_ratios.std()
            
            perf_row['valid_perf'] = s.valid_perf_ratios.mean()
            perf_row['valid_perf_std'] = s.valid_perf_ratios.std()
            
            perf_row.append()
            
            strf_grp = f.createGroup(preproc_param_grp, 'strf')
            strf_grp._v_attrs.high_freq = s.agg_model.high_freq
            strf_grp._v_attrs.low_freq = s.agg_model.low_freq
            f.createArray(strf_grp, 'freqs', s.agg_model.freqs)
            f.createArray(strf_grp, 'time_lags', s.agg_model.timeLags)
            f.createArray(strf_grp, 'weights', s.agg_model.strf)
            f.createArray(strf_grp, 'std', s.agg_model.strfStd)
            
            nl_grp = f.createGroup(preproc_param_grp, 'output_nl')
            f.createArray(nl_grp, 'domain', s.agg_model.output_nl.domain)
            f.createArray(nl_grp, 'range', s.agg_model.output_nl.range)
            f.createArray(nl_grp, 'range_std', s.agg_model.output_nl.range_std)
                
    f.close()            
    
    


def write_best_models_to_file(perf_list_file, output_file, exclude_models=[]):
    f = open(perf_list_file, 'r')
    cr = csv.reader(f, delimiter=',')
    cr.next()
    
    best_models = {}
    
    for row in cr:        
        unit = row[0]
        reg =  row[1]
        preproc =  row[2]
        model =  row[3]
        threshold =  float(row[4])
        score = float(row[5])
    
        if unit not in best_models:
            best_models[unit] = {'score':0.0}
        
        if score > best_models[unit]['score']:
            if exclude_models is not None and model in exclude_models:
                continue
            best_models[unit]['score'] = score
            best_models[unit]['region'] = reg
            best_models[unit]['preproc'] = preproc
            best_models[unit]['model'] = model
            best_models[unit]['threshold'] = threshold
    f.close()
    
    fout = open(output_file, 'w')
    for unit,up in best_models.iteritems():
        fout.write('%s,%s,%s,%s,%0.2f,%0.6f\n' % (unit, up['region'], up['preproc'], up['model'], up['threshold'], up['score']))
    fout.close()

def get_good_units_and_region():

    ulist = []    
    units = get_good_units()
    for u in units:

        unitRegions = [a.name for a in u.recsite.areas]
        if len(unitRegions) > 1:
            if 'L' in unitRegions:
                region = 'L'
            else:
                region = unitRegions[0]
        else:
            region = unitRegions[0]
        
        ulist.append( (u.recsite.old_id, region))
        
    return ulist
