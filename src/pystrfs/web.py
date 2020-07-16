
import os
import csv
import operator
import numpy as np

from mako.template import Template
from mako.runtime import Context

from pystrfs import stimClassToGroup, stimGroupToClass
from pystrfs import PYSTRFS_DIR, get_good_units
from pystrfs.analysis import Struct, create_scoring_cubes_tg
from pystrfs.plot import ModelPlot, UnitResponse

MAKO_DIR = 'pystrfs/mako'
OUTPUT_DIR = os.path.join(PYSTRFS_DIR, 'web')

def render_index():
    
    good_cells = get_good_units()

    region_list = {'mld':[], 'ov':[], 'l':[], 'cm':[]}
    #render each unit
    for u in good_cells[0:28]:
        uregs = [a.name for a in u.recsite.areas]
        reg = 'l'
        if 'MLd' in uregs:
            reg = 'mld'
        elif 'CM' in uregs:
            reg = 'cm'
        elif 'OV' in uregs:
            reg = 'ov'
        region_list[reg].append(u)

    main_file = os.path.join(MAKO_DIR, 'index.html')
    output_file = os.path.join(OUTPUT_DIR, 'index.html')
    tpl = Template(filename=main_file)

    #render main index
    f = open(output_file, 'w')
    ctx = Context(f, region_list=region_list)
    tpl.render_context(ctx)
    f.close()



def render_unit(unit):
    preproc_types = ['rawspectrogram', 'lyons', 'surprise']
    model_types = ['linear', 'sepnl_spline', 'sepnl_dists', 'poisson', 'binomial', 'leglm']
    thresholds = [0.0, 0.25, 0.5, 0.75, 1.0]
    
    unitRegions = [a.name for a in unit.recsite.areas]
    extraPreprocs=[]
    if 'L' in unitRegions:
        preproc_types.append('neurogram')        
        extraPreprocs.append('neurogram')
    
    (max_score, score_values, score_cube, info_train_cube, info_validCube, model_cube) = create_scoring_cubes_tg(unit, 5, extraPreprocs=extraPreprocs)

    score_data = {}
    
    unit_resp = UnitResponse(unit, 'Con')
    
    upath = os.path.join(OUTPUT_DIR, 'units', unit.recsite.old_id)
    if not os.path.isdir(upath):
        os.mkdir(upath)
        
    ipath = os.path.join(upath, 'images')
    if not os.path.isdir(ipath):
        os.mkdir(ipath)
    
    clist = [cp for cp in unit.class_resps if stimClassToGroup[cp.class_name] == 'Con']
    classResp = clist[0]
    classPerf = Struct()
    classPerf.info_mean = classResp.performance.info_mean
    classPerf.info_lower = classResp.performance.info_lower
    classPerf.info_upper = classResp.performance.info_upper
    
    all_model_perfs = []
    for k,ptype in enumerate(preproc_types):
        for p,mtype in enumerate(model_types):
            for q,thresh in enumerate(thresholds):
                scoreVals = score_values[ptype][mtype][thresh]
                if len(scoreVals) > 0:
                    score = scoreVals.mean()
                    score_std = scoreVals.std()
                    t = (score, score_std, ptype, mtype, thresh)
                    all_model_perfs.append(t)
    all_model_perfs.sort(key=operator.itemgetter(0), reverse=True)
    top10 = []
    ntop = min(10, len(all_model_perfs))
    for k in range(ntop):        
        t = all_model_perfs[k]
        s = Struct()
        s.perf_mean = t[0]
        s.perf_std = t[1]
        s.preproc_type = t[2]
        s.model_type = t[3]
        s.threshold = t[4]
        top10.append(s)
    
    unit_template = os.path.join(MAKO_DIR, 'unit.html')
    unit_file = os.path.join(upath, 'index.html')

    render_models_best(unit, score_values)

    for k,ptype in enumerate(preproc_types):
        render_models_by_preproc(unit, ptype)
    
    perf_marginals_thresh = {}
    perf_marginals_model = {}
    for k,ptype in enumerate(preproc_types):
        score_data[ptype] = {}
        perf_marginals_thresh[ptype] = {}
        perf_marginals_model[ptype] = {}
        for p,mtype in enumerate(model_types):
            score_data[ptype][mtype] = {}
            if mtype not in perf_marginals_thresh[ptype]:
                perf_marginals_thresh[ptype][mtype] = []
            for q,thresh in enumerate(thresholds):
                if thresh not in perf_marginals_model[ptype]:
                    perf_marginals_model[ptype][thresh] = []                           
                if thresh not in score_values[ptype][mtype]:
                    (r, g, b) = (0, 0, 0)
                    score = 0.0
                    score_std = 0.0
                    print 'Missing value for preprocType=%s, modelType=%s' % (ptype, mtype)
                    nscore = -1.0
                else:
                    scoreVals = score_values[ptype][mtype][thresh]
                    if len(scoreVals) > 0:
                        score = scoreVals.mean()
                        score_std = scoreVals.std()
                        perf_marginals_thresh[ptype][mtype].append(score)
                        perf_marginals_model[ptype][thresh].append(score)
                        nscore = score / max_score
                        g = 0
                        r = int(nscore * 255)
                        if nscore == 0.0:
                            b = 255
                        else:
                            b = int((1 - nscore)*255)
                    else:
                        (r, g, b) = (0.0, 0.0, 0.0)
                        score = 0.0
                        score_std = 0.0
                        print 'No scores for preprocType=%s, modelType=%s' % (ptype, mtype)
                        nscore = -1.0                        

                fname = 'model_%s_%s_%0.2f.html' % (ptype, mtype, thresh)
                url = os.path.join(upath, fname)

                ms = Struct()
                ms.rgb = (r, g, b)
                ms.score = score
                ms.score_std = score_std
                ms.preproc_type = ptype
                ms.model_type = mtype
                ms.threshold = thresh
                ms.url = url

                score_data[ptype][mtype][thresh] = ms
                
                if nscore > -1.0:
                    render_model(unit, unit_resp, model_cube[k][p][q], ptype, mtype, thresh)
                
                
    for k,ptype in enumerate(preproc_types):        
        for p,mtype in enumerate(model_types):
            scores = np.array(perf_marginals_thresh[ptype][mtype])
            perf_marginals_thresh[ptype][mtype] = scores.mean()              
            
        for q,thresh in enumerate(thresholds):
            scores = np.array(perf_marginals_model[ptype][thresh])
            perf_marginals_model[ptype][thresh] = scores.mean()
    
    
    tpl = Template(filename=unit_template)
    f = open(unit_file, 'w')
    ctx = Context(f, unit=unit, unit_regions=unitRegions, score_data=score_data, preproc_types=preproc_types,\
                  model_types=model_types, thresholds=thresholds, class_perf=classPerf, top10=top10, \
                  perf_marginals_model=perf_marginals_model, perf_marginals_thresh=perf_marginals_thresh)
    tpl.render_context(ctx)
    f.close()

def render_model(unit, unit_response, models, preproc_type, model_type, threshold):
    
    #construct output filename template    
    ofile_template = models[0].directory[0:-4] + '%d.h5'
    (root, fname) = os.path.split(ofile_template)
    print 'Rendering template %s...' % fname
    
    #get model plotter
    mplot = ModelPlot(unit, ofile_template)
    if mplot.failed:
        print '\tModel Failure!!!'
        return
    
    upath = os.path.join(OUTPUT_DIR, 'units', unit.recsite.old_id)
    ipath = os.path.join(upath, 'images')
    
    #generate images
    prefix = '%s_%s_%0.2f'  % (preproc_type, model_type, threshold)
    
    all_data = Struct()
    all_data.preproc_type = preproc_type
    all_data.model_type = model_type
    all_data.threshold = threshold
    
    model_data = []
    for k,m in enumerate(models):
        
        mdata = Struct()        
        
        ofp = os.path.join(ipath, '%s_strf_%d' % (prefix, k))        
        mplot.plot_model_strf(k, ofp)
        mdata.strf_prefix = ofp
        
        ofp = os.path.join(ipath, '%s_nl_%d' % (prefix, k))        
        mplot.plot_model_nl(k, ofp)
        mdata.nl_prefix = ofp
        
        vgrps = mplot.modelData[k].validationGroups        
        vplots = {}
        for vg in vgrps:            
            ofp = os.path.join(ipath, '%s_response_%d_%d' % (prefix, k, vg))
            mplot.plot_response(k, vg, ofp, unit_response=unit_response)
            vplots[vg] = ofp
        mdata.vplots = vplots
        
        try:
            trat = m.training_performance.info_mean / m.training_performance.info_bound_mean
            vrat = m.validation_performance.info_mean / m.validation_performance.info_bound_mean
            mdata.training_perf = trat
            mdata.validation_perf = vrat
        except:
            print 'Problem with model perf: unit=%s, preproc=%s, model=%s, thresh=%0.2f' % \
                  (unit.recsite.old_id, preproc_type, model_type, threshold)
            mdata.training_perf = 0
            mdata.validation_perf = 0
        
        model_data.append(mdata)
    
    tperfs = np.array([m.training_perf for m in model_data])
    vperfs = np.array([m.validation_perf for m in model_data])
    
    all_data.training_performance = (tperfs.mean(), tperfs.std())
    all_data.validation_performance = (vperfs.mean(), vperfs.std())
    all_data.models = model_data
    all_data.unit = unit
    
    ofp = os.path.join(ipath, '%s_strf_mean' % prefix)
    mplot.plot_agg_strf(ofp)
    all_data.strf_prefix = ofp
    
    ofp = os.path.join(ipath, '%s_nl_mean' % prefix)
    mplot.plot_agg_nl(ofp)
    all_data.nl_prefix = ofp
    
    ofp = os.path.join(ipath, '%s_pvals_mean' % prefix)
    mplot.plot_agg_pvals(ofp)
    all_data.pvals_prefix = ofp
    
    model_template = os.path.join(MAKO_DIR, 'model.html')
    fname = 'model_%s_%s_%0.2f.html' % (preproc_type, model_type, threshold)
    model_file = os.path.join(upath, fname)    
    tpl = Template(filename=model_template)
    f = open(model_file, 'w')
    ctx = Context(f, all_data=all_data)
    tpl.render_context(ctx)
    f.close()

def render_models_by_preproc(unit, preproc_type):
    
    upath = os.path.join(OUTPUT_DIR, 'units', unit.recsite.old_id)
    ipath = os.path.join(upath, 'images')
    
    model_types = ['linear', 'sepnl_dists', 'sepnl_spline', 'poisson', 'binomial', 'leglm']
    thresholds = [0.0, 0.25, 0.5, 0.75, 1.0]
    
    all_data = Struct()
    all_data.unit = unit
    all_data.preproc_type = preproc_type
    all_data.model_types = model_types
    all_data.model_types_strfs = ['linear', 'poisson', 'binomial', 'leglm']
    all_data.model_types_nls = model_types
    all_data.thresholds = thresholds
    
    img = {}
    for thresh in thresholds:
        img[thresh] = {}
        for mtype in model_types:
            s = Struct()
            prefix = '%s_%s_%0.2f'  % (preproc_type, mtype, thresh)
            s.strf_prefix = os.path.join(ipath, '%s_strf_mean' % prefix)
            s.nl_prefix = os.path.join(ipath, '%s_nl_mean' % prefix)
            img[thresh][mtype] = s

    all_data.img = img
    
    model_template = os.path.join(MAKO_DIR, 'models_by_preproc.html')
    fname = 'models_by_preproc_%s.html' % preproc_type
    model_file = os.path.join(upath, fname)    
    tpl = Template(filename=model_template)
    f = open(model_file, 'w')
    ctx = Context(f, all_data=all_data)
    tpl.render_context(ctx)
    f.close()
    
def render_models_best(unit, score_values):
    
    upath = os.path.join(OUTPUT_DIR, 'units', unit.recsite.old_id)
    ipath = os.path.join(upath, 'images')
    
    preproc_types = ['rawspectrogram', 'lyons', 'surprise']
    model_types = ['linear', 'sepnl_dists', 'sepnl_spline', 'poisson', 'binomial', 'leglm']
    thresholds = [0.0, 0.25, 0.5, 0.75, 1.0]
    
    all_data = Struct()
    all_data.unit = unit
    all_data.preproc_types = preproc_types
    all_data.model_types = model_types
    all_data.model_types_strfs = ['linear', 'poisson', 'binomial', 'leglm']
    all_data.model_types_nls = model_types
    
    best_models = {}
    for ptype in preproc_types:
        best_models[ptype] = {}
        for mtype in model_types:
            tvals = score_values[ptype][mtype].keys()
            svt = score_values[ptype][mtype].values()
            scores_per_thresh = [s.mean() for s in svt]
            scores_per_thresh = np.array(scores_per_thresh)
            bindx = np.argmax(scores_per_thresh)
            best_thresh = tvals[bindx]
            prefix = '%s_%s_%0.2f'  % (ptype, mtype, best_thresh)
            s = Struct()
            s.strf_prefix = os.path.join(ipath, '%s_strf_mean' % prefix)
            s.nl_prefix = os.path.join(ipath, '%s_nl_mean' % prefix)
            best_models[ptype][mtype] = s

    all_data.img = best_models

    model_template = os.path.join(MAKO_DIR, 'models_best.html')
    fname = 'models_best.html'
    model_file = os.path.join(upath, fname)    
    tpl = Template(filename=model_template)
    f = open(model_file, 'w')
    ctx = Context(f, all_data=all_data)
    tpl.render_context(ctx)
    f.close()

def render_regional_summary(perf_file):
    
    regions = ['MLd', 'OV', 'L', 'CM']
    preprocTypes = ['rawspectrogram', 'lyons', 'surprise']
    modelTypes = ['linear', 'sepnl_spline', 'poisson', 'binomial', 'leglm']
    thresholds = [0.0, 0.25, 0.5, 0.75, 1.0]

    data_per_region = {}
    
    f = open(perf_file, 'r')
    cr = csv.reader(f, delimiter=',')
    cr.next()
        
    for row in cr:        
        
        s = Struct()
        
        s.unit = row[0]
        s.reg =  row[1]
        s.preproc =  row[2]
        s.model =  row[3]
        s.threshold =  float(row[4])
        s.score = float(row[5])
        s.infoBound = float(row[6])
        
        if s.reg not in data_per_region:
            data_per_region[s.reg] = {}
        if s.preproc not in data_per_region[s.reg]:
            data_per_region[s.reg][s.preproc] = {}
        if s.model not in data_per_region[s.reg][s.preproc]:
            data_per_region[s.reg][s.preproc][s.model] = {}
        if s.unit not in data_per_region[s.reg][s.preproc][s.model]:
            data_per_region[s.reg][s.preproc][s.model][s.unit] = []
        
        #print '%s,%s,%s,%0.2f,%0.2f' % (s.unit, s.preproc, s.model, s.threshold, s.score)
        
        data_per_region[s.reg][s.preproc][s.model][s.unit].append(s)                  

    f.close()
    
    
    score_data = {}
    
    for r,reg in enumerate(regions):
        score_data[reg] = {}
        for k,pt in enumerate(preprocTypes):
            score_data[reg][pt] = {}
            for p,mt in enumerate(modelTypes):
                score_data[reg][pt][mt] = []
                for unit,models in data_per_region[reg][pt][mt].iteritems():
                    scores = np.array([s.score for s in models])
                    best_score = scores.max()
                    print '%s,%s,%s,%s: %0.2f' % (reg, pt, mt, unit, best_score)
                    score_data[reg][pt][mt].append(best_score)
                    
    print score_data
   
    
    smax = {}
    smin = {}
    for reg in regions:
        smax[reg] = 0.0
        smin[reg] = 200.0
    for r,reg in enumerate(regions):    
        for k,pt in enumerate(preprocTypes):    
            for p,mt in enumerate(modelTypes):
                scores = np.array(score_data[reg][pt][mt])
                print scores
                s = Struct()
                smean = scores.mean()
                s.score_mean = smean
                s.score_std = scores.std()
                print '(%s,%s,%s) score=%0.2f +/- %0.2f' % (reg, pt, mt, s.score_mean, s.score_std)
                
                if smean > smax[reg]:
                    smax[reg] = smean
                if smean > 0.0 and smean < smin[reg]:
                    smin[reg] = smean
                
                score_data[reg][pt][mt] = s
    
    for r,reg in enumerate(regions):    
        for k,pt in enumerate(preprocTypes):    
            for p,mt in enumerate(modelTypes):
                s = score_data[reg][pt][mt]
                score = s.score_mean
                ns = (score - smin[reg]) / (smax[reg]-smin[reg])
                r = int(ns*255)
                g = 0
                b = int((1.0 - ns)*255)
                s.rgb = (r, g, b)

    model_template = os.path.join(MAKO_DIR, 'regional_summary.html')    
    regsum_file = os.path.join(PYSTRFS_DIR, 'web', 'regional_summary.html')    
    tpl = Template(filename=model_template)
    f = open(regsum_file, 'w')
    ctx = Context(f, score_data=score_data, preprocTypes=preprocTypes, modelTypes=modelTypes, regions=regions)
    tpl.render_context(ctx)
    f.close()
    
    