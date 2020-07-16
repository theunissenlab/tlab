import operator
import shutil
import numpy as np

from browser.model import meta

from browser.model.stims.protocol import Protocol, ProtocolStim
from browser.model.physiology import Recsite, ElectrodePosition, Block
from browser.model.physiology.extracellular import Unit

from pystrfs import *
from pystrfs.preprocess import *
from pystrfs.slurmbot import *
from pystrfs.matlab import MatlabTemplate,MatlabProcess


class ModelRun:

    def __init__(self):
        self.unitFile = None
        self.unitName = None
        self.templateFile = None
        self.trainingGroups = []
        self.validationGroups = []
        self.groupNumber = None
        self.stimClass = None
        self.preprocName = None
        self.preprocFile = None
        self.stimTransforms = []
        self.outputDir = None
        self.numCpus = 2

    def set_unit(self, unitFile, unitName):
        self.unitFile = unitFile
        self.unitName = unitName
        self.unitDir = os.path.join(PYSTRFS_DIR, 'units', unitName)
        self.outputDir = os.path.join(self.unitDir, 'output')


    def set_preproc(self, stimClass, preprocName, preprocFile):

        pstrs = preprocName.split('.')
        preprocShort = pstrs[0]

        self.stimClass = stimClass
        self.preprocName = preprocName
        self.preprocFile = preprocFile

        if preprocShort == 'stft':
            self.stimTransforms = ['log', 'zscore']
        elif preprocShort == 'lyons':
            self.stimTransforms = ['zscore']
        else:
            self.stimTransforms = []

    def set_data_groups(self, trainingGroups, validationGroups, earlyStopping, groupNum):
        self.trainingGroups = trainingGroups
        self.validationGroups = validationGroups
        self.earlyStoppingGroups = earlyStopping
        self.groupNumber = groupNum

    def get_group_str(self, grps):
        gstr = ('%s' % grps.ravel()).strip().replace('\n', '')
        return gstr

    def get_template_params(self):
        return None

    def get_file_desc(self):
        return None

    def get_output_file(self):
        return None

    def get_matlab_run(self):
        mt = MatlabTemplate()
        mt.template_from_file(self.templateFile)
        tparams = self.get_template_params()
        mp = mt.to_matlab_process(tparams)
        return mp


class DirectFitRun(ModelRun):

    def __init__(self):
        ModelRun.__init__(self)
        self.numCpus = 2
        self.templateFile = os.path.join(TEMPLATE_DIR, 'direct_fit+outputnl.m')

    def get_template_params(self):

        tgStr = self.get_group_str(self.trainingGroups)
        vgStr = self.get_group_str(self.validationGroups)

        rfiles = self.get_output_file()
        if self.modelType == 'linear':
            responseFile = rfiles[0]
            distsResponseFile = rfiles[1]
            splineResponseFile = rfiles[2]
        else:
            responseFile = rfiles
            distsResponseFile = '';
            splineResponseFile = '';

        tstrs = ['\'%s\'' % x for x in self.stimTransforms]

        tparams = {'SRC_ROOT_DIR':SRC_ROOT_DIR, 'UNIT_FILE':self.unitFile, 'PREPROC_FILE':self.preprocFile,
                   'STIM_CLASS':self.stimClass, 'TRAINING_GROUPS':tgStr, 'VALIDATION_GROUPS':vgStr,
                   'DF_RESPONSE_FILE':responseFile, 'DF_DISTS_RESPONSE_FILE':distsResponseFile,
                   'DF_SPLINE_RESPONSE_FILE':splineResponseFile, 'TRANSFORMS':','.join(tstrs)}

        return tparams

    def get_file_desc(self):
        fileDesc = 'direct_fit.%s.%s.%s.%d' % (self.unitName, self.preprocName, self.stimClass, self.groupNumber)
        return fileDesc

    def get_output_file(self):
        responseFileName = 'directfit.%s.%s.%d.h5' % (self.stimClass, self.preprocName, self.groupNumber)
        responseFile = os.path.join(self.outputDir, responseFileName)

        distsResponseFileName = 'directfit.nl_dists.%s.%s.%d.h5' % (self.stimClass, self.preprocName, self.groupNumber)
        distsResponseFile = os.path.join(self.outputDir, distsResponseFileName)

        splineResponseFileName = 'directfit.nl_spline.%s.%s.%d.h5' % (self.stimClass, self.preprocName, self.groupNumber)
        splineResponseFile = os.path.join(self.outputDir, splineResponseFileName)

        return [responseFile, distsResponseFile, splineResponseFile]


class ThreshGradDebugRun(ModelRun):

    def __init__(self, threshold):
        ModelRun.__init__(self)
        self.numCpus = 4
        self.threshold = threshold
        self.templateFile = os.path.join(TEMPLATE_DIR, 'threshgrad_debug.m')

    def get_output_file(self):
        responseFileName = 'threshgrad.lin.debug.%s.%s.thresh_%0.2f.%d.mat' % \
                          (self.stimClass, self.preprocName, self.threshold, self.groupNumber)
        responseFile = os.path.join(self.outputDir, responseFileName)

        return responseFile


    def get_template_params(self):

        tgStr = self.get_group_str(self.trainingGroups)
        vgStr = self.get_group_str(self.validationGroups)

        responseFile = self.get_output_file()

        tstrs = ['\'%s\'' % x for x in self.stimTransforms]

        tparams = {'SRC_ROOT_DIR':SRC_ROOT_DIR, 'UNIT_FILE':self.unitFile, 'PREPROC_FILE':self.preprocFile,
                   'STIM_CLASS':self.stimClass, 'TRAINING_GROUPS':tgStr, 'VALIDATION_GROUPS':vgStr,
                   'RESPONSE_FILE':responseFile, 'TRANSFORMS':','.join(tstrs), 'THRESHOLD': '%0.2f' % self.threshold}

        return tparams

    def get_file_desc(self):
        fileDesc = 'threshgrad.debug.%s.%s.%s.%d' % \
                    (self.unitName, self.preprocName, self.stimClass, self.groupNumber)
        return fileDesc


class ThreshGradRun(ModelRun):

    def __init__(self, threshold, modelType):
        ModelRun.__init__(self)
        self.numCpus = 4
        self.version = 3
        self.threshold = threshold
        self.modelType = modelType
        self.computeOutputNL = True
        self.templateFile = os.path.join(TEMPLATE_DIR, 'threshgrad+outputnl.m')

        if self.modelType != 'linear':
            self.computeOutputNL = False


    def get_output_file(self):
        vstr = 'v%d.' % self.version
        if self.version == 1:
            vstr = ''

        responseFileName = 'threshgrad.%s%s.%s.%s.thresh_%0.2f.%d.h5' % \
                             (vstr, self.modelType, self.stimClass, self.preprocName, self.threshold, self.groupNumber)
        responseFile = os.path.join(self.outputDir, responseFileName)

        if self.computeOutputNL:
            distsResponseFileName = 'threshgrad.%s%s.nl_dists.%s.%s.thresh_%0.2f.%d.h5' % \
                                    (vstr, self.modelType, self.stimClass, self.preprocName, self.threshold, self.groupNumber)
            distsResponseFile = os.path.join(self.outputDir, distsResponseFileName)

            splineResponseFileName = 'threshgrad.%s%s.nl_spline.%s.%s.thresh_%0.2f.%d.h5' % \
                                    (vstr, self.modelType, self.stimClass, self.preprocName, self.threshold, self.groupNumber)
            splineResponseFile = os.path.join(self.outputDir, splineResponseFileName)

            return [responseFile, distsResponseFile, splineResponseFile]
        else:
            return [responseFile]


    def get_template_params(self):

        tgStr = self.get_group_str(self.trainingGroups)
        vgStr = self.get_group_str(self.validationGroups)
        egStr = self.get_group_str(self.earlyStoppingGroups)

        rfiles = self.get_output_file()
        responseFile = rfiles[0]

        if self.modelType == 'linear':
            distsResponseFile = rfiles[1]
            splineResponseFile = rfiles[2]
        else:
            distsResponseFile = '';
            splineResponseFile = '';


        compnl = 0
        if self.computeOutputNL:
            compnl = 1

        tstrs = ['\'%s\'' % x for x in self.stimTransforms]

        tparams = {'SRC_ROOT_DIR':SRC_ROOT_DIR, 'UNIT_FILE':self.unitFile, 'PREPROC_FILE':self.preprocFile,
                   'STIM_CLASS':self.stimClass, 'TRAINING_GROUPS':tgStr, 'VALIDATION_GROUPS':vgStr, 'EARLY_STOPPING_GROUPS':egStr,
                   'TG_RESPONSE_FILE':responseFile, 'TG_DISTS_RESPONSE_FILE':distsResponseFile,
                   'TG_SPLINE_RESPONSE_FILE':splineResponseFile, 'TRANSFORMS':','.join(tstrs), 'THRESHOLD': '%0.2f' % self.threshold,
                   'MODEL_TYPE':self.modelType, 'COMPUTE_OUTPUTNL':compnl}

        return tparams

    def get_file_desc(self):
        vstr = 'v%d.' % self.version
        if self.version == 1:
            vstr = ''
        fileDesc = 'threshgrad.%s%s.%s.%s.%s.thresh_%0.2f.%d' % \
                    (vstr, self.modelType, self.unitName, self.preprocName, self.stimClass, self.threshold, self.groupNumber)
        return fileDesc


class ModelRunner:

    def __init__(self):
        self.slurmBot = SlurmBot()

    def create_directfit_runs(self, unit):

        mruns = []

        stimClasses = ['Con']
        unitName =unit.recsite.old_id
        unitFile = map_dcp_unit_to_fs(unit)

        allGroups = get_cv_groups_for_unit(unit, 10)
        allPreprocFiles = get_all_preproc_files(stimClasses)

        for sc in stimClasses:

            cvGrps = allGroups[sc]

            for k,grps in enumerate(cvGrps):

                trainingGroups = grps[0]
                validationGroups = grps[1]

                for ptup in allPreprocFiles[sc]:
                    preprocName = ptup[0]
                    preprocFile = ptup[1]

                    dfRun = DirectFitRun()
                    dfRun.set_unit(unitFile, unitName)
                    dfRun.set_preproc(sc, preprocName, preprocFile)
                    dfRun.set_data_groups(trainingGroups, validationGroups, k)
                    mruns.append(dfRun)

        print 'Created %d runs for unit %s' % (len(mruns), unit.recsite.old_id)
        return mruns

    def create_threshgrad_runs(self, unit, modelType, version = 5, dataGroupsByStim = None, allPreprocFiles=None, thresholds=None):

        mruns = []

        if thresholds is None:
            thresholds = [0.0, 0.25, 0.5, 0.75, 1.0]

        stimClasses = ['Con']
        unitName = unit.recsite.old_id
        unitFile = map_dcp_unit_to_fs(unit)

        if allPreprocFiles is None:
            allPreprocFiles = get_all_preproc_files(stimClasses)
                        
        if dataGroupsByStim is None:
            dataGroupsByStim = get_data_groups(unit, True)

        for thresh in thresholds:
            for sc in stimClasses:

                dg = dataGroupsByStim[sc]

                for ptup in allPreprocFiles[sc]:
                    preprocName = ptup[0]
                    preprocFile = ptup[1]

                    for k,glist in enumerate(dg):
                        tg = glist[0]
                        vg = glist[1]
                        eg = glist[2]

                        tgRun = ThreshGradRun(thresh, modelType)
                        tgRun.set_unit(unitFile, unitName)
                        tgRun.set_preproc(sc, preprocName, preprocFile)
                        tgRun.set_data_groups(tg, vg, eg, k)
                        tgRun.version = version
                        mruns.append(tgRun)

        print 'Created %d runs for unit %s' % (len(mruns), unit.recsite.old_id)
        return mruns

    def create_threshgrad_runs_min(self, unit, version = 5):
        modelTypes = ['linear']
        mruns = []
        dataGroupsByStim = get_data_groups(unit, True)
        for modelType in modelTypes:
            mruns.extend(self.create_threshgrad_runs(unit, modelType, version, dataGroupsByStim, thresholds=[0.75]))
        return mruns


    def create_threshgrad_runs_all(self, unit, version = 5):
        modelTypes = ['linear', 'poisson', 'binomial', 'leglm']
        mruns = []
        dataGroupsByStim = get_data_groups(unit, True)
        for modelType in modelTypes:
            mruns.extend(self.create_threshgrad_runs(unit, modelType, version, dataGroupsByStim))
        return mruns
    
    def create_threshgrad_runs_all_neurogram(self, unit, version = 5):
        modelTypes = ['linear', 'leglm']
        mruns = []
        dataGroupsByStim = get_data_groups(unit, True)
        preprocFiles = {'Con': [('neurogram', os.path.join(PREPROC_DIR, 'neurogram.h5'))]}
        for modelType in modelTypes:
            mruns.extend(self.create_threshgrad_runs(unit, modelType, version, dataGroupsByStim,
                                                     allPreprocFiles=preprocFiles, thresholds=[0.50, 1.00]))
        return mruns


    def add_directfit(self, unit):

        mruns = self.create_directfit_runs(unit)
        self.add_models(mruns)

    def add_models(self, modelRuns):

        for mr in modelRuns:
            mp = mr.get_matlab_run()
            fileDesc = mr.get_file_desc()
            sout = os.path.join(SLURM_DIR, '%s.%%j.out.txt' % fileDesc)
            serr = sout
            scriptName = os.path.join(SLURM_DIR, '%s.sbatch' % fileDesc)

            sparams = {'out': sout,
                       'err': serr,
                       'partition':'all',
                       'cpus': mr.numCpus,
                       'script': scriptName}

            self.slurmBot.add(mp.get_commands(), sparams)

    def run_models(self):
        print 'Running %d models...' % len(self.slurmBot.get_queued_jobs())
        self.slurmBot.run()


    def get_failed_threshgrad_runs(self, xunit, version = 5):
        mruns = self.create_threshgrad_runs_all(unit, version)
        failedRuns = []
        for m in mruns:
            ofiles = m.get_output_file()
            runFailed = False
            for ofile in ofiles:
                if not os.path.exists(ofile):
                    print 'Missing file: %s' % ofile
                    runFailed = True
            failedRuns.append(m)

        return failedRuns
