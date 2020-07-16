from browser.model import meta

from browser.model.stims.protocol import Protocol, ProtocolStim
from browser.model.physiology import Recsite, ElectrodePosition, Block
from browser.model.physiology.extracellular import Unit

from browser.model.analysis.model import Preprocesser
from browser.model.analysis.strflab.spectrogram import RawSpectrogramPreprocessor
from browser.model.analysis.strflab.lyons import LyonsPreprocessor
from browser.model.analysis.strflab.surprise import SurprisePreprocessor
from browser.model.analysis.strflab.neurogram import NeurogramPreprocessor

from pystrfs import *
from pystrfs.slurmbot import *
from pystrfs.matlab import MatlabTemplate,MatlabProcess


def run_coherence_bound_all(individual=False):

    slurmBot = SlurmBot()
    umap = get_all_dcp_units()
    for uId, u in umap.iteritems():
        run_unit_coherence_bound(u, slurmBot, individual)

    slurmBot.run()

def run_unit_coherence_bound(unit, slurmBot = None, individual=False):

    if not individual:
        fdesc = 'coherence_bound.%s' % unit.recsite.old_id
        templateFile = os.path.join(TEMPLATE_DIR, 'coherence_bound.m')
    else:
        fdesc = 'coherence_bound_individual.%s' % unit.recsite.old_id
        templateFile = os.path.join(TEMPLATE_DIR, 'coherence_bound_individual.m')

    unitFile = map_dcp_unit_to_fs(unit)
    mt = MatlabTemplate()
    mt.template_from_file(templateFile)
    tparams = {'SRC_ROOT_DIR':SRC_ROOT_DIR, 'UNIT_FILE':unitFile}
    mp = mt.to_matlab_process(tparams)

    if slurmBot is None:
        mp.run()
    else:
        sout = os.path.join(SLURM_DIR, '%s.%%j.out.txt' % fdesc)
        serr = sout
        scriptName = os.path.join(SLURM_DIR, '%s.sbatch' % fdesc)

        sparams = {'out': sout,
                   'err': serr,
                   'partition':'all',
                   'cpus': 1,
                   'script': scriptName}

        slurmBot.add(mp.get_commands(), sparams)


def get_preprocessor_from_filename(filePath):

    #extract parameters from file
    (rootDir, fileName) = os.path.split(filePath)
    nparts = fileName.split('.')
    preprocParams = {}
    preprocType = nparts[0]
    for npart in nparts[1:]:
        if npart != 'h5':
            pprts = npart.split('_')
            if len(pprts) > 1:
                preprocParams[pprts[0]] = pprts[1]
            else:
                preprocParams['class'] = npart

    #find preproc type in DB
    preproc = None
    
    if preprocType == 'stft':
        pq = meta.Session.query(RawSpectrogramPreprocessor).filter_by(output_file=filePath)
        if pq.count() == 0:
            preproc = RawSpectrogramPreprocessor()
            preproc.zscore = 1
            preproc.nstd = int(preprocParams['nstd'])
            preproc.fband = int(preprocParams['fband'])
            preproc.low_freq = 375.0
            preproc.high_freq = 8000.0
            preproc.output_file = filePath
            preproc.transforms = 'log,zscore'
            meta.Session.add(preproc)
            meta.Session.commit()
        else:
            preproc = pq.one()

    elif preprocType == 'lyons':
        pq = meta.Session.query(LyonsPreprocessor).filter_by(output_file=filePath)
        if pq.count() == 0:
            #kludge...
            stepStr = '.'.join(nparts[-3:-1])
            stepPrts = stepStr.split('_')

            preproc = LyonsPreprocessor()
            preproc.zscore = 1
            preproc.earq = int(preprocParams['earQ'])
            preproc.agc = int(preprocParams['agc'])
            preproc.differ = int(preprocParams['agc'])
            preproc.step = float(stepPrts[1])
            preproc.tau = 3.0
            preproc.low_freq = 100.0
            preproc.high_freq = 15000.0
            preproc.output_file = filePath
            preproc.transforms = 'zscore'
            meta.Session.add(preproc)
            meta.Session.commit()
        else:
            preproc = pq.one()

    elif preprocType == 'surprise':
        pq = meta.Session.query(SurprisePreprocessor).filter_by(output_file=filePath)
        if pq.count() == 0:
            preproc = SurprisePreprocessor()
            preproc.zscore = 0
            preproc.freq_width = preprocParams['dfw']
            preproc.time_width = preprocParams['dtw']
            preproc.gap = preprocParams['dg']
            preproc.stim_class = preprocParams['class']
            preproc.low_freq = 375.0
            preproc.high_freq = 8000.0
            preproc.output_file = filePath
            preproc.transforms = ''
            meta.Session.add(preproc)
            meta.Session.commit()
        else:
            preproc = pq.one()
            
    elif preprocType == 'neurogram':
        pq = meta.Session.query(NeurogramPreprocessor).filter_by(output_file=filePath)
        if pq.count() == 0:
            preproc = NeurogramPreprocessor()
            preproc.zscore = 0
            preproc.output_file = filePath            
            meta.Session.add(preproc)
            meta.Session.commit()
        else:
            preproc = pq.one()
        
    return preproc



def get_all_preproc_files(stimClasses, lyonsParams = None):

    allFiles = {}

    stftFiles = get_stft_preproc_files(stimClasses)
    lyonsFiles = get_lyons_preproc_files(stimClasses, lyonsParams)
    surpriseFiles = get_surprise_preproc_files(stimClasses)

    for sc in stimClasses:
        allFiles[sc] = []

        for pfile in stftFiles[sc]:
            prts = os.path.split(pfile)
            fprts = os.path.splitext(prts[1])
            preprocName = fprts[0]
            allFiles[sc].append( (preprocName, pfile) )
        for pfile in lyonsFiles[sc]:
            prts = os.path.split(pfile)
            fprts = os.path.splitext(prts[1])
            preprocName = fprts[0]
            allFiles[sc].append( (preprocName, pfile) )
        for pfile in surpriseFiles[sc]:
            prts = os.path.split(pfile)
            fprts = os.path.splitext(prts[1])
            preprocName = fprts[0]
            allFiles[sc].append( (preprocName, pfile) )

    return allFiles


def get_stft_preproc_files(stimClasses):

    pf = {}
    for sc in stimClasses:
        pf[sc] = []
        pf[sc].append(os.path.join(PREPROC_DIR, 'stft.nstd_6.fband_125.h5'))
    return pf


def get_lyons_preproc_files(stimClasses, params):

    pf = {}
    for sc in stimClasses:
        pf[sc] = []

    if params is not None:
        steps = params['step']
        earQs = params['earq']
    else:
        #steps = [0.5, 0.25]  #after looking at 12 cells, step=0.5, earq=8 worked the best
        #earQs = [4, 6, 8]
        steps = [0.50]
        earQs = [8]

    for step in steps:
        for earQ in earQs:
            fname = 'lyons.agc_1.earQ_%d.step_%0.2f.h5' % (earQ, step)
            for sc in stimClasses:
                ptup = os.path.join(PREPROC_DIR, fname)
                pf[sc].append(ptup)
    return pf


def get_surprise_preproc_files(stimClasses):
    pf = {}
    for sc in stimClasses:
        pf[sc] = []
        pf[sc].append(os.path.join(PREPROC_DIR, 'surprise.dfw_3.dtw_3.dg_4.%s.h5' % sc))
    return pf


def write_crcns_stimfile(stimdataFile):

    if not os.path.exists(stimdataFile):
        pq = meta.Session.query(Protocol.id).filter(Protocol.name.like('%CRCNS%'))
        plist = [pid[0] for pid in pq]
        psq = meta.Session.query(ProtocolStim).filter(ProtocolStim.protocol_id.in_(plist))

        uniqueFiles = {}
        for pstim in psq:
            md5 = pstim.file.MD5
            filePath = pstim.file.name
            if md5 not in uniqueFiles:
                uniqueFiles[md5] = filePath

        lns = [','.join([md5,filePath]) for md5,filePath in uniqueFiles.iteritems()]
        outStr = '\n'.join(lns)

        f = open(stimdataFile, 'w')
        f.write(outStr)
        f.close()


def preprocess_crcns_spectrograms():

    stimdataFile = os.path.join(PYSTRFS_DIR, 'scratch', 'spectrogram_files.csv')
    write_crcns_stimfile(stimdataFile)

    mt = MatlabTemplate()
    #templatePath = os.path.join(TEMPLATE_DIR, 'preprocess_spectrogram.m')
    templatePath = 'pystrfs/templates/preprocess_spectrogram.m'
    mt.template_from_file(templatePath)

    libraryFile = os.path.join(PYSTRFS_DIR, 'preproc', 'stft.nstd_6.fband_125.h5')
    params = {'STIMDATA_FILE':stimdataFile, 'OUTPUT_FILE':libraryFile, 'SRC_ROOT_DIR':SRC_ROOT_DIR}

    mp = mt.to_matlab_process(params)
    mp.run()

def preprocess_crcns_lyons():

    stimdataFile = os.path.join(PYSTRFS_DIR, 'scratch', 'lyons_files.csv')
    write_crcns_stimfile(stimdataFile)

    mt = MatlabTemplate()

    templatePath = 'pystrfs/templates/preprocess_lyons.m'
    mt.template_from_file(templatePath)

    steps = [0.5, 0.25]
    earQs = [4, 6, 8]
    agcs = [0, 1]

    for step in steps:
        for earQ in earQs:
            for agc in agcs:

                fname = 'lyons.agc_%d.earQ_%d.step_%0.2f.h5' % (agc, earQ, step)
                libraryFile = os.path.join(PYSTRFS_DIR, 'preproc', fname)
                params = {'STIMDATA_FILE':stimdataFile, 'OUTPUT_FILE':libraryFile, 'SRC_ROOT_DIR':SRC_ROOT_DIR,
                          'LYONS_EARQ':earQ, 'LYONS_AGC':agc, 'LYONS_STEP':step}

                mp = mt.to_matlab_process(params)
                mp.run()


def write_crcns_stimfile_surprise(fileNameTemplate):

    stimFiles = []
    stimClasses = []

    pq = meta.Session.query(Protocol.id).filter(Protocol.name.like('%CRCNS%'))
    plist = [pid[0] for pid in pq]
    psq = meta.Session.query(ProtocolStim).filter(ProtocolStim.protocol_id.in_(plist))

    stimsByClass = {}
    for pstim in psq:
        stimClass = pstim.group.name
        if stimClass not in stimsByClass:
            stimsByClass[stimClass] = []
        stimsByClass[stimClass].append(pstim)

    for stimClass,classStims in stimsByClass.iteritems():
        uniqueFiles = {}
        for pstim in classStims:
            md5 = pstim.file.MD5
            filePath = pstim.file.name
            if md5 not in uniqueFiles:
                uniqueFiles[md5] = filePath

        stimFilesName = fileNameTemplate % stimClass
        stimFilesPath = os.path.join(PYSTRFS_DIR, 'scratch', stimFilesName)

        lns = [','.join([md5,filePath]) for md5,filePath in uniqueFiles.iteritems()]
        outStr = '\n'.join(lns)

        if not os.path.exists(stimFilesPath):
            f = open(stimFilesPath, 'w')
            f.write(outStr)
            f.close()

            print 'Wrote file %s' % stimFilesPath
        stimFiles.append(stimFilesPath)
        stimClasses.append(stimClass)

    return zip(stimClasses, stimFiles)


def preproc_crcns_surprise():

    fileNameTemplate = 'surprise_files.%s.csv'
    classData = write_crcns_stimfile_surprise(fileNameTemplate)

    slurmBot = SlurmBot()

    mt = MatlabTemplate()

    templatePath = 'pystrfs/templates/preprocess_surprise.m'
    mt.template_from_file(templatePath)

    domainFrequencyWidth = 3
    domainTimeWidth = 3
    domainGap = 4

    preprocessFile = os.path.join(PYSTRFS_DIR, 'preproc', 'stft.nstd_6.fband_125.h5')

    for cdata in classData:
        stimClass = cdata[0]
        stimFile = cdata[1]
        fileDesc = 'surprise.dfw_%d.dtw_%d.dg_%d.%s' % (domainFrequencyWidth, domainTimeWidth, domainGap, stimClass)
        outFileName = '%s.h5' % fileDesc
        libraryFile = os.path.join(PYSTRFS_DIR, 'preproc', outFileName)
        params = {'STIMDATA_FILE':stimFile, 'OUTPUT_FILE':libraryFile, 'SRC_ROOT_DIR':SRC_ROOT_DIR,
                  'DOMAIN_FREQUENCY_WIDTH':domainFrequencyWidth, 'DOMAIN_TIME_WIDTH':domainTimeWidth,
                  'DOMAIN_GAP':domainGap, 'PREPROCESS_FILE':preprocessFile}

        mp = mt.to_matlab_process(params)

        sout = os.path.join(SLURM_DIR, 'preproc_crcns.%s.%%j.out.txt' % fileDesc)
        serr = sout
        scriptName = os.path.join(SLURM_DIR, 'preproc_crcns.%s.sbatch' % fileDesc)

        sparams = {'out': sout,
                   'err': serr,
                   'partition':'all',
                   'cpus': 2,
                   'script': scriptName}

        slurmBot.add(mp.get_commands(), sparams)


    slurmBot.run()
