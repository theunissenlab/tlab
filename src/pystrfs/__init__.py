import os
import operator
import numpy as np

import h5py

from browser.model import meta, ProtocolResponse

from browser.model.physiology import Recsite, ElectrodePosition, Block
from browser.model.physiology.extracellular import Unit
from browser.model.analysis.spikeinfo import ResponsePerformance

SRC_ROOT_DIR = '/auto/fdata/mschachter'
PYSTRFS_DIR = '/auto/k6/mschachter/pystrfs'
SLURM_DIR = os.path.join(PYSTRFS_DIR, 'slurm')
UNITS_DIR = os.path.join(PYSTRFS_DIR, 'units')
PREPROC_DIR = os.path.join(PYSTRFS_DIR, 'preproc')
TEMPLATE_DIR = 'pystrfs/templates'

RANDOM_SEED = 19382760346

stimClassToGroup = {'zf song (crcns)':'Con', 'songripple zf song (crcns)':'Songrip', 'ml noise (crcns)':'Flatrip'}
stimGroupToClass = dict([(stimClassToGroup[k], k) for k in stimClassToGroup.keys()])

def get_good_units():
    goodfile = os.path.join(UNITS_DIR, 'allgoodcells.csv')
    ulist = []
    f = open(goodfile, 'r')
    for ln in f.readlines():
        ln = ln.strip()
        if len(ln) > 0:
            u = get_unit_by_old_id(ln)
            ulist.append(u)
    f.close()
    return ulist

def get_good_units_by_area(area_name, ulist=None):
    if ulist is None:
        ulist = get_good_units()
    fulist = []
    for u in ulist:
        anames = [a.name for a in u.recsite.areas]
        if area_name in anames:
            fulist.append(u)
    return fulist


def create_response_performance_all():
    umap = get_all_dcp_units(True)

    for uId,u in umap.iteritems():
        create_response_performance(u)


def create_response_performance(unit):
    unitFile = map_dcp_unit_to_fs(unit)
    f = h5py.File(unitFile, 'r')

    blk = unit.recsite.blocks[0]

    for classResp in unit.class_resps:
        stimGroup = stimClassToGroup[classResp.presentation.class_name]
        infoLower = float(f['class_info'][stimGroup]['coherence']['info']['lower'][0])
        infoUpper = float(f['class_info'][stimGroup]['coherence']['info']['upper'][0])
        infoMean = float(f['class_info'][stimGroup]['coherence']['info']['mean'][0])

        #print '%s: %f (%f, %f)' % (stimGroup, infoMean, infoLower, infoUpper)

        if classResp.performance is None:
            rp = ResponsePerformance()
            rp.info_mean = infoMean
            rp.info_upper = infoUpper
            rp.info_lower = infoLower
            rp.response = classResp
            meta.Session.add(rp)
            meta.Session.commit()
        else:
            classResp.performance.info_mean = infoMean
            classResp.performance.info_upper = infoUpper
            classResp.performance.info_lower = infoLower

    f.close()

def create_response_performance_individual(unit):

    unitFile = map_dcp_unit_to_fs(unit)

    f = h5py.File(unitFile, 'r')

    print 'Recording individual performances for %s' % unit.recsite.old_id

    for classResp in unit.class_resps:

        stimGroup = stimClassToGroup[classResp.class_name]

        for k,repResp in enumerate(classResp.repeated):
            presNum = k+1
            respName = '%d' % presNum
            infoGrp = f['class_info'][stimGroup]['responses'][respName]['coherence']['info']
            infoLower = float(infoGrp['lower'][0])
            if np.isnan(infoLower) or np.isinf(infoLower):
                infoLower = 0.0
            infoUpper = float(infoGrp['upper'][0])
            if np.isnan(infoUpper) or np.isinf(infoUpper):
                infoUpper = 0.0
            infoMean = float(infoGrp['mean'][0])
            if np.isnan(infoMean) or np.isinf(infoMean):
                infoMean = 0.0

            if repResp.performance is None:
                rp = ResponsePerformance()
                rp.info_mean = infoMean
                rp.info_upper = infoUpper
                rp.info_lower = infoLower
                rp.response = repResp
                repResp.performance = rp
                meta.Session.add(rp)
                meta.Session.commit()
            else:
                repResp.performance.info_mean = infoMean
                repResp.performance.info_upper = infoUpper
                repResp.performance.info_lower = infoLower
                meta.Session.commit()

            #print 'Created ResponsePerformance for resp %d, info=%0.1f' % (repResp.id, infoMean)

    f.close()


def map_dcp_units_to_fs():
    umap = get_all_crcns_units()
    for uId, u in umap.iteritems():
        map_dcp_unit_to_fs(u)


def get_blacklist():
    blistPath = os.path.join(UNITS_DIR, 'units.blacklist.txt')
    f = open(blistPath, 'r')
    lns = f.readlines()
    f.close()
    lns = [s.strip() for s in lns]
    blist = filter(lambda(x): len(x) > 0, lns)
    return blist


def get_all_dcp_units(useBlacklist = False):

    blackList = []
    if useBlacklist:
        blackList = get_blacklist()

    umap = {}

    prs = meta.Session.query(ProtocolResponse).filter(ProtocolResponse.protocol_name.like('%crcns%'))
    for pr in prs:
        u = pr.unit
        uId = u.recsite.old_id
        if uId not in umap and uId not in blackList:
            umap[uId] = u

    return umap


def get_unit_by_old_id(oldId):

    rq = meta.Session.query(Recsite).filter(Recsite.old_id.like(oldId))
    rs = rq.one()
    u = rs.units[0]

    return u


def get_groups_for_unit(unit):

    groups = {}
    blk = unit.recsite.blocks[0]
    for classPres in blk.class_presns:
        stimGroup = stimClassToGroup[classPres.class_name]
        cnt = len(classPres.repeated)
        nHoldOut = int(cnt / 10)
        nTrain = cnt - nHoldOut
        trainingGroups = np.arange(nTrain) + 1
        vg = []
        for k in range(nHoldOut):
            vg.append(cnt - k)
        validationGroups = np.array(vg)

        groups[stimGroup] = (trainingGroups, validationGroups)

    return groups


def get_holdout_sets(unit, holdoutFraction = 0.10, numSets = 5):

    holdoutsByStim = {}

    for classResp in unit.class_resps:
        stimGroup = stimClassToGroup[classResp.presentation.class_name]
        rperfs = [r.performance.info_mean for r in classResp.repeated]
        rperfs = np.array(rperfs)
        rmean = rperfs.mean()
        rstd = rperfs.std()

        nsamps = int(len(rperfs)*holdoutFraction)
        totalHoldoutsNeeded = nsamps*numSets

        #find a multiplier where enough samples are within the standard deviation
        rmult = 0.5
        while (np.abs(rperfs - rmean) < rmult*rstd).sum() < totalHoldoutsNeeded:
            rmult *= 1.1
        print 'Found %d holdouts that are within %0.2f stds of the mean' % (totalHoldoutsNeeded, rmult)

        hoIndx = (np.abs(rperfs - rmean) < rmult*rstd) * (np.arange(0, len(rperfs)) + 1)
        allPotentialHoldouts = hoIndx[hoIndx > 0] - 1

        #random sample within holdouts
        allHoldouts = []
        while len(allHoldouts) < totalHoldoutsNeeded:
            indx = np.random.randint(0, len(allPotentialHoldouts))
            hoIndx = allPotentialHoldouts[indx]
            if hoIndx not in allHoldouts:
                allHoldouts.append(hoIndx)

        if nsamps == 1:
            holdoutSets = []
            for ho in allHoldouts:
                holdoutSets.append([ho])
        else:
            #organize holdout sets so each one has similar mean info value
            holdoutSets = []
            hoPerfs = [rperfs[ho] for ho in allHoldouts]
            hoList = zip(allHoldouts, hoPerfs)
            hoList.sort(key=operator.itemgetter(1))
            for k in range(int(totalHoldoutsNeeded/2)):   #assumes nsamps is even
                hoSet = [allHoldouts[k], allHoldouts[-(k+1)]]
                holdoutSets.append(hoSet)

        holdoutsByStim[stimGroup] = holdoutSets

    return holdoutsByStim

def get_es_sets(unit, holdoutsByStim, numEsSets = 3):

    esByStim = {}

    for classResp in unit.class_resps:
        stimGroup = stimClassToGroup[classResp.presentation.class_name]
        rperfs = [r.performance.info_mean for r in classResp.repeated]
        rperfs = np.array(rperfs)
        rmean = rperfs.mean()
        rstd = rperfs.std()

        esSets = []
        holdoutSets = holdoutsByStim[stimGroup]
        for hoSet in holdoutSets:
            nonHoIndx = np.array(filter(lambda x: x not in hoSet, range(0, len(rperfs))))
            nonHoPerfs = np.array([rperfs[k] for k in nonHoIndx])

            totalEsNeeded = numEsSets*len(hoSet)

            rmult = 0.5
            while (np.abs(nonHoPerfs - rmean) < rmult*rstd).sum() < totalEsNeeded:
                rmult *= 1.1

            nhindx = np.abs(nonHoPerfs - rmean) < rmult*rstd
            nonHoIndx = nonHoIndx[nhindx]

            esSet = []
            while len(esSet) < totalEsNeeded:
                indx = np.random.randint(0, len(nonHoIndx))
                if nonHoIndx[indx] not in esSet:
                    esSet.append(nonHoIndx[indx])
            esSets.append(esSet)

        esByStim[stimGroup] = esSets

    return esByStim

def get_data_groups(unit, addOne = False):

    np.random.seed(RANDOM_SEED)

    holdoutsByStim = get_holdout_sets(unit)
    esByStim = get_es_sets(unit, holdoutsByStim)

    dataGroupsByStim = {}

    for classResp in unit.class_resps:
        stimGroup = stimClassToGroup[classResp.presentation.class_name]
        hoSets = holdoutsByStim[stimGroup]
        esSets = esByStim[stimGroup]

        dataGroups = []

        for k,hoSet in enumerate(hoSets):
            esSet = esSets[k]
            tgSet = []
            for m in range(0, classResp.repeated.count()):
                if m not in hoSet:
                    tgSet.append(m)
            dg = [np.array(tgSet), np.array(hoSet), np.array(esSet)]
            if addOne:
                for m in range(0, len(dg)):
                    dg[m] = dg[m] + 1
            dataGroups.append(dg)
        dataGroupsByStim[stimGroup] = dataGroups

    return dataGroupsByStim


def read_nvpair_file(filePath):
    nvPairs = {}
    f = open(filePath, 'r')
    lns = f.readlines()
    for ln in lns:
        ln = ln.strip()
        if len(ln) > 0:
            p = ln.split('=')
            if len(p) == 2:
                nvPairs[p[0]] = p[1].strip()

    f.close()
    return nvPairs


def get_stim_classes_for_unit(unit):
    groups = []
    blk = unit.recsite.blocks[0]
    for classPres in blk.class_presns:
        groups.append(stimClassToGroup[classPres.class_name])

    return groups


def get_unit_dir(unit):
    desc = unit.recsite.old_id
    unitDir = os.path.join(PYSTRFS_DIR, 'units', desc)
    return unitDir

def map_dcp_unit_to_fs(unit):

    if not unit.__class__.__name__ == 'DCPUnit':
        print 'No mapping method for unit of type %s' % unit.__class__.__name__
        return

    blk = unit.recsite.blocks[0]
    unitDir = get_unit_dir(unit)

    if not os.path.exists(unitDir):
        print 'Creating directory %s' % unitDir
        os.mkdir(unitDir)
        outputDir = os.path.join(unitDir, 'output')
        os.mkdir(outputDir)

    unitFile = os.path.join(unitDir, 'unit.h5')
    if not os.path.exists(unitFile):
        create_dcp_unit_file(unit, unitFile)

    return unitFile


def create_dcp_unit_file(unit, unitFile):

    blk = unit.recsite.blocks[0]

    f = h5py.File(unitFile, 'w')
    f.attrs['description'] = unit.recsite.old_id
    f.attrs['protocol_id'] = blk.protocol.id
    f.attrs['unit_id'] = unit.id

    print 'Creating unit file for %s' % unit.recsite.old_id

    for classPres in blk.class_presns:

        stimGroup = stimClassToGroup[classPres.class_name]
        classGroupName = '/class_responses/%s' % stimGroup
        cpGroup = f.create_group(classGroupName)
        cpGroup.attrs['num_stims'] = len(classPres.repeated)

        for k,repPres in enumerate(classPres.repeated):
            presNum = k+1
            spikeTrials = unit.spiketrials(repPres.stim_group, presNum)

            stimGroupName = '%s/%d' % (classGroupName, presNum)
            print 'Creating stim group: %s' % stimGroupName

            stimGroup = f.create_group(stimGroupName)
            stimGroup.attrs['stim_md5'] = repPres.stim.file.MD5
            stimGroup.attrs['protocol_stim_id'] = repPres.stim.protocol_stim_id
            stimGroup.attrs['num_trials'] = len(spikeTrials)
            stimGroup.attrs['trial_duration'] = float(repPres.stim.file.duration)

            for trialNum,spikes in enumerate(spikeTrials):
                trialResponseName = '%s/%d' % (stimGroupName, trialNum+1)
                #print '\tCreating trial response with %d spikes: %s' % (len(spikes), trialResponseName)
                trialResponse = f.create_dataset(trialResponseName, (1, len(spikes)), 'f')
                trialResponse[:] = spikes

    print 'Writing unit hdf5 file to %s...' % unitFile
    f.close()

def get_best_lyons_params(unit):
    unitDir = get_unit_dir(unit)
    lyonsFile = os.path.join(unitDir, 'output', 'bestlyons.txt')
    if not os.path.exists(lyonsFile):
        print 'No Lyons best param file exists at: %s' % lyonsFile
        return None
    return read_nvpair_file(lyonsFile)
