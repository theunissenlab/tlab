import os
import crcns

from browser.model import meta
from browser.model.colony import bird as Bird

from browser.model.physiology import session, Penetration, Recsite, ElectrodePosition, Block
from browser.model.physiology.extracellular import DCPUnit as Unit
from browser.model.physiology.spikesort import null as NullSort
from browser.model.physiology.presentation import (Protocol as ProtocolPresentation,\
                                                   ClassPresentation,\
                                                   repeated_stim as RepeatedStimPresentation,\
                                                   trial as TrialPresentation)
from browser.model.anatomy import Area


from browser.model.stims import protocol, protocol_stim, UniqueSoundFile, Sound, soundfile, stim_group

from browser.model.dcp import (Cell as DCPCell, SpikeFile as DCPSpikeFile, Spike as DCPSpike)


def get_or_create_unique_soundfile(wavFile):

   uq = meta.Session.query(UniqueSoundFile).filter_by(MD5=wavFile.md5)
   if uq.count() > 0:
      uniqSound = uq.one()

      sound = uniqSound.unique_sound
      sound.stimtype = stimTypeNames[wavFile.stimType]
      sound.stimtype_detail = wavFile.stimType
      meta.Session.commit()

   else:
      #create unique sound
      print 'creating sound object, type=%s' % stimTypeNames[wavFile.stimType]
      sound = Sound(stimTypeNames[wavFile.stimType], wavFile.stimType)
      meta.Session.add(sound)

      #create unique soundfile
      uniqSound = UniqueSoundFile(wavFile.fullPath)
      uniqSound.unique_sound=sound
      print 'creating unique sound, md5=%s' % wavFile.md5
      meta.Session.add(uniqSound)
      meta.Session.commit()

   #create mapping to filesystem
   sq = meta.Session.query(soundfile).filter_by(name=wavFile.fullPath)
   if sq.count() == 0:   
      #sfile = soundfile(old_stimtype=wavFile.stimType, filename=wavFile.fullPath,
      #                  unique_sound=sound, unique_soundfile=uniqSound)
      sfile = soundfile(wavFile.fullPath)
      sfile.unique_sound_id = sound.id
      print 'creating soundfile, path=%s' % wavFile.fullPath

      meta.Session.add(sfile)
      meta.Session.commit()
      

def insert_unique_stims():

   wavData = get_wavdata()

   for wfile in wavData:
      get_or_create_unique_soundfile(wfile)

def insert_unique_protocols():

   (cellList, protocols) = crcns.identify_unique_protocols()

   stimGroupMap = {'conspecific':'Con', 'flatrip':'Flatrip', 'songrip':'Songrip', 'bengalese':'Bengal'}
      
   for proto in protocols:

      #does a protocol already exist with same name in DB?
      pq = meta.Session.query(protocol).filter_by(name=proto.name)
      if pq.count() == 0:
         print 'Inserting new protocol \'%s\'' % proto.name
         dbProto = protocol(proto.name, 'STRF')
         meta.Session.add(dbProto)
         meta.Session.commit()
      else:
         dbProto = pq.one()
      
      #insert stimuli for protocol         
      for k,md5Name in enumerate(proto.stims):
         stimNum = k+1

         #get stim group
         sgName = stimGroupMap[proto.stimTypes[k]]
         scq = meta.Session.query(stim_group).filter(stim_group.name.like(sgName))
         stimGroup = scq.one()

         #get sound file
         sq = meta.Session.query(soundfile).filter_by(MD5=md5Name)\
                                           .filter(soundfile.name.like('/auto/k6/tlab/crcns/%'))

         #get soundfile corresponding to stim
         if sq.count() == 0:
            print 'No sound file corresponding to md5=%s!!' % md5Name            
         else:
            sfile = sq.one()
            dbProtoStim = protocol_stim(sfile, stimGroup, stimNum, dbProto)
            meta.Session.add(dbProtoStim)

      meta.Session.commit()

def insert_unique_blocks():

   birdList = crcns.get_bird_hierarchy()

   for bird in birdList.values():
      
      for blk in bird.blocks.values():
         oldId = 'CRCNS %s_%s' % (bird.name, blk.name)

         #make sure all cells in block correspond to same protocol
         protoList = []
         for cell in blk.units.values():
            proto = cell.protocol
            if proto.name not in protoList:
               protoList.append(proto.name)
            #print '%s: %s' % (cell.name(), proto.name)
         if len(protoList) > 1:
            print 'Block %s has more than one protocol defined:' % oldId
            for cell in blk.units.values():
               print '\t%s' % cell
            exit()

         #query for already-existing block
         bq = meta.Session.query(Block).filter(Block.old_block_id.like(oldId))
         if bq.count() == 0:
            #get protocol defined for block
            pq = meta.Session.query(protocol).filter(protocol.name.like(blk.protocol.name))
            if pq.count() == 0:
               print 'No protocol found with name %s!' % blk.protocol.name
               exit()
            else:
               proto = pq.one()

            print 'Creating new block, old_id=%s' % oldId
            b = Block(protocol=proto, type='dcp', old_block_id=oldId)
            meta.Session.add(b)
            meta.Session.commit()


def insert_unique_birds():
   birdList = crcns.get_bird_hierarchy()

   for bird in birdList.values():   
      bq = meta.Session.query(Bird).filter(Bird.bandnumber.like(bird.name))
      if bq.count() == 0:
         print 'Creating bird %s...' % bird.name
         b = Bird(bandnumber=bird.name, sex='male')
         meta.Session.add(b)

   meta.Session.commit()

def insert_unique_units():
   
   birdList = crcns.get_bird_hierarchy()
   penMap = {'A':1, 'B':2, 'C':3, 'D':4}

   for bird in birdList.values():   
      bq = meta.Session.query(Bird).filter(Bird.bandnumber.like(bird.name))
      birdObj = bq.one()

      for (blkId, blk) in bird.blocks.iteritems():         
         for (pairId, cell) in blk.units.iteritems():
            penNum = penMap[pairId]
            penId = '%s_%s' % (bird.name, pairId)
            rsOldId = cell.name()
            #look for penetration
            pq = meta.Session.query(Penetration).filter(Penetration.old_id.like(penId))
            if pq.count() == 0:
               print 'Creating penetration %s' % penId
               pen = Penetration(old_id=penId, bird=birdObj, number=penNum, electrode_id=1, refpoint='None')
               meta.Session.add(pen)
            else:
               pen = pq.one()
            
            print 'Creating electrode pos/recsite/unit for %s' % cell.name()
            #create electrode position
            eDepth = int(blkId)
            epos = ElectrodePosition(penetration=pen, depth=eDepth)
            meta.Session.add(epos)
            #create recsite
            rec = Recsite(old_id=rsOldId, position=epos, electrode_site_number=1)
            meta.Session.add(rec)
            #create null sort
            ns = NullSort(rec)               
            meta.Session.add(ns)
            #create unit
            u = Unit(rec, ns, 1)
            meta.Session.add(u)

   meta.Session.commit()


def map_blocks_to_epos():

   birdList = crcns.get_bird_hierarchy()
   penMap = {'A':1, 'B':2, 'C':3, 'D':4}

   for bird in birdList.values():   
      bq = meta.Session.query(Bird).filter(Bird.bandnumber.like(bird.name))
      birdObj = bq.one()

      for (blkId, blk) in bird.blocks.iteritems():         

         for (pairId, cell) in blk.units.iteritems():
            #get recsite by old_id
            rsq = meta.Session.query(Recsite).filter(Recsite.old_id.like(cell.name()))
            if rsq.count() == 0:
               print 'Cannot find recsite with old_id=%s' % cell.name()
            else:
               rec = rsq.one()
               epos = rec.position
               oldBlockId = 'CRCNS %s_%s' % (bird.name, blk.name)
               bq = meta.Session.query(Block).filter(Block.old_block_id.like(oldBlockId))
               if bq.count() == 0:
                  print 'Cannot find block with old_block_id=%s' % oldBlockId
               else:
                  blockObj = bq.one()
                  print 'Adding electrode position %s to block %s' % (cell.name(), oldBlockId)
                  blockObj.positions.append(epos)
   meta.Session.commit()


def insert_dcp_data():

   birdList = crcns.get_bird_hierarchy()

   stimTypeToGroup = {'conspecific':'Con', 'flatrip':'Flatrip', 'songrip':'Songrip'}

   for bird in birdList.values():
      for (blkId, blk) in bird.blocks.iteritems():         
         for (pairId, cell) in blk.units.iteritems():
      
            #get unit corresponding to this cell
            rq = meta.Session.query(Recsite).filter(Recsite.old_id.like(cell.name()))
            if rq.count() == 0:
               print 'Cannot find recsite corresponding to cell %s' % cell.name()
               exit()
            recsiteObj = rq.one()
            unitObj = recsiteObj.units[0]

            cq = meta.Session.query(DCPCell).filter(DCPCell.name.like(cell.name()))
            if cq.count() == 0:
               cellObj = DCPCell()
               cellObj.file_path = cell.path
               cellObj.name = cell.name()
               cellObj.unit = unitObj
               meta.Session.add(cellObj)
            else:
               cellObj = cq.one()

            #insert spike files
            for k,spikeFile in enumerate(cell.spikeFiles):
               stimNumber = cell.protocol.stimNumbers[k]
               stimType = cell.protocol.stimTypes[k]
               
               spikeFileObj = DCPSpikeFile()
               spikeFileObj.cell = cellObj
               spikeFileObj.file_name = spikeFile.name
               spikeFileObj.stim_type = stimTypeToGroup[stimType]
               spikeFileObj.stim_number = stimNumber

               meta.Session.add(spikeFileObj)
               
               #insert spikes
               for trialNumber,spikeTrial in enumerate(spikeFile.trials):
                  
                  for ti in spikeTrial.spikes:
                     spikeObj = DCPSpike()
                     spikeObj.timestamp = ti
                     spikeObj.spike_file = spikeFileObj
                     spikeObj.trial_number = trialNumber 
                     meta.Session.add(spikeObj)

            meta.Session.commit()
            print 'Added DCPCell %s' % cell.name()
  

def get_num_crcns_trials(crcnsBlock, spikeFileIndex):
   """ Makes sure trial counts match up between cells with the same block """
   cells = []
   numTrials = []
   sfNames = []
   stNames = []
   for (pairId, cell) in crcnsBlock.units.iteritems():
      cells.append(cell)
      stNames.append(cell.protocol.stimTypes[spikeFileIndex])
      sfNames.append(cell.spikeFiles[spikeFileIndex].name)
      numTrials.append(len(cell.spikeFiles[spikeFileIndex].trials))

   #ensure # of trials is equal
   isSame = True
   nt = max(numTrials)
   for ont in numTrials:
      if ont != nt:
         isSame = False

   cstrs = [x.name() for x in cells]
   ntstrs = [str(x) for x in numTrials]

   if not isSame:
      print 'Trial counts don\'t match: %s | %s | %s | %s' % \
          (','.join(cstrs), ','.join(sfNames), ','.join(stNames), ','.join(ntstrs))

   return nt
   


def insert_presentations():

   birdList = crcns.get_bird_hierarchy()

   stimGroupToType = {'Con':'conspecific', 'Flatrip':'flatrip', 'Songrip':'songrip'}
   stimGroupToClass = {'Con':'zf song (crcns)', 'Songrip':'songripple zf song (crcns)', 'Flatrip':'ml noise (crcns)'}
   stimTypeToClass = {'conspecific':'zf song (crcns)', 'songrip':'songripple zf song (crcns)', 'flatrip':'ml noise (crcns)'}

   for bird in birdList.values():

      for (blkId, blk) in bird.blocks.iteritems():         
         oldId = 'CRCNS %s_%s' % (bird.name, blk.name)

         # get block
         bq = meta.Session.query(Block).filter(Block.old_block_id.like(oldId))
         if bq.count() == 0:
            print 'No block found! old_block_id=%s' % oldId
            exit()
         blockObj = bq.one()

         # create a ProtocolPresentation instance for this block
         ppq = meta.Session.query(ProtocolPresentation).filter_by(block=blockObj)
         if ppq.count() > 0:
            protocolPres = ppq.one()
         else:
            print 'Creating new protocol presentation for block %s' % oldId
            protocolPres = ProtocolPresentation()
            protocolPres.block = blockObj
            meta.Session.add(protocolPres)

         classPresMap = {}

         for k,protoStimObj in enumerate(blockObj.protocol.stims):
            
            #get or create a ClassPresentation object
            classPresName = stimGroupToClass[protoStimObj.group.name]
            if classPresName not in classPresMap:               
               cpq = meta.Session.query(ClassPresentation).filter_by(block_id=blockObj.id,class_name=classPresName)
               if cpq.count() > 0:
                  classPres = cpq.one()
               else:
                  classPres = create_class_presentation(blockObj, classPresName)
                  print '\tCreating class presentation %s, id=%d' % (classPresName, classPres.id)
               classPresMap[classPresName] = classPres

            else:
               classPres = classPresMap[classPresName]

            #create a repeated presentation object
            rpq = meta.Session.query(RepeatedStimPresentation).filter_by(block=blockObj,group=protoStimObj.group,stim_number=protoStimObj.number)
            if rpq.count() > 0:
               repPresObj = rpq.one()
            else:
               print '\tCreating repeated stim pres...'
               repPresObj = RepeatedStimPresentation(protoStimObj, blockObj)
               meta.Session.add(repPresObj)
               #add repeated presentation to class presentation                                                                            
               cname = stimGroupToClass[protoStimObj.group.name]

            #classPres.repeated.append(repPresObj)

            #create trial presentation objects
            numTrials = get_num_crcns_trials(blk, k)
            print 'Creating %d TrialPresentations...' % numTrials
            for m in range(numTrials):
               repNumber = m+1
               trialPresObj = TrialPresentation(repPresObj, repNumber, blockObj)
               meta.Session.add(trialPresObj)

      meta.Session.commit()


def map_responses():

   birdList = crcns.get_bird_hierarchy()

   for bird in birdList.values():

      for (blkId, blk) in bird.blocks.iteritems():         
         oldId = 'CRCNS %s_%s' % (bird.name, blk.name)

         for (pairId, cell) in blk.units.iteritems():
            rsq = meta.Session.query(Recsite).filter(Recsite.old_id.like(cell.name()))
            if rsq.count() == 0:
               print 'Cannot find recsite with old_id=%s' % cell.name()
               exit()
            else:
               recObj = rsq.one()

            if recObj.units.count() == 0:
               print 'No units for recsite with old_id=%s' % cell.name()
            
            unitObj = recObj.units[0]
            print 'Mapping responses for unit %s' % cell.name()            
            unitObj.map_responses()


def map_regions():

   cellList = crcns.get_all_cells()

   regStrMap = {'mld':'MLd', 'L':'L', 'L1':'L1', 'L2a':'L2a', 'L2b':'L2b', 'L3':'L3', 'CM':'CM', 'OV':'OV'}
   regMap = {}

   for crcnsReg,dbReg in regStrMap.iteritems():
      aq = meta.Session.query(Area).filter_by(name=dbReg)
      if aq.count() == 0:
         print 'No region associated with string %s' % dbReg
         return
      a = aq.one()
      regMap[crcnsReg] = a
      
   for cell in cellList.values():
      rsq = meta.Session.query(Recsite).filter(Recsite.old_id.like(cell.name()))
      if rsq.count() == 0:
         print 'No recsite found for cell %s' % cell.name()

      rs = rsq.one()            
      epos = rs.position

      if cell.region != 'None':
         reg = regMap[cell.region]
         epos.areas.append(reg)

         if cell.region in ['L1','L2a','L2b','L3']:
            epos.areas.append(regMap['L'])

      meta.Session.commit()


def test_responses():

   birdList = crcns.get_bird_hierarchy()

   stimGroupToType = {'Con':'conspecific', 'Flatrip':'flatrip', 'Songrip':'songrip'}
   stimTypeToClass = {'conspecific':'zf song', 'songrip':'songripple zf song', 'flatrip':'ml noise'}

   for bird in birdList.values():   
      for (blkId, blk) in bird.blocks.iteritems():
         for (pairId, cell) in blk.units.iteritems():
            #get recsite by old_id
            rsq = meta.Session.query(Recsite).filter(Recsite.old_id.like(cell.name()))
            if rsq.count() == 0:
               print 'Cannot find recsite with old_id=%s' % cell.name()
               exit()
            else:
               recObj = rsq.one()

            if recObj.units.count() == 0:
               print 'No units for recsite with old_id=%s' % cell.name()
            
            unitObj = recObj.units[0]
            

def create_class_presentation(blockObj, className):
   cp = ClassPresentation(className, blockObj)   
   meta.Session.add(cp)
   meta.Session.commit()

   sql = 'INSERT INTO physiology.class_presentations VALUES (%d,\'%s\')' % (cp.id, className)
   meta.Session.execute(sql)
   meta.Session.commit()
  
   return cp
            
         

def test_stims():

   stimGroupMap = {'conspecific':'Con', 'flatrip':'Flatrip', 'songrip':'Songrip', 'bengalese':'Bengal'}
   birdList = crcns.get_bird_hierarchy()

   for bird in birdList.values():   
      for (blkId, blk) in bird.blocks.iteritems():
         for (pairId, cell) in blk.units.iteritems():
            #get recsite by old_id
            rsq = meta.Session.query(Recsite).filter(Recsite.old_id.like(cell.name()))
            if rsq.count() == 0:
               print 'Cannot find recsite with old_id=%s' % cell.name()
            else:
               rec = rsq.one()
               #get block
               epos = rec.position
               blockObj = epos.blocks[0]
               print 'Cell %s' % cell.name()
               for stype in cell.stimTypes:
                  stimGroup = stimGroupMap[stype]
                  stimDir = os.path.join(cell.path, stype)

                  (stimMd5s, stimTypes) = crcns.get_stim_md5s_from_dir(stimDir, stype)
      
                  psq = blockObj.protocol.stims.filter(protocol_stim.stim_group.like(stimGroup))
                  
                  print '  Stim type %s:' % (stimGroup)
                  

                  if len(stimMd5s) != psq.count():
                     print '\tStimulus count mismatch, # on filesystem=%d, # in DB=%d' \
                           % (len(tsimMd5s), psq.count())

                  for k,pstim in enumerate(psq):                  
                     otherPath = os.path.join(crcns.STIMS_DIR, '%s.wav' % stimMd5s[k])                     
                     #print '%d,%s' % (k, pstim.file.name)

                     if otherPath != pstim.file.name:
                        print '\tFilename mismatch!'
                        print '\t\tDB:   %s' % pstim.file.name
                        print '\t\tReal: %s' % otherPath

