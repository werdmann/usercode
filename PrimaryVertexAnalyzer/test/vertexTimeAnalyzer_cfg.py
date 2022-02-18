import FWCore.ParameterSet.Config as cms

# command-line arguments
import FWCore.ParameterSet.VarParsing as VarParsing
opts = VarParsing.VarParsing('analysis')

opts.register('skipEvents', 0,
              VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'number of events to be skipped')

opts.register('dumpPython', None,
              VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'path to python file with content of cms.Process')

opts.register('numThreads', 1,
              VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'number of threads')

opts.register('numStreams', 0,
              VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'number of streams')

opts.register('reco', 'Phase2_D77',
              VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'keyword to select base configuration file')

opts.register('globalTag', None,
              VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.string,
              'argument of process.GlobalTag')

opts.register('wantSummary', False,
              VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.bool,
              'show cmsRun summary at job completion')

opts.register('verbosity', 0,
              VarParsing.VarParsing.multiplicity.singleton,
              VarParsing.VarParsing.varType.int,
              'level of output verbosity')

# custom defaults
opts.setDefault('maxEvents', 10)
opts.setType('outputFile', VarParsing.VarParsing.varType.string)
opts.setDefault('outputFile', 'tmp.root')

opts.parseArguments()

###
### base configuration file
###
cmsDriverArgsDict = {

  'Phase2_D77': """
 --conditions auto:phase2_realistic_T21
 --era Phase2C11I13M9
 --geometry Extended2026D77
 --filein /store/relval/CMSSW_12_3_0_pre4/RelValTTbar_14TeV/GEN-SIM-RECO/PU_123X_mcRun4_realistic_v3_2026D77PU200-v2/10000/b9af4417-47fc-4a82-bc60-de376040fb83.root
 --customise SLHCUpgradeSimulations/Configuration/aging.customise_aging_1000
 --step RECO
""",
}

if opts.reco in cmsDriverArgsDict:
  cmsDriverCmd = cmsDriverArgsDict[opts.reco]
  # remove newlines and spurious whitespaces
  cmsDriverCmd = ' '.join(cmsDriverCmd.replace('\n','').split())

  if len(opts.inputFiles) > 0:
    cmsDriverCmd += ' --filein '+','.join(opts.inputFiles)
else:
  raise RuntimeError('invalid argument for option "reco": "'+opts.reco+'" (allowed values: '+str(sorted(cmsDriverArgsDict.keys()))+')')

# create, load and remove the configuration file
import glob
import os

cfg_filepath = os.environ['PWD']+'/tmp_cfg.py'

# Process
def EXE(cmd, suspend=True, verbose=False, dry_run=False):
    if verbose: print('\033[1m'+'>'+'\033[0m'+' '+cmd)
    if dry_run: return
    _exitcode = os.system(cmd)
    _exitcode = min(255, _exitcode)
    if _exitcode and suspend:
       raise RuntimeError(_exitcode)
    return _exitcode

EXE('cmsDriver.py '+cmsDriverCmd+' --processName ANALYSIS --number 100 --no_exec --python_filename '+cfg_filepath)
from tmp_cfg import cms, process
for _tmpf in glob.glob(cfg_filepath+'*'):
  os.remove(_tmpf)

# remove OutputModules of base configuration
for _modname in process.outputModules_():
  _mod = getattr(process, _modname)
  if type(_mod) == cms.OutputModule:
    process.__delattr__(_modname)
    if opts.verbosity > 0:
      print('> removed cms.OutputModule:', _modname)

# Input source
if opts.inputFiles != []:
  process.source.fileNames = opts.inputFiles
  process.source.secondaryFileNames = opts.secondaryInputFiles
process.source.skipEvents = cms.untracked.uint32(opts.skipEvents)

# Output file (TFileService)
process.TFileService = cms.Service('TFileService',
    fileName = cms.string(opts.outputFile)
)

# Job configuration parameters (special PSets)
process.maxEvents.input = opts.maxEvents

process.options.numberOfConcurrentLuminosityBlocks = 0
process.options.numberOfConcurrentRuns = 1
process.options.numberOfStreams = max(opts.numStreams, 0)
process.options.numberOfThreads = max(opts.numThreads, 1)
process.options.wantSummary = opts.wantSummary

# GlobalTag
if opts.globalTag is not None:
    from Configuration.AlCa.GlobalTag import GlobalTag
    process.GlobalTag = GlobalTag(process.GlobalTag, opts.globalTag, '')

###
### Schedule
###
process.setSchedule_(cms.Schedule())

# PSet for vertex-time algorithm "VertexTimeAlgorithmFromTracksPID"
# called implicitly via offlinePrimaryVerticesWithTimeFromTracksPID.algorithm.refToPSet_
from usercode.PrimaryVertexAnalyzer.VertexTimeAlgorithmFromTracksPID_cfi import VertexTimeAlgorithmFromTracksPID as _VertexTimeAlgorithmFromTracksPID
process.PSetVertexTimeAlgorithmFromTracksPID = _VertexTimeAlgorithmFromTracksPID.clone(
    trackMTDTimeVMapTag = 'trackExtenderWithMTD:generalTracktmtd',
    trackMTDTimeErrorVMapTag = 'trackExtenderWithMTD:generalTracksigmatmtd',
    trackMTDTimeQualityVMapTag = 'mtdTrackQualityMVA:mtdQualMVA',
    trackMTDMomentumVMapTag = 'trackExtenderWithMTD:generalTrackp',
    trackMTDPathLengthVMapTag = 'trackExtenderWithMTD:generalTrackPathLength',

    minTrackVtxWeight = 0.5,
    minTrackTimeQuality = 0.,
    massPion = 0.139570,
    massKaon = 0.493677,
    massProton = 0.938272,
    probPion = 0.7,
    probKaon = 0.2,
    probProton = 0.1,
    coolingFactor = 0.5,
)

from usercode.PrimaryVertexAnalyzer.vertexTimeProducer_cfi import vertexTimeProducer as _vertexTimeProducer
process.offlinePrimaryVerticesWithTimeFromTracksPID = _vertexTimeProducer.clone(
    vertices = 'offlinePrimaryVertices',
    produceVertices = True,
    algorithm = dict(
      refToPSet_ = cms.string('PSetVertexTimeAlgorithmFromTracksPID'),
    )
)

def histoPSet(_valueMapTag, _nbins, _xmin, _xmax):
  return cms.PSet(valueMapTag = cms.InputTag(_valueMapTag), nbins = cms.uint32(_nbins), xmin = cms.double(_xmin), xmax = cms.double(_xmax))

from usercode.PrimaryVertexAnalyzer.vertexTimeAnalyzer_cfi import vertexTimeAnalyzer as _vertexTimeAnalyzer
process.VertexHistograms_offlinePrimaryVertices = _vertexTimeAnalyzer.clone(
    vertices = 'offlinePrimaryVertices',
    maxNumberOfVertices = -1,
    histogramNamePrefix = 'vertex_',
    histogramPSet = dict(
        t_fromTracksPID = histoPSet('offlinePrimaryVerticesWithTimeFromTracksPID:time', 600, -30, 30),
        tError_fromTracksPID = histoPSet('offlinePrimaryVerticesWithTimeFromTracksPID:timeError', 300, 0, 30)
    )
)

process.vertexAnalysisPath = cms.Path(
    process.offlinePrimaryVerticesWithTimeFromTracksPID
  + process.VertexHistograms_offlinePrimaryVertices
)

process.schedule_().append(process.vertexAnalysisPath)

## PrimaryVertexAnalyzer4PU
parameters = { # can be overwritten by arguments of the same name
  '4D': cms.untracked.bool(True),
  'selNdof': cms.untracked.double(4.0),
  'selNdofWithBS': cms.untracked.double(7.0),
#  'splitmethod': cms.untracked.int32(0),
  'usefit': cms.untracked.bool(False),
  'use_tp': cms.untracked.bool(True),
  'use2file': cms.untracked.bool(False),
#  'ptrunc': cms.untracked.double(1.e-50),
  'fill_track_histos': cms.untracked.bool(True),
  'fplo': cms.untracked.double(0),
  'chi2cutoff': cms.double(2.5),
  'verboseAnalyzer': cms.untracked.bool(False),
  'verboseProducer': cms.untracked.bool(False),
  'verboseClusterizer': cms.untracked.bool(False),
  'verboseClusterizer2D': cms.untracked.bool(False),
  'zdumpcenter': cms.untracked.double(0.0),
  'zdumpwidth': cms.untracked.double(20.0),
  'reco': cms.untracked.bool(False),
  'miniaod': cms.untracked.bool(False),
  'autodump': cms.untracked.int32(0),
  'nDump': cms.untracked.int32(0),
  'nDumpTracks': cms.untracked.int32(0),
#  'mintrkweight': cms.untracked.double(0.5),
  'uniquetrkweight': cms.double(0.8),
  'uniquetrkminp': cms.double(0.0),
  'zmerge': cms.double(1.e-2),
  'coolingFactor': cms.double(0.6),
  'Tmin': cms.double(2.0),
  'Tpurge': cms.double(2.0),
  'Tstop': cms.double(0.5),
  'vertexSize': cms.double(0.006),
  'vertexSizeTime': cms.double(0.008),
  'd0CutOff': cms.double(3.),
  'dzCutOff': cms.double(3.),
  'zrange': cms.double(4),
  'convergence_mode': cms.int32(0),
  'delta_lowT': cms.double(1.e-3),
  'delta_highT': cms.double(1.e-2),
  # track selection
  'maxNormalizedChi2': cms.double(10.0),
  'minPixelLayersWithHits': cms.int32(2),
  'minSiliconLayersWithHits': cms.int32(5),
  'maxD0Significance': cms.double(4.0),
  'minPt': cms.double(0.0),
  'maxEta': cms.double(4.0),
  'trackQuality': cms.string('any'),
  # track selection, experimental
  'maxDzError': cms.double(1.0),
  'maxD0Error': cms.double(1.0),
  # vertex selection
  'minNdof': cms.double(0.0),
  # temp
  'trackTimeQualityThreshold': cms.untracked.double(0.8),
  'purge': cms.untracked.int32(0),
  'rho0mode': cms.untracked.int32(0),  # /nt, as before
  'mergenotc': cms.untracked.bool(False),
  'mergeafterpurge': cms.untracked.bool(False),
  'fillzmerge': cms.untracked.double(0.0060), # default = vertexSize
  'use_hitpattern': cms.untracked.int32(0),
  'use_pt': cms.untracked.int32(0)
}

from RecoVertex.Configuration.RecoVertex_cff import unsortedOfflinePrimaryVertices as _unsortedOfflinePrimaryVertices
tkFilterParameters = _unsortedOfflinePrimaryVertices.TkFilterParameters.clone()
print('original trackFilterParameters (z-clustering)')
print(tkFilterParameters)
for par_name in [
  'maxNormalizedChi2',
  'minPixelLayersWithHits',
  'minSiliconLayersWithHits',
  'maxD0Significance',
  'maxD0Error',
  'maxDzError',
  'minPt',
  'maxEta',
  'trackQuality',
]:
  try:
    default = getattr(tkFilterParameters, par_name)
    if default != parameters[par_name]:
      print('changing tkFilter parameter', par_name, 'from', default, 'to', parameters[par_name])
    setattr(tkFilterParameters, par_name, parameters[par_name])
  except ValueError:
    print('pva_cfg: attribute tkFilterParameters.'+par_name, 'not found')

## PV analyser
process.vertexAnalyser = cms.EDAnalyzer('PrimaryVertexAnalyzer4PU',
  info = cms.untracked.string('test'),
  f4D = parameters['4D'],
  selNdof = parameters['selNdof'],
  selNdofWithBS = parameters['selNdofWithBS'],
  beamSpot = cms.InputTag('offlineBeamSpot'),
  simG4 = cms.InputTag('g4SimHits'),
  outputFile = cms.untracked.string(opts.outputFile[:-5]+'_pvtx.root'),
  verbose = parameters['verboseAnalyzer'],
  veryverbose = cms.untracked.bool(False),
  recoTrackProducer = cms.untracked.string('generalTracks'),
  zmatch = cms.untracked.double(0.05),
  autodump = parameters['autodump'],
  nDump = parameters['nDump'],
  nDumpTracks = parameters['nDumpTracks'],
  RECO = parameters['reco'],
  MINIAOD = parameters['miniaod'],
  use_tp = parameters['use_tp'],
  fill_track_histos = parameters['fill_track_histos'],
  track_timing = cms.untracked.bool(True),
  TkFilterParameters = tkFilterParameters,
  trackingParticleCollection = cms.untracked.InputTag('mix', 'MergedTrackTruth'),
  trackingVertexCollection = cms.untracked.InputTag('mix', 'MergedTrackTruth'),
  trackAssociatorMap = cms.untracked.InputTag('trackingParticleRecoTrackAsssociation'),
  TrackTimesLabel = cms.untracked.InputTag('tofPID:t0safe'),
  TrackTimeResosLabel = cms.untracked.InputTag('tofPID:sigmat0safe'),
#  TrackTimesLabel = cms.untracked.InputTag('tofPID:t0safe'),
#  TrackTimeResosLabel = cms.untracked.InputTag('tofPID:sigmat0safe'),
  TrackTimeQualityMapLabel = cms.untracked.InputTag('mtdTrackQualityMVA:mtdQualMVA'),
  TrackTimeQualityThreshold = cms.untracked.double(parameters['trackTimeQualityThreshold'].value()),
  vertexAssociator = cms.untracked.InputTag('VertexAssociatorByPositionAndTracks'),
  lumiInfoTag = cms.untracked.InputTag('LumiInfo', 'brilcalc'),
  useVertexFilter = cms.untracked.bool(False),
  compareCollections = cms.untracked.int32(0),
  vertexRecoCollections = cms.VInputTag(
    'offlinePrimaryVerticesWithTimeFromTracksPID',
  ),
)

process.vertexAnalysisEndPath = cms.EndPath(process.vertexAnalyser)
process.schedule_().append(process.vertexAnalysisEndPath)

## recoTrack<->TrackingParticle associator
if parameters['use_tp']:
  process.load('Validation.RecoTrack.TrackValidation_cff')
  process.trackingParticleAssociationTask = cms.Task(
    process.tpClusterProducer,
    process.quickTrackAssociatorByHits,
    process.trackingParticleRecoTrackAsssociation
  )
  process.trackingParticleAssociationPath = cms.Path(process.trackingParticleAssociationTask)
  process.schedule_().append(process.trackingParticleAssociationPath)

# prune and dump content of cms.Process to python file
if opts.dumpPython is not None:
    process.prune()
    open(opts.dumpPython, 'w').write(process.dumpPython())

# printouts
if opts.verbosity > 0:
    print('--- vertexTimeAnalyser_cfg.py ---')
    print('')
    print('option: outputFile =', opts.outputFile)
    print('option: dumpPython =', opts.dumpPython)
    print('')
    print('process.GlobalTag =', process.GlobalTag.dumpPython())
    print('process.source =', process.source.dumpPython())
    print('process.maxEvents =', process.maxEvents.dumpPython())
    print('process.options =', process.options.dumpPython())
    print('-------------------------------')
