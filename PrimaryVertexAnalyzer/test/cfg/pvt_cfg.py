
import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList

from Configuration.StandardSequences.Eras import eras

import sys, os

# example
#cmsRun pvt_cfg.py input=/store/relval/CMSSW_11_2_0_pre1/RelValZMM_14/GEN-SIM-RECO/110X_mcRun4_realistic_v3_2026D49noPU-v1/10000/79C28F5A-A01A-8A40-9C5A-95FA87D2E6DD.root events=10 info=RUN4PU_new_A new=True redoPV=True 4D=True era=Phase2 maxEta=4.0 zrange=4. convergence_mode=0 delta_lowT=1e-3 delta_highT=1e-2

DO_VTX_RECO = False

args =""
nevent=-1
info = "test"
debug = False
eventsToProcess = "" #cms.untracked.VEventRange('1:4474')
zdumpcenter = 0
zdumpwidth = 100.
use_tp = False
use_hlt = False
json = "goodList.json"
era = "Run3"  # or Phase2 or Run2_2018,  see eras.__dict__.keys()


inputfile="/store/relval/CMSSW_11_2_0_pre1/RelValZMM_14/GEN-SIM-RECO/110X_mcRun4_realistic_v3_2026D49noPU-v1/10000/79C28F5A-A01A-8A40-9C5A-95FA87D2E6DD.root"
source_files = cms.untracked.vstring('')
vcollections= cms.VInputTag(
    "offlinePrimaryVertices"
#   "offlinePrimaryVertices4D"
#   "offlinePrimaryVerticesWithBS",
)


parameters={  # can be overwritten by arguments of the same name
  "clustering":cms.string("DA_vect"),
  "clustering2D":cms.string("DA2D_vect"),
  "4D": cms.untracked.bool(False),
#  "splitmethod":cms.untracked.int32(0),
  "usefit":cms.untracked.bool(False),
  "ptrunc":cms.untracked.double(1.e-50),
  "fplo":cms.untracked.double(0),
  "chi2cutoff":cms.double(2.5), 
  "verboseAnalyzer":cms.untracked.bool(False),
  "verboseProducer":cms.untracked.bool(False),
  "newClusterizer":cms.untracked.bool(False),
  "verboseClusterizer":cms.untracked.bool(False),
  "verboseClusterizer2D":cms.untracked.bool(False),
  "zdumpcenter" : cms.untracked.double(0.0),
  "zdumpwidth" : cms.untracked.double(20.0),
  "reco":cms.untracked.bool(False),
  "autodump":cms.untracked.int32(0),
  "nDump":cms.untracked.int32(0),
#  "mintrkweight":cms.untracked.double(0.5),
  "uniquetrkweight":cms.double(0.8),
  "uniquetrkminp":cms.double(0.0),
  "zmerge":cms.double(1.e-2),
  "coolingFactor":cms.double(0.6),
  "Tmin": cms.double(0), # 0 = not set, use default
  "Tpurge":cms.double(0),
  "Tstop": cms.double(0),
  "vertexSize" : cms.double(0.006),
  "vertexSizeTime" : cms.double(0.008),
  "d0CutOff" : cms.double(3.),
  "dzCutOff" : cms.double(3.),
  "zrange" : cms.double(4),
  "convergence_mode" : cms.int32(0),
  "delta_lowT" : cms.double(1.e-3),
  "delta_highT" : cms.double(1.e-2),
# track selection
  "maxNormalizedChi2":cms.double(10.0),
  "minPixelLayersWithHits":cms.int32(2),
  "minSiliconLayersWithHits":cms.int32(5),
  "maxD0Significance":cms.double(4.0), 
  "minPt":cms.double(0.0),
  "maxEta":cms.double(2.4),
  "trackQuality":cms.string("any"),
# track selection, experimental
  "maxDzError":cms.double(1.0),
  "maxD0Error":cms.double(1.0),
# vertex selection
  "minNdof": cms.double( 0.0 ),
  "trackTimeQualityThreshold" : cms.double(0.8),
#
  "purge_method" : cms.untracked.int32(0),
  "split_method" : cms.untracked.int32(0)
}

# temporary fix, should not be needed
args=[]
for a in sys.argv:
    args += a.split()
        

for a in args:


    try:
        key, value = a.split("=")
    except ValueError:
        continue


    if key == "nevent" or key == "events":
        nevent = int(value)

    elif key == "input":
        inputfile = value
        if inputfile.endswith(".root"):
            source_files = cms.untracked.vstring(inputfile)
        else:
            files = tuple( l.strip() for l in open(inputfile).readlines() )
            source_files = cms.untracked.vstring(*files)

    elif key == "json":
        json = None if value in ("None","") else value

    elif key == "era":
        era = value

    elif key == "4D":
        parameters["4D"] = cms.untracked.bool(value == "True")
        if (value=="True"):
            vcollections= cms.VInputTag(
                "offlinePrimaryVertices",
                "offlinePrimaryVertices4D"
            )
            era = "Phase2"
    elif key == "redopv" or key == "redoPV":
        DO_VTX_RECO = (value in ("True","yes","Yes") )

    elif key == "zdump":
        debug = True
        try:
            zdumpcenter,zdumpwidth = value.split(":")
            zdumpcenter = float(zdumpcenter)
            zdumpwidth = float(zdumpwidth)
            parameters["zdumpcenter"] = cms.untracked.double(zdumpcenter)
            parameters["zdumpwidth"] = cms.untracked.double(zdumpwidth)
            nevent = 1
        except ValueError:
            nevent = int(value)

    elif key == "debug":
        debug = True
        try:
            run,event,zdumpcenter,zdumpwidth = value.split(":")
            eventsToProcess = '%s:%s'%(run,event)
            zdumpcenter = float(zdumpcenter)
            zdumpwidth = float(zdumpwidth)
            parameters["zdumpcenter"] = cms.untracked.double(zdumpcenter)
            parameters["zdumpwidth"] = cms.untracked.double(zdumpwidth)
            nevent = 1
        except ValueError:
            nevent = int(value)
        
    elif key == "verbose":
        if (value in ("True","yes","Yes") ):
            parameters["verboseClusterizer"] = cms.untracked.bool(True)
            parameters["verboseAnalyzer"] = cms.untracked.bool(True)
            parameters["verboseProducer"] = cms.untracked.bool(True)
    elif key == "new":
        if (value in ("True","yes","Yes") ):
            parameters["newClusterizer"] = cms.untracked.bool(True)
    elif key == "info":
        info  = value
    elif key in parameters.keys():
        typename = parameters[key].configTypeName()
        print key, value, type(parameters[key]),"'%s'"%typename
        if typename == "double":
            parameters[key] = cms.double( float(value) )
        elif type(parameters[key])==type(cms.string("")):
            parameters[key] = cms.string( value )
        elif typename == "untracked double":
            parameters[key] = cms.untracked.double( float(value) )
        elif typename == "int32":
            parameters[key] = cms.int32( int(value) )
        elif typename == "untracked int32":
            parameters[key] = cms.untracked.int32( int(value) )
        elif typename == "bool":
            parameters[key] = cms.bool( value == "True")
        elif typename == "untracked bool":
            parameters[key] = cms.untracked.bool( value == "True")
    else:
        print "!! pvt_cfg.py  :  unknown key ",key


print "pvt_cfg.py"
print "args           ",args
print "events         ",nevent
print "info           ",info
print "json           ",json
print "era            ",era
print "DO_VTX_RECO    ",DO_VTX_RECO
print "newClusterizer ",parameters["newClusterizer"]
print "all parameters"
print parameters
if len(source_files) == 0:
    print "empty source file list!"


if era == "Phase2":
    process = cms.Process("RERECO", eras.Phase2)
    autotag = "auto:phase2_realistic"
elif era == "Run3":
    process = cms.Process("RERECO", eras.Run3)
    autotag = "111X_mcRun3_2021_realistic_v4"
    parameters["4D"] = cms.untracked.bool(False)
else:
    print "unknown era",era
    sys.exit(1)

# edmProvDump root:....
#ESSource: GlobalTag RECO
#  @module_label: string tracked  = 'GlobalTag'
#  globaltag: string tracked  = '111X_mcRun3_2021_realistic_v4'



# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000   # don't be too noisy
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('CommonTools.ParticleFlow.EITopPAG_cff')
process.load('Configuration.StandardSequences.AlCaRecoStreamsMC_cff')
process.load('Configuration.StandardSequences.Validation_cff')
process.load('DQMOffline.Configuration.DQMOfflineMC_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

if era == "Run3":
    process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
elif era== "Phase2":
    process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')  

if True and use_tp:
    print "heyho", "*"*40
    process.load("Validation.RecoTrack.TrackValidation_cff")
    process.theTruth = cms.Sequence(
        process.tpClusterProducer +
        process.quickTrackAssociatorByHits +
        process.trackingParticleRecoTrackAsssociation
    )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(nevent))


# Input source
process.source = cms.Source("PoolSource",
    fileNames = source_files,
    secondaryFileNames = cms.untracked.vstring()
    #eventsToProcess = cms.untracked.VEventRange()
)
if not eventsToProcess == "":
    process.source.eventsToProcess = cms.untracked.VEventRange( eventsToProcess )


if json is not None  and os.path.exists(json) and (len(source_files)==0 or source_files[0].find("SIM")<0):
    process.source.lumisToProcess = LumiList.LumiList(filename = json).getVLuminosityBlockRange()


# no HLT for MC
if len(source_files)>0 and source_files[0].find("SIM") >= 0:
    use_hlt = False



#process.Tracer = cms.Service("Tracer")

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True)
)



from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, autotag, '')

if DO_VTX_RECO:
    process.Timing = cms.Service("Timing", summaryOnly = cms.untracked.bool(True)  )

process.load("RecoLocalCalo.Configuration.hcalLocalReco_cff")
process.load("RecoLocalCalo.HcalRecProducers.HBHEIsolatedNoiseReflagger_cfi")
process.load("RecoVertex.Configuration.RecoVertex_cff")

if DO_VTX_RECO:
    process.caloTowerForTrk.hbheInput = 'reducedHcalRecHits:hbhereco'
    process.caloTowerForTrk.hoInput = 'reducedHcalRecHits:horeco'
    process.caloTowerForTrk.hfInput = 'reducedHcalRecHits:hfreco'
    process.caloTowerForTrk.ecalInputs = ['reducedEcalRecHitsEB', 'reducedEcalRecHitsEE']

process.ak4CaloJetsForTrk.srcPVs = 'unsortedOfflinePrimaryVertices'
process.vertexreco.remove(process.generalV0Candidates)
process.vertexreco.remove(process.inclusiveVertexing)

    
tkFilterParameters = cms.PSet(
    algorithm=cms.string('filter'),
    maxNormalizedChi2 = parameters["maxNormalizedChi2"],
    minPixelLayersWithHits = parameters["minPixelLayersWithHits"],
    minSiliconLayersWithHits = parameters["minSiliconLayersWithHits"],
    maxD0Significance = parameters["maxD0Significance"], 
    maxD0Error = parameters["maxD0Error"], 
    maxDzError = parameters["maxDzError"], 
    minPt = parameters["minPt"],
    maxEta = parameters["maxEta"],
    trackQuality = parameters["trackQuality"]
    )

if DO_VTX_RECO:
    process.PV = cms.Path(process.vertexreco)

    if False:
        print "1d"
        print process.unsortedOfflinePrimaryVertices.TkFilterParameters
        print process.unsortedOfflinePrimaryVertices.TkClusParameters
        print "noPID"
        print process.unsortedOfflinePrimaryVertices4DnoPID.TkFilterParameters
        print process.unsortedOfflinePrimaryVertices4DnoPID.TkClusParameters
        print "PID"
        print process.unsortedOfflinePrimaryVertices4D.TkFilterParameters
        print process.unsortedOfflinePrimaryVertices4D.TkClusParameters

    process.unsortedOfflinePrimaryVertices.verbose = parameters["verboseProducer"]
    #process.unsortedOfflinePrimaryVertices.new = parameters["newClusterizer"]
    process.unsortedOfflinePrimaryVertices.TkClusParameters.TkDAClusParameters.verbose =  parameters["verboseClusterizer"]
    process.unsortedOfflinePrimaryVertices.TkClusParameters.TkDAClusParameters.zrange = parameters["zrange"]
    process.unsortedOfflinePrimaryVertices.TkClusParameters.TkDAClusParameters.convergence_mode = parameters["convergence_mode"]
    process.unsortedOfflinePrimaryVertices.TkClusParameters.TkDAClusParameters.delta_lowT = parameters["delta_lowT"]
    process.unsortedOfflinePrimaryVertices.TkClusParameters.TkDAClusParameters.delta_highT = parameters["delta_highT"]
    process.unsortedOfflinePrimaryVertices.TkClusParameters.TkDAClusParameters.vertexSize = parameters["vertexSize"]
    process.unsortedOfflinePrimaryVertices.TkClusParameters.TkDAClusParameters.uniquetrkminp = parameters["uniquetrkminp"]
    process.unsortedOfflinePrimaryVertices.TkClusParameters.TkDAClusParameters.purge_method = parameters["purge_method"]
    process.unsortedOfflinePrimaryVertices.TkClusParameters.TkDAClusParameters.split_method = parameters["split_method"]
    process.unsortedOfflinePrimaryVertices.TkClusParameters.TkDAClusParameters.zmerge = parameters["zmerge"]

    process.unsortedOfflinePrimaryVertices.TkFilterParameters = tkFilterParameters.clone()
    #
    for p in (process.unsortedOfflinePrimaryVertices4DnoPID, process.unsortedOfflinePrimaryVertices4D):
        p.verbose = parameters["verboseProducer"]
        p.TrackTimeQualityThreshold = parameters["trackTimeQualityThreshold"]
        p.TrackTimeQualityMapLabel = cms.InputTag("mtdTrackQualityMVA:mtdQualMVA")
        #p.new = parameters["newClusterizer"]
        p.TkClusParameters.TkDAClusParameters.verbose =  parameters["verboseClusterizer2D"]
        p.TkClusParameters.TkDAClusParameters.zdumpcenter = parameters["zdumpcenter"]
        p.TkClusParameters.TkDAClusParameters.zdumpwidth = parameters["zdumpwidth"]
        p.TkClusParameters.TkDAClusParameters.zrange = parameters["zrange"]
        p.TkClusParameters.TkDAClusParameters.convergence_mode = parameters["convergence_mode"]
        p.TkClusParameters.TkDAClusParameters.delta_lowT = parameters["delta_lowT"]
        p.TkClusParameters.TkDAClusParameters.delta_highT = parameters["delta_highT"]
        # remember to put any parameter used here into msub.py
        for par_name in ("Tmin", "Tpurge", "zmerge", "vertexSizeTime", "uniquetrkminp"):
            if parameters[par_name] > 0:
                setattr(p.TkClusParameters.TkDAClusParameters, par_name, parameters[par_name])

        p.TkFilterParameters = tkFilterParameters.clone()
    




# analysis
process.oldVertexAnalysis = cms.EDAnalyzer("PrimaryVertexAnalyzer4PU",
    info=  cms.untracked.string(info),
    f4D = parameters["4D"],
    beamSpot = cms.InputTag('offlineBeamSpot'),
    simG4 = cms.InputTag("g4SimHits"),
    outputFile = cms.untracked.string("pv.root"),
    verbose = parameters["verboseAnalyzer"],
    veryverbose = cms.untracked.bool(False),
    recoTrackProducer = cms.untracked.string("generalTracks"),
    zmatch = cms.untracked.double(0.05),
    autodump = parameters["autodump"],
    nDump = parameters["nDump"],
    RECO = parameters["reco"],
    track_timing = cms.untracked.bool(True),
    TkFilterParameters = tkFilterParameters,
    trackingParticleCollection = cms.untracked.InputTag("mix", "MergedTrackTruth"),
    trackingVertexCollection = cms.untracked.InputTag("mix", "MergedTrackTruth"),
    trackAssociatorMap = cms.untracked.InputTag("trackingParticleRecoTrackAsssociation"),
    TrackTimesLabel = cms.untracked.InputTag("tofPID4DnoPID:t0safe"),  # as opposed to "tofPID:t0safe"
    TrackTimeResosLabel = cms.untracked.InputTag("tofPID4DnoPID:sigmat0safe"),
    TrackTimeQualityMapLabel = cms.untracked.InputTag("mtdTrackQualityMVA:mtdQualMVA"),
    TrackTimeQualityThreshold = cms.untracked.double( parameters["trackTimeQualityThreshold"].value()),
    vertexAssociator = cms.untracked.InputTag("VertexAssociatorByPositionAndTracks"),
    useVertexFilter = cms.untracked.bool(False),
    compareCollections = cms.untracked.int32(0),
    vertexRecoCollections = vcollections
)


#---------------------------------------------------------------------------------------------------
# Trigger Selection

process.hltSelection = cms.EDFilter("HLTHighLevel",
                               TriggerResultsTag =
                               cms.InputTag("TriggerResults","","HLT"),
                               HLTPaths = cms.vstring(["HLT_ZeroBias_v*"]),
                               eventSetupPathsKey = cms.string(''),
                               andOr = cms.bool(True),
                               throw = cms.bool(False)  )


if use_tp:
    process.analyze =  cms.Path( process.theTruth * process.oldVertexAnalysis )
elif use_hlt:
    print "using hlt seletion"
    process.analyze =  cms.Path( process.hltSelection * process.oldVertexAnalysis )
else:
    process.analyze =  cms.Path( process.oldVertexAnalysis )


print "============================================="
print process.source
print "============================================="

if False:
    process.content = cms.EDAnalyzer("EventContentAnalyzer")
    process.dump = cms.Path(process.content)
    process.schedule = cms.Schedule(process.dump)

else:

    # Schedule definition
    if DO_VTX_RECO:
        # Schedule definition
        print "running vertex reco and analysis"
        process.schedule = cms.Schedule(process.PV, process.analyze)
    else:
        print "running analysis only"
        process.schedule = cms.Schedule(process.analyze)
