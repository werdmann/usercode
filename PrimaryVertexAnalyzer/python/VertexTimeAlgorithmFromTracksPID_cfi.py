import FWCore.ParameterSet.Config as cms

VertexTimeAlgorithmFromTracksPID = cms.PSet(
    ComponentType = cms.string('VertexTimeAlgorithmFromTracksPID'),

    trackMTDTimeVMapTag = cms.InputTag('trackExtenderWithMTD:generalTracktmtd'),
    trackMTDTimeErrorVMapTag = cms.InputTag('trackExtenderWithMTD:generalTracksigmatmtd'),
    trackMTDTimeQualityVMapTag = cms.InputTag('mtdTrackQualityMVA:mtdQualMVA'),
    trackMTDMomentumVMapTag = cms.InputTag('trackExtenderWithMTD:generalTrackp'),
    trackMTDPathLengthVMapTag = cms.InputTag('trackExtenderWithMTD:generalTrackPathLength'),

    minTrackVtxWeight = cms.double(0.5),
    minTrackTimeQuality = cms.double(0.),

    massPion = cms.double(0.139570),
    massKaon = cms.double(0.493677),
    massProton = cms.double(0.938272),

    probPion = cms.double(0.7),
    probKaon = cms.double(0.2),
    probProton = cms.double(0.1),

    coolingFactor = cms.double(0.5)
)
