#ifndef usercode_PrimaryVertexAnalyzer_VertexTimeAlgorithmFromTracksPID_h
#define usercode_PrimaryVertexAnalyzer_VertexTimeAlgorithmFromTracksPID_h

#include "VertexTimeAlgorithmBase.h"

#include "FWCore/Utilities/interface/EDGetToken.h"
#include "DataFormats/Common/interface/ValueMap.h"

class VertexTimeAlgorithmFromTracksPID : public VertexTimeAlgorithmBase {
 public:
  VertexTimeAlgorithmFromTracksPID(const edm::ParameterSet& conf, edm::ConsumesCollector& iC);
  ~VertexTimeAlgorithmFromTracksPID() override = default;

  static void fillPSetDescription(edm::ParameterSetDescription& iDesc);

  void setEvent(edm::Event& iEvent, edm::EventSetup const& iSetup) override;

  bool vertexTime(float& vtxTime, float& vtxTimeError, reco::Vertex const& vtx) const override;

protected:
  float trackTime(float const mass, float const mtdTime, float const mtdPathLength, float const mtdMomentum) const;

  struct TrackInfo {
    double trkWeight;
    double trkTime;
    double trkTimeError;
    double trkTimeHyp[3];
  };

  edm::EDGetTokenT<edm::ValueMap<float>> const trackMTDTimeToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> const trackMTDTimeErrorToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> const trackMTDTimeQualityToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> const trackMTDMomentumToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> const trackMTDPathLengthToken_;

  double const minTrackVtxWeight_;
  double const minTrackTimeQuality_;
  double const massPion_;
  double const massKaon_;
  double const massProton_;
  double const probPion_;
  double const probKaon_;
  double const probProton_;
  double const coolingFactor_;

  edm::ValueMap<float> trackMTDTimes_;
  edm::ValueMap<float> trackMTDTimeErrors_;
  edm::ValueMap<float> trackMTDTimeQualities_;
  edm::ValueMap<float> trackMTDMomenta_;
  edm::ValueMap<float> trackMTDPathLengths_;
};

#endif
