#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ValidatedPluginMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "VertexTimeAlgorithmFromTracksPID.h"
#include "VertexTimeAlgorithmBaseFactory.h"

#ifdef PVTX_DEBUG
#define LOG edm::LogPrint("VertexTimeAlgorithmFromTracksPID")
#else
#define LOG LogDebug("VertexTimeAlgorithmFromTracksPID")
#endif

VertexTimeAlgorithmFromTracksPID::VertexTimeAlgorithmFromTracksPID(edm::ParameterSet const& iConfig, edm::ConsumesCollector& iCC) :
  VertexTimeAlgorithmBase(iConfig, iCC),
  trackMTDTimeToken_(iCC.consumes(iConfig.getParameter<edm::InputTag>("trackMTDTimeVMapTag"))),
  trackMTDTimeErrorToken_(iCC.consumes(iConfig.getParameter<edm::InputTag>("trackMTDTimeErrorVMapTag"))),
  trackMTDTimeQualityToken_(iCC.consumes(iConfig.getParameter<edm::InputTag>("trackMTDTimeQualityVMapTag"))),
  trackMTDMomentumToken_(iCC.consumes(iConfig.getParameter<edm::InputTag>("trackMTDMomentumVMapTag"))),
  trackMTDPathLengthToken_(iCC.consumes(iConfig.getParameter<edm::InputTag>("trackMTDPathLengthVMapTag"))),
  minTrackVtxWeight_(iConfig.getParameter<double>("minTrackVtxWeight")),
  minTrackTimeQuality_(iConfig.getParameter<double>("minTrackTimeQuality")),
  massPion_(iConfig.getParameter<double>("massPion")),
  massKaon_(iConfig.getParameter<double>("massKaon")),
  massProton_(iConfig.getParameter<double>("massProton")),
  probPion_(iConfig.getParameter<double>("probPion")),
  probKaon_(iConfig.getParameter<double>("probKaon")),
  probProton_(iConfig.getParameter<double>("probProton")),
  coolingFactor_(iConfig.getParameter<double>("coolingFactor")) {}

void VertexTimeAlgorithmFromTracksPID::fillPSetDescription(edm::ParameterSetDescription& iDesc) {
  VertexTimeAlgorithmBase::fillPSetDescription(iDesc);

  iDesc.add<edm::InputTag>("trackMTDTimeVMapTag", edm::InputTag("trackExtenderWithMTD:generalTracktmtd"))->setComment("");
  iDesc.add<edm::InputTag>("trackMTDTimeErrorVMapTag", edm::InputTag("trackExtenderWithMTD:generalTracksigmatmtd"))->setComment("");
  iDesc.add<edm::InputTag>("trackMTDTimeQualityVMapTag", edm::InputTag("mtdTrackQualityMVA:mtdQualMVA"))->setComment("");
  iDesc.add<edm::InputTag>("trackMTDMomentumVMapTag", edm::InputTag("trackExtenderWithMTD:generalTrackp"))->setComment("");
  iDesc.add<edm::InputTag>("trackMTDPathLengthVMapTag", edm::InputTag("trackExtenderWithMTD:generalTrackPathLength"))->setComment("");

  iDesc.add<double>("minTrackVtxWeight", 0.5)->setComment("");
  iDesc.add<double>("minTrackTimeQuality", 0.)->setComment("");

  iDesc.add<double>("massPion", 0.139570)->setComment("");
  iDesc.add<double>("massKaon", 0.493677)->setComment("");
  iDesc.add<double>("massProton", 0.938272)->setComment("");

  iDesc.add<double>("probPion", 0.7)->setComment("");
  iDesc.add<double>("probKaon", 0.2)->setComment("");
  iDesc.add<double>("probProton", 0.1)->setComment("");

  iDesc.add<double>("coolingFactor", 0.5)->setComment("");
}

void VertexTimeAlgorithmFromTracksPID::setEvent(edm::Event& iEvent, edm::EventSetup const&) {
  // additional collections required for vertex-time calculation
  trackMTDTimes_ = iEvent.get(trackMTDTimeToken_);
  trackMTDTimeErrors_ = iEvent.get(trackMTDTimeErrorToken_);
  trackMTDTimeQualities_ = iEvent.get(trackMTDTimeQualityToken_);
  trackMTDMomenta_ = iEvent.get(trackMTDMomentumToken_);
  trackMTDPathLengths_ = iEvent.get(trackMTDPathLengthToken_);
}

float VertexTimeAlgorithmFromTracksPID::trackTime(float const mass, float const mtdTime, float const mtdPathLength, float const mtdMomentum) const {
  // speed of light, c = 29.9792458 cm/ns
  return (mtdTime - mtdPathLength * std::sqrt(1.f + mass*mass/mtdMomentum/mtdMomentum) / 29.9792458f);
}

bool VertexTimeAlgorithmFromTracksPID::vertexTime(float& vtxTime, float& vtxTimeError, reco::Vertex const& vtx) const {
  if(vtx.tracksSize() == 0){
    return false;
  }

  auto const vtxTime_init = vtxTime;
  auto const vtxTimeError_init = vtxTimeError;

  double tsum = 0;
  double wsum = 0;
  double w2sum = 0;

  double const a[3] = {probPion_, probKaon_, probProton_};

  LOG << "vertexTimeFromTracks: vtx x=" << vtx.x() << " y=" << vtx.y() << " z=" << vtx.z() << " t=" << vtx.t();

  std::vector<TrackInfo> v_trackInfo;
  v_trackInfo.reserve(vtx.tracksSize());

  // initial guess
  for (auto trk = vtx.tracks_begin(); trk != vtx.tracks_end(); ++trk) {
    auto const trkWeight = vtx.trackWeight(*trk);
    if (trkWeight > minTrackVtxWeight_) {

      auto const trkTimeQuality = trackMTDTimeQualities_[*trk];

      if (trkTimeQuality >= minTrackTimeQuality_) {

        auto const trkTime = trackMTDTimes_[*trk];
        auto const trkTimeError = trackMTDTimeErrors_[*trk];
        auto const trkPathLength = trackMTDPathLengths_[*trk];
        auto const trkMomentum = trackMTDMomenta_[*trk];

        v_trackInfo.emplace_back();
        auto& trkInfo = v_trackInfo.back();

        trkInfo.trkWeight = trkWeight;
        trkInfo.trkTime = trkTime;
        trkInfo.trkTimeError = trkTimeError;

        if(trkPathLength > 0){
          trkInfo.trkTimeHyp[0] = trackTime(massPion_, trkTime, trkPathLength, trkMomentum);
          trkInfo.trkTimeHyp[1] = trackTime(massKaon_, trkTime, trkPathLength, trkMomentum);
          trkInfo.trkTimeHyp[2] = trackTime(massProton_, trkTime, trkPathLength, trkMomentum);
        } else {
          trkInfo.trkTimeHyp[0] = 0.f;
          trkInfo.trkTimeHyp[1] = 0.f;
          trkInfo.trkTimeHyp[2] = 0.f;
        }

        auto const wgt = trkWeight / (trkTimeError*trkTimeError);
        wsum += wgt;

        for(uint j = 0; j < 3; ++j) {
          tsum += wgt * trkInfo.trkTimeHyp[j] * a[j];
        }

        LOG << "vertexTimeFromTracks:     track"
            << " pt=" << (*trk)->pt() << " eta=" << (*trk)->eta() << " phi=" << (*trk)->phi()
            << " vtxWeight=" << trkWeight << " time=" << trkTime << " timeError=" << trkTimeError
            << " timeQuality=" << trkTimeQuality << " pathLength=" << trkPathLength << " momentum=" << trkMomentum
            << " timeHyp[pion]=" << trkInfo.trkTimeHyp[0] << " timeHyp[kaon]=" << trkInfo.trkTimeHyp[1]
            << " timeHyp[proton]=" << trkInfo.trkTimeHyp[2];
      }
    }
  }

  if (wsum > 0) {

    LOG << "vertexTimeFromTracks:   wsum = " << wsum << " tsum = " << tsum << " t0 = " << (wsum > 0 ? tsum/wsum : 0) << " trec = " << vtx.t();

    auto t0 = tsum / wsum;
    auto beta = 1./256.;
    int nit = 0;
    while ( (nit++) < 100 ) {
      tsum = 0;
      wsum = 0;
      w2sum = 0;

      for (auto const& trkInfo : v_trackInfo) {
        double dt = trkInfo.trkTimeError;
        double e[3] = {0,0,0};
        double Z = exp(-beta * 0.5* 3.* 3.);
        for(unsigned int j = 0; j < 3; j++){
          auto const tpull =  (trkInfo.trkTimeHyp[j] - t0) / dt;
          e[j] = exp(- 0.5 * beta * tpull * tpull);
          Z += a[j] * e[j];
        }

        double wsum_trk = 0;
        for(uint j = 0; j < 3; j++){
          double wt = a[j] * e[j] / Z;
          double w = wt * trkInfo.trkWeight / (dt * dt);
          wsum_trk += w;
          tsum += w * trkInfo.trkTimeHyp[j];
        }

        wsum += wsum_trk;
        w2sum += wsum_trk * wsum_trk * (dt * dt) / trkInfo.trkWeight;
      }

      if (wsum < 1e-10) {
        LOG << "vertexTimeFromTracks:   failed while iterating";
        return false;
      }

      vtxTime = tsum / wsum;

      LOG << "vertexTimeFromTracks:   iteration=" << nit << ", T= " << 1/beta << ", t=" << vtxTime << ", t-t0=" << vtxTime-t0;

      if ((std::abs(vtxTime - t0) < 1e-4 / std::sqrt(beta)) and beta >= 1.) {

        vtxTimeError = std::sqrt(w2sum) / wsum;

        LOG << "vertexTimeFromTracks:   tfit = " << vtxTime << " +/- " << vtxTimeError << " trec = " << vtx.t() << ", iteration=" << nit;

        return true;
      }

      if ((std::abs(vtxTime - t0) < 1e-3) and beta < 1.) {
        beta = std::min(1., beta / coolingFactor_);
      }

      t0 = vtxTime;
    }

    LOG << "vertexTimeFromTracks: failed to converge";
  }
  else {
    LOG << "vertexTimeFromTracks: has no track timing info";
  }

  vtxTime = vtxTime_init;
  vtxTimeError = vtxTimeError_init;

  return false;
}

DEFINE_EDM_VALIDATED_PLUGIN(VertexTimeAlgorithmBaseFactory, VertexTimeAlgorithmFromTracksPID, "VertexTimeAlgorithmFromTracksPID");
