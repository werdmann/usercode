#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include <memory>
#include <utility>

#ifdef PVTX_DEBUG
#define LOG edm::LogPrint("VertexTimeAnalyzer")
#else
#define LOG LogDebug("VertexTimeAnalyzer")
#endif

class VertexTimeProducer : public edm::global::EDProducer<> {
public:
  explicit VertexTimeProducer(edm::ParameterSet const&);
  ~VertexTimeProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

protected:
  void produce(edm::StreamID, edm::Event&, edm::EventSetup const&) const override;

  float trackTime(float const mass, float const mtdTime, float const mtdPathLength, float const mtdMomentum) const;

  bool vertexTimeFromTracks(float& vtxTime, float& vtxTimeError, reco::Vertex const& vtx, edm::ValueMap<float> const& vmapTrackMTDTimes, edm::ValueMap<float> const& vmapTrackMTDTimeErrors, edm::ValueMap<float> const& vmapTrackMTDTimeQualities, edm::ValueMap<float> const& vmapTrackMTDMomenta, edm::ValueMap<float> const& vmapTrackMTDPathLengths) const;

  struct TrackInfo {
    double trkWeight;
    double trkTime;
    double trkTimeError;
    double trkTimeHyp[3];
  };

  edm::EDGetTokenT<reco::VertexCollection> const vtxsToken_;
  bool const produceVertices_;

  edm::EDGetTokenT<edm::ValueMap<float>> const trackMTDTimeToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> const trackMTDTimeErrorToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> const trackMTDTimeQualityToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> const trackMTDMomentumToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> const trackMTDPathLengthToken_;
  float const minTrackVtxWeight_;
  float const minTrackTimeQuality_;
};

VertexTimeProducer::VertexTimeProducer(edm::ParameterSet const& iConfig) :
  vtxsToken_(consumes(iConfig.getParameter<edm::InputTag>("vertices"))),
  produceVertices_(iConfig.getParameter<bool>("produceVertices")),

  trackMTDTimeToken_(consumes(edm::InputTag("trackExtenderWithMTD:generalTracktmtd"))),
  trackMTDTimeErrorToken_(consumes(edm::InputTag("trackExtenderWithMTD:generalTracksigmatmtd"))),
  trackMTDTimeQualityToken_(consumes(edm::InputTag("mtdTrackQualityMVA:mtdQualMVA"))),
  trackMTDMomentumToken_(consumes(edm::InputTag("trackExtenderWithMTD:generalTrackp"))),
  trackMTDPathLengthToken_(consumes(edm::InputTag("trackExtenderWithMTD:generalTrackPathLength"))),
  minTrackVtxWeight_(iConfig.getParameter<double>("minTrackVtxWeight")),
  minTrackTimeQuality_(iConfig.getParameter<double>("minTrackTimeQuality")) {

  produces<edm::ValueMap<float>>("time");
  produces<edm::ValueMap<float>>("timeError");

  if(produceVertices_){
    produces<reco::VertexCollection>();
  }
}

void VertexTimeProducer::produce(edm::StreamID, edm::Event& iEvent, edm::EventSetup const& iSetup) const {

  auto const vtxsHandle = iEvent.getHandle(vtxsToken_);

  // additional collections required for vertex-time calculation
  auto const& trackMTDTimes = iEvent.get(trackMTDTimeToken_);
  auto const& trackMTDTimeErrors = iEvent.get(trackMTDTimeErrorToken_);
  auto const& trackMTDTimeQualities = iEvent.get(trackMTDTimeQualityToken_);
  auto const& trackMTDMomenta = iEvent.get(trackMTDMomentumToken_);
  auto const& trackMTDPathLengths = iEvent.get(trackMTDPathLengthToken_);

  std::vector<float> vtxTimes;
  vtxTimes.reserve(vtxsHandle->size());

  std::vector<float> vtxTimeErrors;
  vtxTimeErrors.reserve(vtxsHandle->size());

  auto newVtxs = std::make_unique<reco::VertexCollection>();

  if(produceVertices_){
    newVtxs->reserve(vtxsHandle->size());
  }

  for (auto const& vtx : *vtxsHandle) {
    auto vtxTime(0.f), vtxTimeError(-1.f);
    vertexTimeFromTracks(vtxTime, vtxTimeError, vtx,
      trackMTDTimes, trackMTDTimeErrors, trackMTDTimeQualities, trackMTDMomenta, trackMTDPathLengths);

    LOG << "Vertex #" << vtxTimes.size() << ": x=" << vtx.x() << ", y=" << vtx.y() << ": z=" << vtx.z() << ", t=" << vtxTime << ", tError=" << vtxTimeError;

    // fill values for time and timeError ValueMaps
    vtxTimes.emplace_back(vtxTime);
    vtxTimeErrors.emplace_back(vtxTimeError);

    // fill collection of output vertices
    if(produceVertices_){
      // https://github.com/cms-sw/cmssw/blob/cbf578c08dd47e201eec8b64f17ae4c646e6b1a6/RecoVertex/VertexPrimitives/src/TransientVertex.cc#L301
      if (vtx.isValid()){
        auto err = vtx.covariance4D();
        err(3, 3) = vtxTimeError*vtxTimeError;
        auto vtx_new = reco::Vertex(vtx.position(), err, vtxTime, vtx.chi2(), vtx.ndof(), vtx.tracksSize());

        // add tracks to vertex
        if (vtx.hasRefittedTracks()) {
          for(auto const& trk : vtx.refittedTracks()){
            auto const& trk_orig = vtx.originalTrack(trk);
            vtx_new.add(trk_orig, trk, vtx.trackWeight(trk_orig));
          }
        } else {
          for(auto trk = vtx.tracks_begin(); trk != vtx.tracks_end(); ++trk){
            vtx_new.add(*trk, vtx.trackWeight(*trk));
          }
        }

        newVtxs->emplace_back(vtx_new);
      }
      else {
        // use default reco::Vertex ctor for invalid vertex
        newVtxs->emplace_back();
      }
    }
  }

  auto vtxsTimeVMap = std::make_unique<edm::ValueMap<float>>();
  edm::ValueMap<float>::Filler vtxsTimeVMapFiller(*vtxsTimeVMap);
  vtxsTimeVMapFiller.insert(vtxsHandle, vtxTimes.begin(), vtxTimes.end());
  vtxsTimeVMapFiller.fill();
  iEvent.put(std::move(vtxsTimeVMap), "time");

  auto vtxsTimeErrorVMap = std::make_unique<edm::ValueMap<float>>();
  edm::ValueMap<float>::Filler vtxsTimeErrorVMapFiller(*vtxsTimeErrorVMap);
  vtxsTimeErrorVMapFiller.insert(vtxsHandle, vtxTimeErrors.begin(), vtxTimeErrors.end());
  vtxsTimeErrorVMapFiller.fill();
  iEvent.put(std::move(vtxsTimeErrorVMap), "timeError");

  if(produceVertices_){
    iEvent.put(std::move(newVtxs));
  }
}

void VertexTimeProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("vertices", edm::InputTag("offlinePrimaryVertices"))->setComment("input collection of vertices");
  desc.add<bool>("produceVertices", false)->setComment("produce a collection of vertices with updated time and timeError");
  desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"))->setComment("input collection of tracks");
  desc.add<double>("minTrackVtxWeight", -1.)->setComment("");
  desc.add<double>("minTrackTimeQuality", -1.)->setComment("");
  descriptions.addWithDefaultLabel(desc);
}

float VertexTimeProducer::trackTime(float const mass, float const mtdTime, float const mtdPathLength, float const mtdMomentum) const {
//!!  // speed of light, c = 29.9792458 cm/ns
//!!  return (mtdTime - mtdPathLength * std::sqrt(1.f + mass*mass/mtdMomentum/mtdMomentum) / 29.9792458f);
  double gammasq= 1. + mtdMomentum * mtdMomentum / (mass * mass);
  double v = 2.99792458e1 * std::sqrt(1. - 1. / gammasq);  // cm / ns
  return mtdTime - mtdPathLength / v;
}

bool VertexTimeProducer::vertexTimeFromTracks(float& t, float& tError, reco::Vertex const& vtx, edm::ValueMap<float> const& vmapTrackMTDTimes, edm::ValueMap<float> const& vmapTrackMTDTimeErrors, edm::ValueMap<float> const& vmapTrackMTDTimeQualities, edm::ValueMap<float> const& vmapTrackMTDMomenta, edm::ValueMap<float> const& vmapTrackMTDPathLengths) const {
  if(vtx.tracksSize() == 0){
    return false;
  }

  auto const t_init = t;
  auto const tError_init = tError;

  double tsum = 0;
  double wsum = 0;
  double w2sum = 0;

  double const a[3] = {0.7, 0.2, 0.1};
  double const cooling_factor = 0.5;
  double const minquality = 0.0;
  double const mintrkweight = 0.5;

  double const massPion = 0.139570;
  double const massKaon = 0.493677;
  double const massProton = 0.938272;

  bool const verbose = true;

  if (verbose) {
    LOG << "vertexTimeFromTracks: vtx x=" << vtx.x() << " y=" << vtx.y() << " z=" << vtx.z() << " t=" << vtx.t();
  }

  std::vector<TrackInfo> v_trackInfo;
  v_trackInfo.reserve(vtx.tracksSize());

  // initial guess
  for (auto trk = vtx.tracks_begin(); trk != vtx.tracks_end(); ++trk) {
    auto const trkWeight = vtx.trackWeight(*trk);
    if (trkWeight > mintrkweight) {

      auto const trkTimeQuality = vmapTrackMTDTimeQualities[*trk];

      if (trkTimeQuality >= minquality) {

        auto const trkTime = vmapTrackMTDTimes[*trk];
        auto const trkTimeError = vmapTrackMTDTimeErrors[*trk];
        auto const trkPathLength = vmapTrackMTDPathLengths[*trk];
        auto const trkMomentum = vmapTrackMTDMomenta[*trk];

        v_trackInfo.emplace_back();
        auto& trkInfo = v_trackInfo.back();

        trkInfo.trkWeight = trkWeight;
        trkInfo.trkTime = trkTime;
        trkInfo.trkTimeError = trkTimeError;

        if(trkPathLength > 0){
          trkInfo.trkTimeHyp[0] = trackTime(massPion, trkTime, trkPathLength, trkMomentum);
          trkInfo.trkTimeHyp[1] = trackTime(massKaon, trkTime, trkPathLength, trkMomentum);
          trkInfo.trkTimeHyp[2] = trackTime(massProton, trkTime, trkPathLength, trkMomentum);
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

        if (verbose) {
          LOG << "vertexTimeFromTracks:     track"
              << " pt=" << (*trk)->pt() << " eta=" << (*trk)->eta() << " phi=" << (*trk)->phi()
              << " vtxWeight=" << trkWeight << " time=" << trkTime << " timeError=" << trkTimeError
              << " timeQuality=" << trkTimeQuality << " pathLength=" << trkPathLength << " momentum=" << trkMomentum
              << " timeHyp[pion]=" << trkInfo.trkTimeHyp[0] << " timeHyp[kaon]=" << trkInfo.trkTimeHyp[1]
              << " timeHyp[proton]=" << trkInfo.trkTimeHyp[2];
        }
      }
    }
  }

  if (wsum > 0) {

    if(verbose) {
      LOG << "vertexTimeFromTracks:   wsum = " << wsum << " tsum = " << tsum << " t0 = " << (wsum > 0 ? tsum/wsum : 0) << " trec = " << vtx.t();
    }

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

      t = tsum / wsum;

      if (verbose) {
        LOG << "vertexTimeFromTracks:   iteration=" << nit << ", T= " << 1/beta << ", t=" << t << ", t-t0=" << t-t0;
      }

      if ((std::abs(t - t0) < 1e-4 / std::sqrt(beta)) and beta >= 1.) {

        tError = std::sqrt(w2sum) / wsum;

        if (verbose) {
          LOG << "vertexTimeFromTracks:   tfit = " << t << " +/- " << tError << " trec = " << vtx.t() << ", iteration=" << nit;
        }

        return true;
      }

      if ((fabs(t - t0) < 1e-3) and beta < 1.) {
        beta = std::min(1., beta / cooling_factor);
      }

      t0 = t;
    }

    LOG << "vertexTimeFromTracks: failed to converge";
  }
  else {
    LOG << "vertexTimeFromTracks: has no track timing info";
  }

  t = t_init;
  tError = tError_init;

  return false;
}

DEFINE_FWK_MODULE(VertexTimeProducer);
