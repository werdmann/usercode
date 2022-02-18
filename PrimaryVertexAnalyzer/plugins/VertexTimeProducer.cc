#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/PluginDescription.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "VertexTimeAlgorithmBaseFactory.h"
#include "VertexTimeAlgorithmBase.h"

#include <memory>
#include <utility>

#ifdef PVTX_DEBUG
#define LOG edm::LogPrint("VertexTimeProducer")
#else
#define LOG LogDebug("VertexTimeProducer")
#endif

class VertexTimeProducer : public edm::global::EDProducer<> {
public:
  explicit VertexTimeProducer(edm::ParameterSet const&);

  static void fillDescriptions(edm::ConfigurationDescriptions&);

protected:
  void produce(edm::StreamID, edm::Event&, edm::EventSetup const&) const override;

  edm::EDGetTokenT<reco::VertexCollection> const vtxsToken_;
  bool const produceVertices_;

  std::unique_ptr<VertexTimeAlgorithmBase> const vertexTimeAlgorithmBase_;
};

namespace {
  std::unique_ptr<VertexTimeAlgorithmBase> createVertexTimeAlgorithmBase(edm::ParameterSet const& pset, edm::ConsumesCollector&& iC) {
    return VertexTimeAlgorithmBaseFactory::get()->create(pset.getParameter<std::string>("ComponentType"), pset, iC);
  }
}

VertexTimeProducer::VertexTimeProducer(edm::ParameterSet const& iConfig) :
  vtxsToken_(consumes(iConfig.getParameter<edm::InputTag>("vertices"))),
  produceVertices_(iConfig.getParameter<bool>("produceVertices")),
  vertexTimeAlgorithmBase_(createVertexTimeAlgorithmBase(iConfig.getParameter<edm::ParameterSet>("algorithm"), consumesCollector())) {

  produces<edm::ValueMap<float>>("time");
  produces<edm::ValueMap<float>>("timeError");

  if(produceVertices_){
    produces<reco::VertexCollection>();
  }
}

void VertexTimeProducer::produce(edm::StreamID, edm::Event& iEvent, edm::EventSetup const& iSetup) const {

  auto const vtxsHandle = iEvent.getHandle(vtxsToken_);

  // initialise algorithm for vertex-time computation
  vertexTimeAlgorithmBase_->setEvent(iEvent, iSetup);

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

    vertexTimeAlgorithmBase_->vertexTime(vtxTime, vtxTimeError, vtx);

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

  edm::ParameterSetDescription psd;
  psd.addNode(edm::PluginDescription<VertexTimeAlgorithmBaseFactory>("ComponentType", true));
  desc.add<edm::ParameterSetDescription>("algorithm", psd);

  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(VertexTimeProducer);
