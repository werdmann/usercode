#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include <memory>
#include <string>
#include <unordered_map>

#include <TH1D.h>

#ifdef PVTX_DEBUG
#define LOG edm::Print("VertexTimeAnalyzer")
#else
#define LOG LogDebug("VertexTimeAnalyzer")
#endif

class VertexTimeAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit VertexTimeAnalyzer(edm::ParameterSet const&);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void analyze(const edm::Event &, const edm::EventSetup &) override;

  edm::InputTag const vtx_tag_;
  edm::EDGetTokenT<edm::View<reco::Vertex>> const vtx_token_;

  int const maxNumberOfVertices_;
  std::string const histogramNamePrefix_;

  edm::ParameterSet const histogramPSet_;

  typedef edm::ValueMap<float> edmVMapType;

  struct vMapEntry {
    std::string name;
    edm::EDGetTokenT<edmVMapType> token;
    edm::Handle<edmVMapType> vmapHandle;
  };

  std::unordered_map<std::string, TH1D*> h1map_;
  std::vector<vMapEntry> vmapVec_;
};

VertexTimeAnalyzer::VertexTimeAnalyzer(const edm::ParameterSet &iConfig)
    : vtx_tag_(iConfig.getParameter<edm::InputTag>("vertices")),
      vtx_token_(consumes(vtx_tag_)),
      maxNumberOfVertices_(iConfig.getParameter<int>("maxNumberOfVertices")),
      histogramNamePrefix_(iConfig.getParameter<std::string>("histogramNamePrefix")),
      histogramPSet_(iConfig.getParameter<edm::ParameterSet>("histogramPSet")) {
  usesResource(TFileService::kSharedResource);

  edm::Service<TFileService> fs;

  if (not fs) {
    throw cms::Exception("Configuration") << "TFileService is not registered in cfg file";
  }

  h1map_.clear();

  h1map_["mult"] = fs->make<TH1D>((histogramNamePrefix_+"mult").c_str(), (histogramNamePrefix_+"mult").c_str(), 300, 0, 300);
  h1map_["x"] = fs->make<TH1D>((histogramNamePrefix_+"x").c_str(), (histogramNamePrefix_+"x").c_str(), 600, -0.1, 0.1);
  h1map_["y"] = fs->make<TH1D>((histogramNamePrefix_+"y").c_str(), (histogramNamePrefix_+"y").c_str(), 600, -0.1, 0.1);
  h1map_["z"] = fs->make<TH1D>((histogramNamePrefix_+"z").c_str(), (histogramNamePrefix_+"z").c_str(), 600, -30, 30);
  h1map_["t"] = fs->make<TH1D>((histogramNamePrefix_+"t").c_str(), (histogramNamePrefix_+"t").c_str(), 600, -30, 30);
  h1map_["tError"] = fs->make<TH1D>((histogramNamePrefix_+"tError").c_str(), (histogramNamePrefix_+"tError").c_str(), 600, 0, 30);
  h1map_["normChi2"] = fs->make<TH1D>((histogramNamePrefix_+"normChi2").c_str(), (histogramNamePrefix_+"normChi2").c_str(), 600, 0, 12);
  h1map_["ndof"] = fs->make<TH1D>((histogramNamePrefix_+"ndof").c_str(), (histogramNamePrefix_+"ndof").c_str(), 120, 0, 480);
  h1map_["nTracks"] = fs->make<TH1D>((histogramNamePrefix_+"nTracks").c_str(), (histogramNamePrefix_+"nTracks").c_str(), 120, 0, 480);

  auto const& histogramPSetNames = histogramPSet_.getParameterNamesForType<edm::ParameterSet>();
  vmapVec_.clear();
  vmapVec_.reserve(histogramPSetNames.size());
  for(auto const& psetName : histogramPSetNames){
    auto const hname = histogramNamePrefix_+psetName;
    auto const& pset = histogramPSet_.getParameter<edm::ParameterSet>(psetName);
    auto const nbins = pset.getParameter<uint>("nbins");
    auto const xmin = pset.getParameter<double>("xmin");
    auto const xmax = pset.getParameter<double>("xmax");
    h1map_[psetName] = fs->make<TH1D>(hname.c_str(), hname.c_str(), nbins, xmin, xmax);
    vmapVec_.emplace_back();
    vmapVec_.back().name = psetName;
    vmapVec_.back().token = consumes(pset.getParameter<edm::InputTag>("valueMapTag"));
  }
}

void VertexTimeAnalyzer::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup) {
  auto const& vertices = iEvent.getHandle(vtx_token_);

  if (vertices.isValid()) {

    for(auto& vmap_i : vmapVec_){
      vmap_i.vmapHandle = iEvent.getHandle(vmap_i.token);
    }

    auto const vtxMax = maxNumberOfVertices_ < 0 ? vertices->size() : std::min(vertices->size(), uint(maxNumberOfVertices_));
    h1map_["mult"]->Fill(vtxMax);

    for (uint idx = 0; idx < vtxMax; ++idx) {
      auto const& vtx = (*vertices)[idx];

      h1map_["x"]->Fill(vtx.x());
      h1map_["y"]->Fill(vtx.y());
      h1map_["z"]->Fill(vtx.z());
      h1map_["t"]->Fill(vtx.t());
      h1map_["tError"]->Fill(vtx.tError());
      h1map_["normChi2"]->Fill(vtx.normalizedChi2());
      h1map_["ndof"]->Fill(vtx.ndof());
      h1map_["nTracks"]->Fill(vtx.nTracks());

      for(auto const& vmap_i : vmapVec_){
        h1map_[vmap_i.name]->Fill((*vmap_i.vmapHandle)[vertices->ptrAt(idx)]);
      }
    }
  } else {
    edm::LogWarning("Input") << "invalid handle to reco::VertexCollection : " << vtx_tag_.encode();
  }
}

void VertexTimeAnalyzer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("vertices", edm::InputTag("offlinePrimaryVertices"))->setComment("edm::InputTag of reco::VertexCollection");
  desc.add<int>("maxNumberOfVertices", -1)->setComment("max of number of vertices used to fill histograms");
  desc.add<std::string>("histogramNamePrefix", "vertex_")->setComment("prefix for names of all output histograms");

  edm::ParameterSetDescription histogramPSetDesc;
  histogramPSetDesc.setAllowAnything();
  desc.add<edm::ParameterSetDescription>("histogramPSet", histogramPSetDesc);

  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(VertexTimeAnalyzer);
