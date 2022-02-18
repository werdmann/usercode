#ifndef usercode_PrimaryVertexAnalyzer_VertexTimeAlgorithmBase_h
#define usercode_PrimaryVertexAnalyzer_VertexTimeAlgorithmBase_h

namespace edm {
  class Event;
  class EventSetup;
  class ParameterSet;
  class ParameterSetDescription;
  class ConsumesCollector;
}

namespace reco {
  class Vertex;
}

class VertexTimeAlgorithmBase {
 public:
  VertexTimeAlgorithmBase(const edm::ParameterSet& conf, edm::ConsumesCollector& iC) {}
  virtual ~VertexTimeAlgorithmBase() = default;
  VertexTimeAlgorithmBase(const VertexTimeAlgorithmBase&) = delete;
  VertexTimeAlgorithmBase& operator=(const VertexTimeAlgorithmBase&) = delete;

  static void fillPSetDescription(edm::ParameterSetDescription& iDesc) {}

  virtual void setEvent(edm::Event& iEvent, edm::EventSetup const& iSetup) = 0;
  virtual bool vertexTime(float& vtxTime, float& vtxTimeError, reco::Vertex const& vtx) const = 0;
};

#endif
