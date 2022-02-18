#ifndef usercode_PrimaryVertexAnalyzer_VertexTimeAlgorithmBaseFactory_h
#define usercode_PrimaryVertexAnalyzer_VertexTimeAlgorithmBaseFactory_h

#include "FWCore/PluginManager/interface/PluginFactory.h"
#include "VertexTimeAlgorithmBase.h"

namespace edm {
  class ParameterSet;
  class ConsumesCollector;
}

typedef edmplugin::PluginFactory<VertexTimeAlgorithmBase*(const edm::ParameterSet&, edm::ConsumesCollector& iC)> VertexTimeAlgorithmBaseFactory;

#endif
