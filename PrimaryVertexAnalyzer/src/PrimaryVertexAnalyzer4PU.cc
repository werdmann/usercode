//#include "Validation/RecoVertex/interface/PrimaryVertexAnalyzer4PU.h"
#include "usercode/PrimaryVertexAnalyzer/interface/PrimaryVertexAnalyzer4PU.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "MagneticField/Engine/interface/MagneticField.h"

// reco track and vertex
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

// simulated vertices,..., add <use name=SimDataFormats/Vertex> and <../Track>
#include <SimDataFormats/Vertex/interface/SimVertex.h>
#include <SimDataFormats/Vertex/interface/SimVertexContainer.h>
#include <SimDataFormats/Track/interface/SimTrack.h>
#include <SimDataFormats/Track/interface/SimTrackContainer.h>

#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// Lumi
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "DataFormats/Luminosity/interface/LumiInfo.h"
// obsolete?
#include "DataFormats/Luminosity/interface/LumiDetails.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"

#include "DataFormats/Common/interface/ConditionsInEdm.h"

// AOD et al
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//generator level + CLHEP
#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"

// TrackingParticle
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
//associator
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"

// fit
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "TrackingTools/TrajectoryParametrization/interface/TrajectoryStateExceptions.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

// Root
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TProfile.h>
//#include <TEfficiency.h>

#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

// cluster stufff
//#include "DataFormats/TrackRecoTrack.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
//#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"  // starting with CMSSW_11
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"

// some of this probably redundant or unused
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
//#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
//#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"

// TCDS / reset orbit info
#include "DataFormats/TCDS/interface/TCDSRecord.h"

#include <sstream>
#include <fstream>

#include <assert.h>


#include<chrono>

using namespace edm;
using namespace reco;
using namespace std;
//
// constants, enums and typedefs
//
//
// static data member definitions
//
namespace {
  template <typename T, size_t N>
  std::array<T, N + 1> makeLogBins(const double min, const double max) {
    const double minLog10 = std::log10(min);
    const double maxLog10 = std::log10(max);
    const double width = (maxLog10 - minLog10) / N;
    std::array<T, N + 1> ret;
    ret[0] = std::pow(10, minLog10);
    const double mult = std::pow(10, width);
    for (size_t i = 1; i <= N; ++i) {
      ret[i] = ret[i - 1] * mult;
    }
    return ret;
  }
}  // namespace

//
// constructors and destructor
//
PrimaryVertexAnalyzer4PU::PrimaryVertexAnalyzer4PU(const ParameterSet& iConfig)
    : verbose_(iConfig.getUntrackedParameter<bool>("verbose", false)),
      //do_pixel_cluster_analysis_( iConfig.getUntrackedParameter<bool>( "doPixelClusterAnalysis", false ) ),
      //do_pixel_cluster_analysis_all_layers_( iConfig.getUntrackedParameter<bool>( "doPixelClusterAnalysisAllLayers", false ) ),
      do_vertex_analysis_(iConfig.getUntrackedParameter<bool>("doVertexAnalysis", true)),
      doMatching_(iConfig.getUntrackedParameter<bool>("matching", false)),
      analyzeLS_(iConfig.getUntrackedParameter<int>("LS", -1)),
      DEBUG_(false),
      eventcounter_(0),
      dumpcounter_(0),
      ndump_(0),
      ndump_tracks_(0),
      run_(0),
      luminosityBlock_(0),
      event_(0),
      bunchCrossing_(0),
      orbitNumber_(0),
      fBfield_(0.),
      simUnit_(1.0),    // cm everywher starting with CMSSW_1_2_x ??
      simtUnit_(1.e9),  // nanoseconds  in rec, do the same for sim
      zmatch_(iConfig.getUntrackedParameter<double>("zmatch", 0.0500)),
      sigmaZoverride_(iConfig.getUntrackedParameter<double>("sigmaZ", 0.0)),
      sigmaZ_(0),
      wxy2_(0.),
      outputFile_(iConfig.getUntrackedParameter<std::string>("outputFile")),
      theTrackFilter(iConfig.getParameter<edm::ParameterSet>("TkFilterParameters")),
      recoTrackProducer_(iConfig.getUntrackedParameter<std::string>("recoTrackProducer")),
      vecPileupSummaryInfoToken_(consumes<std::vector<PileupSummaryInfo>>(edm::InputTag(std::string("addPileupInfo")))),
      genParticleCollection_Token_(consumes<GenParticleCollection>(edm::InputTag(std::string("genParticles")))),
      recoTrackCollectionToken_(consumes<reco::TrackCollection>(
          edm::InputTag(iConfig.getUntrackedParameter<std::string>("recoTrackProducer")))),
      recoBeamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
      edmView_recoTrack_Token_(consumes<edm::View<reco::Track>>(
          edm::InputTag(iConfig.getUntrackedParameter<std::string>("recoTrackProducer")))),
      edmSimVertexContainerToken_(consumes<edm::SimVertexContainer>(iConfig.getParameter<edm::InputTag>("simG4"))),
      edmSimTrackContainerToken_(consumes<edm::SimTrackContainer>(iConfig.getParameter<edm::InputTag>("simG4"))),
      edmHepMCProductToken_(consumes<edm::HepMCProduct>(
          edm::InputTag(std::string("generatorSmeared")))),
      trackingParticleCollectionToken_(consumes<TrackingParticleCollection>(
          iConfig.getUntrackedParameter<edm::InputTag>("trackingParticleCollection"))),
      trackingVertexCollectionToken_(
          consumes<TrackingVertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("trackingVertexCollection"))),
      simToRecoAssociationToken_(
          consumes<reco::SimToRecoCollection>(iConfig.getUntrackedParameter<edm::InputTag>("trackAssociatorMap"))),
      recoToSimAssociationToken_(
          consumes<reco::RecoToSimCollection>(iConfig.getUntrackedParameter<edm::InputTag>("trackAssociatorMap"))),
      vertexAssociatorToken_(consumes<reco::VertexToTrackingVertexAssociator>(
          iConfig.getUntrackedParameter<edm::InputTag>("vertexAssociator"))),
      lumiDetailsToken_(consumes<LumiDetails, edm::InLumi>(edm::InputTag("lumiProducer"))),
      lumiSummaryToken_(consumes<LumiSummary, edm::InLumi>(edm::InputTag("lumiProducer"))),
      lumiInfoToken_(consumes<LumiInfo>(iConfig.getUntrackedParameter<edm::InputTag>("lumiInfoTag")))
{

  // just a test
  auto tkf = iConfig.getParameter<edm::ParameterSet>("TkFilterParameters");
  std::cout << "the track filter parameters " << tkf << std::endl;
  // v34 and higher, avoid biases
  minNumberOfRecTrks_ = 0.;  // FIXME make this configurable or maybe even better obsolete (are these empty BXs?)
  minNumberOfSelTrks_ = 0.;  // FIXME make this configurable

  // definition of visible vertices
  etaMaxVisible_ = 4.0;
  ptMinVisible_ = 0.2;
  numTrkHitsVisible_ = 4;

  fill_track_histos_ = iConfig.getUntrackedParameter<bool>("fill_track_histos", false);
  selNdofNoBS_ = iConfig.getUntrackedParameter<double>("selNdof", 4.);
  std::cout << "selNDof_(noBS) = " << selNdofNoBS_ << std::endl;
  selNdofWithBS_ = iConfig.getUntrackedParameter<double>("selNdofWithBS", 7.);
  std::cout << "selNDofWithBS_ = " << selNdofWithBS_ << std::endl;
  selNdof_ = selNdofNoBS_; // to be changed later according to the collection name
  ndof0trk_ = 0.;

  f4D_ = (iConfig.getUntrackedParameter<bool>("f4D", true));
  if (f4D_) {
    std::cout << "PrimaryVertexAnalyzer4PU  TrackTimeResosLabel "
              << iConfig.getUntrackedParameter<edm::InputTag>("TrackTimeResosLabel") << std::endl;
    trkTimesToken_ = consumes<edm::ValueMap<float>>(iConfig.getUntrackedParameter<edm::InputTag>("TrackTimesLabel"));
    trkTimeResosToken_ =
      consumes<edm::ValueMap<float>>(iConfig.getUntrackedParameter<edm::InputTag>("TrackTimeResosLabel"));  
    trkTimeQualityThreshold_ = iConfig.getUntrackedParameter<double>("TrackTimeQualityThreshold", 0.8);
    trkTimeQualityToken_ =  consumes<edm::ValueMap<float> >(iConfig.getUntrackedParameter<edm::InputTag>("TrackTimeQualityMapLabel"));//mtdTrackQualityMVA:mtdQualMVA"
    //
    //trackExtenderWithMTD:generalTrackPathLength (path length from mtd to beamline)
    //trackExtenderWithMTD:generalTracktmtd  (time at mtd)
    //trackExtenderWithMTD:generalTrackp "Input ValueMap for track momentum magnitude (normally from refit with MTD hits)"
    // see RecoMTD/TimingIDTools/plugins/TOFPIDProducer.cc
    MTD_pathlength_Token_ = consumes<edm::ValueMap<float> >(edm::InputTag("trackExtenderWithMTD:generalTrackPathLength")); 
    MTD_time_Token_ = consumes<edm::ValueMap<float> >(edm::InputTag("trackExtenderWithMTD:generalTracktmtd"));
    MTD_timeerror_Token_ = consumes<edm::ValueMap<float> >(edm::InputTag("trackExtenderWithMTD:generalTracksigmatmtd")); 
    MTD_momentum_Token_ = consumes<edm::ValueMap<float> >(edm::InputTag("trackExtenderWithMTD:generalTrackp")); 
  } else {
    std::cout << "PrimaryVertexAnalyzer4PU: no timing" << std::endl;
  }

  reco_vertex_collections_ = iConfig.getParameter<std::vector<edm::InputTag>>("vertexRecoCollections");

  for (auto const& l : reco_vertex_collections_) {
    std::cout << "vertex collection " << l.label() << std::endl;
    vertexCollectionLabels_.push_back(l.label());

    auto token = edm::EDGetTokenT<reco::VertexCollection>(consumes<reco::VertexCollection>(edm::InputTag(l)));
    vertexCollectionTokens_[l.label()] = token;
    reco_vertex_view_tokens_.push_back(edm::EDGetTokenT<edm::View<reco::Vertex>>(consumes<edm::View<reco::Vertex>>(l)));
  }

  // open output file to store histogram
  rootFile_ = TFile::Open(outputFile_.c_str(), "RECREATE");
  info_ = new TObjString("info");
  info_->SetString(iConfig.getUntrackedParameter<std::string>("info", "").c_str());
  build_ = new TObjString("build");
  cout << "build = " << Form("%s %s", __DATE__, __TIME__) << endl;
  build_->SetString(Form("%s %s", __DATE__, __TIME__));

  veryverbose_ = iConfig.getUntrackedParameter<bool>("veryverbose", false);
  if (verbose_) {
    cout << "veryverbose = " << veryverbose_ << endl;
    cout << "    extra level of verbosity " << endl;
    cout << endl;
  }

  if (!do_vertex_analysis_)
    return;

  if (verbose_) {
    cout << "PrimaryVertexAnalyzer4PU: zmatch=" << zmatch_ << endl;
    cout << "sigmaZ = " << sigmaZoverride_ << endl;
    cout << "     if 0.0 : use the value from the beamspot" << endl;
    cout << endl;
  }

  nPUmin_ = std::abs(iConfig.getUntrackedParameter<int>("PUmin", 0));
  nPUmax_ = std::abs(iConfig.getUntrackedParameter<int>("PUmax", 1000000));
  if (verbose_) {
    cout << "nPUMin = " << nPUmin_ << endl;
    cout << "nPUMax = " << nPUmax_ << endl;
    cout << "     in MC, only analyze events with nPUmin <  N(simvertex) < nPUmax " << endl;
    cout << endl;
  }

  useVertexFilter_ = iConfig.getUntrackedParameter<bool>("useVertexFilter", false);

  nEventSummary_ = iConfig.getUntrackedParameter<int>("eventSummaries", 0);
  ndump_ = iConfig.getUntrackedParameter<int>("nDump", 0);
  ndump_tracks_ = iConfig.getUntrackedParameter<int>("nDumpTracks", 0);

  nCompareCollections_ =
      iConfig.getUntrackedParameter<int>("compareCollections", 0);  //-1= compare all, >0 dump n events

  zmatch_ = iConfig.getUntrackedParameter<double>("zmatch", 0.0500);
  if (verbose_) {
    cout << "zmatch =" << zmatch_ << endl;
    cout << "     cut-off for matching sim to reco by z" << endl;
    cout << endl;
  }

  zWosMatchMax_ = iConfig.getUntrackedParameter<double>("zMatchMax", 1.);
  if (verbose_) {
    cout << "zMatchMax = " << zWosMatchMax_ << endl;
    cout << "     cut-off for insane recvertex <-> simvertex  matches" << endl;
    cout << "     (TrackingParticles, matching by weight/sigma^2) " << endl;
    cout << "     default is 1 cm, configurable for exotic events where" << endl;
    cout << "     all tracks appear far away from the vertex   " << endl;
    // such as LongLivedChi0ToNuLL_MSquark-1000_MChi-148_TuneZ2Star_8TeV-pythia6
    cout << endl;
  }

  RECO_ = iConfig.getUntrackedParameter<bool>("RECO", false);
  if (verbose_) {
    cout << "RECO = " << RECO_ << endl;
    cout << "      use RECO information (pixel hits and some trackextra)" << endl;
    cout << endl;
  }

  autoDumpCounter_ = iConfig.getUntrackedParameter<int>("autodump", 0);
  if (verbose_) {
    cout << "autodump = " << autoDumpCounter_ << endl;
    cout << "      dump detailed information about the first <autodump> events " << endl;
    cout << endl;
  }

  trackAssociatorMin_ = 0.5;

  /* initialize counters */
  eventcounter_ = 0;
  emptyeventcounter_ = 0;
  dumpcounter_ = 0;
  eventSummaryCounter_ = 0;
  nEventNsel_ = 0;

  currentLS_ = -1;

  // profiling
  timers_.clear();
  
}

PrimaryVertexAnalyzer4PU::~PrimaryVertexAnalyzer4PU() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  delete rootFile_;
}

//
// member functions
//

std::map<std::string, TH1*> PrimaryVertexAnalyzer4PU::bookVertexHistograms(TDirectory * dir) {
  std::map<std::string, TH1*> h;  // will be returned


  /* cpu time related */
  dir->mkdir("cputime")->cd();
  addn(h, new TProfile("tcluvsnsel", "clustering time vs nvertex", 300, 0., 300., 0., 1.e8));
  addn(h, new TProfile("tfitvsnsel", "fitting time vs nvertex", 300, 0., 300., 0., 1.e8));
  addn(h, new TProfile("tcluvsLPU", "clustering time vs pu", 300, 0., 300., 0., 1.e8));
  addn(h, new TProfile("tfitvsLPU", "fitting time vs pu", 300, 0., 300., 0., 1.e8));
  addn(h, new TProfile("tcluvsSimPU", "clustering time vs #simvtx", 300, 0., 300., 0., 1.e8));
  addn(h, new TProfile("tfitvsSimPU", "fitting time vs #simvtx", 300, 0., 300., 0., 1.e8));
  dir->cd();

  /* pseudo physics variables, ptvis */
  // 0   |dz| < 0.1 
  // 1   |dz| < 0.2
  // 2   |dz| < 3 sigma
  // 3   nearest vertex
  // 4   weighted

  dir->mkdir("ptvis")->cd();
  for(unsigned int i = 0; i < 5; i++){ 
    addn(h, new TH1F(Form("ptvis%d", i), Form("ptvis%d",i), 100, 0., 2.));
    addn(h, new TProfile(Form("ptvis%dvsz", i), Form("ptvis%d vs z", i), 30, -15., 15., 0., 100.));
    addn(h, new TProfile(Form("ptvis%dvspu", i), Form("ptvis%d vs pu", i), 30, 0., 300., 0., 100.));
    addn(h, new TProfile(Form("ptvis%dvsdensity", i), Form("ptvis%d vs density", i), 30, 0., 30., 0., 100.));
    // note that the vs density plot is a bit tricky to interprete, because at low density it contains more of the high z region

    addn(h, new TH1F(Form("ptmiss_sim%d", i), Form("sim ptmiss%d",i), 100, 0., 200.));
    addn(h, new TH1F(Form("ptmiss%d", i), Form("ptmiss%d",i), 100, 0., 200.));
    addn(h, new TProfile(Form("ptmiss%dvsz", i), Form("ptmiss%d vs z", i), 30, -15., 15., 0., 1000.));
    addn(h, new TProfile(Form("ptmiss%dvspu", i), Form("ptmiss%d vs pu", i), 30, 0., 300., 0., 1000.));
    addn(h, new TProfile(Form("ptmiss%dvsdensity", i), Form("ptmiss%d vs density", i), 30, 0., 30., 0., 1000.));
  }
  dir->cd();

  
  const int nbinzdiffrec = 800.;
  const double npumax = 300.;
  const int npubin2 = 30;
  const int nzbin_wide = 30;
  const int nzbin_wide_fine = 120;
  const float zrange_wide = 30.;

  //const int nzbin_normal = 30; until used
  const int nzbin_normal_fine = 60;
  const float zrange_normal = 15.;

  const int nrecmax = 300;
  const int ntrkmax = 6000.;

  int nvtxbin = 150;
  float nvtxrange = 300.;

  add(h, new TH1F("zofzgap02mm", "z of zgap < 0.2 mm", 300, -15., 15.));
  add(h, new TH1F("zofzgap03mm", "z of zgap < 0.3 mm", 300, -15., 15.));
  add(h, new TH1F("zofzgap05mm", "z of zgap < 0.5 mm", 300, -15., 15.));
  add(h, new TH1F("zofzgap1mm", "z of zgap < 1 mm", 300, -15., 15.));
  add(h, new TH1F("zofselvtx", "z of selected vertices", 300, -15., 15.));  // reference, same binning
  add(h, new TH1F("zofzgap02mm_tail", "z of zgap < 0.2 mm", 300, -15., 15.));
  add(h, new TH1F("zofzgap03mm_tail", "z of zgap < 0.3 mm", 300, -15., 15.));
  add(h, new TH1F("zofzgap05mm_tail", "z of zgap < 0.5 mm", 300, -15., 15.));
  add(h, new TH1F("zofzgap1mm_tail", "z of zgap < 1 mm", 300, -15., 15.));
  add(h, new TH1F("zofselvtx_tail", "z of selected vertices", 300, -15., 15.));  // reference, same binning

  addSP(h, new TH1F("nsim_recvtx_sel", "number of sim vertices in rec vtx", 10, 0, 10));
  addSP(h, new TH1F("nsim_any_recvtx_sel", "number of sim vertices in rec vtx", 10, 0, 10));
  addSP(h, new TH1F("nsim_any_not_matched_recvtx_sel", "number of sim vertices in rec vtx", 10, 0, 10));
  addSP(h, new TH1F("zrec_recvtx_sel", "zrec of selected vertices", 100, -10., 10));
  addSP(h, new TH1F("zrec_recvtx_selnsim0", "zrec of selected vertices with no simvtx", 100, -10., 10));
  addSP(h, new TH1F("zrec_recvtx_selnsim1", "zrec of selected vertices with one simvtx", 100, -10., 10.));
  addSP(h,
        new TH1F("zrec_recvtx_selnsimgt1",
                 "zrec of selected vertices with more than one simvtx",
                 100,
                 -10.,
                 10.));  // for the merge fraction
  addSP(h, new TProfile("nsimvsz_recvtx_sel", "number of sim vertices in rec vtx vs z", 100, -10., 10., 0., 100.));
  // fake rate vs distance
  add(h, new TH1F("countzrecsignal","", 1, 0., 2.)); // count relevant signal vertices, denominator for the next two histos
  add(h, new TH1F("zrecsimsignal_unmatchedPU", "distance of a fake vertex to simulated signal position", 200, -1., 1.));
  add(h, new TH1F("zrecrecsignal_unmatchedPU", "distance of a fake vertex to reconstructed signal position", 200, -1., 1.));
  add(h, new TH1F("zrecsimsignal_unmatchedPU_hr", "distance of a fake vertex to simulated signal position", 200, -0.1, 0.1));
  add(h, new TH1F("zrecrecsignal_unmatchedPU_hr", "distance of a fake vertex to reconstructed signal position", 200, -0.1, 0.1));
  // fake fraction vs distance
  add(h, new TProfile("zrecsimsignal_unmatchedPUfraction", "fake fraction vs distance to simulated signal position",200, -1., 1., 0., 2.));
  add(h, new TProfile("zrecrecsignal_unmatchedPUfraction", "fake fraction vs distance to reconstructed signal position", 200, -1., 1., 0., 2.));
  add(h, new TProfile("zrecsimsignal_unmatchedPUfraction_hr", "fake fraction vs distance to simulated signal position",200, -0.1, 0.1, 0., 2.));
  add(h, new TProfile("zrecrecsignal_unmatchedPUfraction_hr", "fake fraction vs distance to reconstructed signal position", 200, -0.1, 0.1, 0., 2.));

  add(h, new TProfile("effvsdzsimsignal_PU", "PU efficiency vs distance to signal vertex", 200, -1., 1., 0., 2.));
  add(h, new TProfile("effvsdzsimPU_PU", "PU efficiency vs distance to other PU vtxs", 200, -1., 1., 0., 2.));
  add(h, new TProfile("effvsdzsimsignal_PU_hr", "PU efficiency vs distance to signal vertex", 200, -0.1, 0.1, 0., 2.));
  add(h, new TProfile("effvsdzsimPU_PU_hr", "PU efficiency vs distance to other PU vtxs", 200, -0.1, 0.1, 0., 2.));
  add(h, new TProfile("effvsdzsimmatchedPU_PU", "PU efficiency vs distance to other matched PU vtxs", 200, -1., 1., 0., 2.));
  add(h, new TProfile("effvsdzsimunmatchedPU_PU", "PU efficiency vs distance to other unmatched PU vtxs", 200, -1., 1., 0., 2.));
  add(h, new TProfile("effvsdzsimmatchedPU_PU_hr", "PU efficiency vs distance to other matched PU vtxs", 200, -0.1, 0.1, 0., 2.));
  add(h, new TProfile("effvsdzsimunmatchedPU_PU_hr", "PU efficiency vs distance to other unmatched PU vtxs", 200, -0.1, 0.1, 0., 2.));

  add(h, new TProfile("effselvsdzsimsignal_PU", "selected PU efficiency vs distance to signal vertex", 200, -1., 1., 0., 2.));
  add(h, new TProfile("effselvsdzsimPU_PU", "selected PU efficiency vs distance to other PU vtxs", 200, -1., 1., 0., 2.));
  add(h, new TProfile("effselvsdzsimsignal_PU_hr", "selected PU PU efficiency vs distance to signal vertex", 200, -0.1, 0.1, 0., 2.));
  add(h, new TProfile("effselvsdzsimPU_PU_hr", "selected PU  efficiency vs distance to other PU vtxs", 200, -0.1, 0.1, 0., 2.));
  add(h, new TProfile("effselvsdzsimmatchedPU_PU", "selected PU efficiency vs distance to other matched PU vtxs", 200, -1., 1., 0., 2.));
  add(h, new TProfile("effselvsdzsimunmatchedPU_PU", "selected PU efficiency vs distance to other unmatched PU vtxs", 200, -1., 1., 0., 2.));
  add(h, new TProfile("effselvsdzsimmatchedPU_PU_hr", "selected PU efficiency vs distance to other matched PU vtxs", 200, -0.1, 0.1, 0., 2.));
  add(h, new TProfile("effselvsdzsimunmatchedPU_PU_hr", "selected PU efficiency vs distance to other unmatched PU vtxs", 200, -0.1, 0.1, 0., 2.));

  add(h, new TProfile("effallvsdzsimminPU_PU", "PU efficiency vs min distance to other simulated PU", 200, -1., 1., 0., 2.));
  add(h, new TProfile("effselvsdzsimminPU_PU", "selected PU efficiency vs min distance to other simulated PU", 200, -1., 1., 0., 2.));

  add(h, new TH1F("nbtksinvtx", "reconstructed tracks in vertex", 40, -0.5, 39.5));
  add(h, new TH1F("nbtksinvtxPU", "reconstructed tracks in vertex", 40, -0.5, 39.5));
  addSP(h, new TH1F("resx", "residual x", 200, -0.02, 0.02));
  addSP(h, new TH1F("resy", "residual y", 200, -0.02, 0.02));
  addSP(h, new TH1F("resz", "residual z", 200, -0.05, 0.05));
  addSP(h, new TH1F("resz10", "residual z", 200, -0.5, 0.5));
  addSP(h, new TH1F("pullx", "pull x", 100, -25., 25.));
  addSP(h, new TH1F("pully", "pull y", 100, -25., 25.));
  addSP(h, new TH1F("pullz", "pull z", 100, -25., 25.));

  add(h, new TH1F("vtxcxy2", "xy-beam compatibility", 100, 0., 50.));
  add(h, new TH1F("vtxcxy2Matched", "xy-beam compatibility", 100, 0., 50.));
  add(h, new TH1F("vtxprobxyMatched", "xy-beam compatibility", 100, 0., 1.));
  add(h, new TH1F("vtxcxy2Fake", "xy-beam compatibility", 100, 0., 50.));
  add(h, new TH1F("vtxprobxyFake", "xy-beam compatibility", 100, 0., 1.));
  add(h, new TH1F("ntpfake0", "fake vertices", 20, 0., 20.));
  add(h, new TH1F("ntpfake4", "fake vertices with ndof>4", 20, 0., 20.));
  add(h, new TH1F("ndofFake", "ndof of fake vertices", 500, 0., 100.));
  add(h, new TH1F("logndofFake", "log10(ndof) of fake vertices", 500, -1., 10.));
  add(h, new TH1F("ndofMatched", "ndof of matched vertices", 500, 0., 100.));

  add(h, new TH1F("ntpfound0", "found vertices", 300, 0., 300.));
  add(h, new TH1F("ntpfound4", "found vertices with ndof>4", 300, 0., 300.));

  add(h, new TH1F("indexempty", "index of the empty vertex", 200, 0., 200.));
  add(h, new TH1F("trksumpt2rank_signal", "track sumpt2 ranking signal", 200, 0. ,200.));
  add(h, new TH1F("trksumpt2rank_pu", "track sumpt2 ranking pu", 200, 0. ,200.));
  add(h, new TH1F("trksumpt2rank_splitfromsignal", "track sumpt2 ranking split from signal", 200, 0. ,200.));
  add(h, new TH1F("trksumpt2rank_splitfrompu", "track sumpt2 ranking split from pu", 200, 0. ,200.));
  add(h, new TH1F("trksumpt2rank_otherfake", "track sumpt2 ranking other fake", 200, 0. ,200.));

  add(h, new TH1F("index_signal", "index signal", 200, 0. ,200.));
  add(h, new TH1F("index_pu", "index pu", 200, 0. ,200.));
  add(h, new TH1F("index_splitfromsignal", "index split from signal", 200, 0. ,200.));
  add(h, new TH1F("index_splitfrompu", "index  split from pu", 200, 0. ,200.));
  add(h, new TH1F("index_otherfake", "index other fake", 200, 0. ,200.));
  

  add(h, new TH1F("ndof-vtx0", "ndof(vtx0)", 500, 0., 200.));
  add(h, new TH1F("ndof-vtx1", "ndof(vtx1)", 500, 0., 200.));
  add(h, new TH1F("ndof-vtx2", "ndof(vtx2)", 500, 0., 200.));
  add(h, new TH1F("ndof-vtx3", "ndof(vtx3)", 500, 0., 200.));
  add(h, new TH1F("ndof-vtx4", "ndof(vtx4)", 500, 0., 200.));
  add(h, new TH1F("ndof-vtx5", "ndof(vtx5)", 500, 0., 200.));
  add(h, new TH1F("ndof-vtx6", "ndof(vtx6)", 500, 0., 200.));
  add(h, new TH1F("ndof-vtx7", "ndof(vtx7)", 500, 0., 200.));
  add(h, new TH1F("ndof-vtx8", "ndof(vtx8)", 500, 0., 200.));
  add(h, new TH1F("ndof-vtx9", "ndof(vtx9)", 500, 0., 200.));
  add(h, new TH1F("ndof-vtxgt9", "ndof(vtxg9)", 500, 0., 200.));

  add(h, new TH1F("logndof-vtx0", "log ndof(vtx0)", 100, 0., 3.));
  add(h, new TH1F("logndof-vtx1", "log ndof(vtx1)", 100, 0., 3.));
  add(h, new TH1F("logndof-vtx2", "log ndof(vtx2)", 100, 0., 3.));
  add(h, new TH1F("logndof-vtx3", "log ndof(vtx3)", 100, 0., 3.));
  add(h, new TH1F("logndof-vtx4", "log ndof(vtx4)", 100, 0., 3.));
  add(h, new TH1F("logndof-vtx5", "log ndof(vtx5)", 100, 0., 3.));
  add(h, new TH1F("logndof-vtx6", "log ndof(vtx6)", 100, 0., 3.));
  add(h, new TH1F("logndof-vtx7", "log ndof(vtx7)", 100, 0., 3.));
  add(h, new TH1F("logndof-vtx8", "log ndof(vtx8)", 100, 0., 3.));
  add(h, new TH1F("logndof-vtx9", "log ndof(vtx9)", 100, 0., 3.));
  add(h, new TH1F("logndof-vtxgt9", "log ndof(vtxg9)", 500, 0., 3.));

  add(h, new TH1F("zdiffrec-vtx1", "z(vtx1) -z(vtx0)", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("zdiffrec-vtx2", "z(vtx2) -z(vtx0)", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("zdiffrec-vtx3", "z(vtx3) -z(vtx0)", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("zdiffrec-vtx4", "z(vtx4) -z(vtx0)", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("zdiffrec-vtx5", "z(vtx5) -z(vtx0)", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("zdiffrec-vtx6", "z(vtx6) -z(vtx0)", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("zdiffrec-vtx7", "z(vtx7) -z(vtx0)", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("zdiffrec-vtx8", "z(vtx8) -z(vtx0)", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("zdiffrec-vtx9", "z(vtx9) -z(vtx0)", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("zdiffrec-vtxgt9", "z(vtxgt9) -z(vtx0)", nbinzdiffrec, -2., 2.));

  add(h, new TH1F("logsumpt2-vtx0", "logsumpt2(vtx0)", 200, -1, 9.));
  add(h, new TH1F("logsumpt2-vtx1", "logsumpt2(vtx1)", 200, -1, 9.));
  add(h, new TH1F("logsumpt2-vtx2", "logsumpt2(vtx2)", 200, -1, 9.));
  add(h, new TH1F("logsumpt2-vtx3", "logsumpt2(vtx3)", 200, -1, 9.));
  add(h, new TH1F("logsumpt2-vtx4", "logsumpt2(vtx4)", 200, -1, 9.));
  add(h, new TH1F("logsumpt2-vtx5", "logsumpt2(vtx5)", 200, -1, 9.));
  add(h, new TH1F("logsumpt2-vtx6", "logsumpt2(vtx6)", 200, -1, 9.));
  add(h, new TH1F("logsumpt2-vtx7", "logsumpt2(vtx7)", 200, -1, 9.));
  add(h, new TH1F("logsumpt2-vtx8", "logsumpt2(vtx8)", 200, -1, 9.));
  add(h, new TH1F("logsumpt2-vtx9", "logsumpt2(vtx9)", 200, -1, 9.));
  add(h, new TH1F("logsumpt2-vtxgt9", "logsumpt2(vtxgt9)", 200, -1, 9.));

  add(h, new TH1F("zdiffrec10-selsel", "z(vtx1) -z(vtx0)", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("zdiffrec10-signalsel", "z(vtx1) -z(vtx0, signal)", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("zdiffrec10-signalreal", "z(vtx1, real) -z(vtx0, signal)", nbinzdiffrec, -2., 2.));//for bwd compatibility only, obsolete
  add(h, new TH1F("zdiffrec10-signalrealpu", "z(vtx1, real pu) -z(vtx0, signal)", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("zdiffrec10-signalfake", "z(vtx1, fake) -z(vtx0, signal)", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("zdiffrec10-signalsplitfromsignal", "z(vtx1, split from signal) -z(vtx0, signal)", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("zdiffrec10-signalsplitfrompu", "z(vtx1, split from pu) -z(vtx0, signal)", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("zdiffrec10-signalotherfake", "z(vtx1, other fake) -z(vtx0, signal)", nbinzdiffrec, -2., 2.));

  add(h, new TH1F("zdiffrecsignalsel",
               "z-distance between reconstructed signal and selected split-off",
               nbinzdiffrec, -2., 2.));

  add(h, new TH1F("zdiffrecselfound", "#Delta z_rec truth-matched ndof>4 vertices", nbinzdiffrec, -2., 2.));
  add(h,
      new TH1F(
          "zdiffrecselfakereal", "z-distance between truth-matched and fake ndof>4 vertices", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("zdiffrecselfakefake", "z-distance between two fake ndof>4 vertices", nbinzdiffrec, -2., 2.));

  add(h,
      new TH1F("zdiffrecselsignalreal",
               "z-distance between reconstructed signal and other real  ndof>4 vertices",
               nbinzdiffrec,
               0.,
               2.));
  add(h,
      new TH1F("zdiffrecselsignalfake",
               "z-distance between reconstructed signal and fake ndof>4 vertices",
               nbinzdiffrec,
               0.,
               2.));
  add(h,
      new TH1F("zdiffrecselsignalsplit",
               "z-distance between reconstructed signal and split-off selected vertices",
               nbinzdiffrec,
               0.,
               2.));
  add(h,
      new TH1F("zdiffrecselsignalrealSimpv",
               "z-distance between reconstructed signal and other real  ndof>4 vertices",
               nbinzdiffrec,
               0.,
               2.));
  add(h,
      new TH1F("zdiffrecselsignalfakeSimpv",
               "z-distance between reconstructed signal and fake ndof>4 vertices",
               nbinzdiffrec,
               0.,
               2.));
  add(h,
      new TH1F("zdiffrecselpxysignalrealSimpv",
               "z-distance between reconstructed signal and other real  ndof>4 vertices",
               nbinzdiffrec,
               0.,
               2.));
  add(h,
      new TH1F("zdiffrecselpxysignalfakeSimpv",
               "z-distance between reconstructed signal and fake ndof>4 vertices",
               nbinzdiffrec,
               0.,
               2.));
  add(h,
      new TH1F("zdiffrecselsignalsplitSimpv",
               "z-distance between reconstructed signal and split-off ndof>4 vertices",
               nbinzdiffrec,
               0.,
               2.));
  add(h,
      new TH1F("zdiffrecselpxysignalsplitSimpv",
               "z-distance between reconstructed signal and split-off ndof>4 vertices",
               nbinzdiffrec,
               0.,
               2.));
  add(h,
      new TH1F("zdiffrecselPUrealPU",
               "z-distance between reconstructed PU and other real PU ndof>4 vertices",
               nbinzdiffrec,
               0.,
               2.));
  add(h,
      new TH1F(
          "zdiffrecselPUfake", "z-distance between reconstructed PU and fake ndof>4 vertices", nbinzdiffrec, 0., 2.));

  add(h, new TH1F("zdiffrecselsignalfake", "z-distance between reconstructed Signal and fake selected vertices", nbinzdiffrec, 0., 2.));
  add(h, new TH1F("zdiffrecselsignalrealpu", "z-distance between reconstructed Signal and real selected vertices", nbinzdiffrec, 0., 2.));
  add(h, new TH1F("zdiffrecselsignalrealpuidx1", "z-distance between reconstructed Signal and real selected vertices at index 1", nbinzdiffrec, 0., 2.));
  add(h, new TH1F("zdiffrecselsignalsplitfromsignal", "z-distance between reconstructed Signal and fake selected vertices", nbinzdiffrec, 0., 2.));
  add(h, new TH1F("zdiffrecselsignalsplitfrompu", "z-distance between reconstructed Signal and split from pu selected vertices", nbinzdiffrec, 0., 2.));
  add(h, new TH1F("zdiffrecselsignalotherfake", "z-distance between reconstructed Signal and combinatorial fake selected vertices", nbinzdiffrec, 0., 2.));

  add(h, new TH1F("zdiffmin4fakereal", "z-distance between fake and nearest real rec vertex", nbinzdiffrec, 0, 2.));
  add(h, new TH1F("zdiffmin4realreal", "z-distance between real and nearest real rec vertex", nbinzdiffrec, 0, 2.));
  add(h, new TH1F("zdiffmin4realfake", "z-distance between real and nearest fake vertex", nbinzdiffrec, 0, 2.));
  add(h, new TH1F("zdiffmin4signalreal", "z-distance between signal and nearest real rec vertex", nbinzdiffrec, 0, 2.));
  add(h, new TH1F("zdiffmin4signalfake", "z-distance between signal and nearest fake vertex", nbinzdiffrec, 0, 2.));

  add(h, new TH1F("ndofsignalsplit", "ndof of vertices split from the signal vertex", 100, 0., 20.));
  add(h, new TH1F("ndofsignalsplitSimpv", "ndof of vertices split from the signal vertex", 100, 0., 20.));
  add(h, new TH1F("ndofpxysignalsplitSimpv", "ndof of vertices split from the signal vertex", 100, 0., 20.));

  add(h, new TProfile("ntpfake4prof", "fake vertices with ndof>4", 20, 0., 100., 0., 200.));
  add(h, new TProfile("ntpfound4prof", "found vertices with ndof>4", 20, 0., 100., 0., 200.));

  add(h, new TProfile("ntpselvssimPU", "selected vertices", nvtxbin, 0., nvtxrange, 0., nvtxrange * 2.));
  add(h, new TProfile("ntpfakeselvssimPU", "fake selected vertices", nvtxbin, 0., nvtxrange, 0., nvtxrange * 2.));
  add(h, new TProfile("ntprealselvssimPU", "matched selected vertices", nvtxbin, 0., nvtxrange, 0., nvtxrange * 2.));
  add(h, new TProfile("ntpsplitselvssimPU", "split selected vertices", nvtxbin, 0., nvtxrange, 0., nvtxrange * 2.));
  add(h, new TProfile("ntpotherfakeselvssimPU", "other fake selected vertices", nvtxbin, 0., nvtxrange, 0., nvtxrange * 2.));

  add(h, new TH1F("ntpsplitselfromsignal", "split from signal", 10, 0., 10.));
  add(h, new TH1F("ntpsplitselfrompu", "split from PU", 10, 0., 10.));

  add(h, new TH2F("wornk_matchedsel","n70% vs sumwos", 100, -1, 7., 5, 0., 5.));
  add(h, new TH2F("wornk_unmatchedsel","n70% vs sumwos", 100, -1., 7., 5, 0., 5.));
  add(h, new TH2F("wornk_splitsel","n70% vs sumwos", 100, -1., 7., 5, 0., 5.));
  add(h, new TH2F("wornk_othersel","n70% vs sumwos", 100, -1., 7., 5, 0., 5.));
  add(h, new TH2F("wornk_matchedsim","n70% vs sumwos", 100, -1., 7., 5, 0., 5.));
  add(h, new TH2F("wornk_unmatchedsim","n70% vs sumwos", 100, -1., 7., 5, 0., 5.));

  addSP(h, new TH1F("vtxtrkpurity", "track purity in vertex", 100, 0., 1.));
  addSP(h, new TProfile("vtxtrkpurityvspu", "track purity in vertex vs PU", npubin2, 0., npumax, 0., 2.));
  addSP(h, new TProfile("vtxtrkpurityvsz", "vertex purity", nzbin_wide_fine, -zrange_wide, zrange_wide, 0., 1.));
  // now the purity with only matched tracks in the denominator
  addSP(h, new TH1F("vtxtrkpuritym", "track purity in vertex", 100, 0., 1.));
  addSP(h, new TProfile("vtxtrkpuritymvspu", "track purity in vertex vs PU (matched only)", npubin2, 0., npumax, 0., 2.));
  addSP(h, new TProfile("vtxtrkpuritymvsz", "vertex purity (matched only)", nzbin_wide_fine, -zrange_wide, zrange_wide, 0., 1.));
  // assignment efficiency (counts matched tracks only anyway)
  addSP(h, new TH1F("trkAssignmentEfficiency", "conditional track to vertex assignment efficiency", 101, 0., 1.01));
  addSP(h, new TH1F("primtrkAssignmentEfficiency", "conditional track to vertex assignment efficiency", 101, 0., 1.01));
  addSP(h, new TProfile("trkAssignmentEfficiencyvspu", "conditional track to vertex assignment efficiency vs PU", npubin2, 0., npumax, 0., 2.));
  addSP(h, new TProfile("trkAssignmentEfficiencyvsntrk", "conditional track to vertex assignment efficiency vs number of tracks",
			51, 0., 51, 0., 2.));  // last bin is overflow
  addSP(h, new TH1F("utrkAssignmentEfficiency", "track to vertex assignment efficiency", 101, 0., 1.01));
  addSP(h, new TH1F("uprimtrkAssignmentEfficiency", "track to vertex assignment efficiency", 101, 0., 1.01));
  addSP(h, new TProfile("utrkAssignmentEfficiencyvspu", "track to vertex assignment efficiency vs PU", npubin2, 0., npumax, 0., 2.));
  addSP(h, new TProfile("utrkAssignmentEfficiencyvsntrk", "track to vertex assignment efficiency vs number of tracks",
			51, 0., 51, 0., 2.));  // last bin is overflow

 
  addSP(h, new TProfile("effallvsnrectp", "efficiency (all) vs # truth matched rec tracks", 50, 0., 50., 0, 1.));
  addSP(h, new TProfile("effallvsngentp", "efficiency (all) vs # visible tracking particles", 50, 0., 50., 0, 1.));
  addSP(h, new TProfile("effselvsnrectp", "efficiency (sel) vs # truth matched rec tracks", 50, 0., 50., 0, 1.));
  addSP(h, new TProfile("effselvsngentp", "efficiency (sel) vs # visible tracking particles", 50, 0., 50., 0, 1.));
  addSP(h, new TProfile("effallvssumpt", "efficiency (all) vs # pt sum of charged particles", 50, 0., 20., 0, 1.));
  addSP(h, new TProfile("effselvssumpt", "efficiency (sel) vs # pt sum of charged particles", 50, 0., 20., 0, 1.));
  addSP(h, new TProfile("effallvssqrtsumpt2", "efficiency (sel) vs # sqrt(pt2 sum of charged particles)", 100, 0., 100, 0, 1.));
  addSP(h, new TProfile("effselvssqrtsumpt2", "efficiency (sel) vs # sqrt(pt2 sum of charged particles)", 100, 0., 100, 0, 1.));
  addSP(h, new TProfile("effallvstrksqrtsumpt2", "efficiency (sel) vs # sqrt(pt2 sum of matched tracks)", 100, 0., 100, 0, 1.));//FIXME not filled yet
  addSP(h, new TProfile("effselvstrksqrtsumpt2", "efficiency (sel) vs # sqrt(pt2 sum of matched tracks)", 100, 0., 100, 0, 1.));

  add(h, new TH1F("vtxndf_tagged", "degrees of freedom (tagged)", 5000, 0., 1000.));
  add(h, new TH1F("vtxndfc", "expected lower ndof of two", 5000, 0., 1000.));
  add(h, new TH1F("zndof2", "zndof2", 10, -10., 10.));  //dummy

  add(h, new TH1F("vtxndfIso", "degrees of freedom (isolated vertex)", 5000, 0., 1000.));
  add(h, new TH1F("vtxndfoverntk", "ndof / #tracks", 40, 0., 2.));
  add(h, new TH1F("vtxndf2overntk", "(ndof+2) / #tracks", 40, 0., 2.));

  // raw
  add(h, new TH1F("isFake", "fake vertex", 2, -0.5, 1.5));
  add(h, new TH1F("isFake1", "fake vertex or ndof<0", 2, -0.5, 1.5));
  add(h, new TH1F("bunchCrossing", "bunchCrossing", 3565, 0., 3565.));
  add(h, new TH1F("bunchCrossing_PU", "lumi per bunch crossing", 3565, 0., 3565.));

  add(h, new TH1F("highpurityTrackFraction", "fraction of high purity tracks", 20, 0., 1.));
  add(h, new TH1F("trkchi2overndof", "vertices chi2 / ndof", 50, 0., 5.));

  add(h, new TH1F("z0trk", "track z0 (tight, eta<1.5, pt>0.5))", 100., -20., 20.));
  add(h, new TH1F("resx50", "residual x (ndof>50)", 100, -0.04, 0.04));
  add(h, new TH1F("resy50", "residual y (ndof>50)", 100, -0.04, 0.04));
  add(h, new TH1F("resz50", "residual z (ndof>50)", 100, -0.1, 0.1));
  add(h, new TH1F("pullxr", "relative pull x", 100, -25., 25.));
  add(h, new TH1F("pullyr", "relative pull y", 100, -25., 25.));
  add(h, new TH1F("pullzr", "relative pull z", 100, -25., 25.));
  add(h, new TH1F("vtxprob", "chisquared probability", 100, 0., 1.));

  addSP(h, new TH1F("eff0", "efficiency (zmatch)", 2, -0.5, 1.5));
  addSP(h, new TH1F("eff4", "efficiency (zmatch, ndof>4)", 2, -0.5, 1.5));
  addSP(h, new TH1F("effallTP", "efficiency (TPmatch, all)", 2, -0.5, 1.5));
  addSP(h, new TH1F("effselTP", "efficiency (TPmatch, selected)", 2, -0.5, 1.5));

  addSP(h, new TH1F("effsel3sigma", "efficiency (TPmatch, selected, 3 sigma z)", 2, -0.5, 1.5));

  add(h, new TH1F("zdistancenearest", "z-distance between generated nearest rec", 100, -0.1, 0.1));
  add(h, new TH1F("abszdistancenearest", "z-distance between generated and nearest rec", 1000, 0., 1.0));
  add(h, new TH1F("abszdistancenearestcum", "z-distance between generated and nearest rec", 1000, 0., 1.0));
  add(h, new TH1F("indexnearest", "index of nearest rec vertex", 20, 0., 20.));

  add(h, new TProfile("effvsptsq", "efficiency vs ptsq", 20, 0., 10000., 0, 1.));
  add(h, new TProfile("effvspt_hat", "efficiency vs pt_hat", 20, 0., 10000., 0, 1.));
  add(h, new TProfile("effvsnsimtrk", "efficiency vs # simtracks", 150, 0., 300., 0, 1.));
  add(h, new TProfile("effvsnsimtrkacc", "efficiency vs # simtracks", 150, 0., 300., 0, 1.));
  add(h, new TProfile("effvsnrectrk", "efficiency vs # rectracks", 50, 0., 1000., 0, 1.));
  add(h, new TProfile("effvsz", "efficiency vs z", 20, -20., 20., 0, 1.));
  add(h, new TProfile("effvsnsimvtx", "signal efficiency vs #PU", 30, 0., 300., 0., 1.));
  add(h, new TProfile("effvsnrecvtx", "signal efficiency vs #reconstructed vertices", 30, 0., 300., 0., 1.));
  addSP(h, new TProfile("effallvsnsimevt", "efficiency (all) vs #PU", 30, 0., 300., 0., 1.));
  addSP(h, new TProfile("effselvsnsimevt", "efficiency (selected) vs #PU", 30, 0., 300., 0., 1.));
  addSP(h, new TProfile("effsel3sigmavsnsimevt", "efficiency (selected) vs #PU", 30, 0., 300., 0., 1.));

  addSP(h, new TProfile("effSignalvsnsimtrkacc", "efficiency vs # simtracks in acceptance", 200, 0., 200., 0, 1.));
  addSP(
      h,
      new TProfile("effSignal4vsnsimtrkacc", "efficiency vs # simtracks in acceptance (ndof>4)", 200, 0., 200., 0, 1.));
  addSP(h, new TH1F("nsimtrkacc", "# of simulated tracks in acceptance", 200, 0., 200.));
  addSP(h, new TH1F("nsimtrkaccprim", "# of simulated primary tracks in acceptance", 200, 0., 200.));
  addSP(h, new TH1F("nsimtrkaccsec", "# of simulated secondary tracks in acceptance", 200, 0., 200.));

  add(h, new TH1F("effSignal", "signal efficiency", 2, -0.5, 1.5));                                            // simpv
  add(h, new TProfile("effSignalvsNpu0", "signal efficiency vs #PU", npubin2, 0., npumax, 0., 2.));            //simpv
  add(h, new TH1F("effSignal4", "signal efficiency (ndof>4)", 2, -0.5, 1.5));                                  // simpv
  add(h, new TProfile("effSignal4vsNpu0", "signal efficiency (ndof>4) vs #PU", npubin2, 0., npumax, 0., 2.));  //simpv

  add(h, new TProfile("effvsz2", "efficiency vs z (2mm)", 20, -20., 20., 0, 1.));
  add(h, new TProfile("effvsr", "efficiency vs r", 20, 0., 1., 0, 1.));
  add(h, new TProfile("xresvsntrk", "xresolution vs # vertex tracks", 40, 0., 200., 0, 0.01));
  add(h, new TProfile("yresvsntrk", "yresolution vs # vertex tracks", 40, 0., 200., 0, 0.01));
  add(h, new TProfile("zresvsntrk", "zresolution vs # vertex tracks", 40, 0., 200., 0, 0.01));

  add(h, new TH1F("nbtksinvtx2", "reconstructed tracks in vertex", 40, 0., 200.));
  add(h, new TH1F("nbtksinvtxPU2", "reconstructed tracks in vertex", 40, 0., 200.));

  add(h, new TH1F("xrec", "reconstructed x", 100, -0.1, 0.1));
  add(h, new TH1F("yrec", "reconstructed y", 100, -0.1, 0.1));
  add(h, new TH1F("zrec", "reconstructed z", 100, -20., 20.));

  add(h, new TH1F("xrecBeam", "reconstructed x - beam x", 100, -0.1, 0.1));
  add(h, new TH1F("yrecBeam", "reconstructed y - beam y", 100, -0.1, 0.1));
  add(h, new TH1F("zrecBeam", "reconstructed z - beam z", 100, -20., 20.));
  add(h, new TProfile("xrecBeamvszprof", "reconstructed x - beam x vs z-z0", 20, -20., 20., -0.1, 0.1));
  add(h, new TProfile("yrecBeamvszprof", "reconstructed y - beam y vs z-z0", 20, -20., 20., -0.1, 0.1));

  add(h, new TProfile("xrecBeamvsNdofprof", "reconstructed x - beam x vs ndof", 10, 0., 200., -0.1, 0.1));
  add(h, new TProfile("yrecBeamvsNdofprof", "reconstructed y - beam y vs ndof", 10, 0., 200., -0.1, 0.1));

  add(h, new TProfile("resxvsNdofprof", "reconstructed x - simulated x vs ndof", 10, 0., 200., -0.1, 0.1));
  add(h, new TProfile("resyvsNdofprof", "reconstructed y - simulated y vs ndof", 10, 0., 200., -0.1, 0.1));

  add(h, new TProfile("resxvsNdofSpread", "reconstructed x - simulated x vs ndof", 10, 0., 200., -0.1, 0.1, "S"));
  add(h, new TProfile("resyvsNdofSpread", "reconstructed y - simulated y vs ndof", 10, 0., 200., -0.1, 0.1, "S"));

  //add(h, new TH2F("zrecLS","reconstructed z vs LS",1000, 0. , 1000., 40., -20., 20.));
  add(h, new TProfile("dxcorr", "x correction vs event", 1000, 0., 10000.));
  add(h, new TProfile("dycorr", "y correction vs event", 1000, 0., 10000.));
  add(h, new TProfile("dzcorr", "z correction vs event", 1000, 0., 10000.));
  add(h, new TProfile("dxcorrLS", "x correction vs LS", 1000, 0., 1000.));
  add(h, new TProfile("dycorrLS", "y correction vs LS", 1000, 0., 1000.));
  add(h, new TProfile("dzcorrLS", "z correction vs LS", 1000, 0., 1000.));

  double dxbins[] = {0., 0.0005, 0.0010, 0.0015, 0.0020, 0.0030, 0.0040, 0.0060, 0.0080, 0.0100, 0.0200};
  double dx2bins[] = {0.,
                      pow(0.0005, 2),
                      pow(0.0010, 2),
                      pow(0.0015, 2),
                      pow(0.0020, 2),
                      pow(0.0030, 2),
                      pow(0.0040, 2),
                      pow(0.0060, 2),
                      pow(0.0080, 2),
                      pow(0.0100, 2),
                      pow(0.0200, 2)};
  if (DO_BEAMSPOT_ANALYSIS) {
    add(h, new TProfile("dxbinavg", "dx bin average", 10, dxbins));
    add(h, new TProfile("dybinavg", "dy bin average", 10, dxbins));
    add(h, new TProfile("dx2binavg", "dx2 bin average", 10, dx2bins));
    add(h, new TProfile("dy2binavg", "dy2 bin average", 10, dx2bins));
    add(h, new TH2F("xrecBeamvsdx", "reconstructed x - beam x vs resolution", 10, dxbins, 400, -0.1, 0.1));
    add(h, new TH2F("yrecBeamvsdy", "reconstructed y - beam y vs resolution", 10, dxbins, 400, -0.1, 0.1));
    add(h, new TProfile("xrecBeamvsdxprof", "reconstructed x - beam x vs resolution", 10, dxbins, -0.1, 0.1));
    add(h, new TProfile("yrecBeamvsdyprof", "reconstructed y - beam y vs resolution", 10, dxbins, -0.1, 0.1));
    add(h, new TProfile("xrecBeam2vsdx2prof", "(reconstructed x - beam x)^{2} vs resolution^{2}", 10, dx2bins));
    add(h, new TProfile("yrecBeam2vsdy2prof", "(reconstructed y - beam y)^{2} vs resolution^{2}", 10, dx2bins));
  }

  add(h, new TH1F("xrecb", "reconstructed x - beam x", 100, -0.01, 0.01));
  add(h, new TH1F("yrecb", "reconstructed y - beam y", 100, -0.01, 0.01));
  add(h, new TH1F("zrecb", "reconstructed z - beam z", 100, -20., 20.));
  add(h, new TH1F("xrec1", "reconstructed x", 100, -4, 4));
  add(h, new TH1F("yrec1", "reconstructed y", 100, -4, 4));  // should match the sim histos
  add(h, new TH1F("zrec1", "reconstructed z", 100, -80., 80.));
  add(h, new TH1F("xrec2", "reconstructed x", 100, -1, 1));
  add(h, new TH1F("yrec2", "reconstructed y", 100, -1, 1));
  add(h, new TH1F("zrec2", "reconstructed z", 200, -40., 40.));

  add(h, new TH1F("xrec3", "reconstructed x", 100, -0.1, 0.1));
  add(h, new TH1F("yrec3", "reconstructed y", 100, -0.1, 0.1));
  add(h, new TH1F("zrec3", "reconstructed z", 100, -20., 20.));
  add(h, new TH1F("zrec3a", "reconstructed z", 100, -1., 1.));

  add(h, new TH1F("xrecBeamPull", "normalized residuals reconstructed x - beam x", 100, -10, 10));
  add(h, new TH1F("yrecBeamPull", "normalized residuals reconstructed y - beam y", 100, -10, 10));
  add(h, new TH1F("zrecBeamPull", "normalized residuals reconstructed z - beam z", 100, -10, 10));
  add(h, new TH1F("zrecBeamPull0", "normalized residuals reconstructed z - beam z", 100, -10, 10));
  add(h, new TH1F("zrecBeamPull12", "normalized residuals reconstructed z - beam z (ndof>12)", 100, -10, 10));
  add(h, new TH1F("zrecBeamPullndoflt4", "normalized residuals reconstructed z - beam z (ndof<4)", 100, -10, 10));
  add(h, new TH1F("zrecBeamPullndof4-10", "normalized residuals reconstructed z - beam z (ndof=4-10)", 100, -10, 10));

  add(h, new TProfile("zvsls", "z vs ls", 200, 0., 2000., -20., 20.));
  add(h, new TProfile("zbeamvsls", "zbeam vs ls", 200, 0., 2000., -20., 20.));

  add(h, new TH1F("zdiffrec", "z-distance between vertices", nbinzdiffrec, -20., 20.));
  add(h, new TH1F("zdiffrechr", "z-distance between vertices", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("zdiffrec_tagged", "z-distance between tagged and other vertices", nbinzdiffrec, -20., 20.));
  add(h, new TH1F("zdiffrechr_tagged", "z-distance between tagged and other vertices", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("zdiffrecsel", "z-distance between selected vertices", nbinzdiffrec, -20., 20.));
  add(h, new TH1F("zdiffrecselhr", "z-distance between selected vertices", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("zdiffrecsel_tagged", "z-distance between tagged and selected vertices", nbinzdiffrec, -20., 20.));
  add(h, new TH1F("zdiffrecselhr_tagged", "z-distance between tagged and selected vertices", nbinzdiffrec, -2., 2.));
  add(h,
      new TH1F("zdiffrecselhr_tagged_lowndof",
               "z-distance between tagged low ndof and selected vertices",
               nbinzdiffrec,
               -2.,
               2.));
  add(h, new TH1F("zdiffrecselhr_gt1sigma", "z-distance between vertices |z|>1sigma_z", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("zdiffrecselhr_lt1sigma", "z-distance between vertices |z|<1sigma_z", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("zrec_tagged_hpu", "reconstructed z tagged", 200, -40., 40.));
  add(h, new TH1F("zrec_tagged_lowndof_hpu", "reconstructed z tagged low ndof", 200, -40., 40.));
  add(h, new TH1F("zrecsel", "z-zbeam selected vertices", 8 * nbinzdiffrec, -20., 20.));
  add(h, new TProfile("zrecselvsordprof", "z-zbeam selected vertices", 100, 0., 100., -20., 20.));
  add(h, new TH1F("nevtLPU", "event counter for zrec4vsLPU", 10, 0., 100.));
  add(h, new TH2F("zrecselvsLPU", "z-zbeam selected vertices", 10, 0., 100., 100, -20., 20.));
  add(h, new TH1F("zdiffrecselfold", "expected z-distance from z for ndof>4", nbinzdiffrec, -20., 20.));
  add(h, new TH1F("zdiffrecsel", "z-distance between ndof>4 vertices", nbinzdiffrec, -20., 20.));
  add(h, new TH1F("zdiffrecselhr", "z-distance between ndof>4 vertices", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("zdiffrecselhr_zgt1", "z-distance between ndof>4 vertices", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("zdiffrecselhr_tail", "z-distance between ndof>4 vertices", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("zdiffrecselhr_tailzgt1", "z-distance between ndof>4 vertices", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("zdiffrecselhr_hipu", "z-distance between ndof>4 vertices", nbinzdiffrec, -2., 2.));

  add(h, new TH1F("zdiffrecselhrzlt5", "z-distance between ndof>4 vertices for |z|<5 cm", nbinzdiffrec, -2., 2.));
  add(h,
      new TH2F("zdiffrecselhrvsL",
               "z-distance between ndof>4 vertices vs instbxlumi",
               100,
               0.,
               lumiHistoRange_,
               nbinzdiffrec,
               -2.,
               2.));

  add(h, new TH1F("zdiffrecseldeselhr", "z-distance between selected and deselected vertices", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("zdiffrecdeselhr", "z-distance between deselected vertices", nbinzdiffrec, -2., 2.));

  add(h, new TH1F("zdiffrec6", "z-distance between ndof>6 vertices", nbinzdiffrec, -20., 20.));
  add(h, new TH1F("zdiffrec7", "z-distance between ndof>7 vertices", nbinzdiffrec, -20., 20.));
  add(h, new TH1F("zdiffrec8", "z-distance between ndof>8 vertices", nbinzdiffrec, -20., 20.));
  add(h, new TH1F("zdiffrec12", "z-distance between ndof>12 vertices", 2 * nbinzdiffrec, -20., 20.));
  add(h, new TH1F("zdiffrec20", "z-distance between ndof>20 vertices", 2 * nbinzdiffrec, -20., 20.));
  add(h, new TH1F("zdiffrecp", "normalized z-distance between vertices", nbinzdiffrec, -20., 20.));
  if (DO_ZVSZ_ANALYSIS) {
    add(h, new TH2F("zvszrecsel", "z positions of multiple vertices", 100, -20., 20., 100, -20., 20.));
    add(h, new TH1F("zvszrecselref", "z positions", 100, -20., 20.));
    add(h, new TH2F("pzvspzsel", "prob(z) of multiple vertices", 100, 0., 1., 100, 0., 1.));
  }
  add(h, new TH1F("dzreccentral", "z-distance between vertices", 100, 0., 2.));
  add(h, new TProfile("ndofcentral", "<ndof> vs z-distance between vertices", 200, 0., 2., 0., 500.));
  add(h, new TProfile("ndoflocentral", "lower ndof vs z-distance between vertices", 200, 0., 2., 0., 500.));
  add(h, new TProfile("ndofhicentral", "higher ndof vs z-distance between vertices", 200, 0., 2., 0., 500.));
  add(h, new TProfile("ndofsumcentral", "sum of ndof vs z-distance between vertices", 200, 0., 2., 0., 500.));

  add(h, new TH1F("dzreccentral_tagged", "z-distance between vertices", 100, 0., 2.));
  add(h, new TProfile("ndofcentral_tagged", "<ndof> vs z-distance between vertices", 200, 0., 2., 0., 500.));
  add(h, new TProfile("ndoflocentral_tagged", "lower ndof vs z-distance between vertices", 200, 0., 2., 0., 500.));
  add(h, new TProfile("ndofhicentral_tagged", "higher ndof vs z-distance between vertices", 200, 0., 2., 0., 500.));
  add(h, new TProfile("ndofsumcentral_tagged", "sum of ndof vs z-distance between vertices", 200, 0., 2., 0., 500.));
  add(h,
      new TProfile(
          "ndofothercentral_tagged", "ndof of other vertex vs z-distance between vertices", 200, 0., 2., 0., 500.));
  add(h,
      new TProfile("ndofothercentral_tagged_lowndof",
                   "ndof of other vertex vs z-distance between vertices",
                   200,
                   0.,
                   2.,
                   0.,
                   500.));
  add(h,
      new TProfile("ndofother_tagged", "ndof of other vertex vs z-distance between vertices", 200, 0., 2., 0., 500.));
  add(h,
      new TProfile(
          "ndofother_tagged_lowndof", "ndof of other vertex vs z-distance between vertices", 200, 0., 2., 0., 500.));

  add(h, new TH1F("zdiffsimmerge", "z distance of merged or lost simulated vertices", 100, -5., 5.));
  add(h, new TH1F("zdiffsimfound", "z distance of found simulated vertices", 4000, -10., 10));
  add(h, new TH1F("zdiffsimall", "z distance of simulated vertices", 4000, -10., 10));
  add(h, new TH1F("zdiffsimfoundTP", "delta-zsim of found simulated distance with at least 4 tracks", 4000, -10., 10));
  add(h,
      new TH1F(
          "zdiffsimfound4TP", "delta-zsim of found simulated distance with at least 4 tracks, ndof>4", 4000, -10., 10));
  add(h,
      new TH1F(
          "zdiffsimfound4SignalTP", "delta-zsim of found simulated distance Signal-PU, ntrk>3,ndof>4", 4000, -10., 10));
  add(h, new TH1F("zdiffsimfoundTP2", "z distance of found simulated distance (2)", 4000, -10., 10));
  add(h, new TH1F("zdiffsimallTP", "z distance of simulated distance", 4000, -10., 10));
  add(h, new TH2F("zdiffrecvssim", "z distance rec vertices vs simulated distance", 100, -1., 1., 100, -1., 1.));
  add(h, new TH2F("zdiffrecvssimTP", "z distance rec vertices vs simulated distance", 100, -1., 1., 100, -1., 1.));
  add(h,
      new TH2F("zdiffrecselvssim", "z distance rec vertices (nd>4) vs simulated distance", 100, -1., 1., 100, -1., 1.));
  add(h,
      new TH2F("zdiffrecselvssimTP", "z distance rec vertices (nd>4) vs simulated distance", 100, -1., 1., 100, -1., 1.));
  add(h,
      new TH2F("zdiffrec12vssim", "z distance rec vertices (nd>12) vs simulated distance", 100, -1., 1., 100, -1., 1.));
  add(h,
      new TH2F(
          "zdiffrec12vssimTP", "z distance rec vertices (nd>12) vs simulated distance", 100, -1., 1., 100, -1., 1.));
  add(h, new TH2F("zdiffrecvsdsim", "z distance rec vertices vs simulated distance", 100, 0., 1., 100, -0.5, 0.5));
  add(h, new TH2F("zdiffrecvsdsimTP", "z distance rec vertices vs simulated distance", 100, -0.5, 0.5, 100, -0.5, 0.5));
  add(h,
      new TProfile("zdiffrecvsdsimprof", "z distance rec vertices vs simulated distance", 200, -1.0, 1.0, -0.5, 0.5));
  add(h,
      new TProfile("zdiffrecvsdsimTPprof", "z distance rec vertices vs simulated distance", 200, -1.0, 1.0, -0.5, 0.5));
  add(h, new TH1F("zdiffrec4f", "z-distance between ndof>4 vertices", 1000, -10., 10.));
  add(h,
      new TH1F(
          "zdiffsimfound4signalSimpv", "delta-zsim of found simulated distance Signal-PU, ndof>4", 4000, -10., 10));
  add(h,
      new TH1F("zdiffsimfound4signalSimpva",
               "delta-zsim of found simulated distance Signal-PU, ndof>4",
               nbinzdiffrec,
               0.,
               2.));

  add(h, new TH1F("mergerate_denominator_lin", "mergerate denominator", 100, 0., 0.1));
  add(h, new TH1F("mergerate_numerator_lin", "mergerate numerator", 100, 0., 0.1));
  add(h, new TH1F("mergerate_lin", "mergerate", 100, 0., 0.1));

  auto const log_mergez_bins = makeLogBins<float, 16>(1e-4, 1);
  add(h, new TH1F("mergerate_denominator", "mergerate denominator", 16, log_mergez_bins.data()));
  add(h, new TH1F("mergerate_numerator", "mergerate numerator", 16, log_mergez_bins.data()));
  add(h, new TH1F("mergerate", "mergerate", 16, log_mergez_bins.data()));

  add(h, new TH1F("mergerate_denominator_dqm", "mergerate denominator", 16, log_mergez_bins.data()));
  add(h, new TH1F("mergerate_numerator_dqm", "mergerate numerator", 16, log_mergez_bins.data()));
  add(h, new TH1F("mergerate_dqm", "mergerate", 16, log_mergez_bins.data()));

  float max_density = 10.0;  //vertices/cm
  if (DO_DENSITY_ANALYSIS) {
    for (unsigned int bin = 0; bin < dzbins_.size() - 1; bin++) {
      add(h,
          new TProfile(
              Form("dzzbin_recsel%d", bin), Form("%5.3f--%5.3f", dzbins_[bin], dzbins_[bin + 1]), 100, 0., max_density));
      add(h,
          new TProfile(Form("dzzbin_recsel_2_%d", bin),
                       Form("%5.3f--%5.3f", dzbins_[bin], dzbins_[bin + 1]),
                       100,
                       0.,
                       max_density));
      add(h,
          new TProfile(Form("dzzbin2_recsel%d", bin),
                       Form("%5.3f--%5.3f", dzbins_[bin], dzbins_[bin + 1]),
                       100,
                       0.,
                       max_density));
      add(h,
          new TProfile(Form("dzzbin_recreject%d", bin),
                       Form("%5.3f--%5.3f", dzbins_[bin], dzbins_[bin + 1]),
                       100,
                       0.,
                       max_density));
      add(h,
          new TH1F(
              Form("dzzbin_event_recsel%d", bin), Form("%5.3f--%5.3f", dzbins_[bin], dzbins_[bin + 1]), 200, 0., 20.));
      add(h,
          new TH1F(Form("dzzbin_event_recsel_2_%d", bin),
                   Form("%5.3f--%5.3f", dzbins_[bin], dzbins_[bin + 1]),
                   200,
                   0.,
                   20.));
      add(h,
          new TH1F(Form("dzzbin_event_recreject%d", bin),
                   Form("%5.3f--%5.3f", dzbins_[bin], dzbins_[bin + 1]),
                   200,
                   0.,
                   20.));
    }
  }
  add(h, new TProfile("inclusive", "inclusive", 100, 0., max_density));
  add(h, new TProfile("inclusive_2_", "inclusive", 100, 0., max_density));
  add(h, new TProfile("inclusive_ndoflt10", "inclusive ndof<10", 100, 0., max_density));
  add(h, new TProfile("inclusive_ndofgt10", "inclusive ndof>10", 100, 0., max_density));
  add(h, new TProfile("inclusive_ndofgt20", "inclusive ndof>20", 100, 0., max_density));
  add(h, new TH1F("inclusive_event", "inclusive_event", 200, 0., 20.));
  add(h, new TH1F("inclusive_event_2_", "inclusive_event", 200, 0., 20.));
  add(h, new TH1F("inclusive_event_ndoflt10", "inclusive_event_ndoflt10", 200, 0., 20.));
  add(h, new TH1F("inclusive_event_ndofgt10", "inclusive_event_ndofgt10", 200, 0., 20.));
  add(h, new TH1F("inclusive_event_ndofgt20", "inclusive_event_ndofgt20", 200, 0., 20.));

  // history plots
  nls_glb_ = 2;
  add(h, new TProfile("lPULS", "instantaneous BX lumi x 69.1 mb", int(0.001 * nls_glb_), 0., float(nls_glb_)));
  add(h,
      new TProfile(
          "nselvtxoverPUvsLS", "number of rec vtxs per PU", int(0.001 * nls_glb_), 0., float(nls_glb_), 0., 10.));
  add(h,
      new TProfile("nselvtxoverPU00-20vsLS",
                   "number of rec vtxs per PU (pu= 0-20)",
                   int(0.001 * nls_glb_),
                   0.,
                   float(nls_glb_),
                   0.,
                   10.));
  add(h,
      new TProfile("nselvtxoverPU20-40vsLS",
                   "number of rec vtxs per PU (pu= 20-40)",
                   int(0.001 * nls_glb_),
                   0.,
                   float(nls_glb_),
                   0.,
                   10.));
  add(h,
      new TProfile("nselvtxoverPU40-60vsLS",
                   "number of rec vtxs per PU (pu= 40-60)",
                   int(0.001 * nls_glb_),
                   0.,
                   float(nls_glb_),
                   0.,
                   10.));
  add(h,
      new TProfile("nselvtxoverPU60-80vsLS",
                   "number of rec vtxs per PU (pu= 60-80)",
                   int(0.001 * nls_glb_),
                   0.,
                   float(nls_glb_),
                   0.,
                   10.));
  add(h,
      new TProfile(
          "f10vsLS", "fraction of vertices with ndof < 10", int(0.001 * nls_glb_), 0., float(nls_glb_), 0., 10.));
  add(h,
      new TProfile(
          "f20vsLS", "fraction of vertices with ndof < 20", int(0.001 * nls_glb_), 0., float(nls_glb_), 0., 10.));
  add(h, new TProfile("sigmazLS", "beam spot #sigma_{z}", int(0.001 * nls_glb_), 0., float(nls_glb_), 0., 10.));
  add(h, new TProfile("zbeamLS", "beam spot z", int(0.001 * nls_glb_), 0., float(nls_glb_), -10., 10.));
  add(h, new TProfile("vtxndofvsLS", "average vertex ndof", int(0.001 * nls_glb_), 0., float(nls_glb_), 1., 1000.));
  add(h,
      new TProfile(
          "vtxndofzgt1vsLS", "average vertex ndof |z|>sigma_z", int(0.001 * nls_glb_), 0., float(nls_glb_), 1., 1000.));
  add(h, new TProfile("nvtxtailvsLS", "high nvtx tail", int(0.001 * nls_glb_), 0., float(nls_glb_), 0., 10.));

  add(h, new TH1F("nlPU", "PU (instantaneous BX lumi x 69.1 mb)", 200, 0., 200.));
  add(h, new TH1F("lBX", "inst. BX luminosity", 200, 0., lumiHistoRange_));
  add(h, new TH1F("sigmaZ", "sigma Z", 300, 0., 6.));
  add(h, new TH1F("sigmaz_event", "sigma Z event", 300, 0., 6.));
  add(h, new TH1F("sigmaz_event_over_sigmaZ", "sigma Z ratio", 300, 0., 3.));
  add(h, new TProfile("sigmaz_event_vs_beam", "event sigma z vs beam", 300, 0., 6., 0., 20.));
  add(h, new TProfile("sigmaz_event_vs_beam_tail", "event sigma z vs beam", 300, 0., 6., 0., 20.));
  add(h, new TProfile("sigmaz_event_vs_bx", "event sigma z /  beam vs bx", 3655, 0., 3564, 0., 20.));
  add(h, new TProfile("sigmaz_event_vs_bx_tail", "event sigma z /  beam vs bx", 3655, 0., 3564, 0., 20.));

  add(h, new TH1F("nselvtx", "# selected vertices", 200, 0., 200.));
  add(h, new TH1F("nselvtx1sigmaz", "# selected vertices |z|<1 sigma", 200, 0., 200.));
  add(h, new TH1F("nselvtx2sigmaz", "# selected vertices |z|<2 sigma", 200, 0., 200.));
  add(h, new TH1F("LPU", "expected PU", 200, 0., 200.));

  add(h, new TH1F("nselvtgt100", "# selected vertices", nvtxbin, 0., nvtxrange));

  add(h,
      new TProfile("LPUvsLPU",
                   "instBXLumiPU vs instBXLumiPU",
                   nvtxbin,
                   0.,
                   nvtxrange,
                   0.,
                   500.));  // keep track of the average pu within bins
  add(h,
      new TProfile(
          "LPUvsavgLPU", "instBXLumiPU vs avgBXLumiPU", nvtxbin, 0., nvtxrange, 0., 500.));  // weird things in 2018C
  add(h,
      new TProfile(
          "avgLPUvsLPU", "instBXLumiPU vs avgBXLumiPU", nvtxbin, 0., nvtxrange, 0., 500.));  // weird things in 2018C
  add(h,
      new TH2F("nselvtxvsLPU", "# selected vertices vs instBXLumiPU", nvtxbin, 0., nvtxrange, nvtxbin, 0., nvtxrange));
  add(h,
      new TH2F(
          "nselvtxvssimPU", "# selected vertices vs generated PU", nvtxbin, 0., nvtxrange, nvtxbin, 0., nvtxrange));
  add(h,
      new TH2F(
          "nselvtx_2vsLPU", "# selected vertices vs instBXLumiPU", nvtxbin, 0., nvtxrange, nvtxbin, 0., nvtxrange));
  add(h,
      new TH2F(
          "nselvtxvsavgLPU", "# selected vertices vs avg instBXLumiPU", nvtxbin, 0., nvtxrange, nvtxbin, 0., nvtxrange));
  add(h,
      new TH2F("nselvtxvsLPU1sigmaz",
               "# selected vertices vs instBXLumiPU (|z|<1#sigma)",
               nvtxbin,
               0.,
               nvtxrange,
               nvtxbin,
               0.,
               nvtxrange));
  add(h,
      new TH2F("nselvtxvsLPU2sigmaz",
               "# selected vertices vs instBXLumiPU (|z|<2#sigma)",
               nvtxbin,
               0.,
               nvtxrange,
               nvtxbin,
               0.,
               nvtxrange));
  add(h,
      new TH2F("nselvtxvsLPU3sigmaz",
               "# selected vertices vs instBXLumiPU (|z|<3#sigma)",
               nvtxbin,
               0.,
               nvtxrange,
               nvtxbin,
               0.,
               nvtxrange));
  add(h,
      new TH2F("nselvtxvsLPUgt1sigmaz",
               "# selected vertices vs instBXLumiPU (|z|>1#sigma)",
               nvtxbin,
               0.,
               nvtxrange,
               nvtxbin,
               0.,
               nvtxrange));
  add(h,
      new TH2F("nselvtxvsLPUgt2sigmaz",
               "# selected vertices vs instBXLumiPU (|z|>2#sigma)",
               nvtxbin,
               0.,
               nvtxrange,
               nvtxbin,
               0.,
               nvtxrange));
  add(h,
      new TH2F("nselvtxgt1vslt1sigmaz",
               "# selected vertices (|z|>1#sigma) vs (|z|>1#sigma)",
               nvtxbin,
               0.,
               nvtxrange,
               nvtxbin,
               0.,
               nvtxrange));

  add(h,
      new TH2F("nselvtxvsLPUfirstBX",
               "# selected vertices vs instBXLumiPU first BX in train",
               nvtxbin,
               0.,
               nvtxrange,
               nvtxbin,
               0.,
               nvtxrange));
  add(h,
      new TH2F("nselvtxvsLPUnotfirstBX",
               "# selected vertices vs instBXLumiPU without first BX in train",
               nvtxbin,
               0.,
               nvtxrange,
               nvtxbin,
               0.,
               nvtxrange));
  add(h,
      new TH2F("nselvtxvsLPUBX0",
               "# selected vertices vs instBXLumiPU in BX 0",
               nvtxbin,
               0.,
               nvtxrange,
               nvtxbin,
               0.,
               nvtxrange));
  add(h,
      new TH2F(
          "nselvtxvsLPUrlt200", "# selected vertices vs instBXLumiPU", nvtxbin, 0., nvtxrange, nvtxbin, 0., nvtxrange));
  add(h,
      new TH2F(
          "nselvtxvsLPUrgt200", "# selected vertices vs instBXLumiPU", nvtxbin, 0., nvtxrange, nvtxbin, 0., nvtxrange));

  add(h,
      new TH2F("nselvtxvsLPUndof10",
               "# selected vertices vs instBXLumiPU ndof>10",
               nvtxbin,
               0.,
               nvtxrange,
               nvtxbin,
               0.,
               nvtxrange));
  add(h,
      new TH2F("nselvtxvsLPUndoflt10",
               "# selected vertices vs instBXLumiPU 4 <ndof < 10",
               nvtxbin,
               0.,
               nvtxrange,
               nvtxbin,
               0.,
               nvtxrange));
  add(h,
      new TH2F("nselvtxvsLPUndof20",
               "# selected vertices vs instBXLumiPU ndof>20",
               nvtxbin,
               0.,
               nvtxrange,
               nvtxbin,
               0.,
               nvtxrange));
  add(h,
      new TH2F("nselvtxvsLPUpxytight",
               "# selected vertices vs instBXLumiPU pxy>0.1",
               nvtxbin,
               0.,
               nvtxrange,
               nvtxbin,
               0.,
               nvtxrange));
  add(h,
      new TH2F("nselvtxvsLPUaptsumgt2",
               "# selected vertices vs instBXLumiPU aptsum > 2",
               nvtxbin,
               0.,
               nvtxrange,
               nvtxbin,
               0.,
               nvtxrange));
  add(h,
      new TH2F("nselvtxvsLPUaptsumlt2",
               "# selected vertices vs instBXLumiPU aptsum < 2",
               nvtxbin,
               0.,
               nvtxrange,
               nvtxbin,
               0.,
               nvtxrange));
  add(h,
      new TH2F("nselvtxvsLPUptmax2gt04",
               "# selected vertices vs instBXLumiPU ptmax2 > 0.4",
               nvtxbin,
               0.,
               nvtxrange,
               nvtxbin,
               0.,
               nvtxrange));
  add(h,
      new TH2F("nselvtxvsLPUptmax2lt04",
               "# selected vertices vs instBXLumiPU ptmax2 < 0.4",
               nvtxbin,
               0.,
               nvtxrange,
               nvtxbin,
               0.,
               nvtxrange));

  /*
  add(h, new TH2F("ng1vsLPU","# of reconstructed vertices (ndof>4, gap>0.1 cm) vs PU", nvtxbin, 0., nvtxrange, nvtxbin, 0., nvtxrange));
  add(h, new TH2F("ng05vsLPU","# of reconstructed vertices (ndof>4, gap>0.05 cm) vs PU", nvtxbin, 0., nvtxrange, nvtxbin, 0., nvtxrange));
  add(h, new TH2F("ng03vsLPU","# of reconstructed vertices (ndof>4, gap>0.03 cm) vs PU", nvtxbin, 0., nvtxrange, nvtxbin, 0., nvtxrange));
  add(h, new TH2F("ng02vsLPU","# of reconstructed vertices (ndof>4, gap>0.02 cm) vs PU", nvtxbin, 0., nvtxrange, nvtxbin, 0., nvtxrange));

  add(h, new TProfile("ng1vsLPUprof","# of reconstructed vertices (ndof>4, gap>0.1 cm) vs PU", nvtxbin, 0., nvtxrange));
  add(h, new TProfile("ng05vsLPUprof","# of reconstructed vertices (ndof>4, gap>0.05 cm) vs PU", nvtxbin, 0., nvtxrange));
  add(h, new TProfile("ng03vsLPUprof","# of reconstructed vertices (ndof>4, gap>0.03 cm) vs PU", nvtxbin, 0., nvtxrange));
  add(h, new TProfile("ng02vsLPUprof","# of reconstructed vertices (ndof>4, gap>0.02 cm) vs PU", nvtxbin, 0., nvtxrange));
  */

  add(h, new TProfile("nselvtxvsLPUprof", "# selected vertices vs instBXLumiPU", nvtxbin, 0., nvtxrange));
  add(h, new TProfile("nselvtxvssimPUprof", "# selected vertices vs simulated PU", nvtxbin, 0., nvtxrange));
  add(h, new TProfile("LPUnselvtxvsLPUprof", "instBXLumiPU bin center for nselvtxvsLPUprof", nvtxbin, 0., nvtxrange));
  add(h, new TProfile("nselvtxvsavgLPUprof", "# selected vertices vs avg instBXLumiPU", nvtxbin, 0., nvtxrange));

  add(h,
      new TProfile(
          "nselvtxvsLPU1sigmazprof", "# selected vertices vs instBXLumiPU (|z|<1#sigma)", nvtxbin, 0., nvtxrange));
  add(h,
      new TProfile(
          "nselvtxvsLPU2sigmazprof", "# selected vertices vs instBXLumiPU (|z|<2#sigma)", nvtxbin, 0., nvtxrange));
  add(h,
      new TProfile(
          "nselvtxvsLPU3sigmazprof", "# selected vertices vs instBXLumiPU (|z|<3#sigma)", nvtxbin, 0., nvtxrange));
  add(h,
      new TProfile(
          "nselvtxvsLPUgt1sigmazprof", "# selected vertices vs instBXLumiPU (|z|>1#sigma)", nvtxbin, 0., nvtxrange));

  add(h, new TProfile("lPUbybx", "# expected PU by BX", 3565, 0., 3565));
  add(h, new TProfile("nselvtxbybx", "# selected vertices by BX", 3565, 0., 3565));
  add(h, new TProfile("nselvtxbyLS", "# selected vertices by LS", 1000, 0., 1000.));

  zbinmax_ = 10.;

  add(h, new TH1F("vtxz_recsel", "", nzbins_, -zbinmax_, zbinmax_));
  add(h, new TH1F("vtxz_recsel_hipu", "", nzbins_, -zbinmax_, zbinmax_));
  add(h, new TH1F("vtxz_recsel_2", "", nzbins_, -zbinmax_, zbinmax_));
  add(h, new TH1F("vtxz_recsel_2_hipu", "", nzbins_, -zbinmax_, zbinmax_));
  add(h, new TH1F("vtxz_recsel_3", "", nzbins_, -zbinmax_, zbinmax_));
  add(h, new TH1F("vtxz_recsel_3_hipu", "", nzbins_, -zbinmax_, zbinmax_));
  add(h, new TH1F("vtxz_tag", "", nzbins_, -zbinmax_, zbinmax_));
  add(h, new TH1F("vtxz_tag_tail", "", nzbins_, -zbinmax_, zbinmax_));
  add(h, new TH2F("vtxzbybx_recsel", "", nzbins_, -zbinmax_, zbinmax_, 3565, 0., 3565));
  add(h, new TH1F("vtxzhr_recsel", "", nzbins_ * 10, -zbinmax_, zbinmax_));
  add(h, new TH1F("vtxzhr_recsel_hipu", "", nzbins_ * 10, -zbinmax_, zbinmax_));
  add(h, new TH1F("vtxzhr_recsel_2", "", nzbins_ * 10, -zbinmax_, zbinmax_));
  add(h, new TH1F("vtxzhr_recsel_2_hipu", "", nzbins_ * 10, -zbinmax_, zbinmax_));
  add(h, new TH1F("vtxzhr_recsel_3", "", nzbins_ * 10, -zbinmax_, zbinmax_));
  add(h, new TH1F("vtxzhr_recsel_3_hipu", "", nzbins_ * 10, -zbinmax_, zbinmax_));

  /*
  add(h, new TH2F("vtxmult_recsel","", nzbins_, -zbinmax_, zbinmax_, 10, 0., 10));
  add(h, new TH2F("vtxmult_recsel_hipu","", nzbins_, -zbinmax_, zbinmax_, 10, 0., 10));
  add(h, new TH2F("vtxmult_recsel_tail1","", nzbins_, -zbinmax_, zbinmax_, 10, 0., 10));

  add(h, new TProfile("garg1","<g>", nzbins_, -zbinmax_, zbinmax_, 0., 100.));
  add(h, new TProfile("garg2","<g**2>", nzbins_, -zbinmax_, zbinmax_,  0., 100.));
  add(h, new TProfile("garg3","<g**3>", nzbins_, -zbinmax_, zbinmax_, 0., 100.));
  add(h, new TProfile("garg1_hipu","<g> (hi pu)", nzbins_, -zbinmax_, zbinmax_, 0., 100.));
  add(h, new TProfile("garg2_hipu","<g**2> (hi pu)", nzbins_, -zbinmax_, zbinmax_,  0., 100.));
  add(h, new TProfile("garg3_hipu","<g**3> (hi pu)", nzbins_, -zbinmax_, zbinmax_, 0., 100.));
  */

  double* dzbins = &dzbins_[0];  // I love C++
  unsigned int ndzbin = dzbins_.size() - 1;
  /*
  add(h, new TH2F("dzz_recsel","(selected)", ndzbin, dzbins, nzbins_, -zbinmax_, zbinmax_));
  add(h, new TH2F("dzz_recreject","(rejected)", ndzbin, dzbins, nzbins_, -zbinmax_, zbinmax_));
  */

  add(h, new TH1F("zdiffrecsel-dzbin", "z-distance between selected vertices", ndzbin, dzbins));

  // for MC, splitting rate in the same binning , see also zdiffrec4-dzbin
  add(h,
      new TH1F("zdiffrecselrealreal-dzbin",
               "z-distance between reconstructed real and other real ndof>4 vertices",
               ndzbin,
               dzbins));
  add(h,
      new TH1F("zdiffrecselanyfake-dzbin",
               "z-distance between reconstructed reoal or fake and fake ndof>4 vertices",
               ndzbin,
               dzbins));

  add(h,
      new TH1F("zdiffrecselsignalreal-dzbin",
               "z-distance between reconstructed signal and other real  selected vertices",
               ndzbin,
               dzbins));
  add(h,
      new TH1F("zdiffrecselsignalfake-dzbin",
               "z-distance between reconstructed signal and fake ndof>4 vertices",
               ndzbin,
               dzbins));
  add(h,
      new TH1F("zdiffrecselsignalsplit-dzbin",
               "z-distance between reconstructed signal and split-off selected vertices",
               ndzbin,
               dzbins));
  add(h,
      new TH1F("zdiffrecselsignalrealSimpv-dzbin",
               "z-distance between reconstructed signal and other real  ndof>4 vertices",
               ndzbin,
               dzbins));
  add(h,
      new TH1F("zdiffrecselsignalfakeSimpv-dzbin",
               "z-distance between reconstructed signal and fake ndof>4 vertices",
               ndzbin,
               dzbins));
  add(h,
      new TH1F("zdiffrecselpxysignalrealSimpv-dzbin",
               "z-distance between reconstructed signal and other real  ndof>4 vertices",
               ndzbin,
               dzbins));
  add(h,
      new TH1F("zdiffrecselpxysignalfakeSimpv-dzbin",
               "z-distance between reconstructed signal and fake ndof>4 vertices",
               ndzbin,
               dzbins));
  add(h,
      new TH1F("zdiffrecselsignalsplitSimpv-dzbin",
               "z-distance between reconstructed signal and split-off ndof>4 vertices",
               ndzbin,
               dzbins));
  add(h,
      new TH1F("zdiffrecselpxysignalsplitSimpv-dzbin",
               "z-distance between reconstructed signal and split-off ndof>4 vertices",
               ndzbin,
               dzbins));
  add(h,
      new TH1F("zdiffrecselPUrealPU-dzbin",
               "z-distance between reconstructed PU and other real PU ndof>4 vertices",
               ndzbin,
               dzbins));
  add(h,
      new TH1F(
          "zdiffrecselPUfake-dzbin", "z-distance between reconstructed PU and fake ndof>4 vertices", ndzbin, dzbins));

  add(h, new TH1F("zreciso", "zrec-zsim of isolated sim vertices", 500, 0., 1.));
  add(h, new TH1F("nreciso", "number of rec vertices near isolated sim vertices", 10, 0., 10.));
  add(h, new TH1F("zdiffsimisoall", "z distance of simulated distance (isolated pairs)", 500, 0., 10));
  add(h, new TH1F("zdiffsimiso0", "simulated z distance (isolated pairs, 0 rec)", 500, 0., 10));
  add(h, new TH1F("zdiffsimiso1", "simulated z distance (isolated pairs, 1 rec)", 500, 0., 10));
  add(h, new TH1F("zdiffsimiso2", "simulated z distance (isolated pairs, 2 rec)", 500, 0., 10));
  add(h, new TH1F("zdiffsimiso3", "simulated z distance (isolated pairs, 3 rec)", 500, 0., 10));
  add(h, new TH1F("zdiffreciso2", "reconstructed z distance (isolated pairs, 2 rec)", 500, 0., 10));
  add(h,
      new TH2F(
          "dzrecvssimiso2", "reconstructed vs simulated z distance (isolated pairs, 2 rec)", 200, 0., 2, 200, 0., 2.));

  const int nbinzdiff = 400;
  const float zdiffrange = 20.;
  add(h, new TH1F("zrec8r", "reconstructed (z-z0)*sqrt2 (ndof>8)", nbinzdiff, -zdiffrange, zdiffrange));
  add(h, new TH1F("zrec12r", "reconstructed (z-z0)*sqrt2 (ndof>12)", nbinzdiff, -zdiffrange, zdiffrange));
  add(h, new TH1F("zrec12q", "reconstructed (z-z0)/sqrt2 (ndof>12)", nbinzdiff, -zdiffrange, zdiffrange));

  add(h, new TH1F("zbarFakeEnriched", "zbar fake enriched", 100, -20., 20.));
  add(h, new TH1F("zbarFakeEnriched2", "zbar fake enriched (ndof>2)", 100, -20., 20.));
  add(h, new TH1F("zbarFakeEnriched5", "zbar fake enriched (ndof>5)", 100, -20., 20.));
  add(h, new TH1F("zbarFakeDepleted", "zbar fake depleted", 100, -20., 20.));

  if (DO_ZDIFFVSZ_ANALYSIS) {
    add(h, new TH2F("zdiffvsz", "z-distance vs z  (all)", nbinzdiff, -zdiffrange, zdiffrange, 30, -15., 15.));
    add(h, new TH2F("zdiffvszsel", "z-distance vs z (selected)", nbinzdiff, -zdiffrange, zdiffrange, 30, -15., 15.));
    add(h, new TH2F("zdiffvszselhr", "z-distance vs z (selected)", nbinzdiff, -2., 2., 30, -15., 15.));
    add(h, new TH2F("zdiffvszselhr_hipu", "z-distance vs z (ndof>4, hipu)", nbinzdiff, -2., 2., 30, -15., 15.));
    add(h, new TH2F("zdiffvsz6", "z-distance vs z (ndof>6)", nbinzdiff, -zdiffrange, zdiffrange, 30, -15., 15.));
    add(h, new TH2F("zdiffvsz7", "z-distance vs z (ndof>7)", nbinzdiff, -zdiffrange, zdiffrange, 30, -15., 15.));
    add(h, new TH2F("zdiffvsz8", "z-distance vs z (ndof>8)", nbinzdiff, -zdiffrange, zdiffrange, 30, -15., 15.));
    add(h, new TH2F("zdiffvsz12", "z-distance vs z (ndof>12)", nbinzdiff, -zdiffrange, zdiffrange, 30, -15., 15.));

    add(h, new TH2F("zdiffvszp", "z-distance vs z/sigmaZ", nbinzdiff, -zdiffrange, zdiffrange, 30, -5., 5.));
    add(h, new TH2F("zdiffvszpsel", "z-distance vs z/sigmaZ (selected)", nbinzdiff, -zdiffrange, zdiffrange, 30, -5., 5.));

    add(h,
        new TH2F("zdiffvszselNv2", "z-distance vs z (selected,Nv=2)", nbinzdiff, -zdiffrange, zdiffrange, 30, -15., 15.));
  }

  add(h, new TProfile("eff0vsntrec", "efficiency vs # reconstructed tracks", 50, 0., 50., 0, 1.));
  add(h, new TProfile("eff0vsntsel", "efficiency vs # selected tracks", 50, 0., 50., 0, 1.));
  add(h, new TH1F("ndofOverNtk", "ndof / ntk of candidates (ndof>4)", 100, 0., 2.));
  add(h, new TH1F("sumwoverntk", "sumw / ntk of candidates (ndof>4)", 100, 0., 1.));

  add(h, new TH1D("sumwvsz", "sumw (ndof>4)", nzbin_wide, -zrange_wide, zrange_wide));
  add(h, new TH1D("sumntkvsz", "sumntk (ndof>4)", nzbin_wide, -zrange_wide, zrange_wide));
  add(h, new TH1F("sumwoversumntkvsz", "sumw over sumntk (ndof>4)", nzbin_wide, -zrange_wide, zrange_wide));
  add(h,
      new TProfile(
          "sumwoverntkvsz", "sumw / ntk of candidates (ndof>4)", nzbin_wide, -zrange_wide, zrange_wide, 0., 1.));
  add(h,
      new TProfile(
          "sumwoverntkvszlo", "sumw / ntk of candidates (ndof<10)", nzbin_wide, -zrange_wide, zrange_wide, 0., 1.));
  add(h,
      new TProfile(
          "sumwoverntkvszhi", "sumw / ntk of candidates (ndof>20)", nzbin_wide, -zrange_wide, zrange_wide, 0., 1.));
  add(h,
      new TProfile(
          "sumwoverntkvsztp", "sumw / ntk of candidates (ndof>4)", nzbin_wide, -zrange_wide, zrange_wide, 0., 1.));
  add(h,
      new TProfile("sumwoverntkwgt05vsz",
                   "sumw / ntk(w>0.5) of candidates (ndof>4)",
                   nzbin_wide,
                   -zrange_wide,
                   zrange_wide,
                   0.,
                   1.));
  add(h, new TProfile("ntrkvsz", "<ntrk> vs z (ndof>4)", nzbin_wide, -zrange_wide, zrange_wide, 0., 1000.));
  add(h, new TProfile("ntrkpt1vsz", "<ntrk pt>1.0> vs z (ndof>4)", nzbin_wide, -zrange_wide, zrange_wide, 0., 1000.));
  add(h, new TProfile("ndofvsz0", "<ndof> vs z)", nzbin_wide, -zrange_wide, zrange_wide, 0., 1000.));
  add(h, new TProfile("ndofvsz4", "<ndof> vs z (ndof>4)", nzbin_wide, -zrange_wide, zrange_wide, 0., 1000.));
  add(h,
      new TH2F(
          "sumwoverntkvsz4", "sumw / ntk of candidates (ndof>4)", nzbin_wide, -zrange_wide, zrange_wide, 20, 0., 1.));

  add(h, new TH1F("nrecvtx", "# of reconstructed vertices", nrecmax, 0, float(nrecmax)));
  add(h, new TH1F("nrecvtx2", "# of reconstructed vertices with ndof>2", nrecmax, 0., float(nrecmax)));
  add(h, new TH1F("nrecvtx3", "# of reconstructed vertices with ndof>3", nrecmax, 0., float(nrecmax)));
  add(h, new TH1F("nrecvtx4", "# of reconstructed vertices with ndof>4", nrecmax, 0., float(nrecmax)));
  add(h, new TH1F("nrecvtx5", "# of reconstructed vertices with ndof>5", nrecmax, 0., float(nrecmax)));
  add(h, new TH1F("nrecvtx6", "# of reconstructed vertices with ndof>6", nrecmax, 0., float(nrecmax)));
  add(h, new TH1F("nrecvtx7", "# of reconstructed vertices with ndof>7", nrecmax, 0., float(nrecmax)));
  add(h, new TH1F("nrecvtx8", "# of reconstructed vertices with ndof>8", nrecmax, 0., float(nrecmax)));

  add(h, new TH1F("nrectrk", "# of reconstructed tracks", ntrkmax, 0, float(ntrkmax)));
  add(h, new TH1F("nsimtrk", "# of simulated tracks", 200, 0., 200.));
  add(h, new TH1F("nsimtrkSignal", "# of simulated tracks (Signal)", 200, -0.5, 199.5));
  add(h, new TH1F("nsimtrkPU", "# of simulated tracks (PU)", 100, -0.5, 199.5));

  add(h,
      new TH2F(
          "nrecvtxvsLPU", "# of reconstructed vertices vs PU", 100, 0., lumiPUHistoRange_, nrecmax, 0, float(nrecmax)));
  add(h,
      new TH2F("nrecvtx4vsLPU",
               "# of reconstructed vertices (ndof>4) vs PU",
               100,
               0.,
               lumiPUHistoRange_,
               nrecmax,
               0,
               float(nrecmax)));
  add(h,
      new TH2F("nrecvtx4vsavgLPU",
               "# of reconstructed vertices (ndof>4) vs PU",
               100,
               0.,
               lumiPUHistoRange_,
               nrecmax,
               0,
               float(nrecmax)));

  add(h, new TProfile("nrecvtxvsLPUprof", "# of reconstructed vertices vs PU", 100, 0., lumiPUHistoRange_, 0., 200.));
  add(h,
      new TProfile(
          "nrecvtx4vsLPUprof", "# of reconstructed vertices (ndof>4) vs PU", 100, 0., lumiPUHistoRange_, 0, 200.));
  add(h,
      new TProfile(
          "nrecvtx4vsavgLPUprof", "# of reconstructed vertices (ndof>4) vs PU", 100, 0., lumiPUHistoRange_, 0, 200.));
  add(h, new TProfile("LPUnrecvtx4vsavgLPUprof", "LPU vs PU for nrecvtx4", 100, 0., lumiPUHistoRange_, 0, 200.));
  add(h, new TProfile("LPUnrecvtx4vsLPUprof", "LPU vs LPU for nrecvtx4", 100, 0., lumiPUHistoRange_, 0, 200.));

  add(h, new TH2F("nbtksinvtxvsL", "reconstructed tracks in vertex", 100, 0., lumiHistoRange_, 80, 0., 80));
  add(h, new TH2F("nbtksinvtxvsLPU", "reconstructed tracks in vertex", 100, 0., lumiPUHistoRange_, 80, 0., 80));
  add(h, new TH2F("nbtksinvtx4vsL", "reconstructed tracks in vertex(ndof>4)", 100, 0., lumiHistoRange_, 80, 0., 80.));
  add(h,
      new TH2F("nbtksinvtx4vsLPU", "reconstructed tracks in vertex(ndof>4)", 100, 0., lumiPUHistoRange_, 80, 0., 80.));
  add(h, new TH1F("nsimtrkSimpv", "Number of simulated tracks", 1000, 0, 1000));

  add(h,
      new TProfile(
          "nrectrkvsLPUprof", "# of reconstructed tracks vs PU", 100, 0., lumiPUHistoRange_, 0., float(ntrkmax)));
  add(h,
      new TProfile("nseltrkvsLPUprof", "# of selected tracks vs PU", 100, 0., lumiPUHistoRange_, 0., float(ntrkmax)));
  add(h,
      new TProfile("nseltrkptlt04vsLPUprof",
                   "# of selected tracks with pt<0.4 vs PU",
                   100,
                   0.,
                   lumiPUHistoRange_,
                   0.,
                   float(ntrkmax)));
  add(h,
      new TProfile("nseltrkptlt03vsLPUprof",
                   "# of selected tracks with pt<0.3 vs PU",
                   100,
                   0.,
                   lumiPUHistoRange_,
                   0.,
                   float(ntrkmax)));
  add(h,
      new TProfile("nseltrkptlt02vsLPUprof",
                   "# of selected tracks with pt<0.2 vs PU",
                   100,
                   0.,
                   lumiPUHistoRange_,
                   0.,
                   float(ntrkmax)));

  add(h,
      new TH2F("nrecvtx4vsnrectrk",
               "# of reconstructed vertices (ndof>4) vs reconstruted tracks",
               100,
               0.,
               float(ntrkmax),
               nrecmax,
               0,
               float(nrecmax)));
  add(h,
      new TProfile("nrecvtx4vsnrectrkprof",
                   "# of reconstructed vertices (ndof>4) vs reconstruted tracks",
                   100,
                   0.,
                   float(ntrkmax),
                   0,
                   float(nrecmax)));

  add(h,
      new TH2F("nrecvtx4vsnseltrk",
               "# of reconstructed vertices (ndof>4) vs selected tracks",
               100,
               0.,
               float(ntrkmax),
               nrecmax,
               0,
               float(nrecmax)));
  add(h,
      new TH2F("nrecvtx4vsnseltrk_tail",
               "# of reconstructed vertices (ndof>4) vs selected tracks",
               100,
               0.,
               float(ntrkmax),
               nrecmax,
               0,
               float(nrecmax)));
  add(h,
      new TH2F("nrecvtx4vsnseltrk_hipu",
               "# of reconstructed vertices (ndof>4) vs selected tracks",
               100,
               0.,
               float(ntrkmax),
               nrecmax,
               0,
               float(nrecmax)));
  add(h,
      new TProfile("nrecvtx4vsnseltrkprof",
                   "# of reconstructed vertices (ndof>4) vs selected tracks",
                   100,
                   0.,
                   float(ntrkmax),
                   0,
                   float(nrecmax)));
  add(h,
      new TProfile("nrecvtx4vsnseltrkprof_hipu",
                   "# of reconstructed vertices (ndof>4) vs selected tracks",
                   100,
                   0.,
                   float(ntrkmax),
                   0,
                   float(nrecmax)));
  add(h,
      new TProfile("nrecvtx4vsnseltrkprof_tail",
                   "# of reconstructed vertices (ndof>4) vs selected tracks",
                   100,
                   0.,
                   float(ntrkmax),
                   0,
                   float(nrecmax)));

  add(h,
      new TH2F("nvtxgt1vslt1sigmaz",
               "# of selected vertices in z > 1sigma vs z < 1 sigma",
               nrecmax,
               0.,
               float(nrecmax),
               nrecmax,
               0,
               float(nrecmax)));
  add(h,
      new TH2F("nvtxgt1vslt1sigmaz_tail",
               "# of selected vertices in z > 1sigma vs z < 1 sigma",
               nrecmax,
               0.,
               float(nrecmax),
               nrecmax,
               0,
               float(nrecmax)));
  add(h,
      new TH2F("nvtxgt1vslt1sigmaz_tail1",
               "# of selected vertices in z > 1sigma vs z < 1 sigma",
               nrecmax,
               0.,
               float(nrecmax),
               nrecmax,
               0,
               float(nrecmax)));
  add(h,
      new TH2F("nvtxgt1vslt1sigmaz_hipu",
               "# of selected vertices in z > 1sigma vs z < 1 sigma",
               nrecmax,
               0.,
               float(nrecmax),
               nrecmax,
               0,
               float(nrecmax)));

  h["nsimtrk"]->StatOverflows(kTRUE);
  h["nsimtrkPU"]->StatOverflows(kTRUE);
  h["nsimtrkSignal"]->StatOverflows(kTRUE);

  add(h, new TH1F("nrectrk0vtx", "# rec tracks no vertex ", 100, -0.5, 99.5));
  add(h, new TH1F("nseltrk0vtx", "# rec tracks no vertex ", 100, -0.5, 99.5));
  add(h, new TProfile("trackAssEffvsPt", "track association efficiency vs pt", 20, 0., 100., 0, 1.));

  add(h, new TH1F("nseltrk", "# of reconstructed tracks selected for PV", ntrkmax, 0, float(ntrkmax)));
  add(h, new TH1F("zlost1", "z of lost vertices (bad z)", nzbin_wide_fine, -zrange_wide, zrange_wide));

  // properties of fake vertices  (MC only)_
  add(h, new TH1F("fakeVtxZNdofgt05", "z of fake vertices with ndof>0.5", nzbin_wide_fine, -zrange_wide, zrange_wide));
  add(h, new TH1F("fakeVtxZNdofgt2", "z of fake vertices with ndof>2", nzbin_wide_fine, -zrange_wide, zrange_wide));
  add(h, new TH1F("fakeVtxZNdofgt4", "z of fake vertices with ndof>4", nzbin_wide_fine, -zrange_wide, zrange_wide));
  add(h, new TH1F("fakeVtxZNdofgt8", "z of fake vertices with ndof>8", nzbin_wide_fine, -zrange_wide, zrange_wide));
  add(h, new TH1F("fakeVtxZ", "z of fake vertices", nzbin_wide_fine, -zrange_wide, zrange_wide));
  add(h, new TH1F("fakeVtxNdof", "ndof of fake vertices", 500, 0., 100.));
  add(h, new TH1F("fakeVtxNtrk", "number of tracks in fake vertex", 20, -0.5, 19.5));
  add(h, new TH1F("matchedVtxNdof", "ndof of matched vertices", 500, 0., 100.));

  add(h, new TH2I("cores_live", "cores vs LS", 1000, 0., 1000., 2560, 0., 2560.));


  
  // histograms of vertex properties, for fillVertexHistos(*)
  vector<string> vtypes = {"rec",
                           "selected",
                           "tagged",
			   "matchedvtx",
			   "matchedvtxsel",
			   "splitvtxsel",
			   "otherfakevtxsel",
			   "splitvtxselfromsignal",
			   "splitvtxselfrompu",
			   "matchedsignalvtxsel",
			   "matchedpuvtxsel",
                           //"hipu",
                           //"tail",
                           //"tailzgt1",
                           "sel"
                           //"selzgt1"
  };
  unsigned int nvtype = vtypes.size();
  for (unsigned int t = 0; t < nvtype; t++) {
    string st = vtypes[t];
    dir->mkdir(vtypes[t].c_str())->cd();
    //addn :  add a histogram in a subdirectory and use the subdirectory name in the map key
    addn(h, new TH1F("c2xy", "c2xy", 100, 0., 10.));
    addn(h, new TH1F("probxy", "probxy", 100, 0., 1.));
    addn(h, new TH1F("chi2", "chi**2", 100, 0., 100.));
    addn(h, new TH1F("chi2overntk", "chi**2", 100, 0., 4.));
    addn(h, new TH1F("r", "r", 100, 0., 0.1));
    addn(h, new TH1F("index", "index", 200, 0., 200.));
    addn(h, new TH2F("logndofvsindex", "log ndof vs index", 200, 0., 200., 100, 0., 3.));
    addn(h, new TProfile("ndofvsindex", "mean ndof vs index", 200, 0., 200.,  0., 1e9));
    addn(h, new TH1F("zpullbeam", "z / sigmaz(Beam)", 200, -5., 5.));
    addn(h, new TH1F("ndof", "ndof", 500, 0., 100));
    addn(h, new TProfile("ndofvspu", "ndof vs pu", 300, 0., 300, 0., 1000));
    addn(h, new TH1F("logndof", "ndof", 100, 0., 3.));  // log_10 : 1..1000.
    addn(h, new TH1F("trkweight", "track weight in vertex", 256, 0., 1.));
    addn(h, new TH1F("numlowttrk", "number of tracks with weight below 0.5", 20, 0., 20.));
    addn(h, new TH1F("fraclowttrk", "fraction of tracks with weight below 0.5", 100, 0., 0.5));
    addn(h, new TH1F("trkptnorm", "track pt (weighted)", 100, 0., 5.));
    addn(h, new TH1F("pzsum", "signed sum pz", 500, -100., 100));
    addn(h, new TH1F("ptsum", "signed sum pt", 500, -100., 100));
    addn(h, new TH1F("aptsum", "unsigned sum pt", 250, 0., 100));
    addn(h, new TH2F("sumpt2vssumpt", "sqrt(sum p_{t}^{2}) vs sum|pt|", 50, 0., 1000.,  50., 0.,  1000.));
    addn(h, new TH1F("sumpt2oversumpt", "sqrt(sum p_{t}^{2}) / sum|pt|", 200, 0., 1.));
    addn(h, new TH2F("sumpt2oversumptvssumpt2", "sqrt(sum p_{t}^{2}) vs sum|pt|", 100, 0., 1000.,  50., 0.,  1.));
    addn(h, new TH1F("waptsum", "unsigned weighted sum pt", 250, 0., 100));
    addn(h, new TH1F("apzsum", "unsigned sum pz", 500, -100., 100));
    addn(h, new TH1F("dznearest", "dz to nearest neighbour", nbinzdiffrec, -2., 2.));
    addn(h, new TH1F("apz", "pz asymmetry", 300, -1.5, 1.5));
    addn(h, new TH1F("apt", "pt asymmetry", 300, -1.5, 1.5));
    addn(h, new TH1F("wapz", "weighted pz asymmetry", 300, -1.5, 1.5));
    addn(h, new TH1F("wapt", "weighted pt asymmetry", 300, -1.5, 1.5));
    addn(h, new TH2F("azvsdz", "pz asymmetry vs dz", 200, -0.1, 0.1, 200, -1., 1.));
    addn(h, new TH2F("aptvsdz", "pt asymmetry vs dz", 200, -0.1, 0.1, 200, -1., 1.));
    addn(h, new TH2F("aptvsapz", "pt asymmetry vs pz asymmetry", 20, -1., 1., 20, -1., 1.));
    addn(h, new TH1F("sphericity", "sphericity", 200, 0., 1.));
    addn(h, new TH2F("sphericityvsntrk", "sphericity vs nt", 100, 0., 50., 200, 0., 1.));
    addn(h, new TH2F("vtxndfvsntrk", "ndof vs #tracks", 50, 0., 50, 100, 0., 100.));
    addn(h, new TH1F("avweight", "ndof / #tracks", 50, 0., 1.));
    addn(h, new TH2F("avweightvsndof", "ndof / #tracks", 50, 0., 1., 100, 0., 100.));
    addn(h, new TH1F("avweightX", "ndof / #tracks", 50, 0., 1.));
    addn(h, new TH1F("errx", "error x", 100, 0., 0.02)); // 
    addn(h, new TH1F("erry", "error y", 100, 0., 0.02));
    addn(h, new TH1F("errz", "error z", 100, 0., 0.05));
    //addn(h, new TH1F("ntrkvtx", "#tracks", 200, 0., 200.)); obsolete use nseltrkvtx define below

    addnSP(h, new TH1F("sumpt2", "sum pt**2", 200, 0., 1000));
    addnSP(h, new TH1F("logsumpt2", "log sum pt**2", 200, -1., 9.));
    addnSP(h, new TH1F("ptmax2", "second highest pt", 100, 0., 10.));
    addn(h, new TH1F("zvtx", "z", 400, -20., 20.));

    if (f4D_) {
      addn(h, new TH1F("tvtx", "vtx t", 400, -2., 2.));
      addn(h, new TH1F("terrvtx", "vtx t error", 150, 0., 0.15));
      addn(h, new TH1F("tvtx_fromtracks", "vtx t from tracks", 400, -2., 2.));
      addn(h, new TH1F("terrvtx_fromtracks", "vtx t error from tracks", 150, 0., 0.15));
      addn(h, new TH1F("tvtx_withtracks", "vtx t from tracks", 400, -2., 2.));  // reference for fromtracks
      addn(h, new TH1F("terrvtx_withtracks", "vtx t error from tracks", 150, 0., 0.15));
      addn(h, new TH1F("tvtx_fromtracks_pid", "vtx t from tracks (pid)", 400, -2., 2.));
      addn(h, new TH1F("terrvtx_fromtracks_pid", "vtx t error from tracks (pid)", 150, 0., 0.15));
      addn(h, new TH1F("tvtx_withtracks_pid", "vtx t with tracks (pid)", 400, -2., 2.)); // reference for fromtracks_pid
      addn(h, new TH1F("terrvtx_withtracks_pid", "vtx t error with tracks (pid)", 150, 0., 0.15));
    }

    // for matched vertices
    addnSP(h, new TH1F("nseltrkvtx", "#selected tracks", 200, 0., 200.));
    addnSP(h, new TH1F("zrecsim","zrec - zsim", 100, -0.4, 0.4));
    addnSP(h, new TH1F("zrecerr","zrec uncertainty", 100, 0.0, 0.01));
    addnSP(h, new TH1F("zrecsimHR","zrec - zsim", 200, -0.02, 0.02));
    addnSP(h, new TH1F("xrecsim","xrec - xsim", 100, -0.01, 0.01));
    addnSP(h, new TH1F("xrecerr","xrec uncertainty", 100, 0.0, 0.01));
    addnSP(h, new TH1F("yrecsim","yrec - ysim", 100, -0.01, 0.01));
    addnSP(h, new TH1F("yrecerr","yrec uncertainty", 100, 0.0, 0.01));
    addnSP(h, new TH1F("zrecsimpull","(zrec - zsim)/error", 100, -10, 10));
    addnSP(h, new TH1F("xrecsimpull","(xrec - xsim)/error", 100, -10, 10));
    addnSP(h, new TH1F("yrecsimpull","(yrec - ysim)/error", 100, -10, 10));

    if (f4D_){
      // for matched vertices with timing
      addnSP(h, new TH1F("ntimingvtx", "#timing tracks", 200, 0., 200.));
      addnSP(h, new TH1F("ntimingqual05vtx", "#timing tracks above quality threshold 0.5", 200, 0., 200.));
      addnSP(h, new TH1F("ntimingqual08vtx", "#timing tracks above quality threshold 0.8", 200, 0., 200.));
      addnSP(h, new TH1F("trecsim_fromtracks", "vertex time residual from tracks", 200, -0.1, 0.1));
      addnSP(h, new TH1F("trecerr_fromtracks", "vertex time error from tracks", 500, 0., 0.1));
      addnSP(h, new TH1F("trecsimpull_fromtracks", "tpull from tracks", 200, -10., 10.));
      addnSP(h, new TH1F("trecsim_fromtracks_pid", "vertex time residual from tracks (pid)", 200, -0.1, 0.1));
      addnSP(h, new TH1F("trecerr_fromtracks_pid", "vertex time error from tracks (pid)", 500, 0., 0.1));
      addnSP(h, new TH1F("trecsimpull_fromtracks_pid", "tpull from tracks (pid)", 200, -10., 10.));
      addnSP(h, new TH1F("trecsim_fromtracks_qual_pid", "vertex time residual from tracks with quality (pid)", 200, -0.1, 0.1));
      addnSP(h, new TH1F("trecerr_fromtracks_qual_pid", "vertex time error from tracks with quality (pid)", 500, 0., 0.1));
      addnSP(h, new TH1F("trecsimpull_fromtracks_qual_pid", "tpull from tracks with quality (pid)", 200, -10., 10.));
      addnSP(h,
	     new TH1F("trecsim_withtracks", "vertex time residual from vertex ", 200, -0.1, 0.1));  //reference for fromtracks
      addnSP(h, new TH1F("trecerr_withtracks", "vertex time error from vertex", 500, 0., 0.1));
      addnSP(h, new TH1F("trecsimpull_withtracks", "tpull from vertex", 200, -10., 10.));
      addnSP(h,
	     new TH1F("trecsim_withtracks_pid", "vertex time residual from vertex (pid)", 200, -0.1, 0.1));  //reference for fromtracks_pid
      addnSP(h, new TH1F("trecerr_withtracks_pid", "vertex time error from vertex (pid)", 500, 0., 0.1));
      addnSP(h, new TH1F("trecsimpull_withtracks_pid", "tpull from vertex (pid)", 200, -10., 10.));
      addnSP(h,
	     new TH1F("trecsim_withtracks_qual_pid", "vertex time residual from vertex (pid)", 200, -0.1, 0.1));  //reference for fromtracks_pid
      addnSP(h, new TH1F("trecerr_withtracks_qual_pid", "vertex time error from vertex (pid)", 500, 0., 0.1));
      addnSP(h, new TH1F("trecsimpull_withtracks_qual_pid", "tpull from vertex (pid)", 200, -10., 10.));
      
      addnSP(h, new TH1F("trecsim","trec - tsim", 200, -0.1, 0.1));
      addnSP(h, new TH1F("trecsimpull","(trec - tsim)/error", 100, -10, 10));
      addnSP(h, new TH1F("trecsimpullwide","(trec - tsim)/error", 100, -20, 20));
      addnSP(h, new TH1F("trecerr","trec uncertainty", 500, 0.0, 0.1));
      addnSP(h, new TProfile("trecerrvsntrkprof","trec uncertainty vs # of timing tracks", 100, 0., 100, 0.0, 0.1));
      addnSP(h, new TH2F("trecerrvsntrk","trec uncertainty vs # of timing tracks", 100, 0., 100, 200, -0.1, 0.1));

      addnSP(h, new TH1F("trecsim_sigmatlt01", "trec - tsim (sigma_t lt 0.1)", 200, -0.1, 0.1));
      addnSP(h, new TH1F("trecsimpull_sigmatlt01", "tpull from vertex (sigma_t lt 0.1)", 200, -10., 10.));
      addnSP(h, new TH2F("trecsimpullvserr", "tpull from vertex vs terr", 50, 0., 0.05, 200, -10., 10.));
      addnSP(h, new TProfile("trecsimpullsqvserr", "tpull**2 from vertex vs terr", 50, 0., 0.05, 0., 100.));
    }
    dir->cd();
  }// vtypes



  
  //  histograms of track quality for fillTrackHistos and fillVertexHistos
  //  note that track histograms are automatically also created and filled for all vtypes !
  vector<string> ttypes = {"trkall",
                           "trksel",
                           "sellost",
                           "wgt05",
                           "wlt05",
                           "thipu",
                           "ttail",
                           "ttailzgt1",
                           "hiptcentral",
                           "unmatchedVtx",
                           "seltpmatched",
                           "seltpmatchedSignal",
                           "seltpmatchedPU",
                           "seltpunmatched",
                           "seltpmatched_tgt1",
			   "MTDtail",
                           "fakevtxdriver",
                           "realvtxdriver",
			   "splitvtxselfromsignaldriver",
			   "matchedsignalvtxseldriver",
			   "highetadriver",
			   "alletadriver",
			   "missing_inner0",
			   "missing_inner1",
			   "missing_inner2",
			   "nbarrel_lt2",
			   "nbarrel_eq2",
                           "faraway"};
  unsigned int nttype = ttypes.size();

  for (unsigned int t = 0; t < nttype + nvtype; t++) {

    if (t < nttype)
      { // track-only, make a new directory
	dir->mkdir(ttypes[t].c_str())->cd();
      }
    else
      { // track histos for a histogram class, re-use
	dir->cd(vtypes[t-nttype].c_str());
      }
    // string st = (t < nttype) ? ttypes[t] : vtypes[t - nttype];
    addn(h, new TH1F("ztrk", "z0 ", 200, -40., 40.));
    if (f4D_) {
      addn(h, new TH1F("t0trk", "track t0", 400, -2., 2.));
      addn(h, new TH1F("t0errtrk", "track sigma t0", 150, 0., 0.3));
      addn(h, new TH1F("t0qualtrk", "track time quality", 200, -1., 1.));
      addn(h, new TH1F("trestrk", "track time residual", 400, -2., 2.));
      addn(h, new TH1F("tpulltrk", "track--vertex pull", 200, -10., 10.));
    }
    addn(h, new TH1F("phi", "phi", 80, -3.14159, 3.14159));
    addn(h, new TH1F("eta", "eta ", 80, -4., 4.));
    addn(h, new TH1F("pt", "pt ", 100, 0., 5.));
    addn(h, new TH1F("logpt", "logpt ", 60, -1., 5.));
    addn(h, new TH2F("logpteta", "logpt vs eta ", 80, -4., 4., 60, -1., 5.));
    addn(h, new TH2F("z-eta", "z vs eta ", 48, -2.4, 2.4, 60, -15., 15.));
    addn(h, new TH2F("z-pt", "pt vs z", 60, -15., 15., 100, 0., 5.));
    addn(h, new TH2F("z-logpt", "log pt vs z", 60, -15., 15., 60, -1., 5.));
    addn(h, new TH1F("ptfwd", "pt (forward)", 100, 0., 5.));
    addn(h, new TH1F("ptcentral", "pt (central)", 100, 0., 5.));
    addn(h, new TH1F("found", "found hits", 20, 0., 20.));
    addn(h, new TH1F("lost", "lost hits", 20, 0., 20.));
    addn(h, new TH1F("validfraction", "fraction of valid hits", 50, 0., 1.));
    addn(h, new TH1F("nchi2", "normalized track chi2", 100, 0., 20.));
    addn(h, new TProfile("nchi2vsz", "normalized track chi2", 120, -30., 30., 0., 20.));
    addn(h, new TH1F("rstart", "start radius", 100, 0., 20.));
    addn(h, new TH1F("expectedInner", "expected inner hits ", 10, 0., 10.));
    addn(h, new TH1F("expectedOuter", "expected outer hits ", 10, 0., 10.));
    addn(h, new TH1F("logtresxy", "log10(track r-phi resolution/um)", 100, 0., 5.));
    addn(h, new TH1F("logtresz", "log10(track z resolution/um)", 100, 0., 5.));
    addn(h, new TH1F("tpullxy", "track r-phi pull", 100, -10., 10.));
    addn(h, new TProfile("tpullxyvsz", "track r-phi pull", 120, -30., 30., -10., 10.));
    addn(h, new TProfile("tpullxyvseta", "track r-phi pull", 48, -2.4, 2.4, -10., 10.));
    addn(h, new TProfile("tpullzvsz", "track z pull", 120, -30., 30., -10., 10.));
    addn(h, new TProfile("tpullzvseta", "track z pull", 48, 2.4, 2.4, -10., 10.));
    addn(h, new TH2F("lvseta", "cluster length vs eta", 80, -4., 4., 20, 0., 20));
    addn(h, new TH2F("lvstanlambda", "cluster length vs tan lambda", 60, -6., 6., 20, 0., 20));
    addn(h, new TH1F("longestbarrelhit", "longest barrel cluster", 40, 0., 40.));
    // note : histograms involving the relation wrt a rec vertex are usually filled with the tagged vertex
    // for pile-up this is usually the wrong one, expect a peak on top of a more or less flat backgound
    // useful mostly for MC without PU
    addn(h, new TH1D("zrestrk", "z-residuals (track vs vertex)", 200, -2., 2.));
    addn(h,
        new TH2F("zrestrkvsphi", "z-residuals (track - vertex)", 12, -3.14159, 3.14159, 100, -1., 1.));
    addn(h, new TH2F("zrestrkvseta", "z-residuals (track - vertex)", 16, -4., 4., 100, -1., 1.));
    addn(h, new TH2F("zrestrkvsz", "z-residuals (track - vertex) vs z", 100, -20., 20., 100, -1., 1.));
    addn(h,
        new TH2F("zpulltrkvsphi",
                 "normalized z-residuals (track - vertex)",
                 12,
                 -3.14159,
                 3.14159,
                 100,
                 -5.,
                 5.));
    addn(h,
        new TH2F(
            "zpulltrkvseta", "normalized z-residuals (track - vertex)", 12, -3., 3., 100, -5., 5.));
    addn(h,
        new TH2F(
            "zpulltrkvsz", "normalized z-residuals (track - vertex) vs z", 100, -20., 20., 100, -5., 5.));
    addn(h, new TH1D("zpulltrk", "normalized z-residuals (track vs vertex)", 100, -5., 5.));
    addn(h, new TH1D("zerrtrk", "z-resolution (excluding beam)", 100, 0., 1.));
    addn(h, new TH1D("zerrtrk_withbeam", "z-resolution (including beam)", 100, 0., 1.));
    addn(h, new TH1D("nbarrelhits", "number of pixel barrel hits", 10, 0., 10.));
    addn(h, new TH1D("nbarrelLayers", "number of pixel barrel layers", 10, 0., 10.));
    addn(h, new TH1D("nPxLayers", "number of pixel layers (barrel+endcap)", 10, 0., 10.));
    addn(h, new TH1D("nSiLayers", "number of Tracker layers ", 30, 0., 30.));
    addn(h, new TH1D("n3dLayers", "number of 3d Tracker layers ", 20, 0., 20.));
    addn(h, new TH1D("nOT3dLayers", "number of 3d Outer Tracker layers ", 20, 0., 20.));
    addn(h, new TH1D("trackAlgo", "track algorithm ", 30, 0., 30.));
    addn(h, new TH2F("nPxLayersVsPt", "number of pixel layers (barrel+endcap)", 8, 0., 8., 10, 0., 10.));
    addn(h, new TH1D("trackQuality", "track quality ", 7, -1., 6.));
    addn(h, new TH1D("missing_inner", "number of missing pixel hits", 10, 0., 10.)); // actually lost
    addn(h, new TH1D("missing_outer", "number of missing strip hits", 20, 0., 20.));
    addn(h, new TH1D("missing_fraction", "fraction of missing hits", 10, 0., 1.));
    addn(h, new TH1D("npxminmiss", "number of pixel layers minus missing hits", 12, -2., 10.));
    
    // some histograms wrt to MC truth
    addn(h, new TH1D("tkzrecsim", "zrec-zsim", 200, -0.1, 0.1));
    addn(h, new TH1D("tkzrecsimpull", "(zrec-zsim) / #sigma_{z}", 300, -15., 15.));
    addn(h, new TProfile("tkzrecsimpullsqvseta", "zrec-zsim pull**2 vs eta", 100, -4., 4., 0., 1000.));
    addn(h, new TH2F("tkzrecsimpullvseta", "zrec-zsim pull vs eta", 100, -4., 4., 200, -20., 20.));
    addn(h, new TProfile("tkzrecsimvsz",   "zrec-zsim vs z", 100, -15., 15., -1000., 1000.));
    addn(h, new TProfile("tkzrecsimvseta", "zrec-zsim vs eta", 100, -4., 4., -1000., 1000.));
    addn(h, new TProfile("tkzrecsimvsetaz", "zrec-zsim vs eta with z-flip", 100, -4., 4., -1000., 1000.));
    addn(h, new TH2F("tkzrecsimvseta2d", "zrec-zsim vs eta", 40, -2., 2., 200, -4., 4.));

    addn(h, new TH1D("tktrecsim", "trec-tsim", 200, -0.3, 0.3));
    addn(h, new TH1D("tktrecsimpull", "(trec-tsim)/terr", 200, -10., 10.));
    addn(h, new TH1D("tktrecsim_pid", "trec-tsim (pid)", 200, -0.3, 0.3));
    addn(h, new TH1D("tktrecsimpull_pid", "(trec-tsim)/terr (pid)", 200, -10., 10.));
    addn(h, new TH1D("tktrecsimpullwide", "trec-tsim", 200, -20., 20.));
    addn(h, new TH2F("tktrecsimvseta2d", "trec-tsim vs eta", 40, -4., 4., 200, -1.0, 1.0));
    addn(h, new TProfile("tktrecsimpullsqvserr", "((trec-tsim)/terr)**2 vs terr", 50, 0., 0.5, 0., 100.0));
    addn(h, new TH2F("tktrecsimpullvserr", "(trec-tsim)/terr vs terr", 50, 0., 0.5, 50, 0., 100.0));
    addn(h, new TProfile("tkzrecsimvslogpt",   "(zrec-zsim) vs log pt", 100, -1., 5., -1000., 1000.));
    dir->cd();
  }  // track types
  
  // vertex composition histograms
  addSP(h, new TH1F("pt_majority_frac", "min(pt,10.) fraction of majority tracks", 50, 0., 1.));
  addSP(h, new TH1F("pt_minority_frac", "min(pt,10.) fraction of minority tracks", 50, 0., 1.));
  addSP(h, new TH1F("pt_unmatched_frac", "min(pt,10.) fraction of unmatched tracks", 50, 0., 1.));

  addSP(h,
        new TProfile("pt_majority_frac_vsz",
                     "min(pt,10.) fraction of majority tracks",
                     nzbin_normal_fine,
                     -zrange_normal,
                     zrange_normal,
                     0.,
                     1.));
  addSP(h,
        new TProfile("pt_minority_frac_vsz",
                     "min(pt,10.) fraction of minority tracks",
                     nzbin_normal_fine,
                     -zrange_normal,
                     zrange_normal,
                     0.,
                     1.));
  addSP(h,
        new TProfile("pt_unmatched_frac_vsz",
                     "min(pt,10.) fraction of unmatched tracks",
                     nzbin_normal_fine,
                     -zrange_normal,
                     zrange_normal,
                     0.,
                     1.));

  addSP(h, new TH1F("nt_majority_frac", "#track fraction of majority tracks", 50, 0., 1.));
  addSP(h, new TH1F("nt_minority_frac", "#track fraction of minority tracks", 50, 0., 1.));
  addSP(h, new TH1F("nt_unmatched_frac", "#track fraction of unmatched tracks", 50, 0., 1.));

  addSP(h,
        new TProfile("nt_majority_frac_vsz",
                     "#track fraction of majority tracks",
                     nzbin_normal_fine,
                     -zrange_normal,
                     zrange_normal,
                     0.,
                     1.));
  addSP(h,
        new TProfile("nt_minority_frac_vsz",
                     "#track fraction of minority tracks",
                     nzbin_normal_fine,
                     -zrange_normal,
                     zrange_normal,
                     0.,
                     1.));
  addSP(h,
        new TProfile("nt_unmatched_frac_vsz",
                     "#track fraction of unmatched tracks",
                     nzbin_normal_fine,
                     -zrange_normal,
                     zrange_normal,
                     0.,
                     1.));

  addSP(h, new TH1F("wt_majority_frac", "weight fraction of majority tracks", 50, 0., 1.));
  addSP(h, new TH1F("wt_minority_frac", "weight fraction of minority tracks", 50, 0., 1.));
  addSP(h, new TH1F("wt_unmatched_frac", "weight fraction of unmatched tracks", 50, 0., 1.));

  addSP(h, new TH1F("nsimevt", "number of sim events with tracks in rec vertex", 20, 0., 20.));
  addSP(h, new TH1F("nsimevt_nt2", "number of sim events with at least 2 tracks in rec vertex", 20, 0., 20.));
  addSP(h, new TH1F("nsimevtMTD", "number of sim events with MTD track in rec vertex", 20, 0., 20.));
  const float nMTDTDRbins[] = {0, 4, 8, 15, 25, 50, 500};
  add(h, new TH1F("MTDTDR", "number of PU tracks per primary vertex", 6, nMTDTDRbins));

  add(h, new TH1F("wosfrac", "wos fraction of matched vertex", 256, 0., 1.));
  add(h, new TH1F("nwosmatch", "number of wos matched rec vertices per sim-vertex(>1 = splitting)", 5, 0., 5.));
  //
  add(h, new TH1F("trackWt", "track weight in vertex", 256, 0., 1.));
  add(h, new TH1F("allweight", "track weight in vertex", 256, 0., 1.));
  add(h, new TH1F("minorityweight", "minority track weight", 256, 0., 1.));
  add(h, new TH1F("minorityaweight", "minority(a) track weight", 256, 0., 1.));
  add(h, new TH1F("minoritybweight", "minority(b) track weight", 256, 0., 1.));
  add(h, new TH1F("majorityweight", "majority track weight", 256, 0., 1.));
  add(h, new TH1F("unmatchedweight", "unmatched track weight", 256, 0., 1.));
  add(h, new TH1F("unmatchedvtxtrkweight", "unmatched vertex track weight", 256, 0., 1.));

  add(h, new TProfile("ntrkwgt05vsz", "number of w>0.5 tracks", nzbin_wide_fine, -zrange_wide, zrange_wide, 0., 1000.));
  add(h, new TProfile("ftrkwgt05vsz", "fraction of w>0.5 tracks", nzbin_wide_fine, -zrange_wide, zrange_wide, 0., 1.));
  add(h, new TProfile("trackwtvsz", "track weight vs z (all)", nzbin_wide_fine, -zrange_wide, zrange_wide, 0., 1.));
  add(h,
      new TProfile("trackwtgt05vsz", "track weight vs z (w>0.5)", nzbin_wide_fine, -zrange_wide, zrange_wide, 0.5, 1.));
  add(h, new TProfile("allweightvsz", "track weight vs z (all)", nzbin_wide_fine, -zrange_wide, zrange_wide, 0., 1.));
  add(h,
      new TProfile("minorityweightvsz", "minority track weight", nzbin_wide_fine, -zrange_wide, zrange_wide, 0., 1.));
  add(h,
      new TProfile(
          "minorityaweightvsz", "minority(a) track weight", nzbin_wide_fine, -zrange_wide, zrange_wide, 0., 1.));
  add(h,
      new TProfile(
          "minoritybweightvsz", "minority(b) track weight", nzbin_wide_fine, -zrange_wide, zrange_wide, 0., 1.));
  add(h,
      new TProfile("majorityweightvsz", "majority track weight", nzbin_wide_fine, -zrange_wide, zrange_wide, 0., 1.));
  add(h,
      new TProfile("unmatchedweightvsz", "unmatched track weight", nzbin_wide_fine, -zrange_wide, zrange_wide, 0., 1.));
  add(h,
      new TProfile(
          "minorityfractionvsz", "minority track fraction", nzbin_wide_fine, -zrange_wide, zrange_wide, 0., 1.));
  add(h,
      new TProfile(
          "minorityafractionvsz", "minority(a) track fraction", nzbin_wide_fine, -zrange_wide, zrange_wide, 0., 1.));
  add(h,
      new TProfile(
          "minoritybfractionvsz", "minority(b) track fraction", nzbin_wide_fine, -zrange_wide, zrange_wide, 0., 1.));
  add(h,
      new TProfile(
          "majorityfractionvsz", "majority track fraction", nzbin_wide_fine, -zrange_wide, zrange_wide, 0., 1.));
  add(h,
      new TProfile("unmatchedfractionvsz",
                   "unmatched track fraction (in vtx)",
                   nzbin_wide_fine,
                   -zrange_wide,
                   zrange_wide,
                   0.,
                   1.));
  add(h,
      new TProfile("unmatchedvtxtrkfractionvsz",
                   "unmatched vertex track fraction",
                   nzbin_wide_fine,
                   -zrange_wide,
                   zrange_wide,
                   0.,
                   1.));
  add(h,
      new TProfile("unmatchedvtxtrkweightvsz",
                   "unmatched vertex track weight",
                   nzbin_wide_fine,
                   -zrange_wide,
                   zrange_wide,
                   0.,
                   1.));

  add(h,
      new TProfile("matchedselfractionvsz",
                   "matched fraction of selected tracks",
                   nzbin_wide_fine,
                   -zrange_wide,
                   zrange_wide,
                   0.,
                   1.));
  add(h,
      new TProfile("unmatchedselfractionvsz",
                   "unmatched fraction of selected tracks",
                   nzbin_wide_fine,
                   -zrange_wide,
                   zrange_wide,
                   0.,
                   1.));
  add(h,
      new TProfile("matchedallfractionvsz",
                   "matched fraction of reco tracks",
                   nzbin_wide_fine,
                   -zrange_wide,
                   zrange_wide,
                   0.,
                   1.));
  add(h,
      new TProfile("unmatchedallfractionvsz",
                   "unmatched fraction of reco tracks",
                   nzbin_wide_fine,
                   -zrange_wide,
                   zrange_wide,
                   0.,
                   1.));

  h["nrectrk"]->StatOverflows(kTRUE);
  h["nrectrk"]->StatOverflows(kTRUE);
  h["nrectrk0vtx"]->StatOverflows(kTRUE);
  h["nseltrk0vtx"]->StatOverflows(kTRUE);
  h["nseltrk"]->StatOverflows(kTRUE);
  h["nbtksinvtx"]->StatOverflows(kTRUE);
  h["nbtksinvtxPU"]->StatOverflows(kTRUE);
  h["nbtksinvtx2"]->StatOverflows(kTRUE);
  h["nbtksinvtxPU2"]->StatOverflows(kTRUE);

  const int vmax = 250;

  // for analyzeVertexCollectionSimPvNoSimTracks
  add(h, new TH1F("simPU", "number of simulated vertices", vmax, 0., float(vmax)));  // possible duplication -> "npu"
  add(h, new TH1F("dzrecsim", "distance of sim vertex", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("dzrecsimother", "distance of sim vertex not matched to this recvertex", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("dzrecsimptmax2lt04", "distance of sim vertex ptmax2 < 0.4", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("dzrecsimptmax2gt04", "distance of sim vertex ptmax2 > 0.4", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("dzrecsimmin", "distance of nearest sim vertex", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("dzrecsimminptmax2lt04", "distance of nearest sim vertex ptmax2 < 0.4", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("dzrecsimminptmax2gt04", "distance of nearest sim vertex ptmax2 > 0.4", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("zmatched01", "z of matched vertices", 400, -20., 20.));
  add(h, new TH1F("zmatched02", "z of matched vertices (2mm)", 400, -20., 20.));
  add(h, new TH1F("zunmatched01", "z of unmatched vertices", 400, -20., 20.));
  add(h, new TH1F("zunmatched01ptmax2lt04", "z of unmatched vertices with ptmax2 <0.4", 400, -20., 20.));
  add(h, new TH1F("zunmatched02", "z of unmatched vertices", 400, -20., 20.));
  add(h, new TH1F("zunmatched02ptmax2lt04", "z of unmatched vertices with ptmax2 <0.4", 400, -20., 20.));
  add(h, new TH1F("ptmax2matched01", "ptmax2 of matched vertices", 100, 0., 10.));
  add(h, new TH1F("ptmax2matched02", "ptmax2 of matched vertices", 100, 0., 10.));
  add(h, new TH1F("ptmax2unmatched01", "ptmax2 of unmatched vertices", 100, 0., 10.));
  add(h, new TH1F("ptmax2unmatched02", "ptmax2 of unmatched vertices", 100, 0., 10.));

  // pile-up and track assignment related histograms (MC with TP)
  add(h, new TH1F("npu", "Number of simulated vertices", vmax, 0., float(vmax)));
  add(h, new TH1F("npu0", "Number of simulated vertices", vmax, 0., float(vmax)));

  add(h,
      new TH2F("nrecvsnpu",
               "# of reconstructed vertices vs number sim vertices",
               vmax,
               0.,
               float(vmax),
               vmax,
               0.,
               float(vmax)));
  add(h,
      new TH2F("nrec4vsnpu",
               "# of reconstructed vertices vs number sim vertices",
               vmax,
               0.,
               float(vmax),
               vmax,
               0.,
               float(vmax)));
  add(h,
      new TProfile("nrec4vsnpuprof",
                   "# of reconstructed vertices vs number sim vertices",
                   vmax,
                   0.,
                   float(vmax),
                   0.,
                   float(2 * vmax)));

  add(h,
      new TProfile("nsplitvtx0vsnpuprof",
                   "# of split vertices vs number sim vertices",
                   vmax,
                   0.,
                   float(vmax),
                   0.,
                   float(2 * vmax)));
  add(h,
      new TProfile("nsplitvtx4vsnpuprof",
                   "# of split vertices  (ndof>4) vs number sim vertices",
                   vmax,
                   0.,
                   float(vmax),
                   0.,
                   float(2 * vmax)));
  add(h,
      new TProfile(
          "nfakevtx0vsnpuprof", "# of fake vertices vs number sim vertices", vmax, 0., float(vmax), 0., float(2 * vmax)));
  add(h,
      new TProfile("nfakevtx4vsnpuprof",
                   "# of fake vertices  (ndof>4) vs number sim vertices",
                   vmax,
                   0.,
                   float(vmax),
                   0.,
                   float(2 * vmax)));
  add(h,
      new TProfile("nfoundpuvtx0vsnpuprof",
                   "# of found PU vertices vs number sim vertices",
                   vmax,
                   0.,
                   float(vmax),
                   0.,
                   float(2 * vmax)));
  add(h,
      new TProfile("nfoundpuvtx4vsnpuprof",
                   "# of found PU vertices (ndof>4) vs number sim vertices",
                   vmax,
                   0.,
                   float(vmax),
                   0.,
                   float(2 * vmax)));
  add(h,
      new TProfile("nlostpuvtxvsnpuprof",
                   "# of lost PU vertices vs number sim vertices",
                   vmax,
                   0.,
                   float(vmax),
                   0.,
                   float(2 * vmax)));
  // FIXME these histograms don't seem to be filled
  add(h, new TH1F("sumpt2rec", "sumpt2 of reconstructed and matched tracks", 100, 0., 100.));
  add(h, new TH1F("sumpt2", "sumpt2 of simulated tracks", 100, 0., 100.));
  add(h, new TH1F("sumpt2Signal", "sumpt2 of simulated tracks in Signal events", 100, 0., 200.));
  add(h, new TH1F("sumpt2PU", "sumpt2 of simulated tracks in PU events", 100, 0., 200.));
  add(h, new TH1F("sumpt2rec", "sumpt2 of reconstructed and matched tracks", 100, 0., 100.));
  add(h, new TH1F("sumpt2recSignal", "sumpt2 of reconstructed and matched tracks in Signal events", 100, 0., 200.));
  add(h, new TH1F("sumpt2recPU", "sumpt2 of reconstructed and matched tracks in PU events", 100, 0., 200.));
  add(h, new TH1F("nRecTrkInSimVtx", "number of reco tracks matched to sim-vertex", 100, 0., 100));
  add(h, new TH1F("nRecTrkInSimVtxSignal", "number of reco tracks matched to signal sim-vertex", 100, 0., 100));
  add(h, new TH1F("nRecTrkInSimVtxPU", "number of reco tracks matched to PU-vertex", 100, 0., 100));
  add(h, new TH1F("nPrimRecTrkInSimVtx", "number of reco primary tracks matched to sim-vertex", 100, 0., 100));
  add(h,
      new TH1F(
          "nPrimRecTrkInSimVtxSignal", "number of reco primary tracks matched to signal sim-vertex", 100, 0., 100));
  add(h, new TH1F("nPrimRecTrkInSimVtxPU", "number of reco primary tracks matched to PU-vertex", 100, 0., 100));

  add(h, new TH1F("recmatchPurity", "track purity of all vertices", 101, 0., 1.01));
  add(h, new TH1F("recmatchvtxs", "number of sim vertices contributing to a recvtx", 10, 0., 10.));
  add(h, new TH1F("recmatchvtxsFake", "number of sim vertices contributing to a fake recvtx (ndof>4)", 10, 0., 10.));
  add(h, new TH1F("recmatchvtxsReal", "number of sim vertices contributing to a real recvtx (ndof>4)", 10, 0., 10.));
  add(h,
      new TH2F("nsiminrecvsunmatchedFake",
               "number of sim vertices contributing to a fake recvtx (ndof>4)",
               10,
               0.,
               10.,
               10,
               0.,
               10.));
  add(h,
      new TH2F("nsiminrecvsunmatchedReal",
               "number of sim vertices contributing to a real recvtx (ndof>4)",
               10,
               0.,
               10.,
               10,
               0.,
               10.));
  add(h,
      new TH1F("recmatch30vtxs", "number of sim vertices contributing >30% of their tracks to a recvtx", 10, 0., 10.));
  add(h,
      new TH1F("recmatch50vtxs", "number of sim vertices contributing >50% of their tracks to a recvtx", 10, 0., 10.));
  //
 
  add(h, new TH1F("vtxMultiplicity", "number of rec vertices containing tracks from one true vertex", 10, 0., 10.));
  add(h,
      new TH1F(
          "vtxMultiplicitySignal", "number of rec vertices containing tracks from the Signal Vertex", 10, 0., 10.));
  add(h, new TH1F("vtxMultiplicityPU", "number of rec vertices containing tracks from a PU Vertex", 10, 0., 10.));
  add(h,
      new TH1F(
          "vtxMultiplicity50", "number of rec vertices containing >=50% tracks from one true vertex", 10, 0., 10.));
  add(h,
      new TH1F(
          "vtxMultiplicity50Signal", "number of rec vertices containing tracks from the Signal Vertex", 10, 0., 10.));
  add(h, new TH1F("vtxMultiplicity50PU", "number of rec vertices containing tracks from a PU Vertex", 10, 0., 10.));
  add(h,
      new TH1F(
          "vtxMultiplicity30", "number of rec vertices containing >=30% tracks from one true vertex", 10, 0., 10.));
  add(h,
      new TH1F(
          "vtxMultiplicity30Signal", "number of rec vertices containing tracks from the Signal Vertex", 10, 0., 10.));
  add(h, new TH1F("vtxMultiplicity30PU", "number of rec vertices containing tracks from a PU Vertex", 10, 0., 10.));

  add(h,
      new TProfile(
          "vtxFindingEfficiencyVsNtrk", "finding efficiency vs number of associated rec tracks", 100, 0., 100., 0., 1.));
  add(h,
      new TProfile("vtxFindingEfficiencyVsNtrkSignal",
                   "Signal vertex finding efficiency vs number of associated rec tracks",
                   100,
                   0.,
                   100.,
                   0.,
                   1.));
  add(h,
      new TProfile("vtxFindingEfficiencyVsNtrkPU",
                   "PU vertex finding efficiency vs number of associated rec tracks",
                   100,
                   0.,
                   100.,
                   0.,
                   1.));

  add(h, new TH1F("matchVtxFraction", "fraction of sim vertex tracks found in a recvertex", 101, 0, 1.01));
  add(h, new TH1F("matchVtxFractionSignal", "fraction of sim vertex tracks found in a recvertex", 101, 0, 1.01));
  add(h, new TH1F("matchVtxFractionPU", "fraction of sim vertex tracks found in a recvertex", 101, 0, 1.01));
  add(h, new TH1F("matchVtxFractionCum", "fraction of sim vertex tracks found in a recvertex", 101, 0, 1.01));
  add(h, new TH1F("matchVtxFractionCumSignal", "fraction of sim vertexs track found in a recvertex", 101, 0, 1.01));
  add(h, new TH1F("matchVtxFractionCumPU", "fraction of sim vertex tracks found in a recvertex", 101, 0, 1.01));
  add(h, new TH1F("matchVtxEfficiency", "efficiency for finding matching rec vertex (ntsim>0)", 2, -0.5, 1.5));
  add(h, new TH1F("matchVtxEfficiencySignal", "efficiency for finding matching rec vertex (ntsim>0)", 2, -0.5, 1.5));
  add(h, new TH1F("matchVtxEfficiencyPU", "efficiency for finding matching rec vertex (ntsim>0)", 2, -0.5, 1.5));
  add(h, new TH1F("matchVtxEfficiency2", "efficiency for finding matching rec vertex (nt>1)", 2, -0.5, 1.5));
  add(h, new TH1F("matchVtxEfficiency2Signal", "efficiency for finding matching rec vertex (nt>1)", 2, -0.5, 1.5));
  add(h, new TH1F("matchVtxEfficiency2PU", "efficiency for finding matching rec vertex (nt>1)", 2, -0.5, 1.5));
  add(h, new TH1F("matchVtxEfficiency5", "efficiency for finding matching rec vertex (purity>0.5)", 2, -0.5, 1.5));
  add(h,
      new TH1F("matchVtxEfficiency5Signal", "efficiency for finding matching rec vertex (purity>0.5)", 2, -0.5, 1.5));
  add(h, new TH1F("matchVtxEfficiency5PU", "efficiency for finding matching rec vertex (purity>0.5)", 2, -0.5, 1.5));

  add(h, new TH1F("matchVtxEfficiencyZ", "efficiency for finding matching rec vertex within 1 mm", 2, -0.5, 1.5));
  add(h, new TH1F("matchVtxEfficiencyZSignal", "efficiency for finding matching rec vertex within 1 mm", 2, -0.5, 1.5));
  add(h, new TH1F("matchVtxEfficiencyZPU", "efficiency for finding matching rec vertex within 1 mm", 2, -0.5, 1.5));

  add(h,
      new TH1F("matchVtxEfficiencyZ1", "efficiency for finding matching rec vertex within 1 mm (nt>0)", 2, -0.5, 1.5));
  add(h,
      new TH1F(
          "matchVtxEfficiencyZ1Signal", "efficiency for finding matching rec vertex within 1 mm (nt>0)", 2, -0.5, 1.5));
  add(h,
      new TH1F(
          "matchVtxEfficiencyZ1PU", "efficiency for finding matching rec vertex within 1 mm (nt>0)", 2, -0.5, 1.5));

  add(h,
      new TH1F("matchVtxEfficiencyZ2", "efficiency for finding matching rec vertex within 1 mm (nt>1)", 2, -0.5, 1.5));
  add(h,
      new TH1F(
          "matchVtxEfficiencyZ2Signal", "efficiency for finding matching rec vertex within 1 mm (nt>1)", 2, -0.5, 1.5));
  add(h,
      new TH1F(
          "matchVtxEfficiencyZ2PU", "efficiency for finding matching rec vertex within 1 mm (nt>1)", 2, -0.5, 1.5));

  add(h, new TH1F("matchVtxZ", "z distance to matched recvtx", 100, -0.1, 0.1));
  add(h, new TH1F("matchVtxZPU", "z distance to matched recvtx", 100, -0.1, 0.1));
  add(h, new TH1F("matchVtxZSignal", "z distance to matched recvtx", 100, -0.1, 0.1));

  add(h, new TH1F("matchVtxZCum", "z distance to matched recvtx", 1001, 0, 1.01));
  add(h, new TH1F("matchVtxZCumSignal", "z distance to matched recvtx", 1001, 0, 1.01));
  add(h, new TH1F("matchVtxZCumPU", "z distance to matched recvtx", 1001, 0, 1.01));

  add(h, new TH1F("unmatchedVtx", "number of unmatched rec vertices (fakes)", 10, 0., 10.));
  add(h, new TH1F("unmatchedVtx4", "number of unmatched rec vertices ndof>4 (fakes)", 10, 0., 10.));
  add(h, new TH1F("unmatchedVtxW4", "number of unmatched (by weight) rec vertices ndof>4 (fakes)", 10, 0., 10.));
  add(h, new TH1F("unmatchedVtxNtrk", "number of tracks in unmatched vertex", 20, -0.5, 19.5));
  add(h, new TH1F("unmatchedVtxFrac", "fraction of unmatched rec vertices (fakes)", 1000, 0., 1.0));
  add(h, new TH1F("unmatchedVtxZ", "z of unmatched rec  vertices (fakes)", 100, -zrange_wide, zrange_wide));
  add(h,
      new TH1F("unmatchedVtxDeltaZ",
               "Delta z of unmatched rec  vertices (fakes)",
               nzbin_wide_fine,
               -zrange_wide,
               zrange_wide));
  add(h, new TH1F("unmatchedVtxNdof", "ndof of unmatched rec vertices (fakes)", 500, 0., 100.));
  add(h, new TH1F("unmatchedVtxNdof1", "ndof of unmatched rec vertices (fakes, delta z>1cm)", 500, 0., 100.));
  add(h, new TH1F("unmatchedVtxNdof2", "ndof of unmatched rec vertices (fakes, delta z>2cm)", 500, 0., 100.));
  add(h, new TH2F("correctlyassigned", "pt and eta of correctly assigned tracks", 80, -4., 4., 100, 0, 10.));
  add(h, new TH2F("misassigned", "pt and eta of mis assigned tracks", 80, -4., 4., 100, 0, 10.));
  add(h, new TH2F("correctlyassignedlogpt", "log10 pt and eta of correctly assigned tracks", 80, -4., 4., 40, -1., 3.));
  add(h, new TH2F("misassignedlogpt", "log10 pt and eta of mis assigned tracks", 80, -3., 4., 40, -1., 3.));

  add(h, new TH1F("ptcat", "pt of correctly assigned tracks", 100, 0, 10.));
  add(h, new TH1F("etacat", "eta of correctly assigned tracks", 80, -4., 4.));
  add(h, new TH1F("etacatpt2", "eta of correctly assigned tracks pt>2GeV", 80, -4., 4.));
  add(h, new TH1F("phicat", "phi of correctly assigned tracks", 100, -3.14159, 3.14159));
  add(h, new TH1F("dzcat", "dz of correctly assigned tracks", 100, 0., 1.));

  add(h, new TH1F("ptmis", "pt of mis-assigned tracks", 100, 0, 10.));
  add(h, new TH1F("etamis", "eta of mis-assigned tracks", 80, -4., 4.));
  add(h, new TH1F("etamispt2", "eta mis-correctly assigned tracks pt>2GeV", 80, -4., 4.));
  add(h, new TH1F("phimis", "phi of mis-assigned tracks", 100, -3.14159, 3.14159));
  add(h, new TH1F("dzmis", "dz of mis-assigned tracks", 100, 0., 1.));

  add(h, new TH1F("Tc", "Tc computed with Truth matched Reco Tracks", 100, 0., 20.));
  add(h, new TH1F("TcSignal", "Tc of signal vertices computed with Truth matched Reco Tracks", 100, 0., 20.));
  add(h, new TH1F("TcPU", "Tc of PU vertices computed with Truth matched Reco Tracks", 100, 0., 20.));

  add(h, new TH1F("logTc", "log Tc computed with Truth matched Reco Tracks", 100, -2., 8.));
  add(h, new TH1F("logTcSignal", "log Tc of signal vertices computed with Truth matched Reco Tracks", 100, -2., 8.));
  add(h, new TH1F("logTcPU", "log Tc of PU vertices computed with Truth matched Reco Tracks", 100, -2., 8.));

  add(h, new TH1F("xTc", "Tc of merged clusters", 100, 0., 20.));
  add(h, new TH1F("xTcSignal", "Tc of signal vertices merged with PU", 100, 0., 20.));
  add(h, new TH1F("xTcPU", "Tc of merged PU vertices", 100, 0., 20.));

  add(h, new TH1F("logxTc", "log Tc merged vertices", 100, -2., 8.));
  add(h, new TH1F("logxTcSignal", "log Tc of signal vertices merged with PU", 100, -2., 8.));
  add(h, new TH1F("logxTcPU", "log Tc of merged PU vertices ", 100, -2., 8.));

  add(h, new TH1F("logChisq", "Chisq/ntrk computed with Truth matched Reco Tracks", 100, -2., 8.));
  add(h,
      new TH1F(
          "logChisqSignal", "Chisq/ntrk of signal vertices computed with Truth matched Reco Tracks", 100, -2., 8.));
  add(h, new TH1F("logChisqPU", "Chisq/ntrk of PU vertices computed with Truth matched Reco Tracks", 100, -2., 8.));

  add(h, new TH1F("logxChisq", "Chisq/ntrk of merged clusters", 100, -2., 8.));
  add(h, new TH1F("logxChisqSignal", "Chisq/ntrk of signal vertices merged with PU", 100, -2., 8.));
  add(h, new TH1F("logxChisqPU", "Chisq/ntrk of merged PU vertices", 100, -2., 8.));

  add(h, new TH1F("Chisq", "Chisq/ntrk computed with Truth matched Reco Tracks", 100, 0., 20.));
  add(h,
      new TH1F("ChisqSignal", "Chisq/ntrk of signal vertices computed with Truth matched Reco Tracks", 100, 0., 20.));
  add(h, new TH1F("ChisqPU", "Chisq/ntrk of PU vertices computed with Truth matched Reco Tracks", 100, 0., 20.));

  add(h, new TH1F("xChisq", "Chisq/ntrk of merged clusters", 100, 0., 20.));
  add(h, new TH1F("xChisqSignal", "Chisq/ntrk of signal vertices merged with PU", 100, 0., 20.));
  add(h, new TH1F("xChisqPU", "Chisq/ntrk of merged PU vertices", 100, 0., 20.));

  add(h, new TH1F("dzmax", "dzmax computed with Truth matched Reco Tracks", 100, 0., 2.));
  add(h, new TH1F("dzmaxSignal", "dzmax of signal vertices computed with Truth matched Reco Tracks", 100, 0., 2.));
  add(h, new TH1F("dzmaxPU", "dzmax of PU vertices computed with Truth matched Reco Tracks", 100, 0., 2.));

  add(h, new TH1F("xdzmax", "dzmax of merged clusters", 100, 0., 2.));
  add(h, new TH1F("xdzmaxSignal", "dzmax of signal vertices merged with PU", 100, 0., 2.));
  add(h, new TH1F("xdzmaxPU", "dzmax of merged PU vertices", 100, 0., 2.));

  add(h, new TH1F("dztrim", "dzmax computed with Truth matched Reco Tracks", 100, 0., 2.));
  add(h, new TH1F("dztrimSignal", "dzmax of signal vertices computed with Truth matched Reco Tracks", 100, 0., 2.));
  add(h, new TH1F("dztrimPU", "dzmax of PU vertices computed with Truth matched Reco Tracks", 100, 0., 2.));

  add(h, new TH1F("xdztrim", "dzmax of merged clusters", 100, 0., 2.));
  add(h, new TH1F("xdztrimSignal", "dzmax of signal vertices merged with PU", 100, 0., 2.));
  add(h, new TH1F("xdztrimPU", "dzmax of merged PU vertices", 100, 0., 2.));

  add(h, new TH1F("m4m2", "m4m2 computed with Truth matched Reco Tracks", 100, 0., 100.));
  add(h, new TH1F("m4m2Signal", "m4m2 of signal vertices computed with Truth matched Reco Tracks", 100, 0., 100.));
  add(h, new TH1F("m4m2PU", "m4m2 of PU vertices computed with Truth matched Reco Tracks", 100, 0., 100.));

  add(h, new TH1F("xm4m2", "m4m2 of merged clusters", 100, 0., 100.));
  add(h, new TH1F("xm4m2Signal", "m4m2 of signal vertices merged with PU", 100, 0., 100.));
  add(h, new TH1F("xm4m2PU", "m4m2 of merged PU vertices", 100, 0., 100.));

  add(h, new TH1F("trecvtx_selmatched", "reconstructed time matched selected vertices", 100, -1., 1.));
  add(h, new TH1F("tsimvtx_selmatched", "simulated time matched selected vertices", 100, -1., 1.));
  add(h, new TH1F("tresvtx_selmatched", "reconstructed - simulated time matched selected vertices", 100, -1., 1.));
  add(h, new TH1F("tpullvtx_selmatched", "reconstructed - simulated time matched selected vertices", 100, -10., 10.));

  addSP(h, new TH1F("xresvtx_selmatched", "reconstructed x matched selected vertices", 200, -0.0100, 0.0100));
  addSP(h, new TH1F("yresvtx_selmatched", "reconstructed y matched selected vertices", 200, -0.0100, 0.0100));
  addSP(h, new TH1F("zresvtx_selmatched", "reconstructed z matched selected vertices", 200, -0.0400, 0.0400));

  return h;
}


void PrimaryVertexAnalyzer4PU::bookTrackHistograms(const char * directory_name)
// book histograms for truth-matched tracks,
// filled in analyzeTracksTP
{
  rootFile_->cd();
  TDirectory* dir = rootFile_->mkdir(directory_name);
  dir->cd();

  // this prim-only group may not be needed
  add(hTrk, new TH1F("zpullsec", "reconstructed z- generated z / error for non-primary tracks", 200, -10., 10.));

  // these are filled
  add(hTrk, new TH1F("zpulltrk_primselmatched", "reconstructed z- generated z for primary tracks", 200, -10., 10.));
  add(hTrk,
      new TH1F(
          "zpulltrkt_primselmatched", "reconstructed z- generated z for primary tracks with timing", 200, -10., 10.));
  add(hTrk, new TH1F("zpulltprimsel", "reconstructed z- generated z for primary tracks with timing", 200, -10., 10.));
  add(hTrk, new TH1F("zrestrk_primselmatched", "reconstructed z- generated z for primary tracks", 200, -0.2, 0.2));
  add(hTrk, new TH2F("zpulltprimselvseta", "", 48, -2.4, 2.4, 200, -10., 10.));


  //>>>>>>>>>>>>>>>>>
  for(auto bin = 0u; bin < trkdzbin_.size(); bin++){
    add(hTrk, new TH2F(Form("zpulltprimselvseta_%s", trkdzbin_[bin].c_str()), "", 48, -2.4, 2.4, 200, -10., 10.));
    add(hTrk, new TH2F(Form("zpulltprimselbpxlt2vseta_%s", trkdzbin_[bin].c_str()), "", 48, -2.4, 2.4, 200, -10., 10.));
    add(hTrk, new TH2F(Form("zpulltprimselbpxgt2vseta_%s", trkdzbin_[bin].c_str()), "", 48, -2.4, 2.4, 200, -10., 10.));
    add(hTrk, new TH2F(Form("zpulltprimselvslogpt_%s", trkdzbin_[bin].c_str()), "", 40, -1., 3., 200, -10., 10.));
    add(hTrk, new TProfile(Form("ztailtprimselvslogpt_%s",trkdzbin_[bin].c_str()) , "", 40, -1., 3., 0., 2.));
    add(hTrk, new TProfile(Form("ztailtprimselvslogpt_etahi_%s",trkdzbin_[bin].c_str()) , "", 40, -1., 3., 0., 2.));
    add(hTrk, new TProfile(Form("ztailtprimselvslogpt_etalo_%s",trkdzbin_[bin].c_str()) , "", 40, -1., 3., 0., 2.));
    add(hTrk, new TH2F(Form("ztailtprimselvslogpteta_%s", trkdzbin_[bin].c_str()), "", 48, -2.4, 2.4, 40, -1., 3.));
    add(hTrk, new TH2F(Form("tprimselvslogpteta_%s", trkdzbin_[bin].c_str()), "", 48, -2.4, 2.4, 40, -1., 3.));
  }

  add(hTrk, new TH2F("zpulltprimselvslogpt", "", 40, -1., 3., 200, -10., 10.));

  add(hTrk, new TH2F("ztailtprimselvslogpteta", "", 48, -2.4, 2.4, 40, -1., 3.));
  add(hTrk, new TH2F("tprimselvslogpteta", "", 48, -2.4, 2.4, 40, -1., 3.));
  //<<<<<<<<<<<<<<<<<<<<<<<<

  add(hTrk,
      new TH1F("zrestrkt_primselmatched",
               "reconstructed z- generated z / error for primary tracks with timing",
               200,
               -0.2,
               0.2));
  add(hTrk,
      new TH1F("tpulltrk_primselmatched", "reconstructed t- generated t / error for primary tracks", 200, -10., 10.));
  add(hTrk, new TH1F("ttrk_sim_primselmatched", "simulated t for primary tracks", 200, -1., 1.));
  add(hTrk, new TH1F("ttrk_rec_primselmatched", "reconstructed t for primary tracks", 200, -1., 1.));
  add(hTrk, new TH1F("trestrk_primselmatched", "reconstructed t -simulated t for selected matched primary tracks", 200, -1., 1.)); //???

  float tqualmin=-2.0;
  float tqualmax=2.0;
  add(hTrk, new TH1F("ttrk_rec_all_wide", "reconstructed t for primary tracks", 200, -10., 10.));
  add(hTrk, new TH1F("ttrk_rec_all", "reconstructed t for primary tracks", 200, -1., 1.));
  add(hTrk, new TH1F("terrtrk_rec_all", "reconstructed t error for all", 200, 0., 0.2));
  add(hTrk, new TH1F("terrtrk_rec_sel", "reconstructed t error for all selected tracks", 200, 0., 0.2));
  add(hTrk, new TH1F("tqualtrk_rec_all", "time quality for all", 200, tqualmin, tqualmax));
  add(hTrk, new TH1F("tqualtrk_rec_sel", "time quality for all selected tracks", 200, tqualmin, tqualmax));
  add(hTrk, new TH1F("tqualtrk_sigmat01_rec_sel", "time quality for all selected tracks with sigma_t < 0.1", 200, tqualmin, tqualmax));

  add(hTrk, new TH1F("ttrk_rec_sel_wide", "reconstructed t for selected tracks", 200, -10., 10.));
  add(hTrk, new TH1F("ttrk_rec_sel", "reconstructed t for selected tracks", 200, -1., 1.));
  add(hTrk, new TH1F("ttrk_rec_selmatched_wide", "reconstructed t for matched selected tracks", 200, -10., 10.));
  add(hTrk, new TH1F("ttrk_rec_selmatched", "reconstructed t for matched selected tracks", 200, -1., 1.));
  add(hTrk, new TH1F("terrtrk_rec_selmatched", "reconstructed t error for matched selected tracks", 200, 0., 0.1));
  add(hTrk, new TH1F("terrtrk_rec_selmatched_barrel", "reconstructed t error for matched selected tracks (barrel)", 200, 0., 0.1));
  add(hTrk, new TH1F("terrtrk_rec_selmatched_endcap", "reconstructed t error for matched selected tracks (endcap)", 200, 0., 0.1));
  add(hTrk, new TH1F("tqualtrk_rec_selmatched", "time quality for matched selected tracks", 200, tqualmin, tqualmax));
  add(hTrk, new TH1F("tqualtrk_rec_selmatched_barrel", "time quality for matched selected tracks", 200, tqualmin, tqualmax));
  add(hTrk, new TH1F("tqualtrk_rec_selmatched_endcap", "time quality for matched selected tracks", 200, tqualmin, tqualmax));
  add(hTrk, new TH1F("tqualtrk_rec_selmatched_barrel_hipt", "time quality for matched selected tracks", 200, tqualmin, tqualmax));
  add(hTrk, new TH1F("tqualtrk_rec_selmatched_endcap_hipt", "time quality for matched selected tracks", 200, tqualmin, tqualmax));
  add(hTrk, new TH1F("tqualtrk_rec_selmatched_barrel_lopt", "time quality for matched selected tracks", 200, tqualmin, tqualmax));
  add(hTrk, new TH1F("tqualtrk_rec_selmatched_endcap_lopt", "time quality for matched selected tracks", 200, tqualmin, tqualmax));
  add(hTrk, new TH1F("ttrk_sim_selmatched", "simulated t for matched selected tracks", 200, -1., 1.));
  
  add(hTrk, new TH1F("trestrk_selmatched", "reconstructed t - simulated t for matched selected tracks", 200, -1., 1.));
  add(hTrk, new TH1F("trestrk_selmatched_pid", "reconstructed t - simulated t for matched selected tracks", 200, -1., 1.));
  add(hTrk, new TH1F("trestrk_selmatched_pid_ptgt1", "reconstructed t - simulated t for matched selected tracks", 200, -1., 1.));
  add(hTrk, new TH1F("trestrk_selmatched_qgt08", "reconstructed t - simulated t for matched selected tracks, quality > 0.8", 200, -1., 1.));
  add(hTrk, new TH1F("trestrk_selmatched_q0508", "reconstructed t - simulated t for matched selected tracks, quality = 0.5 .. 0.8", 200, -1., 1.));
  add(hTrk, new TH1F("trestrk_selmatched_qlt05", "reconstructed t - simulated t for matched selected tracks, quality < 0.5", 200, -1., 1.));
  add(hTrk, new TH1F("trestrk_selmatched_pid_qlt05", "reconstructed t - simulated t for matched selected tracks, quality < 0.5", 200, -1., 1.));
  add(hTrk, new TH1F("trestrk_selmatched_pid_qgt08", "reconstructed t - simulated t for matched selected tracks, quality > 0.8", 200, -1., 1.));
  add(hTrk, new TH1F("trestrk_selmatched_pid_q0508", "reconstructed t - simulated t for matched selected tracks, quality = 0.5 .. 0.8", 200, -1., 1.));
  
  add(hTrk, new TH1F("tpulltrk_selmatched", "normalized trec-tsim  for matched selected tracks", 200, -10., 10.));
  add(hTrk, new TH1F("tpulltrk_selmatched_pid", "normalized trec-tsim  for matched selected tracks", 200, -10., 10.));
  add(hTrk, new TH1F("tpulltrk_selmatched_pid_ptgt1", "normalized trec-tsim  for matched selected tracks", 200, -10., 10.));
  add(hTrk, new TH1F("tpulltrk_selmatched_qgt08", "normalized trec-tsim  for matched selected tracks, quality > 0.8 ", 200, -10., 10.));
  add(hTrk, new TH1F("tpulltrk_selmatched_q0508", "normalized trec-tsim  for matched selected tracks, quality = 0.5 .. 0.8", 200, -10., 10.));
  add(hTrk, new TH1F("tpulltrk_selmatched_qlt05", "normalized trec-tsim  for matched selected tracks, quality < 0.5 ", 200, -10., 10.));
  add(hTrk, new TH1F("tpulltrk_selmatched_pid_qgt08", "normalized trec-tsim  for matched selected tracks, quality > 0.8 ", 200, -10., 10.));
  add(hTrk, new TH1F("tpulltrk_selmatched_pid_q0508", "normalized trec-tsim  for matched selected tracks, quality = 0.5 .. 0.8", 200, -10., 10.));
  add(hTrk, new TH1F("tpulltrk_selmatched_pid_qlt05", "normalized trec-tsim  for matched selected tracks, quality < 0.5 ", 200, -10., 10.));

  add(hTrk, new TH1F("trestrk_sigmatlt01_selmatched_pid", "reconstructed t - simulated t for matched selected tracks", 200, -1., 1.));
  add(hTrk, new TH1F("tpulltrk_sigmatlt01_selmatched_pid", "normalized trec-tsim  for matched selected tracks", 200, -10., 10.));

  add(hTrk, new TH1F("trestrk_sigmatlt01_selmatched_qgt08", "reconstructed t - simulated t for matched selected tracks, quality > 0.8", 200, -1., 1.));
  add(hTrk, new TH1F("trestrk_sigmatlt01_selmatched_q0508", "reconstructed t - simulated t for matched selected tracks, quality = 0.5 .. 0.8", 200, -1., 1.));
  add(hTrk, new TH1F("trestrk_sigmatlt01_selmatched_qlt05", "reconstructed t - simulated t for matched selected tracks, quality < 0.5", 200, -1., 1.));
  add(hTrk, new TH1F("tpulltrk_sigmatlt01_selmatched_qgt08", "normalized trec-tsim  for matched selected tracks, quality > 0.8 ", 200, -10., 10.));
  add(hTrk, new TH1F("tpulltrk_sigmatlt01_selmatched_q0508", "normalized trec-tsim  for matched selected tracks, quality = 0.5 .. 0.8", 200, -10., 10.));
  add(hTrk, new TH1F("tpulltrk_sigmatlt01_selmatched_qlt05", "normalized trec-tsim  for matched selected tracks, quality < 0.5 ", 200, -10., 10.));

  add(hTrk, new TH1F("trestrk_selmatched_barrel", "reconstructed t - simulated t for matched selected tracks (barrel)", 200, -1., 1.));
  add(hTrk, new TH1F("trestrk_selmatched_endcap", "reconstructed t - simulated t for matched selected tracks (endcap)", 200, -1., 1.));
  add(hTrk, new TH1F("trestrk_selmatched_barrel_hipt", "reconstructed t - simulated t for matched selected tracks (barrel, pt>1)", 200, -1., 1.));
  add(hTrk, new TH1F("trestrk_selmatched_endcap_hipt", "reconstructed t - simulated t for matched selected tracks (endcap, pt>1)", 200, -1., 1.));
  add(hTrk, new TH1F("trestrk_selmatched_barrel_lopt", "reconstructed t - simulated t for matched selected tracks (barrel, pt<1)", 200, -1., 1.));
  add(hTrk, new TH1F("trestrk_selmatched_endcap_lopt", "reconstructed t - simulated t for matched selected tracks (endcap, pt<1)", 200, -1., 1.));
  add(hTrk, new TH1F("trestrk_selmatched_fwd", "reconstructed t - simulated t for matched selected tracks (fwd)", 200, -1., 1.));
  add(hTrk, new TH1F("trestrk_selmatched_bwd", "reconstructed t - simulated t for matched selected tracks (bwd)", 200, -1., 1.));
  add(hTrk, new TH2F("trestrkvszrestrk_selmatched_barrel",
		     "reconstructed t - simulated t for matched selected tracks (barrel)",
		     150, -2.0, 2.0, 100, -0.1, 0.1));  // expect 33 ps/cm = 0.033 ns/cm correlation
  add(hTrk, new TH2F("trestrkvszrestrk_selmatched_endcap",
		     "reconstructed t - simulated t for matched selected tracks (endcap)",
		     150, -2.0, 2.0, 100, -0.1, 0.1));
  
  add(hTrk, new TH1F("tpulltrk_selmatched_endcap", "normalized trec-tsim  for matched selected tracks (endcap)", 200, -10., 10.));
  add(hTrk, new TH1F("tpulltrk_selmatched_barrel", "normalized trec-tsim  for matched selected tracks (barrel)", 200, -10., 10.));
  add(hTrk, new TH1F("tpulltrk_selmatched_endcap_hipt", "normalized trec-tsim  for matched selected tracks (endcap, pt>1.0)", 200, -10., 10.));
  add(hTrk, new TH1F("tpulltrk_selmatched_barrel_hipt", "normalized trec-tsim  for matched selected tracks (barrel, pt>1.0)", 200, -10., 10.));
  add(hTrk, new TH1F("tpulltrk_selmatched_endcap_lopt", "normalized trec-tsim  for matched selected tracks (endcap, pt<1.0)", 200, -10., 10.));
  add(hTrk, new TH1F("tpulltrk_selmatched_barrel_lopt", "normalized trec-tsim  for matched selected tracks (barrel, pt<1.0)", 200, -10., 10.));
  add(hTrk, new TH1F("trestrk_selpion", "reconstructed t - simulated t for matched selected pions", 200, -1., 1.));
  add(hTrk, new TH1F("trestrk_selkaon", "reconstructed t - simulated t for matched selected kaons", 200, -1., 1.));
  add(hTrk, new TH1F("trestrk_selproton", "reconstructed t - simulated t for matched selected protons", 200, -1., 1.));
  add(hTrk, new TH1F("tpulltrk_selpion", "normalized trec-tsim  for matched selected pions", 200, -10., 10.));
  add(hTrk, new TH1F("tpulltrk_selkaon", "normalized trec-tsim  for matched selected kaons", 200, -10., 10.));
  add(hTrk, new TH1F("tpulltrk_selproton", "normalized trec-tsim  for matched selected protons", 200, -10., 10.));
  
  add(hTrk, new TH1F("trestrkh_selpion", "reconstructed t - simulated t for matched selected pions", 200, -1., 1.));
  add(hTrk, new TH1F("trestrkh_selkaon", "reconstructed t - simulated t for matched selected kaons", 200, -1., 1.));
  add(hTrk, new TH1F("trestrkh_selproton", "reconstructed t - simulated t for matched selected protons", 200, -1., 1.));
  add(hTrk, new TH1F("trestrkh_selhyp", "reconstructed t - simulated t for matched particles", 200, -1., 1.));
  add(hTrk, new TH1F("tpulltrkh_selpion", "normalized trec-tsim  for matched selected pions", 200, -10., 10.));
  add(hTrk, new TH1F("tpulltrkh_selkaon", "normalized trec-tsim  for matched selected kaons", 200, -10., 10.));
  add(hTrk, new TH1F("tpulltrkh_selproton", "normalized trec-tsim  for matched selected protons", 200, -10., 10.));
  add(hTrk, new TH1F("tpulltrkh_selhyp", "normalized trec-tsim  for matched particles", 200, -10., 10.));

  add(hTrk, new TH1F("tqualtrk_sigmat01_selmatched", "time quality for matched tracks with sigma_t < 0.1", 200, tqualmin, tqualmax));
  add(hTrk, new TH1F("tpulltrk_sigmat01_selmatched", "normalized trec-tsim  for matched selected tracks", 200, -10., 10.));
  add(hTrk, new TH1F("tpulltrk_sigmat01_selmatched_pid", "normalized trec-tsim  for matched selected tracks", 200, -10., 10.));
  add(hTrk, new TH1F("tpulltrk_sigmat01_selpion", "normalized trec-tsim  for matched selected pions", 200, -10., 10.));
  add(hTrk, new TH1F("tpulltrk_sigmat01_selkaon", "normalized trec-tsim  for matched selected kaons", 200, -10., 10.));
  add(hTrk, new TH1F("tpulltrk_sigmat01_selproton", "normalized trec-tsim  for matched selected protons", 200, -10., 10.));

  add(hTrk, new TH1F("ttrk_rec_selunmatched_wide", "reconstructed t for unmatched selected", 200, -10., 10.));
  add(hTrk, new TH1F("ttrk_rec_selunmatched", "reconstructed t for unmatched selected tracks", 200, -1., 1.));

  add(hTrk, new TH1F("trecunmatchedwide", "reconstructed t for primary tracks", 200, -10., 10.));
  add(hTrk, new TH1F("trecunmatchedselwide", "reconstructed t for primary tracks", 200, -10., 10.));
  add(hTrk, new TH1F("tresprim", "reconstructed t- generated t for primary tracks", 200, -1., 1.));
  add(hTrk, new TH1F("tresprimsel", "reconstructed t- generated t for primary tracks", 200, -1., 1.));
  add(hTrk, new TH1F("d0pullprim", "d0/error for primary tracks", 200, -10., 10.));
  add(hTrk, new TH1F("d0pullsec", "d0/error for non-primary tracks", 200, -10., 10.));
  add(hTrk, new TH2F("zpullvsd0pullprim", "z pull vs d0-pull for primary tracks", 100, 0., 10., 100, 0., 10.));
  add(hTrk, new TH2F("zpullvsd0pullsec", "z pull vs d0-pull for non-primary tracks", 100, 0., 10., 100, 0., 10.));
  add(hTrk,
      new TH1F(
          "zpulltrkt_primselmatched_central", "reconstructed z- generated z / error for primary tracks |eta|<1.5", 200, -10., 10.));
  add(hTrk,
      new TH1F(
          "zpulltrkt_primselmatched_inter", "reconstructed z- generated z / error for primary tracks 1.5<|eta|<2", 200, -10., 10.));
  add(hTrk,
      new TH1F("zpulltrkt_primselmatched_fwd", "reconstructed z- generated z / error for primary tracks |eta|>2", 200, -10., 10.));
  add(hTrk, new TH1F("tpiprim", "prior track weight primary tracks", 100, 0., 1.));
  add(hTrk, new TH1F("tpisec", "prior track weight non-primary tracks", 100, 0., 1.));
  add(hTrk, new TH1F("r200denom", "r denom(z=2000)", 1, 0., 1.));
  add(hTrk, new TH1F("r200", "r(z=2000)", 1000, 0., 50.));
  add(hTrk, new TH1F("r200b", "r(z=2000)", 1000, 0., 50.));
  add(hTrk, new TH1F("r200c", "r(z=2000)", 1000, 0., 50.));
  add(hTrk, new TProfile("dzvseta", "#sigma(z) vs #eta", 20, -5., 5., 0., 10.));
  add(hTrk, new TProfile("dzvsetapt1", "#sigma(z) vs #eta   0.0 GeV < pt < 0.5 GeV", 20, -5., 5., 0., 10.));
  add(hTrk, new TProfile("dzvsetapt2", "#sigma(z) vs #eta   0.5 GeV < pt < 1 GeV", 20, -5., 5., 0., 10.));
  add(hTrk, new TProfile("dzvsetapt3", "#sigma(z) vs #eta   1 GeV < pt < 2 GeV", 20, -5., 5., 0., 10.));
  add(hTrk, new TProfile("dzvsetapt4", "#sigma(z) vs #eta   2 GeV < pt < 4 GeV", 20, -5., 5., 0., 10.));
  add(hTrk, new TProfile("dzvsetapt5", "#sigma(z) vs #eta   4 GeV < pt < 10 GeV", 20, -5., 5., 0., 10.));
  add(hTrk, new TProfile("dzvsetapt6", "#sigma(z) vs #eta  10 GeV < pt < 20 GeV", 20, -5., 5., 0., 10.));
  add(hTrk, new TProfile("dzvsetapt7", "#sigma(z) vs #eta   pt > 20 GeV", 20, -5., 5., 0., 10.));

  // CH histos, simhit maps, FIXME optional / remove
  add(hTrk, new TH2F("hits", "hits", 600, -300., 300., 500, 0., 50.));
  for (unsigned int i = 0; i < hitCollections_.size(); i++) {
    add(hTrk,
        new TH2F(("hits-" + hitCollections_[i]).c_str(),
                 ("hits-" + hitCollections_[i]).c_str(),
                 600,
                 -300.,
                 300.,
                 500,
                 0.,
                 50.));
  }
}

void PrimaryVertexAnalyzer4PU::bookSimPVHistograms(const char * directory_name)
{
  TDirectory* dir = rootFile_->mkdir(directory_name);
  dir->cd();
  add(hsimPV, new TH1F("xsim", "simulated x", 100, -0.1, 0.1));
  add(hsimPV, new TH1F("ysim", "simulated y", 100, -0.1, 0.1));
  add(hsimPV, new TH1F("zsim", "simulated z", 120, -30., 30.));
  add(hsimPV, new TH1F("xsimb", "simulated x", 100, -0.01, 0.01));  // 0.01cm = 100 um
  add(hsimPV, new TH1F("ysimb", "simulated y", 100, -0.01, 0.01));
  add(hsimPV, new TH1F("zsimb", "simulated z", 120, -30., 30.));
  if (f4D_) {
    add(hsimPV, new TH1F("tsim", "simulated t", 400, -2., 2.));
  }

  add(hsimPV, new TH1F("xbeam", "beamspot x", 100, -1., 1.));
  add(hsimPV, new TH1F("ybeam", "beamspot y", 100, -1., 1.));
  add(hsimPV, new TH1F("zbeam", "beamspot z", 100, -5., 5));
  add(hsimPV, new TH1F("wxbeam", "beamspot sigma x", 100, 0., 0.02));
  add(hsimPV, new TH1F("wybeam", "beamspot sigma y", 100, 0., 0.02));
  add(hsimPV, new TH1F("sigmaZbeam", "beamspot sigma z", 100, 0., 10.));
  const int nsimmax = 250;
  add(hsimPV, new TH1F("nsimvtx", "# of simulated vertices", nsimmax, 0., float(nsimmax)));
  add(hsimPV, new TH1F("nbsimtksinvtx", "simulated tracks in vertex", 200, -0.5, 199.5));
}


void PrimaryVertexAnalyzer4PU::bookEventHistograms(const char * directory_name)
{
  TDirectory* dir = rootFile_->mkdir(directory_name);
  dir->cd();
  add(hEvt, new TH1F("lbx", "instantaneous BX lumi (ub s bx)^-1", 100, 0., lumiHistoRange_));
  add(hEvt, new TH1F("bunchCrossing", "bunchCrossing", 4000, 0., 4000.));
}


void PrimaryVertexAnalyzer4PU::get_luminosity_infos(const edm::Event& iEvent) {

  instBXLumi_ = -1.;
  avginstBXLumi_ = -1.;
  lsglb_ = 0;   // for history plots, a global lumisection number, cumulated across all runs

  const LumiInfo& lumi = iEvent.get(lumiInfoToken_);
  /*
  std::cout << "Luminosity for " << iEvent.run() << " LS " << iEvent.luminosityBlock() << " is "
	    << lumi.getTotalInstLumi() << std::endl;
  */
  avginstBXLumi_ = lumi.getTotalInstLumi();
  instBXLumi_ = lumi.getInstLumiBX(iEvent.bunchCrossing()); // not yet

  lumiPU_ = sigma_pp_ * instBXLumi_;
  avglumiPU_ = sigma_pp_ * avginstBXLumi_;
}

void PrimaryVertexAnalyzer4PU::get_particle_data_table(const edm::EventSetup& iSetup) {
  try {
    iSetup.getData(pdt_);
  } catch (const Exception&) {
    std::cout << "Some problem occurred with the particle data table. This may not work !" << std::endl;
  }
}

bool PrimaryVertexAnalyzer4PU::get_beamspot_data(const edm::Event& iEvent) {
  // load beam spot, fill some histograms
  // return false if no usable beamspot data was found

  if (iEvent.getByToken(recoBeamSpotToken_, recoBeamSpotHandle_)) {
    vertexBeamSpot_ = *recoBeamSpotHandle_;
    wxy2_ = pow(vertexBeamSpot_.BeamWidthX(), 2) + pow(vertexBeamSpot_.BeamWidthY(), 2);
    wx_ = vertexBeamSpot_.BeamWidthX();
    wy_ = vertexBeamSpot_.BeamWidthY();
    sigmaZ_ = vertexBeamSpot_.sigmaZ();
    sigmaT_ = sigmaZ_ / 2.998e1;  // c in cm/ns
    if (sigmaZ_ < 0.1) {
      reportEvent("crazy small sigma Z");
      return false;
    }

    if (sigmaZ_ > 9) {
      reportEvent("huge sigma Z");
      return false;
    }

    if (std::abs(vertexBeamSpot_.z0()) > 2.) {
      reportEvent("suspicious beamspot position  %6.3f", vertexBeamSpot_.z0());
      return false;
    }

    Fill(hsimPV, "xbeam", vertexBeamSpot_.x0());
    Fill(hsimPV, "wxbeam", vertexBeamSpot_.BeamWidthX());
    Fill(hsimPV, "ybeam", vertexBeamSpot_.y0());
    Fill(hsimPV, "wybeam", vertexBeamSpot_.BeamWidthY());
    Fill(hsimPV, "zbeam", vertexBeamSpot_.z0());
    Fill(hsimPV, "sigmaZbeam", vertexBeamSpot_.sigmaZ());
    if (verbose_ && (luminosityBlock_ != currentLS_)) {
      cout << "BEAM " << run_ << " : " << std::setw(10) << luminosityBlock_ << " " << std::setw(8) << std::fixed
           << std::setprecision(4) << vertexBeamSpot_.x0() << ", " << vertexBeamSpot_.y0() << ", "
           << vertexBeamSpot_.z0() << "+/-" << vertexBeamSpot_.z0Error() << ",  wx=" << vertexBeamSpot_.BeamWidthX()
           << ",  wy=" << vertexBeamSpot_.BeamWidthY() << " , sigmaZ=" << vertexBeamSpot_.sigmaZ() << "+/-"
           << vertexBeamSpot_.sigmaZ0Error() << endl;
      currentLS_ = luminosityBlock_;
    }

  } else {
    sigmaZ_ = 0;
    reportEvent("no beamspot found, aborting");
    return false;
  }

  // abort without useful beamspot data
  if (sigmaZ_ == 0) {
    return false;
  }

  return true;
}

bool PrimaryVertexAnalyzer4PU::get_reco_and_transient_tracks(const edm::EventSetup& iSetup,
                                                             const edm::Event& iEvent,
                                                             Tracks& tracks) {
  // requires beamspot
  if (!iEvent.getByToken(edmView_recoTrack_Token_, tracks.trackCollectionH)) {
    if (verbose_) {
      cout << " no reco tracks found, bailing out" << endl;
    }
    return false;
  }

  const View<reco::Track>* recTrks = tracks.trackCollectionH.product();

  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", theB_);
  fBfield_ = ((*theB_).field()->inTesla(GlobalPoint(0., 0., 0.))).z();

  edm::Handle<edm::ValueMap<float>> trackTimeQualityH;
  edm::Handle<edm::ValueMap<float>> trackTimesH;
  edm::Handle<edm::ValueMap<float>> trackTimeResosH;

  bool trktime = f4D_
    && iEvent.getByToken(trkTimesToken_, trackTimesH)
    && iEvent.getByToken(trkTimeResosToken_, trackTimeResosH);
  bool have_timequality = iEvent.getByToken(trkTimeQualityToken_, trackTimeQualityH);

  edm::Handle<edm::ValueMap<float>> MTD_pathlength_H, MTD_time_H, MTD_timeerror_H, MTD_momentum_H;
  bool have_MTD_pathlength = iEvent.getByToken(MTD_pathlength_Token_, MTD_pathlength_H);
  bool have_MTD_time = iEvent.getByToken(MTD_time_Token_, MTD_time_H);
  bool have_MTD_timeerror = iEvent.getByToken(MTD_timeerror_Token_, MTD_timeerror_H);
  bool have_MTD_momentum = iEvent.getByToken(MTD_momentum_Token_, MTD_momentum_H);
  
  if (f4D_){
    if (trktime){
      report_counted("PrimaryVertexAnalyzer4PU::get_reco_and_transient_tracks timing found",1);
      //cout << "pathlength " << have_MTD_pathlength << "  mtdtime " << have_MTD_time << "   mtdtimerror" << have_MTD_timeerror << "   mtdmomentum " << have_MTD_momentum << endl;
    }else{
      report_counted("PrimaryVertexAnalyzer4PU::get_reco_and_transient_tracks timing requested but not found",1);
    }
  }

  if (trktime) {
        t_tks_ = (*theB_).build(
                                tracks.trackCollectionH, vertexBeamSpot_, *(trackTimesH.product()), *(trackTimeResosH.product()));
  }else { 
    t_tks_ = (*theB_).build(tracks.trackCollectionH, vertexBeamSpot_);
  }


  /* fill private track container and make it globally available */
  for (View<reco::Track>::size_type i = 0; i < recTrks->size(); ++i) {
    reco::TransientTrack* tt = &(t_tks_[i]);
    unsigned int key1 = t_tks_[i].trackBaseRef().key();
    RefToBase<reco::Track> trb(tracks.trackCollectionH, i);
    unsigned int key2 = trb.key();
    // paranoia is in bloom
    if (key1 != key2) {
      cout << "get_reco_and_transient_tracks : key confusion" << endl;
    }

    auto tk = RecoTrack(i, &recTrks->at(i), tt, key1, trktime);
    if(trktime && have_timequality){
      tk.timeQuality = (*trackTimeQualityH)[tt->trackBaseRef()];
    }
    
    if(have_MTD_time){
      tk.MTD_time = (*MTD_time_H)[tt->trackBaseRef()];
      //cout << "MTD time = " << tk.MTD_time << ", has_timing=" << tk.has_timing << endl;
    }
    if(have_MTD_timeerror){
      tk.MTD_timeerror = (*MTD_timeerror_H)[tt->trackBaseRef()];
      //cout << "MTD time error= " << tk.MTD_timeerror << "    track time error = " << tk.dt <<  endl;
    }else{
      tk.MTD_timeerror = 1.e10;
    }
    if(have_MTD_momentum){
      tk.MTD_momentum = (*MTD_momentum_H)[tt->trackBaseRef()];
      //cout << "MTD momentum = " << tk.MTD_momentum << "  vs " << tk.pt / std::abs(tan(tk.theta)) << endl;
    }else{
      tk.MTD_momentum = tk.pt / std::abs(tan(tk.theta));
    }
    if(have_MTD_pathlength){
      tk.MTD_pathlength = (*MTD_pathlength_H)[tt->trackBaseRef()];
      //cout << "MTD pathlength = " << tk.MTD_pathlength << endl;
    }else{
      tk.MTD_pathlength = -1.;
    }

    if (have_MTD_time && have_MTD_momentum && have_MTD_pathlength && (tk.MTD_pathlength > 0) ){
      tk.th[0] = tk.get_t_pid( 0.139570 );
      tk.th[1] = tk.get_t_pid( 0.493677 );
      tk.th[2] = tk.get_t_pid( 0.938272 );
      //cout << "t(track)=" <<  tk.t << " t(pi) = " << tk.th[0] << " t(K) = " << tk.th[1] << "  t(p) = " << tk.th[2] << endl;
    }else{
      for(unsigned int j = 0; j < 3; j++){ tk.th[j] = 0; }
    }
    
    tk.selected = theTrackFilter(*tt);

    tracks.push_back(tk);
  }

  if (recTrks->size() < minNumberOfRecTrks_) {
    emptyeventcounter_++;
    if (emptyeventcounter_ < 100) {
      cout << "Event without tracks skipped   run= " << run_ << " ls = " << luminosityBlock_ << " BX=" << bunchCrossing_
           << "  instBXLumi_=" << 1.e3 * instBXLumi_ << " expected PU = " << lumiPU_ << "    selected tracks"
           << recTrks->size() << std::endl;
    }
    if (emptyeventcounter_ == 100) {
      cout << "This is the last message of this kind" << endl;
    }
    return false;
  }

  return true;
}

bool PrimaryVertexAnalyzer4PU::get_MC_truth(const edm::Event& iEvent,
                                            Tracks& tracks,
                                            bool bPuInfo,
                                            PileupSummaryInfo& puInfo,
                                            std::vector<SimEvent>& simEvt,
                                            std::vector<simPrimaryVertex>& simpv,
                                            std::vector<SimPart>& tsim)
{
  std::vector<simPrimaryVertex> simpvHepMC;  //  another list of primary MC vertices
  std::string mcproduct = "generator";       // starting with 3_1_0 pre something
  tracking_truth_available_ = false;

  // genParticles for AOD et al:  https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookGenParticleCandidate
  Handle<reco::GenParticleCollection> genParticlesH;
  bool bgenParticles = iEvent.getByToken(genParticleCollection_Token_, genParticlesH);
  if (bgenParticles) report_counted("found genparticles", 1);

  Handle<HepMCProduct> evtMC;
  Handle<SimVertexContainer> simVtxs;
  Handle<SimTrackContainer> simTrks;
  bool bSimVtxs = false;
  bool bSimTrks = false;
  try {
    bSimVtxs = iEvent.getByToken(edmSimVertexContainerToken_, simVtxs);
    bSimTrks = iEvent.getByToken(edmSimTrackContainerToken_, simTrks);
  } catch (...) {
    cout << "no sim tracks found" << endl;
    MC_ = false;
  }
  
  if (bSimTrks) report_counted("found simTracks", 1);
  if (bSimVtxs) report_counted("found simVtxs", 1);

  edm::Handle<TrackingParticleCollection> TPCollectionH;
  edm::Handle<TrackingVertexCollection> TVCollectionH;
  bool gotTP = iEvent.getByToken(trackingParticleCollectionToken_, TPCollectionH);
  bool gotTV = iEvent.getByToken(trackingVertexCollectionToken_, TVCollectionH);

  if(gotTP)
    {
      report_counted("found tracking particles", 1);
    }
  else
    {
      report_counted("no tracking particles found", 1);
    }      
  if(gotTV) report_counted("found tracking vertices", 1);

  MC_ |= gotTP;
  trkidx2tp_.clear();
  trkidx2simevt_.clear();  // deprecated
  trkkey2simevt_.clear();  // deprecated

  
  // always fill simpvs from HepEvt, when possible
  if (iEvent.getByToken(edmHepMCProductToken_, evtMC)) {
    MC_ = true;
    simpvHepMC = getSimPVs(evtMC);  // the signal vertex in HepMC becomes the first entry in simpv
    tsim = PrimaryVertexAnalyzer4PU::getSimTrkParameters(simTrks, simVtxs, simUnit_);

    if (bPuInfo) {
      report_counted("filling simpvs from hepMCProduct + puInfo", 1);
      for (unsigned int i = 0; i < puInfo.getPU_zpositions().size(); i++) {
        double t = 0;
        if (puInfo.has_times()) {
          t = puInfo.getPU_times()[i];  // appears to be in ns already
        }
        if (veryverbose_) {
          cout << setw(4) << i << "  z=" << setw(8) << fixed << setprecision(4) << puInfo.getPU_zpositions()[i]
               << "  t=" << setw(7) << fixed << setprecision(3) << t << "  pthat= " << scientific
               << puInfo.getPU_pT_hats()[i] << endl;
        }

        PrimaryVertexAnalyzer4PU::simPrimaryVertex v(
            vertexBeamSpot_.x0(), vertexBeamSpot_.y0(), puInfo.getPU_zpositions()[i], t);

        v.type = 2;  // partial
        // nlo cut-off is 0.1 GeV  (includes the high pt tracks)
        // hi cut-off is 0.5 GeV
        v.nGenTrk = 0;  //puInfo.getPU_ntrks_lowpT()[i] ;
        v.sumpT = 0;    //puInfo.getPU_sumpT_lowpT()[i];
        v.pt_hat = puInfo.getPU_pT_hats()[i];
        v.is_visible = (v.pt_hat > 1.);
        simpvHepMC.push_back(v);
      }
    }

    // Now compute the closest distance in z between all simulated vertex  (from slimmed)
    // first initialize
    auto prev_z = simpvHepMC.back().z;
    for (simPrimaryVertex& vsim : simpvHepMC) {
      vsim.closest_vertex_distance_z = std::abs(vsim.z - prev_z);
      prev_z = vsim.z;
    }
    // then calculate
    for (std::vector<simPrimaryVertex>::iterator vsim = simpvHepMC.begin(); vsim != simpvHepMC.end(); vsim++) {
      std::vector<simPrimaryVertex>::iterator vsim2 = vsim;
      vsim2++;
      for (; vsim2 != simpvHepMC.end(); vsim2++) {
        double distance = std::abs(vsim->z - vsim2->z);
        // need both to be complete
        vsim->closest_vertex_distance_z = std::min(vsim->closest_vertex_distance_z, distance);
        vsim2->closest_vertex_distance_z = std::min(vsim2->closest_vertex_distance_z, distance);
      }
    }
  }

  //if (gotTP && gotTV) {
  if (gotTP)
    {
      // fill simEvt with Trackingparticles and associator output
      // FIXME can we do this without TV?
      
      edm::Handle<reco::SimToRecoCollection> simToRecoH;
      iEvent.getByToken(simToRecoAssociationToken_, simToRecoH);
      
      edm::Handle<reco::RecoToSimCollection> recoToSimH;
      iEvent.getByToken(recoToSimAssociationToken_, recoToSimH);
      
      //    s2r_ = simToRecoH.product();
      r2s_ = recoToSimH.product();
      tracking_truth_available_ = true;
      
      /*
      // Vertex associator a la *Slimmed
      edm::Handle<reco::VertexToTrackingVertexAssociator> vertexAssociatorH;
      iEvent.getByToken(vertexAssociatorToken_, vertexAssociatorH);
      const reco::VertexToTrackingVertexAssociator& vertexAssociator = *(vertexAssociatorH.product());
      */
      
      simEvt = getSimEvents(TPCollectionH, tracks);
      analyzeTracksTP(tracks, simEvt);
    }
  else if( !bPuInfo  || (puInfo.getPU_zpositions().size() == 0) )
    {
      // no tracking particles, but MC info and no PU
      simEvt = getSimEvents_notp(simTrks, simVtxs, tracks);
    }

  
  if (gotTV) {
    // fill simpv from Tracking Vertices, if we have them

    MC_ = true;
    if (verbose_) {
      cout << "Found Tracking Vertices " << endl;
    }
    simpv = getSimPVs(TVCollectionH, simpvHepMC);
    if (verbose_) {
      cout << "simpv.size = " << simpv.size() << endl;
    }

  }
  // no tracking vertices from here on
  else if (iEvent.getByToken(edmHepMCProductToken_, evtMC)) {
    simpv = simpvHepMC;

  } else if (iEvent.getByToken(edmHepMCProductToken_, evtMC)) { // FIXME, this is the same if as above
    // fill simpv from HepEvt

    MC_ = true;
    simpv = getSimPVs(evtMC);  // the signal vertex in HepMC becomes the first entry in simpv
    tsim = PrimaryVertexAnalyzer4PU::getSimTrkParameters(simTrks, simVtxs, simUnit_);

    if (bPuInfo) {
      report_counted("filling simpvs from hepMCProduct + puInfo (2)", 1);
      if (verbose_) {
        cout << "filling simpvs from hepMCProduct + puInfo (2)" << endl;
      }
      for (unsigned int i = 0; i < puInfo.getPU_zpositions().size(); i++) {
        double t = 0;
        if (puInfo.has_times()) {
          t = puInfo.getPU_times()[i];
        }

        PrimaryVertexAnalyzer4PU::simPrimaryVertex v(
            vertexBeamSpot_.x0(), vertexBeamSpot_.y0(), puInfo.getPU_zpositions()[i], t);
        v.type = 2;  // partial
        // nlo cut-off is 0.1 GeV  (includes the high pt tracks)
        // hi cut-off is 0.5 GeV
        v.nGenTrk = 0;  //puInfo.getPU_ntrks_lowpT()[i];
        v.sumpT = 0;    //puInfo.getPU_sumpT_lowpT()[i];
        v.pt_hat = puInfo.getPU_pT_hats()[i];
        v.is_visible = (v.pt_hat > 1.);
        simpv.push_back(v);
      }

      cout << "peng " << simpv.size() << endl;
    } else {
      if (verbose_) {
        cout << "filling simpvs from hepMCProduct, no puInfo" << endl;
      }
    }

  } else if (bSimTrks && bSimVtxs) {
    report_counted("filling simpvs from simVtxs and simTrks",1);
    simpv = getSimPVs(simVtxs, simTrks);
    tsim = getSimTrkParameters(simTrks, simVtxs, simUnit_);
    MC_ = true;

  } else if (bgenParticles) {
    simpv = getSimPVs(genParticlesH);
    tsim = getSimTrkParameters(genParticlesH);
    MC_ = true;
    if (verbose_) {
      cout << "Signal vertex  z=" << simpv[0].z << "  n=" << simpv[0].nGenTrk << endl;
    }
    if (bPuInfo) {
      simPU_ = puInfo.getPU_NumInteractions();  // well, +1 if we are not looking at neutrino gun events
      if (verbose_) {
        cout << "PileupSummaryInfo  nPU=" << puInfo.getPU_NumInteractions() << endl;
      }

      for (unsigned int i = 0; i < puInfo.getPU_zpositions().size(); i++) {
        double t = 0;
        if (puInfo.has_times()) {
          t = puInfo.getPU_times()[i];
        }

        PrimaryVertexAnalyzer4PU::simPrimaryVertex v(
            vertexBeamSpot_.x0(), vertexBeamSpot_.y0(), puInfo.getPU_zpositions()[i], t);
        v.type = 2;  // partial
        // nlo cut-off is 0.1 GeV  (includes the high pt tracks)
        // hi cut-off is 0.5 GeV
        v.nGenTrk = 0;  //puInfo.getPU_ntrks_lowpT()[i];
        v.sumpT = 0;    // puInfo.getPU_sumpT_lowpT()[i];
        //v.eventId=puInfo.getPU_EventID()[i];
        v.pt_hat = puInfo.getPU_pT_hats()[i];
        v.is_visible = (v.pt_hat > 1.);
        simpv.push_back(v);
      }
    }

  } else {
    MC_ = false;
    return false;
  }
  return true;
}



void PrimaryVertexAnalyzer4PU::fill_simvtx_histos(std::vector<simPrimaryVertex>& simpv)
// fill some basic sim vertex quantities not using any reco stuff
// booked in bookSimPVHistograms, get their own subdirectory "simpvs"
// called by analyze, information in simpv is filled by get_MC_truth
{
  hsimPV["nsimvtx"]->Fill(simpv.size());
  for (std::vector<simPrimaryVertex>::iterator vsim = simpv.begin(); vsim != simpv.end(); vsim++) {
    if (vsim->type == 1) {
      hsimPV["nbsimtksinvtx"]->Fill(vsim->nGenTrk);
      hsimPV["xsim"]->Fill(vsim->x * simUnit_);
      hsimPV["ysim"]->Fill(vsim->y * simUnit_);
      hsimPV["zsim"]->Fill(vsim->z * simUnit_);
      hsimPV["xsimb"]->Fill(vsim->x * simUnit_ - vertexBeamSpot_.x0());
      hsimPV["ysimb"]->Fill(vsim->y * simUnit_ - vertexBeamSpot_.y0());
      hsimPV["zsimb"]->Fill(vsim->z * simUnit_ - vertexBeamSpot_.z0());

      if (f4D_) {
        hsimPV["tsim"]->Fill(vsim->t);
      }
    }  //type==1
  }
}




void PrimaryVertexAnalyzer4PU::beginJob() {
  nihmatch_ = 1;
  matchsummaries_ = 0;

  MC_ = false;
  dxb_ = 0.0000;
  dyb_ = 0.0000;
  dzb_ = 0.0;
  wxb_ = 0.1;
  wyb_ = 0.1;
  wzb_ = 5.;

  lumiHistoRange_ = 1.;  //0.001; in nb/bx  : must multiply instbxlumi_ by 1e3
  lumiPUHistoRange_ = 250.;
  sigma_pp_ = 69.1e3;  // [ub]  pp-cross section for PU

  reports_.clear();

  /* histogram booking */
  rootFile_->cd();
  
  // vertex collections
  for (std::vector<std::string>::const_iterator vCollection = vertexCollectionLabels_.begin();
       vCollection != vertexCollectionLabels_.end();
       vCollection++) {
    cout << "PrimaryVertexAnalyzer4PU: booking histograms for collection " << *vCollection << endl;
    std::string sdir = *vCollection;
    if (sdir == "offlineSlimmedPrimaryVertices") {
      sdir = "offlinePrimaryVertices";
      cout << "histogram directory changed from offlineSlimmedPrimaryVertices to offlinePrimaryVertices " << endl;
    }
    TDirectory* dir = rootFile_->mkdir(sdir.c_str());
    dir->cd();
    std::string s = *vCollection;
    
    histograms_[s] = bookVertexHistograms(dir);
  }

  // others
  rootFile_->cd();
  bookTrackHistograms("tracks");
  bookSimPVHistograms("simpvs");
  bookEventHistograms("event");

  rootFile_->cd();
}




void PrimaryVertexAnalyzer4PU::endJob() {
  std::cout << flush;
  std::cout << "this is void PrimaryVertexAnalyzer4PU::endJob() " << std::endl;


  // do some normalization / cumultion etc before saving
  for (std::vector<std::string>::const_iterator vCollection = vertexCollectionLabels_.begin();
       vCollection != vertexCollectionLabels_.end();
       vCollection++)
    {
      std::map<std::string, TH1*> h = histograms_[*vCollection];
      
      if (false) {
	// normalize some histograms to the eventcounter
	if (nEventNsel_ > 0) {
	  h["nnsel"]->Scale(1. / nEventNsel_);
	  h["n11sel"]->Scale(1. / nEventNsel_);
	  h["n11rej"]->Scale(1. / nEventNsel_);
        std::cout << " successfully normalized histograms, nEventNsel_ = " << nEventNsel_ << std::endl;
	} else {
	  std::cout << " histograms not normalized, nEventNsel_ = " << nEventNsel_ << std::endl;
	}
      } else {
	std::cout << " histogram normalization not requested " << nEventNsel_ << std::endl;
    }
      
      if ((h.count("sumntkvsz") > 0) && (h.count("sumwoversumntkvsz") > 0) && (h.count("sumwvsz") > 0)) {
	for (int i = 1; i < 101; i++) {
	  if (h["sumntkvsz"]->GetBinContent(i) > 0) {
	    h["sumwoversumntkvsz"]->SetBinContent(i, h["sumwvsz"]->GetBinContent(i) / h["sumntkvsz"]->GetBinContent(i));
        }
	}
      }

      double sum = 0;
      if ((h.count("matchVtxFractionSignal") > 0) && (h.count("matchVtxFractionCumSignal") > 0)) {
	for (int i = 101; i > 0; i--) {
	  sum += h["matchVtxFractionSignal"]->GetBinContent(i) / h["matchVtxFractionSignal"]->Integral();
	  h["matchVtxFractionCumSignal"]->SetBinContent(i, sum);
	}
    }
      
      sum = 0;
      if ((h.count("abszdistancenearest") > 0) && (h.count("abszdistancenearestcum") > 0)) {
	for (int i = 1; i < 1001; i++) {
	  sum += h["abszdistancenearest"]->GetBinContent(i);
	  h["abszdistancenearestcum"]->SetBinContent(i, sum / float(h["abszdistancenearest"]->GetEntries()));
	}
    }
      
      Cumulate(h, "matchVtxZCum");
      Cumulate(h, "matchVtxZCumSignal");
      Cumulate(h, "matchVtxZCumPU");
      
      if ((h.count("vtxndf") > 0) && (h.count("vtxndfc") > 0)) {
      double p;
      unsigned int nbin = h["vtxndf"]->GetNbinsX();
      for (unsigned int i = 1; i <= nbin; i++) {
        if (h["vtxndf"]->GetEntries() > 0) {
          p = h["vtxndf"]->Integral(i, nbin + 1) / h["vtxndf"]->GetEntries();
          h["vtxndfc"]->SetBinContent(i, p * h["vtxndf"]->GetBinContent(i));
        }
      }
      }

      for (int i = 1; i <= h["mergerate"]->GetNbinsX(); i++) {
	int denom = h["mergerate_denominator"]->GetBinContent(i);
	if (denom > 0) {
	  h["mergerate"]->SetBinContent(i, float(h["mergerate_numerator"]->GetBinContent(i)) / denom);
	}
      }
      for (int i = 1; i <= h["mergerate_lin"]->GetNbinsX(); i++) {
	int denom = h["mergerate_denominator_lin"]->GetBinContent(i);
	if (denom > 0) {
	  h["mergerate_lin"]->SetBinContent(i, float(h["mergerate_numerator_lin"]->GetBinContent(i)) / denom);
	}
      }

      for (int i = 1; i <= h["mergerate_dqm"]->GetNbinsX(); i++) {
	int denom = h["mergerate_denominator_dqm"]->GetBinContent(i);
	if (denom > 0) {
	  h["mergerate_dqm"]->SetBinContent(i, float(h["mergerate_numerator_dqm"]->GetBinContent(i)) / denom);
	}
      }
    }
  
  /* now write the rootfile */
  
  rootFile_->cd();

  std::cout << "Info=" << info_->String() << std::endl;
  std::cout << flush;

  if (info_->GetString().Length() > 0) {
    info_->Write("Info");
  }
  if (build_->GetString().Length() > 0) {
    build_->Write("Build");
  }

  /*
  for (std::map<std::string, TH1*>::const_iterator hist = hEvt.begin(); hist != hEvt.end(); hist++) {
    hist->second->Write();
  }

  for (std::map<std::string, TH1*>::const_iterator hist = hsimPV.begin(); hist != hsimPV.end(); hist++) {
    hist->second->Write();
  }
  */
  
  cout << "empty histograms:" << endl;
  int nempty = 0;
  for (std::map<std::string, TH1*>::const_iterator hist = hTrk.begin(); hist != hTrk.end(); hist++) {
    if (hist->second->GetEntries() == 0) {
      cout << hist->first << " ";
      nempty++;
    }
  }
  cout << " found " << nempty << endl;

  for (std::map<std::string, TH1*>::const_iterator hist = hTrk.begin(); hist != hTrk.end(); hist++) {
    //hist->second->Write();
    int nbin = hist->second->GetNbinsX();
    double I = hist->second->GetBinContent(0) + hist->second->GetBinContent(nbin + 1) + hist->second->Integral();
    if (hist->second->GetEntries() > I) {
      cout << "Warning !  possible bin content overflow in '" << hist->first << "' entries=" << hist->second->GetEntries()
           << "   integral=" << I << endl;
    }
  }
  rootFile_->Write();

  /* print out collected reports at the end of the logfile */

  std::cout << "*******************************************************" << std::endl;
  std::cout << "reports " << reports_.size() << std::endl;
  for (auto r : reports_) {
    std::cout << r << endl;
  }
  std::cout << "*******************************************************" << std::endl;

  /* print out counted messages */
  report_counted("summary", 0);

  /* poor mans profiling */
  for (auto it : timers_ ){
    std::cout << setw(50) << it.first   << std::setw(15) << std::fixed << std::setprecision(3) << it.second * 1e-6<< " s" <<std::endl;
  }

    
  std::cout << flush;
  std::cerr << "PrimaryVertexAnalyzer4PU::endJob: done" << std::endl;
  std::cout << flush;
}

// helper functions

bool PrimaryVertexAnalyzer4PU::getPuInfo(const Event& iEvent, PileupSummaryInfo& puInfo) {
  /*
    loads the puInfo for MC events and  fills pseudo "lumi" info
   */

  try {
    //get the pileup information
    Handle<vector<PileupSummaryInfo>> puInfoH;
    if (iEvent.getByToken(vecPileupSummaryInfoToken_, puInfoH)) {
      // ignore out-of-time pile-up
      for (unsigned int puInfo_ite = 0; puInfo_ite < (*puInfoH).size(); ++puInfo_ite) {
        if ((*puInfoH)[puInfo_ite].getBunchCrossing() == 0) {
          puInfo = (*puInfoH)[puInfo_ite];
          break;
        }
      }

      lumiPU_ = 1. + puInfo.getTrueNumInteractions();  // this is the expected number of interactions
      avglumiPU_ = lumiPU_;
      instBXLumi_ = lumiPU_ / sigma_pp_;
      avginstBXLumi_ = lumiPU_ / sigma_pp_;
      simPU_ = puInfo.getPU_NumInteractions();  // this is the actual number of interactions

      if (verbose_) {
        cout << "  puInfo size = " << puInfo.getPU_NumInteractions() << endl;
        cout << "  puInfo true = " << puInfo.getTrueNumInteractions() << endl;
        if (veryverbose_) {
          for (unsigned int i = 0; i < puInfo.getPU_zpositions().size(); i++) {
            cout << "pile-up " << setw(3) << i << ")"
                 << " z= " << setw(8) << setprecision(4) << fixed << puInfo.getPU_zpositions()[i];
            if (puInfo.has_times()) {
              cout << " t = " << setw(7) << setprecision(3) << fixed << puInfo.getPU_times()[i];
            }

            cout << "  pt_hat= " << setw(6) << setprecision(1) << fixed
                 << puInfo.getPU_pT_hats()[i]
                 //<<" ntrks_lowPt=" << setw(6) << setprecision(1) << fixed << puInfo.getPU_ntrks_lowpT()[i]
                 //<< " sumpT_lowPT " << setw(6) << setprecision(1) << fixed << puInfo.getPU_sumpT_lowpT()[i]
                 << endl;
          }
        }
      }

      return true;

    } else {
      return false;
    }

  } catch (edm::Exception const& ex) {
    if (!ex.alreadyPrinted()) {
      std::cout << ex.what() << " Maybe data?" << std::endl;
    }
  }
  return false;
}

std::vector<PrimaryVertexAnalyzer4PU::SimPart> PrimaryVertexAnalyzer4PU::getSimTrkParameters(
    edm::Handle<edm::SimTrackContainer>& simTrks, edm::Handle<edm::SimVertexContainer>& simVtcs, double simUnit) {
  std::vector<SimPart> tsim;
  if (simVtcs->begin() == simVtcs->end()) {
    if (verbose_) {
      cout << "  PrimaryVertexAnalyzer4PU::getSimTrkParameters  no simvtcs" << endl;
    }
    return tsim;
  }
  if (verbose_) {
    cout << "  PrimaryVertexAnalyzer4PU::getSimTrkParameters simVtcs n=" << simVtcs->size() << endl;
    cout << "  PrimaryVertexAnalyzer4PU::getSimTrkParameters 1st position" << setw(8) << setprecision(4)
         << simVtcs->begin()->position() << endl;
  }
  double t0 = simVtcs->begin()->position().e();

  for (edm::SimTrackContainer::const_iterator t = simTrks->begin(); t != simTrks->end(); ++t) {
    if (t->noVertex()) {
      std::cout << "simtrk has no vertex" << std::endl;
    } else {
      // get the vertex position
      //HepLorentzVector v=(*simVtcs)[t->vertIndex()].position();
      math::XYZTLorentzVectorD v((*simVtcs)[t->vertIndex()].position().x(),
                                 (*simVtcs)[t->vertIndex()].position().y(),
                                 (*simVtcs)[t->vertIndex()].position().z(),
                                 (*simVtcs)[t->vertIndex()].position().e());
      int pdgCode = t->type();

      if (pdgCode == -99) {
        // such entries cause crashes, no idea what they are
        std::cout << "funny particle skipped  , code=" << pdgCode << std::endl;
      } else {
        double Q = 0;  //double Q=HepPDT::theTable().getParticleData(pdgCode)->charge();
        if ((pdgCode == 11) || (pdgCode == 13) || (pdgCode == 15) || (pdgCode == -211) || (pdgCode == -2212) ||
            (pdgCode == -321) || (pdgCode == -3222)) {
          Q = -1;
        } else if ((pdgCode == -11) || (pdgCode == -13) || (pdgCode == -15) || (pdgCode == 211) || (pdgCode == 2212) ||
                   (pdgCode == 321) || (pdgCode == 3222)) {
          Q = 1;
        } else {
          //std::cout << pdgCode << " " <<std::endl;
        }
        math::XYZTLorentzVectorD p(t->momentum().x(), t->momentum().y(), t->momentum().z(), t->momentum().e());
        if ((Q != 0) && (p.pt() > 0.1) && (fabs(t->momentum().eta()) < 5.0) && fabs(v.z() * simUnit < 20) &&
            (sqrt(v.x() * v.x() + v.y() * v.y()) < 10.)) {
          double x0 = v.x() * simUnit;
          double y0 = v.y() * simUnit;
          double z0 = v.z() * simUnit;
          double kappa = -Q * 0.002998 * fBfield_ / p.pt();
          double D0 = x0 * sin(p.phi()) - y0 * cos(p.phi()) - 0.5 * kappa * (x0 * x0 + y0 * y0);
          double q = sqrt(1. - 2. * kappa * D0);
          double s0 = (x0 * cos(p.phi()) + y0 * sin(p.phi())) / q;
          double s1;
          if (fabs(kappa * s0) > 0.001) {
            s1 = asin(kappa * s0) / kappa;
          } else {
            double ks02 = (kappa * s0) * (kappa * s0);
            s1 = s0 * (1. + ks02 / 6. + 3. / 40. * ks02 * ks02 + 5. / 112. * pow(ks02, 3));
          }
          SimPart sp;  //ParameterVector par;
          sp.par[reco::TrackBase::i_qoverp] = Q / p.P();
          sp.par[reco::TrackBase::i_lambda] = M_PI / 2. - p.theta();
          sp.par[reco::TrackBase::i_phi] = p.phi() - asin(kappa * s0);
          sp.par[reco::TrackBase::i_dxy] = -2. * D0 / (1. + q);
          sp.par[reco::TrackBase::i_dsz] = z0 * sin(p.theta()) - s1 * cos(p.theta());

          sp.pdg = pdgCode;
          if (v.t() - t0 < 1e-15) {
            sp.type = 0;  // primary
          } else {
            sp.type = 1;  //secondary
          }

          // now get zpca  (get perigee wrt beam)
          //double x1 = x0; double y1 = y0;
          double x1 = x0 - vertexBeamSpot_.x(z0);
          double y1 = y0 - vertexBeamSpot_.y(z0);

          D0 = x1 * sin(p.phi()) - y1 * cos(p.phi()) - 0.5 * kappa * (x1 * x1 + y1 * y1);
          q = sqrt(1. - 2. * kappa * D0);
          s0 = (x1 * cos(p.phi()) + y1 * sin(p.phi())) / q;
          if (fabs(kappa * s0) > 0.001) {
            s1 = asin(kappa * s0) / kappa;
          } else {
            double ks02 = (kappa * s0) * (kappa * s0);
            s1 = s0 * (1. + ks02 / 6. + 3. / 40. * ks02 * ks02 + 5. / 112. * pow(ks02, 3));
          }
          sp.ddcap = -2. * D0 / (1. + q);
          sp.zdcap = z0 - s1 / tan(p.theta());
          sp.zvtx = z0;
          sp.xvtx = x0;
          sp.yvtx = y0;

          sp.phi = p.phi();
          sp.eta = p.eta();

          tsim.push_back(sp);
        }
      }
    }  // has vertex
  }    //for loop
  return tsim;
}

std::vector<PrimaryVertexAnalyzer4PU::SimPart> PrimaryVertexAnalyzer4PU::getSimTrkParameters(
    const Handle<reco::GenParticleCollection> genParticles) {
  std::vector<SimPart> tsim;
  double xp = 0, yp = 0, zp = -99;
  for (size_t i = 0; i < genParticles->size(); ++i) {
    const GenParticle& gp = (*genParticles)[i];
    int pdgCode = gp.pdgId();
    int st = gp.status();

    if ((st == 1) && (xp == 0) && (yp == 0) && (zp == -99)) {
      xp = gp.vx();
      yp = gp.vy();
      zp = gp.vz();
    }
    if (pdgCode == -99) {
      // such entries cause crashes, no idea what they are
      std::cout << "funny particle skipped  , code=" << pdgCode << std::endl;
    } else {
      double Q = gp.charge();
      if ((st == 1) && (Q != 0) && (gp.pt() > 0.1) && (fabs(gp.eta()) < 5.0) && fabs(gp.vz() < 20) &&
          (sqrt(gp.vx() * gp.vx() + gp.vy() * gp.vy()) < 10.)) {
        double x0 = gp.vx();
        double y0 = gp.vy();
        double z0 = gp.vz();
        double kappa = -Q * 0.002998 * fBfield_ / gp.pt();
        double D0 = x0 * sin(gp.phi()) - y0 * cos(gp.phi()) - 0.5 * kappa * (x0 * x0 + y0 * y0);
        double q = sqrt(1. - 2. * kappa * D0);
        double s0 = (x0 * cos(gp.phi()) + y0 * sin(gp.phi())) / q;
        double s1;
        if (fabs(kappa * s0) > 0.001) {
          s1 = asin(kappa * s0) / kappa;
        } else {
          double ks02 = (kappa * s0) * (kappa * s0);
          s1 = s0 * (1. + ks02 / 6. + 3. / 40. * ks02 * ks02 + 5. / 112. * pow(ks02, 3));
        }
        SimPart sp;  //ParameterVector par;
        sp.par[reco::TrackBase::i_qoverp] = Q / gp.p();
        sp.par[reco::TrackBase::i_lambda] = M_PI / 2. - gp.theta();
        sp.par[reco::TrackBase::i_phi] = gp.phi() - asin(kappa * s0);
        sp.par[reco::TrackBase::i_dxy] = -2. * D0 / (1. + q);
        sp.par[reco::TrackBase::i_dsz] = z0 * sin(gp.theta()) - s1 * cos(gp.theta());

        sp.pdg = pdgCode;
        double t = sqrt(pow(gp.vx() - xp, 2) + pow(gp.vy() - yp, 2) + pow(gp.vz() - zp, 2));
        if (t < 1e-6) {
          sp.type = 0;  // primary
        } else {
          sp.type = 1;  //secondary
        }

        // now get zpca  (get perigee wrt beam)
        double x1 = x0;
        double y1 = y0;
        D0 = x1 * sin(gp.phi()) - y1 * cos(gp.phi()) - 0.5 * kappa * (x1 * x1 + y1 * y1);
        q = sqrt(1. - 2. * kappa * D0);
        s0 = (x1 * cos(gp.phi()) + y1 * sin(gp.phi())) / q;
        if (fabs(kappa * s0) > 0.001) {
          s1 = asin(kappa * s0) / kappa;
        } else {
          double ks02 = (kappa * s0) * (kappa * s0);
          s1 = s0 * (1. + ks02 / 6. + 3. / 40. * ks02 * ks02 + 5. / 112. * pow(ks02, 3));
        }
        sp.ddcap = -2. * D0 / (1. + q);
        sp.zdcap = z0 - s1 / tan(gp.theta());
        sp.zvtx = z0;
        sp.xvtx = x0;
        sp.yvtx = y0;

        sp.phi = gp.phi();
        sp.eta = gp.eta();

        tsim.push_back(sp);
      }
    }
  }  //for loop
  return tsim;
}

/*************************************************************************************************************/
int* PrimaryVertexAnalyzer4PU::supf(std::vector<SimPart>& simtrks, const View<Track>& trks) {
  // track rec to sim matching for hepMC simtrks
  bool debug_supf = false;

  unsigned int nrec = trks.size();
  int* rectosim =
      new int[nrec];  // pointer to associated simtrk, will be returned, caller must delete when no longer used
  for (unsigned int i0 = 0; i0 < nrec; i0++)
    rectosim[i0] = -1;

  unsigned int nsim = simtrks.size();
  if (nsim == 0)
    return rectosim;
  for (unsigned int j = 0; j < nsim; j++)
    simtrks[j].rec = -1;

  // get tracks ordered by volume in z-eta-phi space // other choices det(V), # hits
  std::vector<std::pair<double, unsigned int>> vol;
  for (unsigned int i = 0; i < nrec; i++) {
    auto t = &trks.at(i);
    vol.push_back(std::make_pair(t->etaError() * t->phiError() * t->dszError(), i));
  }
  std::stable_sort(vol.begin(), vol.end(), std::less<std::pair<double, unsigned int>>());

  unsigned int nmatch = 0;
  double chi2max = 100.;  // to be increased for subsequent passes
  for (unsigned int pass = 0; pass < 2; pass++) {
    if (pass > 0)
      chi2max *= 10;
    for (unsigned int k = 0; k < nrec; k++) {
      unsigned int i = vol[k].second;
      auto t = &trks.at(i);
      ParameterVector par = t->parameters();
      reco::TrackBase::CovarianceMatrix S = t->covariance();
      // don't believe in 99.99% correlations
      for (unsigned int a = 0; a < 4; a++) {
        for (unsigned int b = a + 1; b < 5; b++) {
          S(a, b) *= 0.5;
        }
      }

      if (!S.Invert()) {
        cout << "covariance matrix inversion failed ";
        break;
      }

      // try the remaining sim tracks
      int jmatch = -1;
      int jnextbest = -1;
      double cmatch = 1e10;
      double cnextbest = 1e10;

      for (unsigned int j = 0; j < nsim; j++) {
        if (simtrks[j].rec >= 0)
          continue;
        SimPart s = simtrks[j];

        double d[5];  // the parameter difference
        for (int k = 0; k < 5; k++) {
          d[k] = par[k] - s.par[k];
        }

        // phi wrap-around
        if (d[2] > M_PI) {
          d[2] -= M_2_PI;
        } else if (d[2] < -M_PI) {
          d[2] += M_2_PI;
        }

        // matching chi**2
        double c = 0;
        for (int k = 0; k < 5; k++) {
          for (int l = 0; l < 5; l++) {
            c += d[k] * d[l] * S(k, l);
          }
        }

        if (c <= 0) {
          reportEvent("supf: non-positive-definite covariance matrix encountered ", false);
        }

        if (debug_supf) {
          double c0 = pow((par[0] - s.par[0]) / t->qoverpError(), 2) * 0.1  //!
                      + pow((par[1] - s.par[1]) / t->lambdaError(), 2) + pow((par[2] - s.par[2]) / t->phiError(), 2) +
                      pow((par[3] - s.par[3]) / t->dxyError(), 2) * 0.1;  //!
          //	+pow((par[4]-s.par[4])/t->dszError(),2)*0.1;
          //pij[i][j]=exp(-0.5*c0);
          if ((c0 < 100000000) || (c < 100000000)) {
            cout << "pass " << pass << endl;
            cout << setw(3) << i << " rec " << scientific << setw(8) << par << "    " << fixed << setprecision(3)
                 << t->phi() << " " << t->eta() << endl;
            cout << setw(3) << j << " sim " << scientific << setw(8) << s.par << "    " << fixed << setprecision(3)
                 << simtrks[j].phi << " " << simtrks[j].eta << " ---> C=" << scientific << c << endl;
            cout << "       " << setw(8) << d[0] << "," << d[1] << "," << d[2] << "," << d[3] << "," << d[4]
                 << " match=" << match(simtrks[j].par, trks.at(i).parameters()) << endl;
            cout << "       " << setw(8) << (par[0] - s.par[0]) / t->qoverpError() << ","
                 << (par[1] - s.par[1]) / t->lambdaError() << "," << (par[2] - s.par[2]) / t->phiError() << ","
                 << (par[3] - s.par[3]) / t->dxyError() << "," << (par[4] - s.par[4]) / t->dszError() << " c0=" << c0
                 << endl
                 << endl;
          }
        }

        if ((jmatch < 0) || (c < cmatch)) {
          jnextbest = jmatch;
          cnextbest = cmatch;
          jmatch = j;
          cmatch = c;
        } else if ((jnextbest < 0) || (c < cnextbest)) {
          jnextbest = j;
          cnextbest = c;
        }
      }

      if (debug_supf) {
        cout << setw(3) << i << " rec --> " << jmatch;
        if (jmatch >= 0) {
          cout << " cmatch = " << cmatch;
        }
        cout << endl << endl << endl;
      }

      if ((jmatch >= 0) && (cmatch < chi2max)) {
        nmatch++;
        rectosim[i] = jmatch;
        simtrks[jmatch].rec = i;
      }

    }  // rec trk
  }    // pass

  if (verbose_) {
    std::cout << "simtracks (pt>1) without a matching rec track: " << std::endl;
    int nunmatched = 0;
    for (unsigned int j = 0; j < nsim; j++) {
      if (simtrks[j].rec < 0) {
        nunmatched++;
        double pt = 1. / simtrks[j].par[0] / tan(simtrks[j].par[1]);
        if ((fabs(pt)) > 1.) {
          std::cout << setw(3) << j << setw(8) << simtrks[j].pdg << setw(8) << setprecision(4) << "  ("
                    << simtrks[j].xvtx << "," << simtrks[j].yvtx << "," << simtrks[j].zvtx << ")"
                    << " pt= " << pt << " phi=" << simtrks[j].par[2]
                    << " eta= " << -log(tan(0.5 * (M_PI / 2 - simtrks[j].par[1]))) << std::endl;
        }
      }
    }
    cout << " total unmatched " << nunmatched << " out of " << nsim << " sim tracks" << endl;
    cout << " total matched   " << nmatch << " out of " << nrec << " rec tracks" << endl;
  }
  return rectosim;  // caller must delete it !!! delete [] rectosim
}  //supf

/*************************************************************************************************************/

/*************************************************************************************************************/
int* PrimaryVertexAnalyzer4PU::supfv(std::vector<SimPart>& simtrks, const vector<Track>& trks) {
  // track rec to sim matching for hepMC simtrks
  bool debug_supf = false;

  unsigned int nrec = trks.size();
  int* rectosim =
      new int[nrec];  // pointer to associated simtrk, will be returned, caller must delete when no longer used
  for (unsigned int i0 = 0; i0 < nrec; i0++)
    rectosim[i0] = -1;

  unsigned int nsim = simtrks.size();
  if (nsim == 0)
    return rectosim;
  for (unsigned int j = 0; j < nsim; j++)
    simtrks[j].rec = -1;

  // get tracks ordered by volume in z-eta-phi space // other choices det(V), # hits
  std::vector<std::pair<double, unsigned int>> vol;
  for (unsigned int i = 0; i < nrec; i++) {
    auto t = &trks.at(i);
    vol.push_back(std::make_pair(t->etaError() * t->phiError() * t->dszError(), i));
  }
  std::stable_sort(vol.begin(), vol.end(), std::less<std::pair<double, unsigned int>>());

  unsigned int nmatch = 0;
  double chi2max = 100.;  // to be increased for subsequent passes
  for (unsigned int pass = 0; pass < 2; pass++) {
    if (pass > 0)
      chi2max *= 10;
    for (unsigned int k = 0; k < nrec; k++) {
      unsigned int i = vol[k].second;
      auto t = &trks.at(i);
      ParameterVector par = t->parameters();
      reco::TrackBase::CovarianceMatrix S = t->covariance();
      // don't believe in 99.99% correlations
      for (unsigned int a = 0; a < 4; a++) {
        for (unsigned int b = a + 1; b < 5; b++) {
          S(a, b) *= 0.5;
        }
      }

      if (!S.Invert()) {
        cout << "covariance matrix inversion failed ";
        break;
      }

      // try the remaining sim tracks
      int jmatch = -1;
      int jnextbest = -1;
      double cmatch = 1e10;
      double cnextbest = 1e10;

      for (unsigned int j = 0; j < nsim; j++) {
        if (simtrks[j].rec >= 0)
          continue;
        SimPart s = simtrks[j];

        double d[5];  // the parameter difference
        for (int k = 0; k < 5; k++) {
          d[k] = par[k] - s.par[k];
        }

        // phi wrap-around
        if (d[2] > M_PI) {
          d[2] -= M_2_PI;
        } else if (d[2] < -M_PI) {
          d[2] += M_2_PI;
        }

        // matching chi**2
        double c = 0;
        for (int k = 0; k < 5; k++) {
          for (int l = 0; l < 5; l++) {
            c += d[k] * d[l] * S(k, l);
          }
        }

        if (c <= 0) {
          reportEvent("supf: non-positive-definite covariance matrix encountered ", false);
        }

        if (debug_supf) {
          double c0 = pow((par[0] - s.par[0]) / t->qoverpError(), 2) * 0.1  //!
                      + pow((par[1] - s.par[1]) / t->lambdaError(), 2) + pow((par[2] - s.par[2]) / t->phiError(), 2) +
                      pow((par[3] - s.par[3]) / t->dxyError(), 2) * 0.1;  //!
          //	+pow((par[4]-s.par[4])/t->dszError(),2)*0.1;
          //pij[i][j]=exp(-0.5*c0);
          if ((c0 < 100000000) || (c < 100000000)) {
            cout << "pass " << pass << endl;
            cout << setw(3) << i << " rec " << scientific << setw(8) << par << "    " << fixed << setprecision(3)
                 << t->phi() << " " << t->eta() << endl;
            cout << setw(3) << j << " sim " << scientific << setw(8) << s.par << "    " << fixed << setprecision(3)
                 << simtrks[j].phi << " " << simtrks[j].eta << " ---> C=" << scientific << c << endl;
            cout << "       " << setw(8) << d[0] << "," << d[1] << "," << d[2] << "," << d[3] << "," << d[4]
                 << " match=" << match(simtrks[j].par, trks.at(i).parameters()) << endl;
            cout << "       " << setw(8) << (par[0] - s.par[0]) / t->qoverpError() << ","
                 << (par[1] - s.par[1]) / t->lambdaError() << "," << (par[2] - s.par[2]) / t->phiError() << ","
                 << (par[3] - s.par[3]) / t->dxyError() << "," << (par[4] - s.par[4]) / t->dszError() << " c0=" << c0
                 << endl
                 << endl;
          }
        }

        if ((jmatch < 0) || (c < cmatch)) {
          jnextbest = jmatch;
          cnextbest = cmatch;
          jmatch = j;
          cmatch = c;
        } else if ((jnextbest < 0) || (c < cnextbest)) {
          jnextbest = j;
          cnextbest = c;
        }
      }

      if (debug_supf) {
        cout << setw(3) << i << " rec --> " << jmatch;
        if (jmatch >= 0) {
          cout << " cmatch = " << cmatch;
        }
        cout << endl << endl << endl;
      }

      if ((jmatch >= 0) && (cmatch < chi2max)) {
        nmatch++;
        rectosim[i] = jmatch;
        simtrks[jmatch].rec = i;
      }

    }  // rec trk
  }    // pass

  if (verbose_) {
    std::cout << "simtracks (pt>1) without a matching rec track: " << std::endl;
    int nunmatched = 0;
    for (unsigned int j = 0; j < nsim; j++) {
      if (simtrks[j].rec < 0) {
        nunmatched++;
        double pt = 1. / simtrks[j].par[0] / tan(simtrks[j].par[1]);
        if ((fabs(pt)) > 1.) {
          std::cout << setw(3) << j << setw(8) << simtrks[j].pdg << setw(8) << setprecision(4) << "  ("
                    << simtrks[j].xvtx << "," << simtrks[j].yvtx << "," << simtrks[j].zvtx << ")"
                    << " pt= " << pt << " phi=" << simtrks[j].par[2]
                    << " eta= " << -log(tan(0.5 * (M_PI / 2 - simtrks[j].par[1]))) << std::endl;
        }
      }
    }
    cout << " total unmatched " << nunmatched << " out of " << nsim << " sim tracks" << endl;
    cout << " total matched   " << nmatch << " out of " << nrec << " rec tracks" << endl;
  }
  return rectosim;  // caller must delete it !!! delete [] rectosim
}  //supf

/*************************************************************************************************************/

bool PrimaryVertexAnalyzer4PU::match(const ParameterVector& a, const ParameterVector& b) {
  double dqoverp = a(0) - b(0);
  double dlambda = a(1) - b(1);
  double dphi = a(2) - b(2);
  double dsz = a(4) - b(4);
  if (dphi > M_PI) {
    dphi -= M_2_PI;
  } else if (dphi < -M_PI) {
    dphi += M_2_PI;
  }
  //  return ( (fabs(dqoverp)<0.2) && (fabs(dlambda)<0.02) && (fabs(dphi)<0.04) && (fabs(dsz)<0.1) );
  return ((fabs(dqoverp) < 0.2) && (fabs(dlambda) < 0.02) && (fabs(dphi) < 0.04) && (fabs(dsz) < 1.0));
}

bool PrimaryVertexAnalyzer4PU::matchVertex(const simPrimaryVertex& vsim, const reco::Vertex& vrec) {
  return (fabs(vsim.z * simUnit_ - vrec.z()) < zmatch_);
}

bool PrimaryVertexAnalyzer4PU::isResonance(const HepMC::GenParticle* p) {
  double ctau = (pdt_->particle(std::abs(p->pdg_id())))->lifetime();
  //std::cout << "isResonance   " << p->pdg_id() << " " << ctau << std::endl;
  return (ctau > 0) && (ctau < 1e-6);
}

bool PrimaryVertexAnalyzer4PU::isFinalstateParticle(const HepMC::GenParticle* p) {
  return (!p->end_vertex() && p->status() == 1);
}

bool PrimaryVertexAnalyzer4PU::isCharged(const HepMC::GenParticle* p) {
  const ParticleData* part = pdt_->particle(p->pdg_id());
  if (part) {
    return part->charge() != 0;
  } else {
    // the new/improved particle table doesn't know anti-particles
    return pdt_->particle(-p->pdg_id()) != 0;
  }
}

double PrimaryVertexAnalyzer4PU::vertex_pxy(const reco::Vertex& v) {
  double z = v.z();
  double vxx = v.covariance(iX, iX) + pow(vertexBeamSpot_.BeamWidthX(), 2);
  double vyy = v.covariance(iY, iY) + pow(vertexBeamSpot_.BeamWidthY(), 2);
  ;
  double vxy = v.covariance(iX, iY);
  double dx = v.x() - vertexBeamSpot_.x(z) - dxb_;
  double dy = v.y() - vertexBeamSpot_.y(z) - dyb_;
  double D = vxx * vyy - vxy * vxy;
  double c2xy = pow(dx, 2) * vyy / D + pow(dy, 2) * vxx / D - 2 * dx * dy * vxy / D;
  return TMath::Prob(c2xy, 2);
}

double PrimaryVertexAnalyzer4PU::vertex_aptsum(const reco::Vertex& v) {
  double aptsum = 0.;
  double waptsum = 0.;
  double wsum = 0.;
  for (trackit_t t = v.tracks_begin(); t != v.tracks_end(); t++) {
    aptsum += t->get()->pt();
    double w = v.trackWeight(*t);
    wsum += w;
    waptsum += t->get()->pt() * w;
  }
  return aptsum;
}

bool PrimaryVertexAnalyzer4PU::vertex_time_from_tracks(const reco::Vertex& v,
                                                       Tracks& tracks,
						       double minquality,
                                                       double& t,
                                                       double& tError) {
  double tsum = 0;
  double wsum = 0;
  double w2sum = 0;
  for (auto tk = v.tracks_begin(); tk != v.tracks_end(); tk++) {
    if (v.trackWeight(*tk) > 0.5) {
      auto trk = tracks.from_ref(*tk);
      if (trk.has_timing && (trk.timeQuality >= minquality)) {
        double w = v.trackWeight(*tk) / (trk.dt * trk.dt);
        wsum += w;
        tsum += w * trk.t;
        //w2sum += w * w * trk.dt * trk.dt;
      }
    }
  }

  if (wsum > 0) {
    double t0 = tsum / wsum;
    // robustify
    int nit = 0;
    while ( (nit++) < 50)
      {
	tsum = 0;
	wsum = 0;
	w2sum = 0;
	for (auto tk = v.tracks_begin(); tk != v.tracks_end(); tk++)
	  {
	    if (v.trackWeight(*tk) > 0.5)
	      {
		auto trk = tracks.from_ref(*tk);
		if (trk.has_timing && (trk.timeQuality >= minquality))
		  {
		    double tpull = (trk.t - t0) / trk.dt;
		    if (fabs(tpull) < 5.)
		      {
			double wt = 1./(1.+ exp(0.5 * tpull * tpull - 0.5 * 9.));
			double w = wt * v.trackWeight(*tk) / (trk.dt * trk.dt);
			wsum += w;
			tsum += w * trk.t;
			w2sum += w * w * trk.dt * trk.dt;
		      }
		  }
	      }
	  }

    	if (wsum > 0)
	  {
	    t = tsum / wsum;
	    if (fabs(t - t0) < 1e-3)
	      {
		tError = sqrt(w2sum) / wsum;
		return true;
	      }
	  }
	t0 = t;
      }
  }
  t = 0;
  tError = 1.e10;
  return false;
}

bool PrimaryVertexAnalyzer4PU::vertex_time_from_tracks_pid(const reco::Vertex& v,
                                                       Tracks& tracks,
        					   double minquality,
                                                       double& t,
                                                       double& tError) {
  double tsum = 0;
  double wsum = 0;
  double w2sum = 0;
  //double a[3] = {0.333333333,0.33333333,0.333333333};
  double a[3] = {0.7,0.2,0.1};
  constexpr double cooling_factor = sqrt(0.5);
  // initial guess
  for (auto tk = v.tracks_begin(); tk != v.tracks_end(); tk++) {
    if (v.trackWeight(*tk) > 0.5) {
      auto trk = tracks.from_ref(*tk);
      if (trk.has_timing && (trk.timeQuality >= minquality)) {
	double w = v.trackWeight(*tk);
        wsum += w;
	for(unsigned int j=0; j < 3; j++){
	  tsum += w * trk.th[j] * a[j];
	}
      }
    }
  }
  
  if (wsum > 0) {
    double t0 = tsum / wsum;
    int nit = 0;
    double beta = 1./32.;
    while ( (nit++) < 100)
      {
	tsum = 0;
	wsum = 0;
	w2sum = 0;
	for (auto tk = v.tracks_begin(); tk != v.tracks_end(); tk++)
	  {
	    if (v.trackWeight(*tk) > 0.5)
	      {
		auto trk = tracks.from_ref(*tk);
		if (trk.has_timing && (trk.timeQuality >= minquality))
		  {
		    double dt = trk.MTD_timeerror;
		    double e[3] = {0,0,0};
		    double Z = exp(-beta * 0.5* 3.* 3.);
		    for(unsigned int j = 0; j < 3; j++){
		      double tpull =  (trk.th[j] - t0) / dt;
		      e[j] = exp(- 0.5 * beta * tpull * tpull);
		      Z += a[j] * e[j];
		    }

		    double wsum_trk = 0;
		    for(unsigned int j = 0; j < 3; j++){
		      double wt = a[j] * e[j] / Z;
		      double w = wt * v.trackWeight(*tk) / (dt * dt);
		      wsum_trk += w;
		      tsum += w * trk.th[j]; 
		    }
		    wsum += wsum_trk;
		    w2sum += wsum_trk * wsum_trk * (dt * dt); // 100 % correlation among th[j]
		  }
	      }
	  }
	if (wsum < 1e-10)
	  {
	    report_counted("PrimaryVertexAnalyzer4PU::vertex_time_from_tracks_pid failed while iterating", 10);
	    return false;
	  }

	t = tsum / wsum;
	//cout << " XX " << nit << " T= " << 1/beta << "  t=" << t <<  "    t-t0=" <<  t-t0 << endl;
	if ((fabs(t - t0) < 1e-4) && (beta >= 1.))
	  {
	    tError = sqrt(w2sum) / wsum;
	    return true;
	  }
	
	if ((fabs(t - t0) < 1e-3) && ((beta < 1.)))
	  { 
	    beta = std::min(1., beta / cooling_factor);
	  }

	t0 = t;

      }
  }
  report_counted("PrimaryVertexAnalyzer4PU::vertex_time_from_tracks_pid failed to converge", 10);
  t = 0;
  tError = 1.e10;
  return false;
}


double PrimaryVertexAnalyzer4PU::vertex_ptmax2(const reco::Vertex& v) {
  double ptmax1 = 0;
  double ptmax2 = 0;

  for (trackit_t t = v.tracks_begin(); t != v.tracks_end(); t++) {
    if (v.trackWeight(*t) > 0.5) {
      double pt = t->get()->pt();
      if (pt > ptmax1) {
        ptmax2 = ptmax1;
        ptmax1 = pt;
      } else if (pt > ptmax2) {
        ptmax2 = pt;
      }
    }
  }
  return ptmax2;
}


double PrimaryVertexAnalyzer4PU::vertex_sumpt2(const reco::Vertex& v) {
  double sumpt2 = 0.;

  for (trackit_t t = v.tracks_begin(); t != v.tracks_end(); t++) {
    // no cut on the track weight here!
    double pt = t->get()->pt();
    sumpt2 += pt * pt;
  }
  return sqrt(sumpt2);
}

double PrimaryVertexAnalyzer4PU::vertex_sumpt(const reco::Vertex& v) {
  double sumpt = 0.;

  for (trackit_t t = v.tracks_begin(); t != v.tracks_end(); t++) {
    // no cut on the track weight here!
    sumpt += t->get()->pt();
  }
  return sumpt;
}

double PrimaryVertexAnalyzer4PU::vertex_r(const reco::Vertex& v) {
  double z = v.z();
  double dx = v.x() - vertexBeamSpot_.x(z) - dxb_;
  double dy = v.y() - vertexBeamSpot_.y(z) - dyb_;
  return sqrt(dx * dx + dy * dy);
}

bool PrimaryVertexAnalyzer4PU::select(const reco::Vertex& v, int level) {
  /* level
     0  !isFake  && ndof>4  (default)
     1  !isFake  && ndof>4 && prob > 0.01 
     2  !isFake  && ndof>4 && prob > 0.01 && ptmax2 > 0.4
  */
  if (v.isFake())
    return false;
  if ((level == 0) && (v.ndof() > selNdof_))
    return true;
  if ((level == 1) && (v.ndof() > selNdof_) && (vertex_pxy(v) > 0.01))
    return true;
  if ((level == 2) && (v.ndof() > selNdof_) && (vertex_pxy(v) > 0.01) && (vertex_ptmax2(v) > 0.4))
    return true;
  if ((level == 3) && (v.ndof() > selNdof_) && (vertex_ptmax2(v) < 0.4))
    return true;
  return false;
}





void PrimaryVertexAnalyzer4PU::fillVertexHistosNoTracks(std::map<std::string, TH1*>& h,
                                                        const std::string& vtype,
                                                        const reco::Vertex* v,
							const int index,
                                                        const double deltaz,
                                                        const bool verbose)
// fill vertex quantities available from the reco vertex itself
{
  timer_start("fillVertexHistosNoTracks");
  // delta z = z_this - z_other
  double z = v->z();
  double vxx = v->covariance(iX, iX) + pow(vertexBeamSpot_.BeamWidthX(), 2);
  double vyy = v->covariance(iY, iY) + pow(vertexBeamSpot_.BeamWidthY(), 2);
  ;
  double vxy = v->covariance(iX, iY);
  double dx = v->x() - vertexBeamSpot_.x(z) - dxb_;
  double dy = v->y() - vertexBeamSpot_.y(z) - dyb_;
  double D = vxx * vyy - vxy * vxy;
  double c2xy = pow(dx, 2) * vyy / D + pow(dy, 2) * vxx / D - 2 * dx * dy * vxy / D;

  if (index >= 0 ) {
    Fill(h, vtype + "/index", float(index));
    Fill(h, vtype + "/logndofvsindex", float(index), log(v->ndof()) / log(10.));
    Fill(h, vtype + "/ndofvsindex", float(index), v->ndof());
  }
  Fill(h, vtype + "/c2xy", c2xy);
  Fill(h, vtype + "/chi2", v->chi2());
  if (v->tracksSize() > 0)
    {
      Fill(h, vtype + "/chi2overntk", v->chi2() / v->tracksSize());
    }
  Fill(h, vtype + "/probxy", TMath::Prob(c2xy, 2));
  Fill(h, vtype + "/r", sqrt(dx * dx + dy * dy));
  Fill(h, vtype + "/zpullbeam", (z - vertexBeamSpot_.z0() - dzb_) / sigmaZ_);
  Fill(h, vtype + "/ndof", v->ndof());
  Fill(h, vtype + "/ndofvspu", lumiPU_, v->ndof());

  if (v->ndof() > 0) {
    Fill(h, vtype + "/logndof", log(v->ndof()) / log(10.));
  }
  Fill(h, vtype + "/vtxndfvsntrk", v->tracksSize(), v->ndof());
  Fill(h, vtype + "/avweight", (v->ndof() - ndof0trk_) / (2. * v->tracksSize()));
  Fill(h, vtype + "/avweightvsndof", (v->ndof() - ndof0trk_) / (2. * v->tracksSize()), v->ndof());
  Fill(h, vtype + "/avweightX", (v->ndof() - ndof0trk_) / (2. * v->tracksSize()), v->ndof()); // avweight weighted with ndof
  Fill(h, vtype + "/errx", v->xError());
  Fill(h, vtype + "/erry", v->yError());
  Fill(h, vtype + "/errz", v->zError());

  Fill(h, vtype + "/zvtx", z);
  if (f4D_) {
    // contains vertices that don't actually have timing information
    Fill(h, vtype + "/terrvtx", v->tError());
    if(v->tError() < 0.1){
      Fill(h, vtype + "/tvtx", v->t());
    }
  }

  if (deltaz != 0) {
    Fill(h, vtype + "/dznearest", deltaz);

    if (verbose) {
      cout << "fillVertexHistos  z=" << setw(10) << setprecision(4) << fixed << v->z() << "    deltaz = " << deltaz
           << endl;
    }
  }
  timer_stop("fillVertexHistosNoTracks");
}



void PrimaryVertexAnalyzer4PU::fillVertexHistos(std::map<std::string, TH1*>& h,
                                                const std::string& vtype,
                                                const reco::Vertex* v,
                                                Tracks& tracks,
						const int index,
                                                const double deltaz,
                                                const bool verbose) {
  // delta z = z_this - z_other
  timer_start("fillVertexHistos");
  fillVertexHistosNoTracks(h, vtype, v, index, deltaz, verbose);
  Fill(h, vtype + "/ptmax2", vertex_ptmax2(*v));
  double sumpt2 = vertex_sumpt2(*v);
  double sumpt = vertex_sumpt(*v);  // not weighted, unlike sumapt
  Fill(h, vtype + "/sumpt2", sumpt2);
  Fill(h, vtype + "/logsumpt2", log(sumpt2)/log(10.));
  Fill(h, vtype + "/sumpt2vssumpt", sumpt, sumpt2);
  if(sumpt >0){
    Fill(h, vtype + "/sumpt2oversumpt", sumpt2/sumpt);
    Fill(h, vtype + "/sumpt2oversumptvssumpt2", sumpt2, sumpt2/sumpt);
  }
  Fill(h, vtype + "/nseltrkvtx", float(v->tracksSize()));



  int nlowt = 0;
  for (trackit_t t = v->tracks_begin(); t != v->tracks_end(); t++) {
    auto tk = tracks.from_ref(*t);
    fillTrackHistos(h, vtype, tk, v);   // select based on weight here?
    Fill(h, vtype + "/trkweight", v->trackWeight(*t));
    Fill(h, vtype + "/trkptnorm", t->get()->pt(), 1./  v->tracksSize());
    if ( v->trackWeight(*t) < 0.6) nlowt += 1;
  }
  Fill(h, vtype + "/numlowttrk", float(nlowt) );
  if (v->tracksSize() > 0) Fill(h, vtype + "/fraclowttrk", float(nlowt) / v->tracksSize());

  /* consider dropping this
  if (f4D_) {
    double t, tError;
    if (vertex_time_from_tracks(*v, tracks, 0, t, tError)) {
      Fill(h, vtype + "/tvtx_fromtracks", t);
      Fill(h, vtype + "/terrvtx_fromtracks", tError);
      Fill(h, vtype + "/tvtx_withtracks", v->t());
      Fill(h, vtype + "/terrvtx_withtracks", v->tError());
    }
  }
  */
  
  if (deltaz != 0) {
    double pzsum = 0;
    double ptsum = 0;
    double apzsum = 0;
    double aptsum = 0;
    double wpzsum = 0, wptsum = 0, wapzsum = 0, waptsum = 0;

    for (trackit_t t = v->tracks_begin(); t != v->tracks_end(); t++) {
      int zsign = (t->get()->pz() * deltaz) > 0 ? 1 : -1;
      double w = v->trackWeight(*t) / pow(t->get()->dzError(), 2);
      wptsum += t->get()->pt() * zsign * w;
      wpzsum += t->get()->pz() * w;
      waptsum += t->get()->pt() * w;
      wapzsum += std::abs(t->get()->pz()) * v->trackWeight(*t);

      if (v->trackWeight(*t) > 0.5) {
        pzsum += t->get()->pz();
        ptsum += t->get()->pt() * zsign;
        apzsum += std::abs(t->get()->pz());
        aptsum += t->get()->pt();
      }

    }

    if (deltaz < 0) {  
      pzsum = -pzsum;
      wpzsum = -wpzsum;
    }

    Fill(h, vtype + "/pzsum", pzsum);
    Fill(h, vtype + "/ptsum", ptsum);
    Fill(h, vtype + "/apzsum", apzsum);
    Fill(h, vtype + "/aptsum", aptsum);
    Fill(h, vtype + "/waptsum", waptsum);

    if (aptsum > 0) {
      Fill(h, vtype + "/apt", ptsum / aptsum);
      Fill(h, vtype + "/aptvsdz", deltaz, ptsum / aptsum);
      Fill(h, vtype + "/wapt", wptsum / waptsum);
    }

    if (apzsum > 0) {
      Fill(h, vtype + "/apz", pzsum / apzsum);
      Fill(h, vtype + "/azvsdz", deltaz, pzsum / apzsum);
      Fill(h, vtype + "/wapz", wpzsum / wapzsum);
    }

    if ((apzsum > 0) && (aptsum > 0)) {
      Fill(h, vtype + "/aptvsapz", pzsum / apzsum, ptsum / aptsum);
    }
  }

  // sphericity
  double Sxx = 0, Syy = 0, Sxy = 0, Spt = 0, Sw = 0;
  for (trackit_t t = v->tracks_begin(); t != v->tracks_end(); t++) {
    double pt = t->get()->pt();
    double px = pt * cos(t->get()->phi());
    double py = pt * sin(t->get()->phi());
    double w = v->trackWeight(*t);
    Sw += w;
    Spt += w * pt;
    Sxx += w * px * px / pt;
    Syy += w * py * py / pt;
    Sxy += w * px * py / pt;
  }
  double Ds = sqrt(pow(Sxx - Syy, 2) + 4 * Sxy);
  double l1 = (Sxx + Syy + Ds) / (2 * Spt);
  double l2 = (Sxx + Syy - Ds) / (2 * Spt);
  double SLT = 2 * l2 / (l1 + l2);
  Fill(h, vtype + "/sphericity", SLT);
  Fill(h, vtype + "/sphericityvsntrk", Sw, SLT);
  timer_stop("fillVertexHistos");
}


void PrimaryVertexAnalyzer4PU::fillVertexHistosMatched(std::map<std::string, TH1*>& h,
						       const std::string& vtype,
						       const reco::Vertex* v,
						       Tracks& tracks,
						       const RSmatch& rs,
						       const std::vector<SimEvent>& simEvt,
						       const int index,
						       const double deltaz,
						       const bool verbose) {
  fillVertexHistos(h, vtype, v, tracks, index, deltaz, verbose);


  // now fill properties wrt the matched simvertex
  assert(rs.sim != NOT_MATCHED);
  const SimEvent& simevt = simEvt[rs.sim];
  
  unsigned int ntiming = 0;
  unsigned int ntiming_qual05 = 0;
  unsigned int ntiming_qual08 = 0;
  unsigned int nseltrk = 0;
  for (trackit_t tv = v->tracks_begin(); tv != v->tracks_end(); tv++) {
    nseltrk++;
    if(f4D_){
      auto tk = tracks.from_ref(*tv);
      if (tk.has_timing) {
	ntiming++;
	if (tk.timeQuality > 0.5) ntiming_qual05++;
	if (tk.timeQuality > 0.8) ntiming_qual08++;
      }
    }
  }
  // these histograms are also filled for unmatched vertices, the false flag prevents double bookings
  Fill(h, vtype + "/nseltrkvtx", float(nseltrk), simevt.is_signal(), false);
  Fill(h, vtype + "/ptmax2", vertex_ptmax2(*v), simevt.is_signal(), false);
  /* OBSOLETE
  if (simevt.is_signal()){ // not using the automatic version because "ptmax2" (no suffix) is filled in fillVertexHistos
    Fill(h, vtype + "/ptmax2Signal", vertex_ptmax2(*v));
    Fill(h, vtype + "/sumpt2Signal", vertex_sumpt2(*v));
  }else{
    Fill(h, vtype + "/ptmax2PU", vertex_ptmax2(*v));
    Fill(h, vtype + "/sumpt2PU", vertex_sumpt2(*v));
  }
  */
  
  // note that this is redundant with some histos filled in analyzevetexcollectiontp
  
  double xsim = simevt.x;
  double ysim = simevt.y;
  double zsim = simevt.z;
  Fill(h, vtype + "/xrecsim", v->x() - xsim, simevt.is_signal());
  Fill(h, vtype + "/yrecsim", v->y() - ysim, simevt.is_signal());
  Fill(h, vtype + "/zrecsim", v->z() - zsim, simevt.is_signal());
  Fill(h, vtype + "/xrecerr", v->xError(), simevt.is_signal());  // like xerr but separated by signal/pu
  Fill(h, vtype + "/yrecerr", v->yError(), simevt.is_signal());
  Fill(h, vtype + "/zrecerr", v->zError(), simevt.is_signal());
  Fill(h, vtype + "/zrecsimHR", v->z() - zsim, simevt.is_signal());
  Fill(h, vtype + "/xrecsimpull", (v->x() - xsim) / v->xError(), simevt.is_signal());
  Fill(h, vtype + "/yrecsimpull", (v->y() - ysim) / v->yError(), simevt.is_signal());
  Fill(h, vtype + "/zrecsimpull", (v->z() - zsim) / v->zError(), simevt.is_signal());

  if (f4D_)
    {
      Fill(h, vtype + "/ntimingvtx", float(ntiming), simevt.is_signal());
      Fill(h, vtype + "/ntimingqual05vtx", float(ntiming_qual05), simevt.is_signal());
      Fill(h, vtype + "/ntimingqual08vtx", float(ntiming_qual08), simevt.is_signal());
      
      double t_fromtracks, tError_fromtracks;
      bool timing_from_tracks = vertex_time_from_tracks(*v, tracks, 0, t_fromtracks, tError_fromtracks);
      
      double t_fromtracks_pid, tError_fromtracks_pid;
      bool timing_from_tracks_pid = vertex_time_from_tracks_pid(*v, tracks, 0, t_fromtracks_pid, tError_fromtracks_pid);
      
      double t_fromtracks_qual_pid, tError_fromtracks_qual_pid;
      bool timing_from_tracks_qual_pid = vertex_time_from_tracks_pid(*v, tracks, 0.8, t_fromtracks_qual_pid, tError_fromtracks_qual_pid);

      double tsim = simevt.t;
      Fill(h, vtype + "/trecsim", v->t() - tsim, simevt.is_signal());
      Fill(h, vtype + "/trecerr", v->tError(), simevt.is_signal());
      Fill(h, vtype + "/trecsimpull", (v->t() - tsim) / v->tError(), simevt.is_signal());
      Fill(h, vtype + "/trecsimpullwide", (v->t() - tsim) / v->tError(), simevt.is_signal());
      Fill(h, vtype + "/trecerrvsntrkprof", float(ntiming), v->tError(), simevt.is_signal());
      Fill(h, vtype + "/trecerrvsntrk", float(ntiming), v->tError(), simevt.is_signal());
      if( v->tError() < 0.1){
	Fill(h, vtype + "/trecsim_sigmatlt01", v->t() - tsim, simevt.is_signal());
	Fill(h, vtype + "/trecsimpull_sigmatlt01", (v->t() - tsim) / v->tError(), simevt.is_signal());
      }
      
      if (timing_from_tracks)
	{
	  // timing from tracks
	  Fill(h, vtype + "/trecsim_fromtracks" , t_fromtracks - simevt.t, simevt.is_signal());
	  Fill(h, vtype + "/trecerr_fromtracks", tError_fromtracks, simevt.is_signal());
	  Fill(h, vtype + "/trecsimpull_fromtracks", (t_fromtracks - simevt.t) / tError_fromtracks, simevt.is_signal());

	  // also fill histos with the default values for the same list of vertices for comparison
	  Fill(h, vtype + "/trecsim_withtracks", v->t() - tsim, simevt.is_signal());
	  Fill(h, vtype + "/trecerr_withtracks", v->tError(), simevt.is_signal());
	  Fill(h, vtype + "/trecsimpull_withtracks", (v->t() - tsim) / v->tError(), simevt.is_signal());
	}
      else
	{
	  report_counted("fillVertexHistosMatched: timing from tracks failed",1);
	  //reportVertex(*v, Form("timing from tracks failed,  with tracks: %f", v->tError()), false);
	}

      if (timing_from_tracks_pid)
	{
	  // timing from tracks
	  Fill(h, vtype + "/trecsim_fromtracks_pid" , t_fromtracks_pid- simevt.t, simevt.is_signal());
	  Fill(h, vtype + "/trecerr_fromtracks_pid", tError_fromtracks_pid, simevt.is_signal());
	  Fill(h, vtype + "/trecsimpull_fromtracks_pid", (t_fromtracks_pid - simevt.t) / tError_fromtracks_pid, simevt.is_signal());

	  // also fill histos with the default values for the same list of vertices for comparison
	  Fill(h, vtype + "/trecsim_withtracks_pid", v->t() - tsim, simevt.is_signal());
	  Fill(h, vtype + "/trecerr_withtracks_pid", v->tError(), simevt.is_signal());
	  Fill(h, vtype + "/trecsimpull_withtracks_pid", (v->t() - tsim) / v->tError(), simevt.is_signal());
	}
      else
	{
	  report_counted("fillVertexHistosMatched: timing from tracks failed",1);
	  //reportVertex(*v, Form("timing from tracks failed,  with tracks: %f", v->tError()), false);
	}

      
      if (timing_from_tracks_qual_pid)
	{
	  // timing from tracks
	  Fill(h, vtype + "/trecsim_fromtracks_qual_pid" , t_fromtracks_qual_pid - simevt.t, simevt.is_signal());
	  Fill(h, vtype + "/trecerr_fromtracks_qual_pid", tError_fromtracks_qual_pid, simevt.is_signal());
	  Fill(h, vtype + "/trecsimpull_fromtracks_qual_pid", (t_fromtracks_qual_pid - simevt.t) / tError_fromtracks_qual_pid, simevt.is_signal());

	  // also fill histos with the default values for the same list of vertices for comparison
	  Fill(h, vtype + "/trecsim_withtracks_qual_pid", v->t() - tsim, simevt.is_signal());
	  Fill(h, vtype + "/trecerr_withtracks_qual_pid", v->tError(), simevt.is_signal());
	  Fill(h, vtype + "/trecsimpull_withtracks_qual_pid", (v->t() - tsim) / v->tError(), simevt.is_signal());
	}
      else
	{
	  report_counted("fillVertexHistosMatched: timing from tracks failed",1);
	  //reportVertex(*v, Form("timing from tracks failed,  with tracks: %f", v->tError()), false);
	}

      
    }
}


void PrimaryVertexAnalyzer4PU::fillTrackHistos(std::map<std::string, TH1*>& h,
                                               const std::string& ttype,
                                               RecoTrack& tk,
                                               const reco::Vertex* v) {
  if(! fill_track_histos_) return;
  timer_start("fillTrackHistos");
  fillRecoTrackHistos(h, ttype, *(tk.trk));
  
  if (f4D_ && tk.has_timing)
    {
      Fill(h, ttype + "/t0trk", tk.t);
      Fill(h, ttype + "/t0errtrk", tk.dt);
      Fill(h, ttype + "/t0qualtrk", tk.timeQuality);
      if (! (v==NULL) && (v->isValid()) && (v->ndof() > 10.) ){
	double deltat = tk.t - v->t();
	Fill(h, ttype + "/trestrk", deltat);
	Fill(h, ttype + "/tpulltrk", deltat / sqrt(pow(tk.dt,2) + pow(v->tError(), 2)));
      }
    }

  
  if ((!(v == NULL)) && (v->isValid()) && (v->ndof() > 10.)) {
    double zvtx = v->position().z();
    double dz2 = tk.dz * tk.dz + (pow(wx_ * cos(tk.phi), 2) + pow(wy_ * sin(tk.phi), 2)) / pow(tan(tk.theta), 2); // really?
    Fill(h, ttype + "/zrestrk", tk.z - zvtx);
    Fill(h, ttype + "/zrestrkvsphi", tk.phi, tk.z - zvtx);
    Fill(h, ttype + "/zrestrkvseta", tk.eta, tk.z - zvtx);
    Fill(h, ttype + "/zpulltrkvsphi", tk.phi, (tk.z - zvtx) / sqrt(dz2));
    Fill(h, ttype + "/zrestrkvsz", zvtx, tk.z - zvtx);
    Fill(h, ttype + "/zpulltrkvsz", zvtx, (tk.z - zvtx) / sqrt(dz2));
    Fill(h, ttype + "/zpulltrkvseta", tk.eta, (tk.z - zvtx) / sqrt(dz2));
    Fill(h, ttype + "/zpulltrkvsz", tk.z, pow(tk.z - zvtx, 2) / (pow(v->zError(), 2) + dz2));
    Fill(h, ttype + "/zpulltrk", (tk.z - zvtx) / sqrt(pow(v->zError(), 2) + dz2));
  }
      
  fillTrackHistosMatched(h, ttype, tk);
  
  timer_stop("fillTrackHistos");

}


// basically obsolete, all of this can be done with RecoTracks
void PrimaryVertexAnalyzer4PU::fillTransientTrackHistos(std::map<std::string, TH1*>& h,
                                                        const std::string& ttype,
                                                        const TransientTrack* tt,
                                                        const reco::Vertex* v) {

  if ((!(v == NULL)) && (v->isValid()) && (v->ndof() > 10.)) {
    double z = (tt->stateAtBeamLine().trackStateAtPCA()).position().z();
    double zvtx = v->position().z();
    double tantheta = tan((tt->stateAtBeamLine().trackStateAtPCA()).momentum().theta());
    double phi = (tt->stateAtBeamLine().trackStateAtPCA()).momentum().phi();
    double eta = (tt->stateAtBeamLine().trackStateAtPCA()).momentum().eta();
    double dz2 = pow(tt->track().dzError(), 2) + (pow(wx_ * cos(phi), 2) + pow(wy_ * sin(phi), 2)) / pow(tantheta, 2);

    // note that these are relative to a reconstructed vertex, (passed as an argument), not wrt MC truth
    Fill(h, ttype + "/zrestrk", z - zvtx);
    Fill(h, ttype + "/zrestrkvsphi", phi, z - zvtx);
    Fill(h, ttype + "/zrestrkvseta", eta, z - zvtx);
    Fill(h, ttype + "/zpulltrkvsphi", phi, (z - zvtx) / sqrt(dz2));
    Fill(h, ttype + "/zrestrkvsz", zvtx, z - zvtx);
    Fill(h, ttype + "/zpulltrkvsz", zvtx, (z - zvtx) / sqrt(dz2));
    Fill(h, ttype + "/zpulltrkvseta", eta, (z - zvtx) / sqrt(dz2));
    Fill(h, ttype + "/zpulltrkvsz", z, pow(z - zvtx, 2) / (pow(v->zError(), 2) + dz2));
    Fill(h, ttype + "/zpulltrk", (z - zvtx) / sqrt(pow(v->zError(), 2) + dz2));

    if (f4D_ && (tt->dtErrorExt() < 0.3) && (tt->dtErrorExt() > 0.0)) {
      double deltat = tt->timeExt() - v->t();
      Fill(h, ttype + "/trestrk", deltat);
      Fill(h, ttype + "/tpulltrk", deltat / sqrt(pow(tt->dtErrorExt(), 2) + pow(v->tError(), 2)));
    }

  }
}




void PrimaryVertexAnalyzer4PU::fillRecoTrackHistos(std::map<std::string, TH1*>& h,
                                                   const std::string& ttype,
                                                   const reco::Track& t) {
  /* reco::track version, no transient tracks here */
  Fill(h, ttype + "/eta", t.eta());
  Fill(h, ttype + "/ztrk", t.vz());
  Fill(h, ttype + "/zerrtrk", t.dzError());
  Fill(h, ttype + "/phi", t.phi());
  Fill(h, ttype + "/pt", t.pt());
  if (fabs(t.eta()) > 2.0)
    Fill(h, ttype + "/ptfwd", t.pt());
  if (fabs(t.eta()) < 1.0)
    Fill(h, ttype + "/ptcentral", t.pt());
  Fill(h, ttype + "/logpt", log(t.pt()) / log(10.));
  Fill(h, ttype + "/logpteta", t.eta(), log(t.pt()) / log(10.));
  Fill(h, ttype + "/z-eta", t.eta(), t.vz());
  Fill(h, ttype + "/z-pt", t.vz(), t.pt());
  Fill(h, ttype + "/z-logpt", t.vz(), log(t.pt())/log(10.));
  Fill(h, ttype + "/found", t.found());
  Fill(h, ttype + "/lost", t.lost());
  Fill(h, ttype + "/validfraction", t.validFraction());
  Fill(h, ttype + "/nchi2", t.normalizedChi2());
  Fill(h, ttype + "/nchi2vsz", t.vz(), t.normalizedChi2());
  if (RECO_)
    Fill(h, ttype + "/rstart", (t.innerPosition()).Rho());  // innerPosition need TrackExtra

  double d0Error = t.d0Error();
  double d0 = t.dxy(vertexBeamSpot_.position());
  if (d0Error > 0) {
    Fill(h, ttype + "/logtresxy", log(d0Error / 0.0001) / log(10.));
    Fill(h, ttype + "/tpullxy", d0 / d0Error);
    Fill(h, ttype + "/tpullxyvsz", t.vz(), pow(d0 / d0Error, 2));
    Fill(h, ttype + "/tpullxyvseta", t.eta(), pow(d0 / d0Error, 2));
  }

  double dzError = t.dzError();
  if (dzError > 0) {
    Fill(h, ttype + "/logtresz", log(dzError / 0.0001) / log(10.)); // in microns
  }

  //
  Fill(h, ttype + "/zerrtrk_withbeam", sqrt(pow(t.dzError(), 2) + wxy2_ / pow(tan(t.theta()), 2)));

  // collect some info on hits and clusters
  Fill(h, ttype + "/nbarrelLayers", static_cast<double>(t.hitPattern().pixelBarrelLayersWithMeasurement()));
  Fill(h, ttype + "/nPxLayers", static_cast<double>(t.hitPattern().pixelLayersWithMeasurement()));
  if (fabs(t.eta() < 2))
    Fill(h, ttype + "/nPxLayersVsPt", t.pt(), static_cast<double>(t.hitPattern().pixelLayersWithMeasurement()));
  Fill(h, ttype + "/nSiLayers", static_cast<double>(t.hitPattern().trackerLayersWithMeasurement()));
  Fill(h, ttype + "/n3dLayers", static_cast<double>(t.hitPattern().numberOfValidStripLayersWithMonoAndStereo()));
  Fill(h, ttype + "/nOT3dLayers",static_cast<double>(t.hitPattern().numberOfValidStripLayersWithMonoAndStereo())
       - static_cast<double>(t.hitPattern().pixelLayersWithMeasurement())  );
  //Fill(h, ttype + "/expectedInner_"+ttype,static_cast<double>(t.trackerExpectedHitsInner().numberOfAllHits())); FIXME
  //Fill(h, ttype + "/expectedOuter_"+ttype,static_cast<double>(t.trackerExpectedHitsOuter().numberOfAllHits())); FIXME
  double mi = t.hitPattern().numberOfLostHits(HitPattern::MISSING_INNER_HITS); // changed from numberOfAllHits
  double mo = t.hitPattern().numberOfLostHits(HitPattern::MISSING_OUTER_HITS);
  Fill(h, ttype + "/missing_inner", mi);
  Fill(h, ttype + "/missing_outer", mo);
  Fill(h, ttype + "/npxminmiss", static_cast<double>(t.hitPattern().pixelLayersWithMeasurement())-mi);
  Fill(h, ttype + "/missing_fraction",
       float(mi + mo) / float(mi + mo + t.hitPattern().numberOfAllHits(HitPattern::TRACK_HITS)));
  Fill(h, ttype + "/trackAlgo", static_cast<double>(t.algo()));
  Fill(h, ttype + "/trackQuality", static_cast<double>(t.qualityMask()));

  //-------------------------------------------------------------------
  if (RECO_) {
    //
    int longesthit = 0, nbarrel = 0;
    for (trackingRecHit_iterator hit = t.recHitsBegin(); hit != t.recHitsEnd(); hit++) {
      if ((**hit).isValid() && (**hit).geographicalId().det() == DetId::Tracker) {
        bool barrel = (**hit).geographicalId().subdetId() == static_cast<int>(PixelSubdetector::PixelBarrel);
        if (barrel) {
          const SiPixelRecHit* pixhit = dynamic_cast<const SiPixelRecHit*>(&(**hit));
          edm::Ref<edmNew::DetSetVector<SiPixelCluster>, SiPixelCluster> const& clust = (*pixhit).cluster();
          if (clust.isNonnull()) {
            nbarrel++;
            if (clust->sizeY() > longesthit)
              longesthit = clust->sizeY();
            if (clust->sizeY() > 20.) {
              Fill(h, ttype + "/lvseta", t.eta(), 19.9);
              Fill(h, ttype + "/lvstanlambda", tan(t.lambda()), 19.9);
            } else {
              Fill(h, ttype + "/lvseta", t.eta(), float(clust->sizeY()));
              Fill(h, ttype + "/lvstanlambda", tan(t.lambda()), float(clust->sizeY()));
            }
          }
        }
      }
    }

    if (nbarrel > 0) {
      Fill(h, ttype + "/longestbarrelhit", float(longesthit));
    }
    Fill(h, ttype + "/nbarrelhits", float(nbarrel));
  }
  //-------------------------------------------------------------------
}


void PrimaryVertexAnalyzer4PU::fillTrackHistosMatched(std::map<std::string, TH1*>& h,
                                               const std::string& ttype,
                                               RecoTrack& tk) 
{
  if(!(tk.matched)) return;
  Fill(h, ttype + "/tkzrecsim", tk.z - tk.zsim);
  Fill(h, ttype + "/tkzrecsimpull", (tk.z - tk.zsim)/tk.dz);
  Fill(h, ttype + "/tkzrecsimpullvseta", tk.eta, (tk.z - tk.zsim)/tk.dz);
  Fill(h, ttype + "/tkzrecsimpullsqvseta", tk.eta,    pow((tk.z - tk.zsim)/tk.dz, 2));
  Fill(h, ttype + "/tkzrecsimvseta", tk.eta, tk.z - tk.zsim);
  Fill(h, ttype + "/tkzrecsimvseta2d", tk.eta, tk.z - tk.zsim);
  Fill(h, ttype + "/tkzrecsimvsz", tk.zsim, tk.z - tk.zsim);
  double flip = tk.eta < 0 ? -1. : 1.;
  Fill(h, ttype + "/tkzrecsimvslogpt", log(tk.pt)/log(10.), (tk.z - tk.zsim) * flip );
  flip = tk.z < 0 ? -1 : 1.;
  Fill(h, ttype + "/tkzrecsimvsetaz", tk.eta, (tk.z - tk.zsim) * flip);
  
  if (f4D_ && tk.has_timing){
    Fill(h, ttype + "/tktrecsim", tk.t - tk.tsim);
    Fill(h, ttype + "/tktrecsimpull", (tk.t - tk.tsim)/tk.dt);
    Fill(h, ttype + "/tktrecsimpullwide", (tk.t - tk.tsim)/tk.dt);
    Fill(h, ttype + "/tktrecsimvseta2d", tk.eta, tk.t - tk.tsim);
    Fill(h, ttype + "/tktrecsimpullsqvserr", tk.dt, pow((tk.t - tk.tsim) / tk.dt, 2));
    Fill(h, ttype + "/tktrecsimpullvserr", tk.dt, (tk.t - tk.tsim) / tk.dt);
  }
}





void PrimaryVertexAnalyzer4PU::fillTrackClusterHistos(std::map<std::string, TH1*>& h,
                                                      const std::string& ttype,
                                                      const reco::Track& t,
                                                      const reco::Vertex* v) {
  if (!RECO_)
    return;

  for (trackingRecHit_iterator hit = t.recHitsBegin(); hit != t.recHitsEnd(); hit++) {
    const DetId& hit_detId = (**hit).geographicalId();

    if ((**hit).isValid() && hit_detId.det() == DetId::Tracker) {
      bool barrel = (**hit).geographicalId().subdetId() == static_cast<int>(PixelSubdetector::PixelBarrel);

      if (barrel) {
        //PXBDetId bid = (PXBDetId) (**hit).geographicalId();//DataFormats/SiPixelDetId/interface/PXBDetId.h
        const SiPixelRecHit* pixhit = dynamic_cast<const SiPixelRecHit*>(&(**hit));
        edm::Ref<edmNew::DetSetVector<SiPixelCluster>, SiPixelCluster> const& clust = (*pixhit).cluster();

        int layer = tTopo_->pxbLayer(hit_detId);
        int ladder = tTopo_->pxbLadder(hit_detId);
        int zindex = tTopo_->pxbModule(hit_detId);
        //int side = 0;
        //if(zindex<5) side=1; else side=2;

        int roc = (clust->minPixelRow() < 80) ? int(clust->minPixelCol() / 52) : (15 - int(clust->minPixelCol() / 52));

        int core = -1;
        if (layer == 1) {
          core = ((ladder - 1) * 8 + (zindex - 1)) * 4 + int(roc / 4);
        } else if (layer == 2) {
          core = 96 * 4 + ((ladder - 1) * 8 + (zindex - 1)) * 2 + int(roc / 8);
        } else if (layer == 3) {
          core = 96 * 4 + 2 * 8 * 28 + ((ladder - 1) * 8 + (zindex - 1)) * 2 + int(roc / 8);
        } else if (layer == 4) {
          core = 96 * 4 + 2 * 8 * 28 + 2 * 8 * 44 + ((ladder - 1) * 8 + (zindex - 1)) * 2 + int(roc / 8);
        }

        /*
	if ((core<0) || (core>2559) || (roc<0) || (roc>15) ){
	  cout << " bad trafo " << core << " "  << roc << endl;
	}
	*/

        if (clust.isNonnull()) {
          Fill(h, "cores_live", float(luminosityBlock_), float(core));
          /*
	  cout << "  min row " << clust->minPixelRow() << "  minCol " << clust->minPixelCol() 
	    //<< "  layer=" << bid.layer() << "  ladder = " << bid.ladder()  << "  module = " << bid.module() 
	       << " topo :  layer=" << layer << "   ladder =  " << ladder << "   zindex = " << zindex  << "   core= " << core << "  roc = " << roc
	       << endl; 
	  */
        }
      }
    }
  }
}

void PrimaryVertexAnalyzer4PU::printRecVtxs(const reco::VertexCollection* recVtxs, std::string title) {
  int ivtx = 0;
  std::cout << std::endl;
  std::cout << std::endl << title << "  nv=" << recVtxs->size() << "" << std::endl;

  for (reco::VertexCollection::const_iterator v = recVtxs->begin(); v != recVtxs->end(); ++v) {
    string vtype = " recvtx  ";
    if (v->isFake()) {
      vtype = " fake   ";
    } else if (v->ndof() == -3) {
      vtype = " event   ";
    }
    std::cout << "vtx " << std::setw(3) << std::setfill(' ') << ivtx++ << vtype << " #trk " << std::fixed
              << std::setprecision(4) << std::setw(3) << v->tracksSize() << " chi2 " << std::fixed << std::setw(5)
              << std::setprecision(1) << v->chi2() << " ndof " << std::fixed << std::setw(6) << std::setprecision(2)
              << v->ndof() << "  x=" << std::setw(7) << std::fixed << std::setprecision(4) << v->x() << " +/-"
              << std::setw(6) << v->xError()                                                // <<  std::endl
              << "  y=" << std::setw(7) << v->y() << " +/-" << std::setw(6) << v->yError()  //<< std::endl
              << "  z=" << std::setw(8) << v->z() << " +/-" << std::setw(6) << v->zError();
    if(f4D_){
      std::cout << " t=" << std::setw(7) << v->t() << " +/-" << std::setw(5) << v->tError(); 
    }
    std::cout << "  dxy = ("
              << std::setw(7) << std::fixed << std::setprecision(4) << v->x() - vertexBeamSpot_.x(v->z()) << ","
              << std::setw(7) << std::fixed << std::setprecision(4) << v->y() - vertexBeamSpot_.y(v->z()) << ")"
              << " pxy= " << std::setw(7) << std::fixed << std::setprecision(4) << vertex_pxy(*v)
              << "  r= " << std::setw(7) << std::fixed << std::setprecision(4) << vertex_r(*v)
              << " ptsum= " << std::setw(6)  << std::fixed << std::setprecision(1) << vertex_aptsum(*v)
              << " sumpt2= " << std::setw(6)  << std::fixed << std::setprecision(1) << vertex_sumpt2(*v)
              << " ptmax2= " << std::setw(5) << std::fixed << std::setprecision(2) << vertex_ptmax2(*v) << std::endl;
  }
  std::cout << std::endl;
}

void PrimaryVertexAnalyzer4PU::printSimVtxs(const Handle<SimVertexContainer> simVtxs) {
  int i = 0;
  for (SimVertexContainer::const_iterator vsim = simVtxs->begin(); vsim != simVtxs->end(); ++vsim) {
    if (vsim->position().x() * vsim->position().x() + vsim->position().y() * vsim->position().y() < 1.) {
      std::cout << i++ << ")" << std::scientific << " evtid=" << vsim->eventId().event() << ","
                << vsim->eventId().bunchCrossing() << " sim x=" << vsim->position().x() * simUnit_
                << " sim y=" << vsim->position().y() * simUnit_ << " sim z=" << vsim->position().z() * simUnit_
                << " sim t=" << vsim->position().t() << " parent=" << vsim->parentIndex() << std::endl;
    }
  }
}

void PrimaryVertexAnalyzer4PU::printSimTrks(const Handle<SimTrackContainer> simTrks) {
  std::cout << " simTrks   type, (momentum), vertIndex, genpartIndex" << std::endl;
  int i = 1;
  for (SimTrackContainer::const_iterator t = simTrks->begin(); t != simTrks->end(); ++t) {
    std::cout << i++ << ")" << t->eventId().event() << "," << t->eventId().bunchCrossing() << (*t)
              << " index=" << (*t).genpartIndex();
    std::cout << std::endl;
  }
}

void PrimaryVertexAnalyzer4PU::printRecTrks(const View<reco::Track>& recTrks) {
  cout << "printRecTrks" << endl;
  int i = 1;
  for (auto t = recTrks.begin(); t != recTrks.end(); ++t) {
    //    reco::TrackBase::ParameterVector  par = t->parameters();

    cout << endl;
    cout << "Track " << i << " ";
    i++;
    //enum TrackQuality { undefQuality=-1, loose=0, tight=1, highPurity=2, confirmed=3, goodIterative=4, qualitySize=5};
    if (t->quality(reco::TrackBase::loose)) {
      cout << "loose ";
    };
    if (t->quality(reco::TrackBase::tight)) {
      cout << "tight ";
    };
    if (t->quality(reco::TrackBase::highPurity)) {
      cout << "highPurity ";
    };
    if (t->quality(reco::TrackBase::confirmed)) {
      cout << "confirmed  ";
    };
    if (t->quality(reco::TrackBase::goodIterative)) {
      cout << "goodIterative  ";
    };
    cout << endl;

    TransientTrack tk = theB_->build(&(*t));
    tk.setBeamSpot(vertexBeamSpot_);
    double ipsig = 0;
    if (tk.stateAtBeamLine().isValid()) {
      ipsig = tk.stateAtBeamLine().transverseImpactParameter().significance();
    } else {
      ipsig = -1;
    }

    cout << Form("pt=%8.3f phi=%6.3f eta=%6.3f z=%8.4f  dz=%6.4f, ipsig=%6.1f",
                 t->pt(),
                 t->phi(),
                 t->eta(),
                 t->vz(),
                 t->dzError(),
                 ipsig)
         << endl;

    cout << Form(" found=%6d  lost=%6d   chi2/ndof=%10.3f ", t->found(), t->lost(), t->normalizedChi2()) << endl;
    //const reco::HitPattern & p= t->hitPattern();
    //cout << "subdet layers valid lost" << endl;
    //    cout << Form(" barrel  %2d  %2d  %2d",p.pixelBarrelLayersWithMeasurement(),p.numberOfValidPixelBarrelHits(), p.numberOfLostPixelBarrelHits()) << endl;
    //    cout << Form(" fwd     %2d  %2d  %2d",p.pixelEndcapLayersWithMeasurement(),p.numberOfValidPixelEndcapHits(), p.numberOfLostPixelEndcapHits()) << endl;
    //cout << Form(" pixel   %2d  %2d  %2d",p.pixelLayersWithMeasurement(), p.numberOfValidPixelHits(), p.numberOfLostPixelHits()) << endl;
    //cout << Form(" tracker %2d  %2d  %2d",p.trackerLayersWithMeasurement(), p.numberOfValidTrackerHits(), p.numberOfLostTrackerHits()) << endl;
    cout << endl;
    //const reco::HitPattern & pinner= t->trackerExpectedHitsInner();
    //const reco::HitPattern & pouter= t->trackerExpectedHitsOuter();
    //int ninner = 0; //pinner.numberOfAllHits();
    //int nouter = 0; //pouter.numberOfAllHits();

    //
    for (trackingRecHit_iterator hit = t->recHitsBegin(); hit != t->recHitsEnd(); hit++) {
      if ((**hit).isValid() && (**hit).geographicalId().det() == DetId::Tracker) {
        bool barrel = (**hit).geographicalId().subdetId() == static_cast<int>(PixelSubdetector::PixelBarrel);
        bool endcap = (**hit).geographicalId().subdetId() == static_cast<int>(PixelSubdetector::PixelEndcap);
        //       bool barrel = DetId::DetId((**hit).geographicalId()).subdetId() == static_cast<int>(PixelSubdetector::PixelBarrel);reco::HitPattern::numberOfLostPixelHits(reco::HitPattern::HitCategory)
        //       bool endcap = DetId::DetId((**hit).geographicalId()).subdetId() == static_cast<int>(PixelSubdetector::PixelEndcap);
        if (barrel) {
          const SiPixelRecHit* pixhit = dynamic_cast<const SiPixelRecHit*>(&(**hit));
          edm::Ref<edmNew::DetSetVector<SiPixelCluster>, SiPixelCluster> const& clust = (*pixhit).cluster();
          if (clust.isNonnull()) {
            cout << Form(" barrel cluster size=%2d   charge=%6.1f wx=%2d  wy=%2d, expected=%3.1f",
                         clust->size(),
                         1.0 * clust->charge(),
                         clust->sizeX(),
                         clust->sizeY(),
                         1. + 2. / fabs(tan(t->theta())))
                 << endl;
          }
        } else if (endcap) {
          const SiPixelRecHit* pixhit = dynamic_cast<const SiPixelRecHit*>(&(**hit));
          edm::Ref<edmNew::DetSetVector<SiPixelCluster>, SiPixelCluster> const& clust = (*pixhit).cluster();
          if (clust.isNonnull()) {
            cout << Form(" endcap cluster size=%2d   charge=%6.1f wx=%2d  wy=%2d",
                         clust->size(),
                         1.0 * clust->charge(),
                         clust->sizeX(),
                         clust->sizeY())
                 << endl;
          }
        }
      }
    }
  }
}

/********************************************************************************************************/
// helpers for z-sorting
namespace {
  bool recTrackLessZ(const std::pair<reco::TransientTrack, unsigned int>& tk1,
                     const std::pair<reco::TransientTrack, unsigned int>& tk2) {
    if (tk1.first.stateAtBeamLine().isValid() && tk2.first.stateAtBeamLine().isValid()) {
      return tk1.first.stateAtBeamLine().trackStateAtPCA().position().z() <
             tk2.first.stateAtBeamLine().trackStateAtPCA().position().z();
    } else {
      return false;
    }
  }
}  // namespace

namespace {
  bool lt(const std::pair<double, unsigned int>& a, const std::pair<double, unsigned int>& b) {
    return a.first < b.first;
  }
}  // namespace
/********************************************************************************************************/

/********************************************************************************************************/
std::vector<bool> PrimaryVertexAnalyzer4PU::trackClass(const reco::Track& t) {
  const int nclass = 5;
  std::vector<bool> tc(nclass);
  tc[0] = t.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::PixelBarrel, 1);
  tc[1] = t.hitPattern().numberOfAllHits(HitPattern::MISSING_INNER_HITS) <
          0.5 * t.hitPattern().trackerLayersWithMeasurement();
  tc[2] = t.hitPattern().numberOfAllHits(HitPattern::MISSING_INNER_HITS) >
          0.5 * t.hitPattern().trackerLayersWithMeasurement();
  tc[3] = (t.pt() < 0.3) && (std::fabs(t.eta()) > 2);
  tc[4] = (t.pt() > 0.5) && (t.hitPattern().pixelBarrelLayersWithMeasurement() == 4);
  return tc;
}
/********************************************************************************************************/




/********************************************************************************************************/
void PrimaryVertexAnalyzer4PU::analyzeTracksTP(Tracks& tracks, std::vector<SimEvent>& simEvt) {
  /********************************************************************************************************/

  //const double vertexSize = 0.1;
  const double d0CutOff = 3.;

  // now histogram all tracks
  for (unsigned int i = 0; i < tracks.size(); i++) {
    RecoTrack tk = tracks(i);

    double t_pi = 1. / (1. + exp(std::pow(tk.ip / tk.dip, 2) - std::pow(d0CutOff, 2)));

    if (tk.has_timing) {
      Fill(hTrk, "ttrk_rec_all_wide", tk.t);
      Fill(hTrk, "ttrk_rec_all", tk.t);
      Fill(hTrk, "terrtrk_rec_all", tk.dt);
      Fill(hTrk, "tqualtrk_rec_all", tk.timeQuality);
    }

    if (tk.selected && tk.has_timing) {
      Fill(hTrk, "ttrk_rec_sel_wide", tk.t);
      Fill(hTrk, "ttrk_rec_sel", tk.t);
      Fill(hTrk, "terrtrk_rec_sel", tk.dt);
      Fill(hTrk, "tqualtrk_rec_sel", tk.timeQuality);
      if(tk.dt < 0.1)  Fill(hTrk, "tqualtrk_sigmat01_rec_sel", tk.timeQuality);

      if (tk.matched) {
	double t_pid = tk.get_t_pid();
	if (std::abs(t_pid) < 10.){
	  double tres = t_pid - tk.tsim;
	  double tpull = tres / tk.MTD_timeerror;
	  Fill(hTrk, "trestrk_selmatched_pid", tres );
	  Fill(hTrk, "tpulltrk_selmatched_pid", tpull );
	  if(tk.MTD_timeerror < 0.1){
	    Fill(hTrk, "trestrk_sigmatlt01_selmatched_pid", tres );
	    Fill(hTrk, "tpulltrk_sigmatlt01_selmatched_pid", tpull );
	  }
	  if(tk.pt > 1.0){
	    Fill(hTrk, "trestrk_selmatched_pid_ptgt1", tres );
	    Fill(hTrk, "tpulltrk_selmatched_pid_ptgt1", tpull );
	  }
	  if (tk.timeQuality < 0.5){
	    Fill(hTrk, "trestrk_selmatched_pid_qlt05", tres );
	    Fill(hTrk, "tpulltrk_selmatched_pid_qlt05", tpull );
	  }else if(tk.timeQuality < 0.8){
	    Fill(hTrk, "trestrk_selmatched_pid_q0508", tres );
	    Fill(hTrk, "tpulltrk_selmatched_pid_q0508", tpull );
	  }else{
	    Fill(hTrk, "trestrk_selmatched_pid_qgt08", tres );
	    Fill(hTrk, "tpulltrk_selmatched_pid_qgt08", tpull );
	  }
	}

        Fill(hTrk, "ttrk_rec_selmatched_wide", tk.t);
        Fill(hTrk, "ttrk_rec_selmatched", tk.t);
        Fill(hTrk, "terrtrk_rec_selmatched", tk.dt);
        Fill(hTrk, "ttrk_sim_selmatched", tk.tsim);
        Fill(hTrk, "trestrk_selmatched", tk.t - tk.tsim);
        Fill(hTrk, "tpulltrk_selmatched", (tk.t - tk.tsim) / tk.dt);
	Fill(hTrk, "tqualtrk_rec_selmatched", tk.timeQuality);

	if(tk.dt < 0.1)
	  {
	    Fill(hTrk, "tpulltrk_sigmat01_selmatched", (tk.t - tk.tsim) / tk.dt);
	    Fill(hTrk, "tqualtrk_sigmat01_selmatched", tk.timeQuality);
	    if(fabs(tk.eta) <  1.4)
	      {// barrel timing layer
		Fill(hTrk, "trestrk_selmatched_barrel", tk.t - tk.tsim);
		Fill(hTrk, "terrtrk_rec_selmatched_barrel", tk.dt );
		Fill(hTrk, "tpulltrk_selmatched_barrel", (tk.t - tk.tsim) /  tk.dt );
		Fill(hTrk, "trestrkvszrestrk_selmatched_barrel", tk.z - tk.zsim, tk.t - tk.tsim);
		Fill(hTrk, "tqualtrk_rec_selmatched_barrel", tk.timeQuality);
		if(tk.pt > 1.0){
		  Fill(hTrk, "tpulltrk_selmatched_barrel_hipt", (tk.t - tk.tsim) /  tk.dt );
		  Fill(hTrk, "trestrk_selmatched_barrel_hipt", tk.t - tk.tsim);
		  Fill(hTrk, "tqualtrk_rec_selmatched_barrel_hipt", tk.timeQuality);
		}else{
		  Fill(hTrk, "trestrk_selmatched_barrel_lopt", tk.t - tk.tsim);
		  Fill(hTrk, "tpulltrk_selmatched_barrel_lopt", (tk.t - tk.tsim) /  tk.dt );
		  Fill(hTrk, "tqualtrk_rec_selmatched_barrel_lopt", tk.timeQuality);
		}
		
	      }
	    else if (fabs(tk.eta) > 1.7)
	      {// endcap
		Fill(hTrk, "trestrk_selmatched_endcap", tk.t - tk.tsim);
		Fill(hTrk, "terrtrk_rec_selmatched_endcap", tk.dt );
		Fill(hTrk, "tpulltrk_selmatched_endcap", (tk.t - tk.tsim) /  tk.dt );
		Fill(hTrk, "tqualtrk_rec_selmatched_endcap", tk.timeQuality);
		if(tk.pt > 1.0){
		  Fill(hTrk, "trestrk_selmatched_endcap_hipt", tk.t - tk.tsim);
		  Fill(hTrk, "tpulltrk_selmatched_endcap_hipt", (tk.t - tk.tsim) /  tk.dt );
		  Fill(hTrk, "tqualtrk_rec_selmatched_endcap_hipt", tk.timeQuality);
		}else{
		  Fill(hTrk, "trestrk_selmatched_endcap_lopt", tk.t - tk.tsim);
		  Fill(hTrk, "tpulltrk_selmatched_endcap_lopt", (tk.t - tk.tsim) /  tk.dt );
		  Fill(hTrk, "tqualtrk_rec_selmatched_endcap_lopt", tk.timeQuality);
		}
		Fill(hTrk, "trestrkvszrestrk_selmatched_endcap", tk.z - tk.zsim, tk.t - tk.tsim);
		if(tk.eta > 1.7)
		  {
		    Fill(hTrk, "trestrk_selmatched_fwd", tk.t - tk.tsim);
		  }
		else
		  {
		    Fill(hTrk, "trestrk_selmatched_bwd", tk.t - tk.tsim);
		  }
	      }
	  }

	// fill residuals end pulls for three regions of the quality variable
	if (tk.timeQuality > 0.8){
	  Fill(hTrk, "trestrk_selmatched_qgt08", tk.t - tk.tsim);
	  Fill(hTrk, "tpulltrk_selmatched_qgt08", (tk.t - tk.tsim) / tk.dt );
	  if (tk.dt < 0.1){
	    Fill(hTrk, "trestrk_sigmatlt01_selmatched_qgt08", tk.t - tk.tsim);
	    Fill(hTrk, "tpulltrk_sigmatlt01_selmatched_qgt08", (tk.t - tk.tsim) / tk.dt );
	  }
	}else if(tk.timeQuality < 0.5){
	  Fill(hTrk, "trestrk_selmatched_qlt05", tk.t - tk.tsim);
	  Fill(hTrk, "tpulltrk_selmatched_qlt05", (tk.t - tk.tsim) / tk.dt );
	  if (tk.dt < 0.1){
	    Fill(hTrk, "trestrk_sigmatlt01_selmatched_qlt05", tk.t - tk.tsim);
	    Fill(hTrk, "tpulltrk_sigmatlt01_selmatched_qlt05", (tk.t - tk.tsim) / tk.dt );
	  }
	}else if((tk.timeQuality >= 0.5) && (tk.timeQuality <= 0.8)){
	  Fill(hTrk, "trestrk_sigmatlt01_selmatched_q0508", tk.t - tk.tsim);
	  Fill(hTrk, "tpulltrk_sigmatlt01_selmatched_q0508", (tk.t - tk.tsim) / tk.dt );
	  if (tk.dt < 0.1){
	    Fill(hTrk, "trestrk_selmatched_q0508", tk.t - tk.tsim);
	    Fill(hTrk, "tpulltrk_selmatched_q0508", (tk.t - tk.tsim) / tk.dt );
	  }
	}
	  
        double mtd_dt = tk.MTD_timeerror;
	if (tk.is_pion()){
	  Fill(hTrk, "trestrk_selpion", tk.t - tk.tsim);
	  Fill(hTrk, "tpulltrk_selpion", (tk.t - tk.tsim) / tk.dt);
	  Fill(hTrk, "trestrkh_selpion", tk.th[0] - tk.tsim);
	  Fill(hTrk, "tpulltrkh_selpion", (tk.th[0] - tk.tsim) / mtd_dt);
	  Fill(hTrk, "trestrkh_selhyp", tk.th[0] - tk.tsim);
	  Fill(hTrk, "tpulltrkh_selhyp", (tk.th[0] - tk.tsim) / mtd_dt);
	  if(tk.dt < 0.1) Fill(hTrk, "tpulltrk_sigmat01_selpion", (tk.t - tk.tsim) / tk.dt);
	}else if (tk.is_kaon()){
	  Fill(hTrk, "trestrk_selkaon", tk.t - tk.tsim);
	  Fill(hTrk, "tpulltrk_selkaon", (tk.t - tk.tsim) / tk.dt);
	  Fill(hTrk, "trestrkh_selkaon", tk.th[1] - tk.tsim);
	  Fill(hTrk, "tpulltrkh_selkaon", (tk.th[1] - tk.tsim) / mtd_dt);
	  Fill(hTrk, "trestrkh_selhyp", tk.th[1] - tk.tsim);
	  Fill(hTrk, "tpulltrkh_selhyp", (tk.th[1] - tk.tsim) / mtd_dt);
	  if(tk.dt < 0.1) Fill(hTrk, "tpulltrk_sigmat01_selkaon", (tk.t - tk.tsim) / tk.dt);
	}else if (tk.is_proton()){
	  Fill(hTrk, "trestrk_selproton", tk.t - tk.tsim);
	  Fill(hTrk, "tpulltrk_selproton", (tk.t - tk.tsim) / tk.dt);
	  Fill(hTrk, "trestrkh_selproton", tk.th[2] - tk.tsim);
	  Fill(hTrk, "tpulltrkh_selproton", (tk.th[2] - tk.tsim) / mtd_dt);
	  Fill(hTrk, "trestrkh_selhyp", tk.th[2] - tk.tsim);
	  Fill(hTrk, "tpulltrkh_selhyp", (tk.th[2] - tk.tsim) / mtd_dt);
	  if(tk.dt < 0.1) Fill(hTrk, "tpulltrk_sigmat01_selproton", (tk.t - tk.tsim) / tk.dt);
	}	  
	
      } else {
        Fill(hTrk, "ttrk_rec_selunmatched_wide", tk.t);
        Fill(hTrk, "ttrk_rec_selunmatched", tk.t);
      }
    }

    // plots for primary, selected tracks
    if (tk.selected && tk.is_primary) {
      double zpull = (tk.z - tk.zsim) / tk.dz;
      Fill(hTrk, "zpulltrk_primselmatched", zpull);
      Fill(hTrk, "zpulltprimselvseta", tk.eta, zpull);

      double logpt = log(tk.pt) / log(10.);
      unsigned int bin = 4;
      if (tk.dz <= 0.01) {bin =0;}
      else if (tk.dz < 0.02) { bin = 1; }
      else if (tk.dz < 0.05) { bin = 2; }
      else if (tk.dz < 0.10) { bin = 3; }
      else { bin = 4;}
      const char * sbin = trkdzbin_[bin].c_str();
      

      Fill(hTrk, "zpulltprimselvslogpt", logpt, zpull);
      Fill(hTrk, Form("zpulltprimselvslogpt_%s", sbin), logpt, zpull);

      auto hitPattern = tk.tt->track().hitPattern();
      auto nbarrel = hitPattern.pixelBarrelLayersWithMeasurement();
      if ( nbarrel < 2){
	Fill(hTrk, Form("zpulltprimselbpxlt2vseta_%s", sbin), tk.eta, zpull);
      }else if (nbarrel > 2){
	Fill(hTrk, Form("zpulltprimselbpxgt2vseta_%s", sbin), tk.eta, zpull);
      }

      Fill(hTrk, Form("ztailtprimselvslogpt_%s", sbin), logpt, (std::abs(zpull) > 3. ? 1. : 0.));

      if (std::abs(tk.eta) > 1.2){
	Fill(hTrk, Form("ztailtprimselvslogpt_etahi_%s", sbin), logpt, (std::abs(zpull) > 3. ? 1. : 0.));
      }else{
	Fill(hTrk, Form("ztailtprimselvslogpt_etalo_%s", sbin), logpt, (std::abs(zpull) > 3. ? 1. : 0.));
      }

      Fill(hTrk, "tprimselvslogpteta", tk.eta, logpt);
      Fill(hTrk, Form("tprimselvslogpteta_%s", sbin), tk.eta, logpt);
      if(std::abs(zpull) > 3.){
	Fill(hTrk, "ztailtprimselvslogpteta", tk.eta, logpt);
	Fill(hTrk, Form("ztailtprimselvslogpteta_%s",sbin), tk.eta, logpt);
      }


      
      Fill(hTrk, "zrestrk_primselmatched", tk.z - tk.zsim);
      if (tk.has_timing) {
        Fill(hTrk, "ttrk_rec_primselmatched", tk.t);
        Fill(hTrk, "ttrk_sim_primselmatched", tk.tsim);
        Fill(hTrk, "trestrk_primselmatched", tk.t - tk.tsim);
        Fill(hTrk, "tpulltrk_primselmatched", (tk.t - tk.tsim) / tk.dt);
        Fill(hTrk, "zpulltrkt_primselmatched", (tk.z - tk.zsim) / tk.dz);
        Fill(hTrk, "zrestrkt_primselmatched", tk.z - tk.zsim);
      }

      Fill(hTrk, "d0pullprim", tk.ip / tk.dip);
      Fill(hTrk, "zpullvsd0pullprim", std::fabs(tk.ip / tk.dip), std::fabs((tk.z - tk.zsim) / tk.dz));
      if (std::fabs(tk.eta) < 1.5) {
        Fill(hTrk, "zpulltrkt_primselmatched_central", (tk.z - tk.zsim) / tk.dz);
      } else if (std::fabs(tk.eta) < 2.0) {
        Fill(hTrk, "zpulltrkt_primselmatched_inter", (tk.z - tk.zsim) / tk.dz);
      } else {
        Fill(hTrk, "zpulltrkt_primselmatched_fwd", (tk.z - tk.zsim) / tk.dz);
      }
      Fill(hTrk, "tpiprim", t_pi);

    } else {
      Fill(hTrk, "zpullsec", (tk.z - tk.zsim) / tk.dz);
      Fill(hTrk, "d0pullsec", tk.ip / tk.dip);
      Fill(hTrk, "zpullvsd0pullsec", std::fabs(tk.ip / tk.dip), std::fabs((tk.z - tk.zsim) / tk.dz));
      Fill(hTrk, "tpisec", t_pi);
    }

  }  // track loop
}

/********************************************************************************************************/
void PrimaryVertexAnalyzer4PU::printPVTrksZT(edm::Handle<edm::View<Track>> trackCollectionH,
                                             const reco::VertexCollection* recVtxs,
                                             std::vector<SimPart>& tsim,
                                             std::vector<simPrimaryVertex>& simpv,
                                             std::vector<SimEvent>& simEvt) {
  // make a printout of the tracks selected for PV reconstructions, show matching MC tracks, too
  /********************************************************************************************************/

  vector<std::pair<TransientTrack, unsigned int>> selTrks;

  bool show_timing = false;
  unsigned int nTimingTracks = 0;
  unsigned int nSelectedTracks = 0;
  //  for(vector<TransientTrack>::iterator tt = t_tks.begin(); tt != t_tks.end(); tt++){
  for (unsigned int i = 0; i < t_tks_.size(); i++) {
    auto tt = t_tks_.at(i);
    selTrks.push_back(make_pair(tt, i));
    if (tt.dtErrorExt() < 0.3) {
      nTimingTracks++;
    }
    if (theTrackFilter(tt)) {
      nSelectedTracks += 1;
    }
  }
  std::cout << "printPVTracksZT : " << nTimingTracks << " tracks out of " << t_tks_.size() << " have timing info"
            << std::endl;
  show_timing = (nTimingTracks > 0);
  cout << " selected " << nSelectedTracks << endl;

  const View<Track> tC = *(trackCollectionH.product());


  if (selTrks.size() == 0)
    return;  // nothing do print

  stable_sort(selTrks.begin(), selTrks.end(), recTrackLessZ);

  // select tracks like for PV reconstruction and match them to sim tracks
  // rectosim is used when we have no trackingparticles
  /* oh well, how can I fill my own View?? Do I really need to keep two versions of supf? */
  vector<Track> selRecTrks;
  for (unsigned int i = 0; i < selTrks.size(); i++) {
    selRecTrks.push_back(selTrks[i].first.track());
  }
  int* rectosim = NULL;
  if (MC_)
    rectosim = supfv(tsim, selRecTrks);

  // now dump in a format similar to the clusterizer
  cout << endl << "printPVTrksZT      run : event  =   " << run_ << " : " << event_ << endl;
  cout << "----                  z +/- dz";
  if (show_timing) {
    cout << "          t +/- dt    ";
  }
  cout << "    q1bfet-ilo:alg CSQ     ip +/-dip         pt   phi   eta";
  if ((tsim.size() > 0) || (simEvt.size() > 0)) {
    cout << "  type      pdg    zvtx    zpull      match";
  }
  cout << endl;

  cout.precision(4);
  int isel = 0;
  double tz0 = -10000;

  for (unsigned int i = 0; i < selTrks.size(); i++) {
    TransientTrack& tt = selTrks[i].first;
    reco::Track t = tt.track();
    unsigned int i0 = selTrks[i].second;  // original index in the track collection before sorting

    // is this track in the tracklist of a recvtx ?
    int vmatch = -1;
    float wmatch = 0;
    int iv = 0;
    for (reco::VertexCollection::const_iterator v = recVtxs->begin(); v != recVtxs->end(); ++v) {
      for (trackit_t tv = v->tracks_begin(); tv != v->tracks_end(); tv++) {
        const reco::Track& RTv = *(tv->get());
        if (t.vz() == RTv.vz()) {
          vmatch = iv;
          wmatch = v->trackWeight(*tv);
        }
      }
      iv++;
    }

    double tz = (tt.stateAtBeamLine().trackStateAtPCA()).position().z();
    double tantheta = tan((tt.stateAtBeamLine().trackStateAtPCA()).momentum().theta());
    double tdz0 = t.dzError();
    double phi = (tt.stateAtBeamLine().trackStateAtPCA()).momentum().phi();
    double tdz2 = pow(t.dzError(), 2) + (pow(wx_ * cos(phi), 2) + pow(wy_ * sin(phi), 2)) / pow(tantheta, 2);

    // print vertices in between tracks
    int iiv = 0;
    for (reco::VertexCollection::const_iterator v = recVtxs->begin(); v != recVtxs->end(); ++v) {
      if ((v->ndof() > 0) && (v->position().z() < tz) && (v->z() > tz0)) {
        cout << "rec [" << setw(4) << iiv << "]     " << setw(8) << fixed << setprecision(4) << v->z() << " +/-"
             << setw(6) << v->zError();
        if (show_timing) {
          if (v->t() == 0) {
            cout << "        +/-     ";
          } else {
            cout << setw(7) << fixed << setprecision(3) << v->t() << " +/-" << setw(5) << v->tError();
          }
        }
        cout << "  ndof=" << v->ndof() << "  ptmax2=" << vertex_ptmax2(*v);
        cout << "  sumpt2=" << vertex_sumpt2(*v);
        if ((v->ndof() > 4) && (vertex_r(*v) > 0.0200)) {
          cout << "    rvtx= " << setw(8) << setprecision(4) << vertex_r(*v);
        }
        cout << " " << endl;
      }
      iiv++;
    }

    // print MC vertices in between tracks
    if (simEvt.size() > 0) {
      // use simEvt if available
      for (unsigned int event = 0; event < simEvt.size(); event++) {
        if ((simEvt[event].z < tz) && (simEvt[event].z > tz0)) {
          if (simEvt[event].type == 1) {
            cout << "sim (" << setw(4) << simEvt[event].eventId.event() << ")     " << setw(8) << fixed
                 << setprecision(4) << simEvt[event].z << "          ";
            if (show_timing) {
              cout << setw(7) << fixed << setprecision(3) << simEvt[event].t << "          ";
            }
            cout << "          "
                 << "  ntrk=" << simEvt[event].rtk.size() << " " << endl;
          } else {
            cout << "sim (PUx)     " << setw(7) << fixed << setprecision(3) << simEvt[event].z << "          "
                 << " " << endl;
          }
        }
      }

    } else if (simpv.size() > 0) {
      // otherwise use simpvs
      for (auto pv = simpv.begin(); pv != simpv.end(); pv++) {
        if ((pv->z < tz) && (pv->z > tz0)) {
          if (pv->type == 1) {
            cout << "sim (MC)     " << setw(8) << fixed << setprecision(3) << pv->z;
            if (show_timing) {
              cout << setw(7) << fixed << setprecision(3) << pv->t;
            }
            cout << "          " << setw(4) << setprecision(1) << sqrt(pv->ptsq) << endl;
          } else {
            cout << "sim (PU)     " << setw(8) << fixed << setprecision(4) << pv->z;
            if (show_timing) {
              cout << setw(7) << fixed << setprecision(3) << pv->t;
            }
            cout << "          "
                 << "pthat = " << setw(4) << setprecision(1) << pv->pt_hat << endl;
          }
        }
      }
    }

    tz0 = tz;

    // for selected tracks print a selected track index
    if (theTrackFilter(tt)) {
      cout << setw(4) << isel;
      isel++;
    } else if (!tt.stateAtBeamLine().isValid()) {
      cout << "XXXX";
    } else {
      cout << "    ";
    }

    if (vmatch > -1) {
      cout << "[" << setw(4) << vmatch << "]" << setw(4) << setprecision(2) << fixed << wmatch << " ";
    } else {
      cout << "           ";
    }

    if (fabs(tz) < 100) {
      cout << setw(8) << fixed << setprecision(4) << tz << " +/-" << setw(6) << sqrt(tdz2);
    } else {
      cout << setw(8) << fixed << setprecision(4) << 99.99 * tz / fabs(tz) << " +/-" << setw(6) << sqrt(tdz2);
    }

    if (show_timing) {
      if (tt.dtErrorExt() < 0.3) {
        cout << setw(7) << fixed << setprecision(3) << tt.timeExt() << " +/-" << setw(5) << tt.dtErrorExt();
      } else {
        cout << "                ";
      }
    }

    // track quality and hit information, see DataFormats/TrackReco/interface/HitPattern.h
    if (t.quality(reco::TrackBase::highPurity)) {
      cout << " *";
    } else {
      cout << "  ";
    }
    if (t.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::PixelBarrel, 1)) {
      cout << "+";
    } else {
      cout << "-";
    }
    // pixel layers with hits
    cout << setw(1) << t.hitPattern().pixelBarrelLayersWithMeasurement();
    cout << setw(1) << t.hitPattern().pixelEndcapLayersWithMeasurement();
    auto tepx = t.hitPattern().pixelLayersWithMeasurement() -
                (t.hitPattern().pixelEndcapLayersWithMeasurement() +
                 t.hitPattern().pixelBarrelLayersWithMeasurement());  // just guessin'
    cout << setw(1) << tepx;
    // outer tracker layers with hits
    int mm = t.hitPattern().trackerLayersWithMeasurement() - t.hitPattern().pixelLayersWithMeasurement();
    if (mm >= 0) {
      cout << setw(1) << hex << mm << dec;
    } else {
      cout << "X";
    }

    // - missing hits  -[inner][outer],  2020-08-12 changed from numberOfAllHits to numberOfLostHits, seer PR #20938
    cout << "-" 
	 << setw(1) << hex  << t.hitPattern().numberOfLostHits(HitPattern::MISSING_INNER_HITS)
         << setw(1) << hex <<  t.lost()
         << setw(1) << hex << t.hitPattern().numberOfLostHits(HitPattern::MISSING_OUTER_HITS) << dec;

    if( t.isLooper()){
      cout << "l";
    }else{
      cout << " ";
    }

    cout << ":" << setw(2) << dec << t.algo();

    // longest barrel hit (if applicable)
    if (RECO_) {
      int longesthit = 0, nbarrel = 0;
      for (trackingRecHit_iterator hit = t.recHitsBegin(); hit != t.recHitsEnd(); hit++) {
        if ((**hit).isValid() && (**hit).geographicalId().det() == DetId::Tracker) {
          bool barrel = (**hit).geographicalId().subdetId() == static_cast<int>(PixelSubdetector::PixelBarrel);
          if (barrel) {
            const SiPixelRecHit* pixhit = dynamic_cast<const SiPixelRecHit*>(&(**hit));
            edm::Ref<edmNew::DetSetVector<SiPixelCluster>, SiPixelCluster> const& clust = (*pixhit).cluster();
            if (clust.isNonnull()) {
              nbarrel++;
              if (clust->sizeY() > longesthit)
                longesthit = clust->sizeY();
            }
          }
        }
      }
      if (nbarrel > 0) {
        cout << setw(3) << setprecision(0) << fixed << longesthit;
      } else {
        cout << "   ";
      }
    } else {  // no RECO
      cout << "   ";
    }

    // track fit chi**2 / ndf
    cout << setw(4) << setprecision(1) << fixed << t.normalizedChi2();

    // xy impact paramter and error
    Measurement1D IP = tt.stateAtBeamLine().transverseImpactParameter();
    cout << setw(8) << setprecision(4) << IP.value() << "+/-" << setw(6) << IP.error();

    // pt
    double relative_pterror = t.ptError() / t.pt();
    if (relative_pterror < 0.1) {
      cout << " " << setw(7) << setprecision(2) << t.pt() * t.charge();
    } else if (relative_pterror < 0.5) {
      cout << " " << setw(6) << setprecision(1) << t.pt() * t.charge() << "-";
    } else {
      cout << " " << setw(6) << setprecision(1) << t.pt() * t.charge() << "*";
    }
    // phi and eta
    cout << " " << setw(5) << setprecision(2) << t.phi() << " " << setw(5) << setprecision(2) << t.eta();

    // print MC info, if available
    if (MC_) {
      if ((simEvt.size() > 0) && (tracking_truth_available_)){
        RefToBase<Track> trtb(trackCollectionH, i0);
        TrackingParticleRef tpr;
        bool matched = truthMatchedTrack(trtb, tpr);

        if (!matched) {
          cout << " not matched              ";
          printTruthMatchValues(trtb);

        } else if (tpr.isNull()) {
          cout << " invalid TrackingParticleRef ??";
          //FIXME, this looks like a bug
          printTruthMatchValues(trtb);

        } else {
          if (tpr->parentVertex().isNull()) {
            cout << " null parent vertex ";

          } else {
            const TrackingVertex* parentVertex = tpr->parentVertex().get();

            double ez1 = 0;
            for (unsigned int j = 0; j < simpv.size(); j++) {
              if (simpv[j].eventId == tpr->eventId()) {
                ez1 = simpv[j].z;
                break;
              }
            }
            if (parentVertex->sourceTracks().size() == 0) {
              if (fabs(parentVertex->position().z() - ez1) < 1e-4) {
                cout << " prim";
              } else {
                cout << " ?   ";  // I don't know what these TrackingParticles are
                // their parentvertex has no source tracks, it is displaced significantly in z from
                // the event vertex, still the tracks point well to the beam spot
              }
            } else {
              cout << " sec ";
            }

            cout << "(" << setw(3) << tpr->eventId().event();
            cout << ")" << setw(5) << tpr->pdgId();
            double vz = parentVertex->position().z();
            cout << " " << setw(9) << setprecision(4) << vz;
            double pz = (tz - vz) / sqrt(tdz2);
            cout << " " << setw(5) << setprecision(1) << pz;
            cout << " ";
            printTruthMatchValues(trtb);
            // mark particles matched to multiple tracks (loopers/double tracks)
            int nl = 0;
            for (View<Track>::size_type j = 0; j < tC.size(); ++j) {
              if (trkidx2tp_[j] == tpr) {
                if (nl > 0)
                  cout << "*";
                nl++;
              }
            }
          }
          /*
	      }catch (...){
	      cout << " not matched1";
	      }
	    */
        }  //
      } else {
        // no tracking particles
        if (rectosim[i] >= 0) {
          if (tsim[rectosim[i]].type == 0) {
            cout << " prim ";
          } else {
            cout << " sec  ";
          }
          cout << " " << setw(5) << tsim[rectosim[i]].pdg;
          cout << " " << setw(8) << setprecision(4) << tsim[rectosim[i]].zvtx;
          cout << " " << setw(8) << setprecision(4) << tsim[rectosim[i]].zdcap;
          cout << " " << setw(8) << setprecision(4) << tsim[rectosim[i]].ddcap;
          double zvtxpull = (tz - tsim[rectosim[i]].zvtx) / sqrt(tdz2);
          cout << setw(6) << setprecision(1) << zvtxpull;
          double zdcapull = (tz - tsim[rectosim[i]].zdcap) / tdz0;
          cout << setw(6) << setprecision(1) << zdcapull;
          double dszpull = (t.dsz() - tsim[rectosim[i]].par[4]) / t.dszError();
          cout << setw(6) << setprecision(1) << dszpull;
        }
      }
    }
    cout << endl;
  }
  cout << "printPVTrks.end ---------------- " << endl;
  if (MC_)
    delete[] rectosim;
}

/********************************************************************************************************/

void PrimaryVertexAnalyzer4PU::getTc(const std::vector<reco::TransientTrack>& tracks,
                                     double& Tc,
                                     double& chsq,
                                     double& dzmax,
                                     double& dztrim,
                                     double& m4m2) {
  if (tracks.size() < 2) {
    Tc = -1;
    chsq = -1;
    dzmax = -1;
    dztrim = -1;
    m4m2 = -1;
    return;
  }

  double sumw = 0, sumwz = 0, sumww = 0, sumwwz = 0, sumwwzz = 0;
  double zmin = 1e10, zmin1 = 1e10, zmax1 = -1e10, zmax = -1e10;
  double m4 = 0, m3 = 0, m2 = 0, m1 = 0, m0 = 0;
  for (vector<reco::TransientTrack>::const_iterator it = tracks.begin(); it != tracks.end(); it++) {
    double tantheta = tan(((*it).stateAtBeamLine().trackStateAtPCA()).momentum().theta());
    reco::BeamSpot beamspot = (it->stateAtBeamLine()).beamSpot();
    double z = ((*it).stateAtBeamLine().trackStateAtPCA()).position().z();
    double dz2 = pow((*it).track().dzError(), 2) + pow(beamspot.BeamWidthX() / tantheta, 2);
    double w = 1. / dz2;  // take p=1
    sumw += w;
    sumwz += w * z;
    sumwwz += w * w * z;
    ;
    sumwwzz += w * w * z * z;
    sumww += w * w;
    m0 += w;
    m1 += w * z;
    m2 += w * z * z;
    m3 += w * z * z * z;
    m4 += w * z * z * z * z;
    if (dz2 < pow(0.1, 2)) {
      if (z < zmin1) {
        zmin1 = z;
      }
      if (z < zmin) {
        zmin1 = zmin;
        zmin = z;
      }
      if (z > zmax1) {
        zmax1 = z;
      }
      if (z > zmax) {
        zmax1 = zmax;
        zmax = z;
      }
    }
  }
  double z = sumwz / sumw;
  double a = sumwwzz - 2 * z * sumwwz + z * z * sumww;
  double b = sumw;
  if (tracks.size() > 1) {
    chsq = (m2 - m0 * z * z) / (tracks.size() - 1);
    Tc = 2. * a / b;
    m4m2 = sqrt((m4 - 4 * m3 * z + 6 * m2 * z * z - 3 * m1 * z * z * z + m0 * z * z * z * z) /
                (m2 - 2 * m1 * z + z * z * m0));
  } else {
    chsq = 0;
    Tc = 0;
    m4m2 = 0;
  }
  dzmax = zmax - zmin;
  dztrim = zmax1 - zmin1;  // truncated
}
/********************************************************************************************************/
void PrimaryVertexAnalyzer4PU::getTc(
    const std::vector<RecoTrack>& tracks, double& Tc, double& chsq, double& dzmax, double& dztrim, double& m4m2) {
  if (tracks.size() < 2) {
    Tc = -1;
    chsq = -1;
    dzmax = -1;
    dztrim = -1;
    m4m2 = -1;
    return;
  }

  double sumw = 0, sumwz = 0, sumww = 0, sumwwz = 0, sumwwzz = 0;
  double zmin = 1e10, zmin1 = 1e10, zmax1 = -1e10, zmax = -1e10;
  double m4 = 0, m3 = 0, m2 = 0, m1 = 0, m0 = 0;

  for (auto t = tracks.begin(); t != tracks.end(); t++) {
    double tantheta = tan(t->theta);
    double z = t->z;
    double dz2 = pow(t->dz, 2) + pow(vertexBeamSpot_.BeamWidthX() / tantheta, 2);
    double w = 1. / dz2;  // take p=1
    sumw += w;
    sumwz += w * z;
    sumwwz += w * w * z;
    ;
    sumwwzz += w * w * z * z;
    sumww += w * w;
    m0 += w;
    m1 += w * z;
    m2 += w * z * z;
    m3 += w * z * z * z;
    m4 += w * z * z * z * z;
    if (dz2 < pow(0.1, 2)) {
      if (z < zmin1) {
        zmin1 = z;
      }
      if (z < zmin) {
        zmin1 = zmin;
        zmin = z;
      }
      if (z > zmax1) {
        zmax1 = z;
      }
      if (z > zmax) {
        zmax1 = zmax;
        zmax = z;
      }
    }
  }
  double z = sumwz / sumw;
  double a = sumwwzz - 2 * z * sumwwz + z * z * sumww;
  double b = sumw;
  if (tracks.size() > 1) {
    chsq = (m2 - m0 * z * z) / (tracks.size() - 1);
    Tc = 2. * a / b;
    m4m2 = sqrt((m4 - 4 * m3 * z + 6 * m2 * z * z - 3 * m1 * z * z * z + m0 * z * z * z * z) /
                (m2 - 2 * m1 * z + z * z * m0));
  } else {
    chsq = 0;
    Tc = 0;
    m4m2 = 0;
  }
  dzmax = zmax - zmin;
  dztrim = zmax1 - zmin1;  // truncated
}
/********************************************************************************************************/

/********************************************************************************************************/
bool PrimaryVertexAnalyzer4PU::truthMatchedTrack(edm::RefToBase<reco::Track> track, TrackingParticleRef& tpr)

/********************************************************************************************************/
// for a reco track select the matching tracking particle, always use this function to make sure we
// are consistent
// after calling truthMatchedTrack, tpr may have changed its value
// to get the TrackingParticle from the TrackingParticleRef, use ->get();
{
  if (r2s_->find(track) == r2s_->end()) {
    return false;
  } else {
    double f = -1e10;
    TrackingParticleRef tf;
    std::vector<std::pair<TrackingParticleRef, double>> tp = (*r2s_)[track];
    for (auto it = tp.begin(); it != tp.end(); ++it) {
      if (it->second > f) {
        tf = it->first;
        f = it->second;
      }
    }
    if (f > trackAssociatorMin_) {
      tpr = tf;
      return true;
    }
  }
  return false;
}
/********************************************************************************************************/

/********************************************************************************************************/
void PrimaryVertexAnalyzer4PU::printTruthMatchValues(edm::RefToBase<reco::Track> track)

/********************************************************************************************************/
// print the tracking truth assocation values for a reco track
{
  if (r2s_->find(track) == r2s_->end()) {
    return;
  }

  cout << "[" << setw(4) << setprecision(2) << fixed;
  std::vector<std::pair<TrackingParticleRef, double>> tp = (*r2s_)[track];
  unsigned int n = 0;
  for (auto it = tp.begin(); it != tp.end(); ++it) {
    if (n > 0)
      cout << ",";
    cout << it->second;
    if (++n > 5)
      break;
  }
  cout << "]";

  return;
}
/********************************************************************************************************/

/********************************************************************************************************/
std::vector<edm::RefToBase<reco::Track>> PrimaryVertexAnalyzer4PU::getTruthMatchedVertexTracks(const reco::Vertex& v,
                                                                                               double min_weight)
// for rec vertex v get a list of tracks for which truth matching is available
/********************************************************************************************************/
{
  std::vector<edm::RefToBase<reco::Track>> b;
  TrackingParticleRef tpr;

  for (trackit_t tv = v.tracks_begin(); tv != v.tracks_end(); tv++) {
    if (v.trackWeight(*tv) >= min_weight) {
      if (truthMatchedTrack(*tv, tpr)) {
        b.push_back(*tv);
      }
    }
  }
  return b;
}
/********************************************************************************************************/


/********************************************************************************************************/
std::vector<PrimaryVertexAnalyzer4PU::SimEvent> PrimaryVertexAnalyzer4PU::getSimEvents_notp
(
 const Handle<SimTrackContainer> simTrks,
 const Handle<SimVertexContainer> simVtxs,
 Tracks& tracks
 )
{
  report_counted("PrimaryVertexAnalyzer4PU::getSimEvents without tracking particles",1);
  SimVertexContainer::const_iterator vsim = simVtxs->begin(); 
  vector<SimEvent> simEvt;
  SimEvent e(0);
  e.type = 1; // partial
  e.eventId = vsim->eventId();
  e.nChTP = 0;
  e.ptvis = 0;
  e.sumpt = 0;
  e.sumpt2 = 0;
  e.pxvis = 0;
  e.pyvis = 0;
  e.sumpt2 = 0;
  e.x = vsim->position().x() * simUnit_;
  e.y = vsim->position().y() * simUnit_;
  e.z = vsim->position().z() * simUnit_;
  e.t = vsim->position().t() * simtUnit_;
  
  simEvt.push_back(e);
  
  for (unsigned int i = 0; i < tracks.size(); i++) {
    RecoTrack& t = tracks(i);
    simEvt[0].trkidx.push_back(i);
    if (t.selected) {
      simEvt[0].rtk.push_back(t);
    }
    
    // fill the matching info for the RecoTrack
    t.matched = true;
    t.simEvt = &simEvt[0];
    t.zsim = simEvt[0].z;
    if (f4D_) {
      t.tsim = simEvt[0].t;
    }
    //if (ipdist < 10)  // originated less than 10 um from the interaction point
    //  {
        t.is_primary = true;
    //  }
  }
  
  return simEvt;
}
/********************************************************************************************************/



/********************************************************************************************************/
std::vector<PrimaryVertexAnalyzer4PU::SimEvent> PrimaryVertexAnalyzer4PU::getSimEvents(
    edm::Handle<TrackingParticleCollection> TPCollectionH,
    Tracks& tracks)
// collect information about simulated collisions from tracking particles / vertices
/********************************************************************************************************/
{

  const TrackingParticleCollection* simTracks = TPCollectionH.product();

  vector<SimEvent> simEvt;
  map<EncodedEventId, unsigned int> eventIdToEventMap;
  map<EncodedEventId, unsigned int>::iterator id;

  bool dumpTP = false;
  for (TrackingParticleCollection::const_iterator it = simTracks->begin(); it != simTracks->end(); it++) {
    if (fabs(it->parentVertex().get()->position().z()) > 100.)
      continue;  // skip funny entries @ -13900

    unsigned int event = 0;  //note, this is no longer the same as eventId().event()
    id = eventIdToEventMap.find(it->eventId());
    // skip out of time pile-up, irrelevant for tracking
    if (it->eventId().bunchCrossing() != 0)
      continue;
    //
    if (id == eventIdToEventMap.end()) {
      // new event here
      event = simEvt.size();
      SimEvent e(event);
      e.type = 1;  //full
      e.eventId = it->eventId();
      e.nChTP = 0;
      e.sumpt = 0;
      e.sumpt2 = 0;
      const TrackingVertex* parentVertex = it->parentVertex().get();
      if (DEBUG_) {
        cout << "getSimEvents: ";
        cout << it->eventId().bunchCrossing() << "," << it->eventId().event() << " z=" << it->vz() << " "
             << parentVertex->eventId().bunchCrossing() << "," << parentVertex->eventId().event()
             << " z=" << parentVertex->position().z() << endl;
      }
      if (it->eventId() == parentVertex->eventId()) {
        e.x = parentVertex->position().x();
        e.y = parentVertex->position().y();
        e.z = parentVertex->position().z();
        e.t = parentVertex->position().t() * simtUnit_;
      } else {
        e.x = it->vx();
        e.y = it->vy();
        e.z = it->vz();
        e.t = 0;  // FIXME timing
      }
      simEvt.push_back(e);
      eventIdToEventMap[e.eventId] = event;
    } else {
      event = id->second;
    }

    simEvt[event].tp.push_back(&(*it));
    if ((fabs(it->eta()) < etaMaxVisible_) && (it->charge() != 0) && (it->numberOfTrackerHits() > numTrkHitsVisible_)) {
      if (it->pt() > ptMinVisible_) {
        simEvt[event].nChTP++;
	simEvt[event].ptvis += it->pt();
	simEvt[event].pxvis += it->pt() * cos(it->phi());
	simEvt[event].pyvis += it->pt() * sin(it->phi());
      }

      simEvt[event].sumpt2 += pow(it->pt(), 2);  // should keep track of decays ?
      simEvt[event].sumpt += it->pt();
    }
  }

  if (dumpTP) {
    for (auto it = simTracks->begin(); it != simTracks->end(); it++) {
      std::cout << *it << std::endl;
    }
  }

  // collect reco tracks  and truth matched tracks for each simevent

  for (unsigned int i = 0; i < tracks.size(); i++) {
    RecoTrack& tk = tracks(i);

    if (truthMatchedTrack(tracks.ref(i), tk.tpr)) {
      if (eventIdToEventMap.find(tk.tpr->eventId()) == eventIdToEventMap.end()) {
        if (tk.tpr->eventId().bunchCrossing() == 0) {
          cout << "Bug in getSimEvents, missing event ? " << tk.tpr->eventId().bunchCrossing() << ","
               << tk.tpr->eventId().event() << endl;
        } else if (DEBUG_) {
          cout << "track without SimEvent " << tk.tpr->eventId().bunchCrossing() << "," << tk.tpr->eventId().event()
               << endl;
        }
        //break;
        continue;
      }

      unsigned int event = eventIdToEventMap[tk.tpr->eventId()];
      simEvt[event].trkidx.push_back(i);
      const TrackingVertex* parentVertex = tk.tpr->parentVertex().get();
      double vx = parentVertex->position().x();  // problems with tpr->vz()
      double vy = parentVertex->position().y();
      double vz = parentVertex->position().z();
      double ipdist = sqrt(pow(simEvt[event].x - vx, 2) + pow(simEvt[event].y - vy, 2) + pow(simEvt[event].z - vz, 2)) *
                      1.e4;  // in um

      // fill the matching info for the RecoTrack
      tk.matched = true;
      tk.simEvt = &simEvt[event];
      tk.zsim = simEvt[event].z;
      if (f4D_) {
        tk.tsim = simEvt[event].t;
      }
      if (ipdist < 10)  // originated less than 10 um from the interaction point
      {
        tk.is_primary = true;
      }

      if (tk.selected) {
	/* FIXME obsolete
        simEvt[event].trkref.push_back(tracks.ref(tk.index));  // deprecated
        simEvt[event].tk.push_back(*tk.tt);                    // deprecated
        if (ipdist < 5)
          simEvt[event].tkprim.push_back(*tk.tt);  // deprecated
        if ((tk.ip / tk.dip) < 5)
          simEvt[event].tkprimsel.push_back(*tk.tt);  // deprecated
	*/
	
        simEvt[event].rtk.push_back(tk);
        if (tk.is_primary)
          simEvt[event].rtkprim.push_back(tk);
        if ((tk.ip / tk.dip) < 4) // should this cut be configurable? or at least a global constant?
          simEvt[event].rtkprimsel.push_back(tk);
      }
      
    } else {
      // track not truth matched
      tk.matched = false;
      tk.simEvt = NULL;
      tk.zsim = 0;
      tk.tsim = 0;
      tk.is_primary = false;
    }



  }  // end of track loop

  for (unsigned int i = 0; i < simEvt.size(); i++) {
    if (simEvt[i].rtkprim.size() > 0) {
      getTc(simEvt[i].rtkprimsel, simEvt[i].Tc, simEvt[i].chisq, simEvt[i].dzmax, simEvt[i].dztrim, simEvt[i].m4m2);
      simEvt[i].zfit = -99;
    } else {
      simEvt[i].Tc = 0;
      simEvt[i].chisq = 0;
      simEvt[i].dzmax = 0;
      simEvt[i].dztrim = 0;
      simEvt[i].m4m2 = 0;
      simEvt[i].zfit = -99;
      simEvt[i].tfit = -99;
    }

    if (DEBUG_) {
      cout << setw(3) << i << " )   nTP=" << setw(4) << simEvt[i].tp.size() << "   z=" << setw(8) << setprecision(4)
           << fixed << simEvt[i].z << "    recTrks=" << setw(3) << simEvt[i].rtk.size() << "    recTrksPrim=" << setw(3)
           << simEvt[i].rtkprim.size() << "    allTrks=" << setw(3) << simEvt[i].trkidx.size() << endl;
    }
  }

  return simEvt;
}

/********************************************************************************************************/





/********************************************************************************************************/
std::vector<PrimaryVertexAnalyzer4PU::simPrimaryVertex> PrimaryVertexAnalyzer4PU::getSimPVs(
    const Handle<SimVertexContainer> simVtxs, const Handle<SimTrackContainer> simTrks)
/********************************************************************************************************/
{
  if (verbose_) {
    std::cout << "getSimPVs from simVtxs/simTrks " << std::endl;
  }

  std::vector<PrimaryVertexAnalyzer4PU::simPrimaryVertex> simpv;
  SimVertexContainer::const_iterator vsim = simVtxs->begin();
  {  // take only the first one as a primary
    simPrimaryVertex sv(vsim->position().x() * simUnit_,
                        vsim->position().y() * simUnit_,
                        vsim->position().z() * simUnit_,
                        vsim->position().t() * simtUnit_);
    sv.type = 1;
    for (edm::SimTrackContainer::const_iterator t = simTrks->begin(); t != simTrks->end(); ++t) {
      int pdgCode = std::abs(t->type());
      if ((pdgCode == 11) || (pdgCode == 13) || (pdgCode == 15) || (pdgCode == -211) || (pdgCode == -2212) ||
          (pdgCode == -321) || (pdgCode == -3222)) {
        if ((t->momentum().Pt() > 0.1) && (fabs(t->momentum().Eta()) < 2.5)) {
          sv.nGenTrk++;
        }
        // count visible tracks inside acceptance
        if ((t->momentum().Pt() > 0.2) && (fabs(t->momentum().Eta()) < 2.4)) {
          sv.nTrk++;
          double xstart = simVtxs->at(t->vertIndex()).position().x() * simUnit_ - sv.x;
          double ystart = simVtxs->at(t->vertIndex()).position().y() * simUnit_ - sv.y;
          if (xstart * xstart + ystart * ystart < 20e-4 * 20e-4) {
            sv.nTrkPrim++;
          } else {
            sv.nTrkSec++;
          }
        }
      }
    }
    simpv.push_back(sv);
  }

  return simpv;
}
/********************************************************************************************************/

/********************************************************************************************************/
std::vector<PrimaryVertexAnalyzer4PU::simPrimaryVertex> PrimaryVertexAnalyzer4PU::getSimPVs(
    const Handle<reco::GenParticleCollection> genParticles)
/********************************************************************************************************/
{
  if (verbose_) {
    std::cout << "getSimPVs from genParticles " << std::endl;
  }

  std::vector<PrimaryVertexAnalyzer4PU::simPrimaryVertex> simpv;
  double x = 0, y = 0, z = 0, t = 0;
  for (size_t i = 0; i < genParticles->size(); ++i) {
    const GenParticle& p = (*genParticles)[i];
    int st = p.status();
    if (st == 1) {
      // FIXME : how dow I get the time information from a genparticle ?
      x = p.vx();
      y = p.vy();
      z = p.vz();
      break;
    }
  }
  simPrimaryVertex sv(x, y, z, t);
  sv.type = 1;
  for (size_t i = 0; i < genParticles->size(); ++i) {
    const GenParticle& p = (*genParticles)[i];
    int pdgCode = std::abs(p.pdgId());
    int st = p.status();
    if ((st == 1) && ((pdgCode == 11) || (pdgCode == 13) || (pdgCode == 15) || (pdgCode == 211) || (pdgCode == 2212) ||
                      (pdgCode == 321) || (pdgCode == 3222))) {
      if ((p.pt() > 0.2) && (fabs(p.eta()) < 2.5)) {
        sv.nGenTrk++;
      }
    }
  }
  simpv.push_back(sv);
  return simpv;
}
/********************************************************************************************************/

/********************************************************************************************************/
std::vector<PrimaryVertexAnalyzer4PU::simPrimaryVertex> PrimaryVertexAnalyzer4PU::getSimPVs(
    const Handle<HepMCProduct> evtMC)
/********************************************************************************************************/
{
  if (verbose_) {
    std::cout << "getSimPVs (HepMCProduct) " << std::endl;
  }

  std::vector<PrimaryVertexAnalyzer4PU::simPrimaryVertex> simpv;
  const HepMC::GenEvent* evt = evtMC->GetEvent();
  if (evt) {
    for (HepMC::GenEvent::vertex_const_iterator vitr = evt->vertices_begin(); vitr != evt->vertices_end();
         ++vitr) {  // loop for vertex ...

      HepMC::FourVector pos = (*vitr)->position();

      if (fabs(pos.z()) > 1000)
        continue;  // skip funny junk vertices

      bool hasMotherVertex = false;
      //std::cout << "mothers" << std::endl;
      for (HepMC::GenVertex::particle_iterator mother = (*vitr)->particles_begin(HepMC::parents);
           mother != (*vitr)->particles_end(HepMC::parents);
           ++mother) {
        HepMC::GenVertex* mv = (*mother)->production_vertex();
        if (mv) {
          hasMotherVertex = true;
        }
        //std::cout << "\t"; (*mother)->print();
      }

      if (hasMotherVertex) {
        continue;
      }

      // could be a new vertex, check  all primaries found so far to avoid multiple entries
      const double mm = 0.1;
      const double s = 1e9;
      simPrimaryVertex sv(pos.x() * mm, pos.y() * mm, pos.z() * mm, pos.t() * s);
      sv.type = 1;
      simPrimaryVertex* vp = NULL;  // will become non-NULL if a vertex is found and then point to it
      for (std::vector<simPrimaryVertex>::iterator v0 = simpv.begin(); v0 != simpv.end(); v0++) {
        if ((fabs(sv.x - v0->x) < 1e-5) && (fabs(sv.y - v0->y) < 1e-5) && (fabs(sv.z - v0->z) < 1e-5)) {
          vp = &(*v0);
          break;
        }
      }
      if (!vp) {
        // this is a new vertex
        simpv.push_back(sv);
        vp = &simpv.back();
      }

      // store the gen vertex barcode with this simpv
      vp->genVertex.push_back((*vitr)->barcode());

      // collect final state descendants and sum up momenta etc
      for (HepMC::GenVertex::particle_iterator daughter = (*vitr)->particles_begin(HepMC::descendants);
           daughter != (*vitr)->particles_end(HepMC::descendants);
           ++daughter) {
        //std::cout << "checking daughter  type " << (*daughter)->pdg_id() << " final :" <<isFinalstateParticle(*daughter) << std::endl;
        if (isFinalstateParticle(*daughter)) {
          if (find(vp->finalstateParticles.begin(), vp->finalstateParticles.end(), (*daughter)->barcode()) ==
              vp->finalstateParticles.end()) {
            vp->finalstateParticles.push_back((*daughter)->barcode());
            HepMC::FourVector m = (*daughter)->momentum();
            vp->ptot.setPx(vp->ptot.px() + m.px());
            vp->ptot.setPy(vp->ptot.py() + m.py());
            vp->ptot.setPz(vp->ptot.pz() + m.pz());
            vp->ptot.setE(vp->ptot.e() + m.e());
            vp->ptsq += (m.perp()) * (m.perp());
            vp->pt_hat = sqrt(vp->ptsq);

            // count relevant particles
            if ((m.perp() > 0.1) && (fabs(m.pseudoRapidity()) < 2.5) && isCharged(*daughter)) {
              vp->nGenTrk++;
            }

            // count visible tracks inside acceptance
            if (isCharged(*daughter) && !isResonance(*daughter) && (m.perp() > 0.2) &&
                (fabs(m.pseudoRapidity()) < 2.4)) {
              vp->nTrk++;
              double xstart = (*daughter)->production_vertex()->point3d().x() * mm - sv.x;
              double ystart = (*daughter)->production_vertex()->point3d().y() * mm - sv.y;
              if (xstart * xstart + ystart * ystart < 20e-4 * 20e-4) {
                vp->nTrkPrim++;
              } else {
                vp->nTrkSec++;
              }
            }

          }  //new final state particle for this vertex
        }    //daughter is a final state particle
      }
    }
  } else {
    if (DEBUG_)
      cout << "PrimaryVertexAnalyzer4PU::getSimPVs   no evt" << endl;
  }

  if (DEBUG_ || verbose_) {
    cout << "------- PrimaryVertexAnalyzer4PU simPVs -------  n=" << simpv.size() << endl;
    for (std::vector<simPrimaryVertex>::iterator v0 = simpv.begin(); v0 != simpv.end(); v0++) {
      cout << "  x=" << setw(10) << setprecision(4) << v0->x << "  y=" << setw(10) << setprecision(4) << v0->y
           << "  z=" << setw(10) << setprecision(4) << v0->z << "  px=" << v0->ptot.px() << "  py=" << v0->ptot.py()
           << "  pz=" << v0->ptot.pz() << "  pt2=" << v0->ptsq << "  nTrkPrim=" << v0->nTrkPrim << "  nTrkSec="
           << v0->nTrkSec
           //	   << "  closest=" << scientific << v0->closest_vertex_distance_z
           << endl;
    }
    cout << "-----------------------------------------------" << endl;
  }
  return simpv;
}
/********************************************************************************************************/

/********************************************************************************************************/
/* get sim pv from TrackingParticles/TrackingVertex */
std::vector<PrimaryVertexAnalyzer4PU::simPrimaryVertex> PrimaryVertexAnalyzer4PU::getSimPVs(
    const edm::Handle<TrackingVertexCollection> tVC)

/********************************************************************************************************/
{
  bool SIMPVDEBUG = DEBUG_ & false;
  std::vector<PrimaryVertexAnalyzer4PU::simPrimaryVertex> simpv;

  if (verbose_) {
    std::cout << "getSimPVs from TrackingVertexCollection " << std::endl;
  }

  for (TrackingVertexCollection::const_iterator v = tVC->begin(); v != tVC->end(); ++v) {
    if ((v->position().t() < -20) || (v->position().t() > 20))
      continue;  // out of time pu

    if (SIMPVDEBUG) {
      std::cout << (v->eventId()).event() << v->position() << v->g4Vertices().size() << " " << v->genVertices().size()
                << "   t=" << v->position().t() * 1.e12 << "    ==0:" << (v->position().t() > 0) << std::endl;
      for (TrackingVertex::g4v_iterator gv = v->g4Vertices_begin(); gv != v->g4Vertices_end(); gv++) {
        std::cout << *gv << std::endl;
      }
      std::cout << "----------" << std::endl;
    }

    if ((unsigned int)v->eventId().event() < simpv.size())
      continue;
    if (fabs(v->position().z()) > 1000)
      continue;  // skip funny junk vertices

    // could be a new vertex, check  all primaries found so far to avoid multiple entries
    const double mm = 1.0;  // for tracking vertices
    const double s = 1.e9;
    simPrimaryVertex sv(v->position().x() * mm, v->position().y() * mm, v->position().z() * mm, v->position().t() * s);

    sv.eventId = v->eventId();
    sv.type = 1;
    // funny loop, why is it here?
    for (auto iTrack = v->daughterTracks_begin(); iTrack != v->daughterTracks_end(); ++iTrack) {
      sv.eventId = (**iTrack).eventId();  // an iterator of Refs, dereference twice
    }

    simPrimaryVertex* vp = NULL;  // will become non-NULL if a vertex is found and then point to it
    for (auto v0 = simpv.begin(); v0 != simpv.end(); v0++) {
      if ((sv.eventId == v0->eventId) ||
          ((fabs(sv.x - v0->x) < 1e-5) && (fabs(sv.y - v0->y) < 1e-5) && (fabs(sv.z - v0->z) < 1e-5))) {
        vp = &(*v0);
        break;
      }
    }

    if (!vp) {
      // this is a new vertex
      if (SIMPVDEBUG) {
        std::cout << "this is a new vertex " << sv.eventId.event() << "   " << sv.x << " " << sv.y << " " << sv.z
                  << std::endl;
      }
      // Loop over daughter tracks
      for (TrackingVertex::tp_iterator iTP = v->daughterTracks_begin(); iTP != v->daughterTracks_end(); ++iTP) {
        //double pt=(*(*iTP)).momentum().perp2();
        const TrackingParticle& tp = (*(*iTP));
        if (!(tp.charge() == 0) && (fabs(tp.momentum().eta()) < 2.5)) {
          if (tp.momentum().perp2() > 0.1 * 0.1) {
            sv.nGenTrk++;
            sv.sumpT += sqrt(tp.momentum().perp2());
          }
        }
      }
      sv.is_visible = (sv.nGenTrk > 2);
      simpv.push_back(sv);
      vp = &simpv.back();
    } else {
      if (SIMPVDEBUG) {
        std::cout << "this is not a new vertex" << sv.x << " " << sv.y << " " << sv.z << std::endl;
      }
    }

    if (SIMPVDEBUG) {
      for (TrackingVertex::tp_iterator iTP = v->daughterTracks_begin(); iTP != v->daughterTracks_end(); ++iTP) {
        std::cout << "  Daughter momentum:      " << (*(*iTP)).momentum();
        std::cout << "  Daughter type     " << (*(*iTP)).pdgId();
        std::cout << std::endl;
      }
    }
  }

  if (DEBUG_) {
    cout << "------- PrimaryVertexAnalyzer4PU simPVs from TrackingVertices -------" << endl;
    int idx = 0;
    for (auto v0 = simpv.begin(); v0 != simpv.end(); v0++) {
      cout << setw(3) << idx << ") "
           << "z=" << setw(8) << setprecision(4) << fixed << v0->z << "  event=" << setw(3) << v0->eventId.event()
           << endl;
      idx++;
    }
    cout << "-----------------------------------------------" << endl;
  }
  return simpv;
}
/********************************************************************************************************/

/********************************************************************************************************/
/* get sim pv from TrackingParticles/TrackingVertex */
std::vector<PrimaryVertexAnalyzer4PU::simPrimaryVertex> PrimaryVertexAnalyzer4PU::getSimPVs(
    const edm::Handle<TrackingVertexCollection> tVC, std::vector<PrimaryVertexAnalyzer4PU::simPrimaryVertex>& pu_simpv)
/********************************************************************************************************/
{
  bool SIMPVDEBUG = DEBUG_ & false;

  std::vector<PrimaryVertexAnalyzer4PU::simPrimaryVertex> simpv;

  if (verbose_) {
    std::cout << "getSimPVs from TrackingVertexCollection and other simpvs" << std::endl;
  }

  std::vector<int> eventIds;

  for (TrackingVertexCollection::const_iterator v = tVC->begin(); v != tVC->end(); ++v) {
    if ((v->eventId()).bunchCrossing() != 0)
      continue;
    if (std::find(eventIds.begin(), eventIds.end(), v->eventId().event()) != eventIds.end())
      continue;
    eventIds.push_back(v->eventId().event());

    if (SIMPVDEBUG) {
      std::cout << "----------" << std::endl;
      std::cout << (v->eventId()).event() << " " << (v->eventId()).bunchCrossing() << " " << v->position()
                << v->position().t() * simtUnit_ << v->g4Vertices().size() << " " << v->genVertices().size()
                << "    ==0:" << (v->position().t() > 0) << std::endl;
      /*
	for( TrackingVertex::g4v_iterator gv=v->g4Vertices_begin(); gv!=v->g4Vertices_end(); gv++)
	  {
	    std::cout << "     " << *gv << std::endl;
	  }
	*/
      std::cout << "----------" << std::endl;
    }

    //if ((unsigned int) v->eventId().event() < simpv.size()) continue;
    if (fabs(v->position().z()) > 1000)
      continue;  // skip funny junk vertices

    // could be a new vertex, check  all primaries found so far to avoid multiple entries
    const double mm = 1.0;  // for tracking vertices
    const double s = 1.e9;
    simPrimaryVertex sv(v->position().x() * mm, v->position().y() * mm, v->position().z() * mm, v->position().t() * s);

    sv.eventId = v->eventId();
    sv.type = 1;
    // funny loop, why is it here?
    for (auto iTrack = v->daughterTracks_begin(); iTrack != v->daughterTracks_end(); ++iTrack) {
      sv.eventId = (**iTrack).eventId();  // an iterator of Refs, dereference twice
    }

    simPrimaryVertex* vp = NULL;  // will become non-NULL if a vertex is found and then point to it
    for (auto v0 = simpv.begin(); v0 != simpv.end(); v0++) {
      if ((sv.eventId == v0->eventId) ||
          ((fabs(sv.x - v0->x) < 1e-5) && (fabs(sv.y - v0->y) < 1e-5) && (fabs(sv.z - v0->z) < 1e-5))) {
        vp = &(*v0);
        break;
      }
    }

    if (!vp) {
      // this is a new vertex
      if (SIMPVDEBUG) {
        std::cout << "this is a new vertex " << sv.eventId.event() << "   " << sv.x << " " << sv.y << " " << sv.z
                  << std::endl;
      }
      // Loop over daughter tracks
      for (TrackingVertex::tp_iterator iTP = v->daughterTracks_begin(); iTP != v->daughterTracks_end(); ++iTP) {
        //double pt=(*(*iTP)).momentum().perp2();
        const TrackingParticle& tp = (*(*iTP));
        if (!(tp.charge() == 0) && (fabs(tp.momentum().eta()) < etaMaxVisible_)) {
          if (tp.momentum().perp2() > (ptMinVisible_ * ptMinVisible_)) {
            sv.nGenTrk++;
            sv.sumpT += sqrt(tp.momentum().perp2());
          }
        }
      }

      sv.is_visible = (sv.nGenTrk > 2);
      simpv.push_back(sv);
      vp = &simpv.back();
    } else {
      if (SIMPVDEBUG) {
        std::cout << "this is not a new vertex " << sv.eventId.event() << " " << sv.x << " " << sv.y << " " << sv.z
                  << std::endl;
      }
    }

    if (SIMPVDEBUG && false) {
      for (TrackingVertex::tp_iterator iTP = v->daughterTracks_begin(); iTP != v->daughterTracks_end(); ++iTP) {
        std::cout << "  Daughter momentum:      " << (*(*iTP)).momentum();
        std::cout << "  Daughter type     " << (*(*iTP)).pdgId();
        std::cout << std::endl;
      }
    }
  }

  // compare the other collection
  if (SIMPVDEBUG) {
    cout << "getSimPVs(TV +) :" << simpv.size() << " " << pu_simpv.size() << endl;
    int idx = 0;
    int match = 0;
    for (auto v0 = simpv.begin(); v0 != simpv.end(); v0++) {
      for (auto v1 = pu_simpv.begin(); v1 != pu_simpv.end(); v1++) {
        if (abs(v0->z - v1->z) < 0.0001) {
          match += 1;
          if (verbose_) {
            cout << setw(3) << idx << ") "
                 << "z=" << setw(8) << setprecision(4) << fixed << v0->z << "  ptHat=" << setw(6) << v1->pt_hat << " "
                 << sqrt(v0->sumpT) << endl;
          }
          v0->pt_hat = v1->pt_hat;
        }
      }
      idx++;
    }
    cout << "getSimPVs(TV +) : " << match << "  matched" << endl;
  }

  if (DEBUG_ || SIMPVDEBUG) {
    cout << "------- PrimaryVertexAnalyzer4PU simPVs from TrackingVertices + -------" << endl;
    int idx = 0;
    for (auto v0 = simpv.begin(); v0 != simpv.end(); v0++) {
      cout << setw(3) << idx << ") "
           << "z=" << setw(8) << setprecision(4) << fixed << v0->z << "  event=" << setw(3) << v0->eventId.event()
           << endl;
      idx++;
    }
    cout << "-----------------------------------------------" << endl;
  }

  return simpv;
}
/********************************************************************************************************/




/********************************************************************************************************/
reco::VertexCollection* PrimaryVertexAnalyzer4PU::vertexFilter(Handle<reco::VertexCollection> pvs, bool filter)
/********************************************************************************************************/
{
  reco::VertexCollection* pv = new reco::VertexCollection;
  if (filter) {
    // ptmin filter
    for (reco::VertexCollection::const_iterator ipv = pvs->begin(); ipv != pvs->end(); ipv++) {
      double ptmin = 0;
      for (trackit_t tv = ipv->tracks_begin(); tv != ipv->tracks_end(); tv++) {
        double pt = tv->get()->pt();
        if (pt > ptmin) {
          ptmin = pt;
        }
      }
      if (ptmin > 0.5) {
        pv->push_back(*ipv);
      }
    }
  } else {
    for (reco::VertexCollection::const_iterator ipv = pvs->begin(); ipv != pvs->end(); ipv++) {
      pv->push_back(*ipv);
    }
  }
  return pv;
}
/********************************************************************************************************/







/********************************************************************************************************/
void PrimaryVertexAnalyzer4PU::analyze(const Event& iEvent, const EventSetup& iSetup)
/********************************************************************************************************/

{
  eventcounter_++;
  run_ = iEvent.id().run();
  luminosityBlock_ = iEvent.luminosityBlock();
  event_ = iEvent.id().event();
  bunchCrossing_ = iEvent.bunchCrossing();
  orbitNumber_ = iEvent.orbitNumber();
  if (sigmaZoverride_ > 0)
    sigmaZ_ = sigmaZoverride_;
  MC_ = false;
  dumpThisEvent_ = false;
  forceDump_ = false;   // use with caution
  lsglb_ = 0;

  // in case we wanted to analyze a specific lumi block only
  if ((analyzeLS_ >= 0) && !(luminosityBlock_ == analyzeLS_))
    return;
  
  // get lumi info
  if(run_ > 1) get_luminosity_infos(iEvent);

  if (verbose_) {
    std::cout << endl
              << "PrimaryVertexAnalyzer4PU::analyze   event counter=" << eventcounter_ << " Run=" << run_
              << "  LumiBlock " << luminosityBlock_ << "  event  " << event_ << " bx=" << bunchCrossing_
              << "  orbit=" << orbitNumber_ << std::endl;
  }

  //Retrieve tracker topology from geometry, need for track cluster infos
  edm::ESHandle<TrackerTopology> tTopoH;
  iSetup.get<TrackerTopologyRcd>().get(tTopoH);
  tTopo_ = tTopoH.product();

  if (!do_vertex_analysis_)
    return;

  // load pile-up info and select events in [nPUMin_, nPUmax_] if requested
  PileupSummaryInfo puInfo;
  bool bPuInfo = getPuInfo(iEvent, puInfo);
  if (bPuInfo) {
    report_counted("found PU info", 1);
  }else{
    report_counted("no PU info found", 1);
  }    
  if (bPuInfo && ((puInfo.getPU_zpositions().size() < nPUmin_) || (puInfo.getPU_zpositions().size() > nPUmax_))) {
    if (verbose_) {
      cout << "skipping event, out of pu window  npu=" << puInfo.getPU_zpositions().size() << endl;
    }
    return;
  }

  if (!get_beamspot_data(iEvent)) {
    return;
  }


  // load vertex collections
  for (std::vector<std::string>::const_iterator vCollection = vertexCollectionLabels_.begin();
       vCollection != vertexCollectionLabels_.end();
       vCollection++) {

    if (vCollection->find("WithBS") !=std::string::npos){
      selNdof_ = selNdofWithBS_;
      ndof0trk_ = -1.;
    }else{
      selNdof_ = selNdofNoBS_;
      ndof0trk_ = -3.;
    }

    Handle<reco::VertexCollection> recVtxsH;
    if (iEvent.getByToken(vertexCollectionTokens_[*vCollection], recVtxsH)) {
      recVtxs_[*vCollection] = vertexFilter(recVtxsH, useVertexFilter_);
      analyzeVertexRecoCPUTime(histograms_[*vCollection], recVtxs_[*vCollection], *vCollection);
      analyzeVertexCollectionRecoNoTracks(histograms_[*vCollection], recVtxs_[*vCollection], *vCollection);
    } else {
      recVtxs_[*vCollection] = NULL;
      cout << "collection " << *vCollection << " not found " << endl;
    }
  }

  Tracks tracks;

  if (!get_reco_and_transient_tracks(iSetup, iEvent, tracks)) {
    // clean up
    for (auto vCollection = vertexCollectionLabels_.begin(); vCollection != vertexCollectionLabels_.end();
         vCollection++) {
      delete recVtxs_[*vCollection];
    }
    return;  //FIXME some things can be done without tracks
  }

  // collect MC information, fills

  std::vector<SimEvent> simEvt;         // full tracking particle information
  std::vector<simPrimaryVertex> simpv;  //  a list of primary MC vertices, filled from whatever we get
  std::vector<SimPart> tsim;            //  simulated tracks filled from whatever we get

  if (run_ < 1000) {
    get_particle_data_table(iSetup);

    MC_ = get_MC_truth(iEvent, tracks, bPuInfo, puInfo, simEvt, simpv, tsim);
    
    if (MC_){
      fill_simvtx_histos(simpv);
      // dump time variables
      /*
      for(auto tk : tracks.recotracks){
	if (tk.has_timing && tk.matched){
	  cout << tk.tpr->pdgId() << " tsim=" << tk.tsim  << " trec=" << tk.t  << " trec(pid)=" << tk.get_t_pid()  << "    pion" << tk.th[0] << "  kaon" << tk.th[1] << " proton" << tk.th[2] << endl;
	}
      }
      */
    }
  }

  // analyze the vertex collections
  for (std::vector<std::string>::const_iterator vCollection = vertexCollectionLabels_.begin();
       vCollection != vertexCollectionLabels_.end();
       vCollection++) {

    if (recVtxs_[*vCollection] == NULL)
      continue;

    if (verbose_) {
      cout << "**** analyzing " << *vCollection << "  with size " << recVtxs_[*vCollection]->size() << endl;
    }

    // replace by set_ndof_globals(*vCollection);
    if (vCollection->find("WithBS") !=std::string::npos){
      selNdof_ = selNdofWithBS_;
      ndof0trk_ = -1.;
    }else{
      selNdof_ = selNdofNoBS_;
      ndof0trk_ = -3.;
    }

    // set up the track-> vertex map
    trkkey2recvtx_.clear();
    unsigned int iv = 0;
    for (auto v = recVtxs_[*vCollection]->begin(); v != recVtxs_[*vCollection]->end(); ++v) {
      for (trackit_t t = v->tracks_begin(); t != v->tracks_end(); t++) {
        trkkey2recvtx_[t->key()] = iv;
        tracks.from_key(t->key()).recv[*vCollection] = iv;
      }
      iv++;
    }

    timer_start("analyzeVertexCollectionReco");
    analyzeVertexCollectionReco(histograms_[*vCollection], recVtxs_[*vCollection], tracks, *vCollection);
    timer_stop("analyzeVertexCollectionReco");
    
    if (MC_) {
      timer_start("analyzeVertexCollectionSimPvNoSimTracks");
      analyzeVertexCollectionSimPvNoSimTracks(
          histograms_[*vCollection], recVtxs_[*vCollection], tracks, simpv, *vCollection);
      timer_stop("analyzeVertexCollectionSimPvNoSimTracks");

    /* don't when there are tracking particles, needs revision anyway
      analyzeVertexCollectionSimPv(
          histograms_[*vCollection], recVtxs_[*vCollection], tracks, simpv, tsim, *vCollection);
      */
      
      //analyzeVertexCollectionDQMMC(histograms_[*vCollection], recVtxs_[*vCollection],  trackCollectionH, recTrks, simpvHepMC, tsim, *vCollection);

      
      timer_start("matching");
      recvmatch_[*vCollection] =
          tpmatch(recVtxs_[*vCollection], simEvt, tracks);
      wos_match(recvmatch_[*vCollection], recVtxs_[*vCollection], simEvt, tracks);
      timer_stop("matching");

      if (matchsummaries_-- > 0)
        printMatchingSummary(recVtxs_[*vCollection], simEvt, recvmatch_[*vCollection], *vCollection);

      
      timer_start("analyzeMergeRateTP");
      analyzeVertexMergeRateTP(
          histograms_[*vCollection], recVtxs_[*vCollection], simEvt, recvmatch_[*vCollection], *vCollection);
      timer_stop("analyzeMergeRateTP");

      
      timer_start("analyzeVertexCollectionTP");
      analyzeVertexCollectionTP(
          histograms_[*vCollection], recVtxs_[*vCollection], tracks, simEvt, recvmatch_[*vCollection], *vCollection);
      timer_stop("analyzeVertexCollectionTP");

      analyzeVertexCollectionPtvis(
				   histograms_[*vCollection], recVtxs_[*vCollection], tracks, simEvt, recvmatch_[*vCollection], *vCollection);
    }
  }

  if (do_vertex_analysis_ && ((nCompareCollections_ < 0) || (eventcounter_ <= nCompareCollections_))) {
    cout << " comparing collections " << endl;
    compareCollections(simEvt, simpv);
  }

  // print summary info
  if ((dumpThisEvent_ && (dumpcounter_ < ndump_)) || (verbose_ && (eventcounter_ < ndump_)) ||
      (autoDumpCounter_-- > 0) || (forceDump_)) {
    cout << endl
         << "Event dump " << dumpcounter_ << endl
         << "event counter=" << eventcounter_ << " Run=" << run_ << "  LumiBlock " << luminosityBlock_ << "  event  "
         << event_ << " bx=" << bunchCrossing_ << "(" << firstBXinTrain_ << ")"
         << " orbit=" << orbitNumber_ << " PU= " << lumiPU_ << " sigma_z= " << sigmaZ_
         << " zbeam = " << vertexBeamSpot_.z0() << std::endl;

    cout << "  dumpThisEvent_ = " << dumpThisEvent_ << "  dumpcounter_ = " << dumpcounter_ << "  ndump_ = " << ndump_
         << "  verbose_ = " << verbose_ << " eventcounter_ = " << eventcounter_ << "  autoDumpCounter_"
         << autoDumpCounter_ << "  forceDump_ = " << forceDump_ << endl;

    dumpcounter_++;
    forceDump_ = false;

    bool trksdumped = (--ndump_tracks_ < 0);

    for (auto vCollection = vertexCollectionLabels_.begin(); vCollection != vertexCollectionLabels_.end();
         vCollection++) {
      if (recVtxs_[*vCollection] != NULL) {
	cout << " dumping collection " << *vCollection << endl;
        printRecVtxs(recVtxs_[*vCollection]);
	// redo the matching for the dump (invalid if another collection has been analyzed)
	recvmatch_[*vCollection] =
          tpmatch(recVtxs_[*vCollection], simEvt, tracks);
	wos_match(recvmatch_[*vCollection], recVtxs_[*vCollection], simEvt, tracks);

	printMatchingSummary(recVtxs_[*vCollection], simEvt, recvmatch_[*vCollection], *vCollection);
        if (!trksdumped) {
          cout << "--- dumping " << *vCollection << endl;
          printPVTrksZT(tracks.trackCollectionH, recVtxs_[*vCollection], tsim, simpv, simEvt);

          cout << "---" << endl;
          trksdumped = true;  // only dump one track list per event
        }
      } else {
        cout << "Vertex collection not found !!!  This should not happen!" << endl;
      }
    }

    if (dumpcounter_ < 2) {
      cout << "beamspot " << vertexBeamSpot_ << endl;
    }
    if (verbose_)
      cout << endl << endl;
  }

  // clean up
  for (auto vCollection = vertexCollectionLabels_.begin(); vCollection != vertexCollectionLabels_.end();
       vCollection++) {
    delete recVtxs_[*vCollection];
  }

  if (verbose_) {
    std::cout << std::endl << " End of PrimaryVertexAnalyzer4PU " << std::endl;
  }
}
// end of "analyze"

/***************************************************************************************/

/***************************************************************************************/
void PrimaryVertexAnalyzer4PU::reportEvent(const char* message, bool dump) {
  stringstream s;
  s << "run:event:ls:bx = " << run_ << ":" << setw(10) << setfill('0') << event_ << ":" << setw(4) << setfill('0')
    << luminosityBlock_ << ":" << setw(4) << setfill('0') << bunchCrossing_ << " PU=" << setfill(' ') << setw(5)
    << setprecision(1) << fixed << lumiPU_ << " avg PU=" << setw(5) << setprecision(1) << avglumiPU_
    << " sigma_z=" << setw(6) << setprecision(2) << fixed << sigmaZ_ << "    ###   " << message;
  cout << " report  " << s.str() << endl;
  reports_.push_back("[report] " + s.str());
  dumpThisEvent_ |= dump;
}
/***************************************************************************************/


/***************************************************************************************/
void PrimaryVertexAnalyzer4PU::report_counted(const std::string msg, const int max_count)
{
  if (msg == "summary")
    {
      for(auto m : counted_messages_)
	{
	  std::cout << std::setw(8) << m.second.first << " " << m.first << std::endl;
	}
      return;
    }

  report_counted(msg, "", max_count);
  /*
  if ( counted_messages_.find(msg) == counted_messages_.end() )
    {
      counted_messages_[msg] = std::make_pair<unsigned int, unsigned int>(1, max_count);
    }
  else
    {
      counted_messages_[msg].first++;
    }    

  if (counted_messages_[msg].first <= counted_messages_[msg].second)
    {
    std::cout << msg << std::endl;
    }
  */
}
/***************************************************************************************/



/***************************************************************************************/
void PrimaryVertexAnalyzer4PU::report_counted(const std::string msg, const std::string msg2, const int max_count)
{
 
  if ( counted_messages_.find(msg) == counted_messages_.end() )
    {
      counted_messages_[msg] = std::make_pair<unsigned int, unsigned int>(1, max_count);
    }
  else
    {
      counted_messages_[msg].first++;
    }    

  if (counted_messages_[msg].first <= counted_messages_[msg].second)
    {
      if (msg2 == ""){
	std::cout << msg << std::endl;
      }else{
	std::cout << msg << " " << msg2 << std::endl;
      }
    }
}
/***************************************************************************************/


/***************************************************************************************/
void PrimaryVertexAnalyzer4PU::reportVertex(const Vertex& v, const char* message, bool dump) {
  reportEvent(message, dump);
  double dx = v.x() - vertexBeamSpot_.x(v.z());
  double dy = v.y() - vertexBeamSpot_.y(v.z());
  std::cout << " #trk " << std::fixed << std::setprecision(4) << std::setw(3) << v.tracksSize() << " chi2 "
            << std::fixed << std::setw(5) << std::setprecision(1) << v.chi2() << " ndof " << std::fixed << std::setw(6)
            << std::setprecision(2) << v.ndof() 
	    << " x "  << std::setw(8) << std::fixed << std::setprecision(4) << v.x() << " dx " << std::setw(6) << v.xError()
            << " y "  << std::setw(8) << v.y() << " dy " << std::setw(6) << v.yError()
            << " z "  << std::setw(8) << v.z() << " dz " << std::setw(6) << v.zError() ;

  if(f4D_){
    std::cout << " t " << std::setw(7)<< std::setprecision(3) << v.t() << " dt " << std::setw(5) << v.tError() ;
  }
  std::cout  << " r  " << std::setw(8)<< std::setprecision(4) << sqrt(dx * dx + dy * dy);
  cout << endl;
}
/***************************************************************************************/

/***************************************************************************************/
void PrimaryVertexAnalyzer4PU::printEventSummary(std::map<std::string, TH1*>& h,
                                                 const reco::VertexCollection* recVtxs,
                                                 Tracks& tracks,
                                                 std::vector<simPrimaryVertex>& simpv,
                                                 const string message) {
  // make a readable summary using simpv (no TrackingParticles, use simparticles or genparticles etc)
  if (simpv.size() == 0)
    return;

  cout << endl << "printEventSummary  (simpv)" << endl;
  cout << "---------------------------" << endl;
  cout << "event counter = " << eventcounter_ << "   " << message << endl;
  cout << "run " << run_ << "    event " << event_ << endl;
  cout << "---------------------------" << endl;

  vector<pair<double, unsigned int>> zrecv;
  for (unsigned int idx = 0; idx < recVtxs->size(); idx++) {
    //if (recVtxs->at(idx).isFake()) continue;
    zrecv.push_back(make_pair(recVtxs->at(idx).z(), idx));
  }
  stable_sort(zrecv.begin(), zrecv.end(), lt);

  // same for simulated vertices
  vector<pair<double, unsigned int>> zsimv;
  for (unsigned int idx = 0; idx < simpv.size(); idx++) {
    zsimv.push_back(make_pair(simpv[idx].z, idx));
  }
  stable_sort(zsimv.begin(), zsimv.end(), lt);

  unsigned int idxrec = 0;
  unsigned int idxsim = 0;
  double zmatch = 0.05;
  cout.precision(4);

  cout << "  rec z"
       << "             "
       << " sim z"
       << "    pT                                                  #trk" << endl;
  while ((idxrec < recVtxs->size()) || (idxsim < simpv.size())) {
    string signal = " ";
    string tag = " ";
    if ((idxsim < simpv.size()) && (zsimv[idxsim].second == 0)) {
      signal = "*";
    }
    if ((idxrec < recVtxs->size()) && (zrecv[idxrec].second == 0)) {
      tag = "*";
    }

    double ndof = 0, pxy = 0, ptmax2 = 0;
    if (idxrec < recVtxs->size()) {
      ndof = recVtxs->at(zrecv[idxrec].second).ndof();
      pxy = vertex_pxy(recVtxs->at(zrecv[idxrec].second));
      ptmax2 = vertex_ptmax2(recVtxs->at(zrecv[idxrec].second));
    }

    if ((idxrec < recVtxs->size()) && (idxsim < simpv.size()) &&
        (std::abs(zrecv[idxrec].first - zsimv[idxsim].first) <
         std::max(zmatch, 3 * recVtxs->at(zrecv[idxrec].second).zError())) &&
        (((idxsim + 1) == simpv.size()) ||
         (fabs(zrecv[idxrec].first - zsimv[idxsim].first) < fabs(zrecv[idxrec].first - zsimv[idxsim + 1].first))) &&
        (((idxrec + 1) == recVtxs->size()) ||
         (fabs(zrecv[idxrec].first - zsimv[idxsim].first) < fabs(zrecv[idxrec + 1].first - zsimv[idxsim].first)))) {
      cout << setw(8) << setprecision(4) << fixed << zrecv[idxrec].first << tag << "   <->    " << setw(8) << fixed
           << zsimv[idxsim].first << signal << setw(5) << setprecision(1) << fixed << simpv[zsimv[idxsim].second].pt_hat
           << " (ndof=" << fixed << setw(6) << setprecision(1) << ndof << ", p=" << setw(6) << setprecision(4) << pxy
           << ", ptmax2=" << setw(4) << setprecision(1) << ptmax2 << ")";
      // cout  << setw(5) << zrecv[idxrec].second; // debugging only

      if (simpv[zsimv[idxsim].second].nGenTrk > 0) {
        cout << "            (" << fixed << setw(3) << simpv[zsimv[idxsim].second].nGenTrk << ")";
        //cout << "            (" << fixed << setw(3) << simpv[zsimv[idxsim].second].pt_hat << ")";
      }

      if (zsimv[idxsim].second == 0) {
        if (tag == " ") {
          cout << "  signal vertex not tagged" << endl;
        } else {
          cout << "  signal vertex found and tagged" << endl;
        }
      } else {
        cout << endl;
      }

      idxrec++;
      idxsim++;

    } else if (((idxrec < recVtxs->size()) && (idxsim < simpv.size()) && (zrecv[idxrec].first < zsimv[idxsim].first)) ||
               ((idxrec < recVtxs->size()) && (idxsim == simpv.size()))) {
      cout << setw(8) << setprecision(4) << fixed << zrecv[idxrec].first << tag << "                       "
           << "  (ndof=" << fixed << setw(6) << setprecision(1) << ndof << ", p=" << setw(6) << setprecision(4) << pxy
           << ", ptmax2=" << setw(4) << setprecision(1) << ptmax2 << ")"
           << "   fake " << endl;
      idxrec++;

    } else if (((idxrec < recVtxs->size()) && (idxsim < simpv.size()) && (zrecv[idxrec].first > zsimv[idxsim].first)) ||
               ((idxrec == recVtxs->size()) && (idxsim < simpv.size()))) {
      cout << "         "
           << "   <->    " << setw(8) << setprecision(4) << fixed << zsimv[idxsim].first << signal << setw(5)
           << setprecision(1) << simpv[zsimv[idxsim].second].pt_hat;
      if (simpv[zsimv[idxsim].second].is_visible) {
        if (zsimv[idxsim].second == 0) {
          cout << "                                        lost signal vertex" << endl;
        } else {
          cout << "                                        lost PU  ";
          if (simpv[zsimv[idxsim].second].nGenTrk > 0) {
            cout << "(" << setw(3) << simpv[zsimv[idxsim].second].nGenTrk << ")";
            //   << "(" << setw(5) << setprecision(3) << simpv[zsimv[idxsim].second].pt_hat << ")"
          }
          cout << endl;
        }
      } else {
        cout << endl;
      }
      idxsim++;
    } else {
      cout << "what else?" << endl;
      break;
    }
  }
  std::cout << "---over and out" << std::endl;
}
/***************************************************************************************/




/***************************************************************************************/
void PrimaryVertexAnalyzer4PU::compareCollections(vector<SimEvent>& simEvt, std::vector<simPrimaryVertex>& simpv) {
  /***************************************************************************************
prints a side-by-side comparison of the selected vertex collections
simulated vertices are shown if they are available
***************************************************************************************/

  bool debug = false;

  unsigned int ncoll = vertexCollectionLabels_.size();

  cout << "----compareCollections------------------------" << endl;
  cout << "event counter = " << eventcounter_ << "   ";
  cout << "run " << run_ << "    event " << event_ << endl;
  cout << "----------------------------------------------" << endl;

  cout << setw(15) << " simulated [ ch, tk]";

  for (unsigned int icol = 0; icol < ncoll; icol++) {
    string h = vertexCollectionLabels_[icol];
    if (h.size() > 17) {
      h = h.substr(h.size() - 17, h.size());
    } else {
      h.resize(17);
    }
    cout << setw(18) << h;
  }
  cout << endl;

  int differences = 0;

  // build the table

  vector<pair<double, unsigned int*>> row;

  // if we don't have simEvt with truth matching (TrackingParticles)
  // put in sim vertices
  if (simEvt.size() == 0) {
    for (unsigned int idx = 0; idx < simpv.size(); idx++) {
      unsigned int* v = new unsigned int[ncoll + 1];
      for (unsigned int j = 0; j < ncoll + 1; j++) {
        v[j] = -1;
      }
      v[ncoll] = idx;
      row.push_back(make_pair(simpv[idx].z, v));
    }
  }

  // unmatched rec vertices from all collections
  for (unsigned int icol = 0; icol < ncoll; icol++) {
    reco::VertexCollection* recVtxs = recVtxs_[vertexCollectionLabels_[icol]];
    if (recVtxs) {
      for (unsigned int idx = 0; idx < recVtxs->size(); idx++) {
	if (recVtxs->at(idx).isFake()) continue;
        if (recvmatch_[vertexCollectionLabels_[icol]][idx].sim == NOT_MATCHED) {
          unsigned int* v = new unsigned int[ncoll + 1];
          for (unsigned int j = 0; j < ncoll + 1; j++) {
            v[j] = NOT_MATCHED;
          }
          v[icol] = idx;
          row.push_back(make_pair(recVtxs->at(idx).z(), v));
        }
      }
    }
  }

  if (row.size() > 1) {
    if (debug) {
      cout << "dump    size=" << row.size() << endl;
      for (unsigned int irow = 0; irow < row.size(); irow++) {
        cout << setw(2) << irow << ")";
        cout << setw(8) << setprecision(4) << fixed << row[irow].first;
        if (row[irow].second == NULL)
          continue;
        for (unsigned int i = 0; i < ncoll + 1; i++) {
          cout << setw(6) << row[irow].second[i];
        }
        cout << endl;
      }
      cout << endl;
    }

    // join rows
    int join = 0;
    while (join >= 0) {
      if (row.size() > 1) {
        stable_sort(row.begin(), row.end());
      }

      if (debug) {
        cout << "dump before joining  size=" << row.size() << endl;
        for (unsigned int irow = 0; irow < row.size(); irow++) {
          cout << setw(2) << irow << ")";
          cout << setw(8) << setprecision(4) << fixed << row[irow].first;
          if (!(row[irow].second == NULL)) {
            for (unsigned int i = 0; i < ncoll + 1; i++) {
              cout << setw(6) << row[irow].second[i];
            }
          }
          cout << endl;
        }
        cout << endl;
      }

      double dzmin = 0.1;
      join = -1;
      for (unsigned int irow = 0; irow < row.size() - 1; irow++) {
        if ((row[irow].second == NULL) || (row[irow + 1].second == NULL))
          continue;
        if ((row[irow + 1].first - row[irow].first) < dzmin) {
          bool joinable = true;
          for (unsigned int i = 0; i < ncoll + 1; i++) {
            if ((row[irow].second[i] != NOT_MATCHED) && (row[irow + 1].second[i] != NOT_MATCHED))
              joinable = false;
          }
          if (joinable) {
            join = irow;
            dzmin = fabs(row[irow + 1].first - row[irow].first);
          }
        }
      }

      if (join >= 0) {
        if (debug)
          cout << "joining " << join << endl;
        if ((row[join].second == NULL) || (row[join + 1].second == NULL)) {
          cout << " bad join=" << join << endl;
        }
        if (join >= int(row.size())) {
          cout << " row pointers screwed up " << join << "   " << row.size() << endl;
        }
        //join
        for (unsigned int i = 0; i < ncoll + 1; i++) {
          if ((row[join].second[i] == NOT_MATCHED) && (row[join + 1].second[i] != NOT_MATCHED)) {
            row[join].second[i] = row[join + 1].second[i];
            if (i == ncoll)
              row[join].first = row[join + 1].first;
          }
        }

        //row z
        if (row[join].second[ncoll] == NOT_MATCHED) {
          double zrow = 0;
          int nv = 0;
          for (unsigned int i = 0; i < ncoll; i++) {
            int iv = row[join].second[i];
            if (iv > int(recVtxs_[vertexCollectionLabels_[i]]->size())) {
              cout << "illegal vertex index " << iv << "    join=" << join << endl;
            }
            if (iv >= 0) {
              reco::VertexCollection* recVtxs = recVtxs_[vertexCollectionLabels_[i]];
              zrow += recVtxs->at(iv).z();
              nv++;
            }
          }
          if (nv > 0) {
            row[join].first = zrow / nv;
          } else {
            // hmmm
          }
        }
        //delete swallowed row
        if (debug)
          cout << "deleting row " << join + 1 << "  row size= " << row.size() << "  ncoll= " << ncoll << endl;
        delete[] row[join + 1].second;
        row[join + 1].second = NULL;
        row.erase(row.begin() + (join + 1));

        if (debug) {
          cout << "dump after joining  " << join << " with " << join + 1 << endl;
          for (unsigned int irow = 0; irow < row.size(); irow++) {
            cout << setw(2) << irow << ")";
            cout << setw(8) << setprecision(4) << fixed << row[irow].first;
            if (!(row[irow].second == NULL)) {
              for (unsigned int i = 0; i < ncoll + 1; i++) {
                cout << setw(6) << row[irow].second[i];
              }
            }
            cout << endl;
          }
          cout << endl;
        }
      }
    }

    if (debug) {
      cout << "dump after joining  size=" << row.size() << endl;
      for (unsigned int irow = 0; irow < row.size(); irow++) {
        cout << setw(2) << irow << ")";
        cout << setw(8) << setprecision(4) << fixed << row[irow].first;
        if (!(row[irow].second == NULL)) {
          for (unsigned int i = 0; i < ncoll + 1; i++) {
            cout << setw(6) << row[irow].second[i];
          }
        }
        cout << endl;
      }
      cout << endl;
    }

  }  // handle un-matched vertices and simpv's

  // fill in sim vertices and matched rec vertices
  unsigned int suppressed = 0;
  for (unsigned int idx = 0; idx < simEvt.size(); idx++) {
    if (simEvt[idx].nChTP > 0) {
      unsigned int* v = new unsigned int[ncoll + 1];
      for (unsigned int j = 0; j < ncoll + 1; j++) {
        v[j] = NOT_MATCHED;
      }
      for (unsigned int jcol = 0; jcol < ncoll; jcol++) {
        reco::VertexCollection* recVtxs = recVtxs_[vertexCollectionLabels_[jcol]];
        if (!(recVtxs == NULL)) {
          int i = NOT_MATCHED;
          for (unsigned int j = 0; j < recVtxs->size(); j++) {
            if (recvmatch_[vertexCollectionLabels_[jcol]][j].sim == idx) {
              i = j;
            }
          }
          v[jcol] = i;
        }
      }
      v[ncoll] = idx;  // the sim vertex
      row.push_back(make_pair(simEvt[idx].z, v));
    } else {
      suppressed++;
    }
  }

  if (row.size() > 1) {
    stable_sort(row.begin(), row.end());
  }

  if (debug) {
    cout << "final dump  size=" << row.size() << endl;
    for (unsigned int irow = 0; irow < row.size(); irow++) {
      cout << setw(2) << irow << ")";
      cout << setw(8) << setprecision(4) << fixed << row[irow].first;
      if (!(row[irow].second == NULL)) {
        for (unsigned int i = 0; i < ncoll + 1; i++) {
          cout << setw(6) << row[irow].second[i];
        }
      }
      cout << endl;
    }
    cout << endl;
    cout << "done" << endl;
  }

  // readable dump
  for (unsigned int irow = 0; irow < row.size(); irow++) {
    if (row[irow].second == NULL)
      continue;

    double z = row[irow].first;

    unsigned int* v = row[irow].second;
    unsigned int idx0 = v[ncoll];

    if (idx0 == NOT_MATCHED) {
      cout << "%                    ";
    } else {
      if (simEvt.size() > 0) {
        // sim vertex
        cout << fixed << setw(10) << setprecision(4) << z << " [" << setw(3) << simEvt[idx0].nChTP << "," << setw(3)
             << simEvt[idx0].rtk.size() << "]";
        if (idx0 == 0) {
          cout << "*";
        } else {
          cout << " ";
        }
      } else {
        cout << fixed << setw(10) << setprecision(4) << z << " [" << setw(3) << simpv[idx0].nGenTrk << "]    ";
        if (idx0 == 0) {
          cout << "*";
        } else {
          cout << " ";
        }
      }
    }

    // count collections that  have a rec vertex for this sim vertex (for reporting)
    unsigned int nrec = 0;
    for (unsigned int jcol = 0; jcol < ncoll; jcol++) {
      if (v[jcol] != NOT_MATCHED) {
        nrec++;
      }
    }
    if ((nrec > 0) && (nrec < ncoll)) {
      differences++;
    }

    for (unsigned int jcol = 0; jcol < ncoll; jcol++) {
      unsigned int idx = v[jcol];
      if (idx != NOT_MATCHED) {
        reco::VertexCollection* recVtxs = recVtxs_[vertexCollectionLabels_[jcol]];
        if (!(recVtxs == NULL)) {
          cout << fixed << setw(10) << setprecision(4) << recVtxs->at(idx).z() << " (" << setw(5) << setprecision(1)
               << recVtxs->at(idx).ndof() << ")";
        } else {
          cout << "                  ";
        }
      } else {
        if (idx0 == NOT_MATCHED) {
          // no sim vertex, "fake not not found" (found by others)
          cout << "       +          ";
        } else {
          if (nrec == 0) {
            cout << "       -          ";
          } else {
            cout << "      ---         ";  // missed one (found by others)
          }
        }
      }
    }
    cout << endl;
  }

  for (unsigned int irow = 0; irow < row.size(); irow++) {
    if (!(row[irow].second == NULL)) {
      delete[] row[irow].second;
      row[irow].second = NULL;
    }
  }

  if (differences > 0) {
    cout << "the collections differ,  " << differences << "  differences " << endl;
  }
  if (suppressed > 0) {
    cout << suppressed << "  sim vertices without tracker hits suppressed " << endl;
  }
}
/***************************************************************************************/

std::string PrimaryVertexAnalyzer4PU::formatMatchList(std::map<unsigned int, double>& valuemap,
                                                      unsigned int nfield,
                                                      bool sim) {
  // sort
  std::vector<std::pair<double, unsigned int>> values;
  for (auto it : valuemap) {
    values.push_back(make_pair(it.second, it.first));
  }
  stable_sort(values.rbegin(), values.rend());  // reverse

  std::ostringstream s;
  unsigned int n = 0;
  while ((n < nfield) && (n < values.size())) {
    if (n > 0) {
      s << ",";
    }

    s << setw(5) << setprecision(1) << fixed << values[n].first;
    if (sim) {
      s << "(" << setw(3) << setprecision(0) << fixed << values[n].second << ")";
    } else {
      s << "[" << setw(3) << setprecision(0) << fixed << values[n].second << "]";
    }
    n++;
  }

  if (n < values.size()) {
    s << "..";
  } else {
    while (n < nfield) {
      s << "           ";
      n++;
    }
    s << "  ";
  }

  return s.str();
}

/***************************************************************************************/
void PrimaryVertexAnalyzer4PU::printMatchingSummary(const reco::VertexCollection* recVtxs,
                                                    vector<SimEvent>& simEvt,
                                                    std::vector<RSmatch>& recvmatch,
                                                    const string message)
// dump details about the matching, tracking particles only!
{
  if (simEvt.size() == 0) {
    return;
  }

  cout << endl << "printMatchingSummary  (simEvt)" << endl;
  cout << "---------------------------" << endl;
  cout << "event counter = " << eventcounter_ << "   " << message << endl;
  cout << "run " << run_ << "    event " << event_ << endl;

  cout << "     z        t   [sim]ntrk   wnt(rec)                             z        t   (rec) ndof   wnt [sim]" << endl;
  //  dump matched and unmatched vertices here, then sort in z
  vector<pair<double, std::string>> row;
  vector<unsigned int> dumped_rec;
  
  // unmatched rec
  for (unsigned int idxrec = 0; idxrec < recVtxs->size(); idxrec++) {
    if (recVtxs->at(idxrec).isFake()) continue;
    if (recvmatch[idxrec].sim == NOT_MATCHED) {
      std::ostringstream s;
      s << "%" << setw(61) << " " << setw(9) << setprecision(4) << fixed << recVtxs->at(idxrec).z() << "," << setw(7)
        << setprecision(3) << fixed << recVtxs->at(idxrec).t() << " " << setw(1) << "(" << setw(3) << idxrec << ") "
	<< setprecision(1) << fixed << setw(5) << recVtxs->at(idxrec).ndof() << " "
        << formatMatchList(recvmatch[idxrec].wnt, 3, false);
      row.push_back(make_pair(recVtxs->at(idxrec).z(), s.str()));
      dumped_rec.push_back(idxrec);
    }
  }

  // unmatched sim
  for (unsigned int idxsim = 0; idxsim < simEvt.size(); idxsim++) {
    if ((simEvt[idxsim].nChTP > 0) && (simEvt[idxsim].rec == NOT_MATCHED)) {
      std::ostringstream s;
      s << setw(9) << setprecision(4) << fixed << simEvt[idxsim].z << "," << setw(7) << setprecision(3) << fixed
        << simEvt[idxsim].t << " " << setw(1) << "[" << setw(3) << idxsim << "]"
	<< setw(4) << simEvt[idxsim].rtk.size() << " "
        << formatMatchList(simEvt[idxsim].wnt, 3, true)
	<< "%";
      row.push_back(make_pair(simEvt[idxsim].z, s.str()));
    }
  }

  // matched rec + sim
  for (unsigned int idxsim = 0; idxsim < simEvt.size(); idxsim++) {
    if (simEvt[idxsim].rec != NOT_MATCHED) {
      unsigned int idxrec = simEvt[idxsim].rec;
      if( recvmatch[idxrec].sim != idxsim ){ cout << "crash and burn :  " << idxrec << "," << recvmatch[idxrec].sim << "," << idxsim << endl;}
      assert(recvmatch[idxrec].sim == idxsim);

      std::ostringstream s;
      s << setw(9) << setprecision(4) << fixed << simEvt[idxsim].z << "," << setw(7) << setprecision(3) << fixed
        << simEvt[idxsim].t << " " << setw(1) << "[" << setw(3) << idxsim << "]"
	<< setw(4) << simEvt[idxsim].rtk.size() << " "
        << formatMatchList(simEvt[idxsim].wnt, 3, true) << setw(9) << setprecision(4) << fixed
        << recVtxs->at(idxrec).z() << "," << setw(7) << setprecision(3) << fixed << recVtxs->at(idxrec).t() << setw(1)
        << "(" << setw(4) << idxrec << ") "
	<< setprecision(1) << fixed << setw(5) << recVtxs->at(idxrec).ndof() << " "
	<< formatMatchList(recvmatch[idxrec].wnt, 3, false);

      dumped_rec.push_back(idxrec);
      row.push_back(make_pair(recVtxs->at(idxrec).z(), s.str()));
    }
  }

  // what have we missed?
  // unmatched sim
  for (unsigned int idxrec = 0; idxrec < recVtxs->size(); idxrec++) {
    if (recVtxs->at(idxrec).isFake()) continue;
    if (std::find(dumped_rec.begin(), dumped_rec.end(), idxrec) == dumped_rec.end()){
      unsigned int idxsim = recvmatch[idxrec].sim;
      
      cout << "missed in printMatchingSummary  "  << idxrec <<  "  " << recVtxs->at(idxrec).z() << "  match=" << recvmatch[idxrec].sim;
      std::ostringstream s;
      s << "@" << setw(61) << " " << setw(9) << setprecision(4) << fixed << recVtxs->at(idxrec).z() << "," << setw(7)
        << setprecision(3) << fixed << recVtxs->at(idxrec).t() << " " << setw(1) << "(" << setw(3) << idxrec << ") "
        << formatMatchList(recvmatch[idxrec].wnt, 3, false);
      row.push_back(make_pair(recVtxs->at(idxrec).z(), s.str()));
      if(idxsim != NOT_MATCHED){
	std::ostringstream s;
	s << "@" <<setw(8) << setprecision(4) << fixed << simEvt[idxsim].z << "," << setw(7) << setprecision(3) << fixed
        << simEvt[idxsim].t << " " << setw(1) << "[" << setw(3) << idxsim << "]"
	  << setw(4) << simEvt[idxsim].rtk.size() << " "
        << formatMatchList(simEvt[idxsim].wnt, 3, true);
	row.push_back(make_pair(simEvt[idxsim].z, s.str()));
	cout << s.str();
	cout << "  matched to  " << simEvt[idxsim].rec;
      }
      cout <<endl;
    }
  } 
  //cout << "printMatchingSummary rows=" << row.size() << endl;
  stable_sort(row.begin(), row.end());
  for (auto r : row) {
    cout << r.second << endl;
  }
}

/***************************************************************************************/







/***************************************************************************************/
void PrimaryVertexAnalyzer4PU::printEventSummary(std::map<std::string, TH1*>& h,
                                                 const reco::VertexCollection* recVtxs,
                                                 Tracks& tracks,
                                                 vector<SimEvent>& simEvt,
                                                 std::vector<RSmatch>& recvmatch,
                                                 const string message) {
  // make a readable summary of the vertex finding if the TrackingParticles are availabe
  if (simEvt.size() == 0) {
    return;
  }

  if (eventSummaryCounter_++ > nEventSummary_){
    report_counted(message + " count limit exceeded", 1);
    return;
  }
  // sort vertices in z ... for nicer printout

  vector<pair<double, unsigned int>> zrecv;
  for (unsigned int idx = 0; idx < recVtxs->size(); idx++) {
    if (recVtxs->at(idx).isFake()) continue;
    zrecv.push_back(make_pair(recVtxs->at(idx).z(), idx));
  }
  stable_sort(zrecv.begin(), zrecv.end(), lt);

  // same for simulated vertices
  vector<pair<double, unsigned int>> zsimv;
  for (unsigned int idx = 0; idx < simEvt.size(); idx++) {
    zsimv.push_back(make_pair(simEvt[idx].z, idx));
  }
  stable_sort(zsimv.begin(), zsimv.end(), lt);

  double zmin = -50;
  double zmax = 50;
  bool fillHistograms = true;

  if (message == "signal split") {
    zmin = simEvt[0].z - 1.0;
    zmax = simEvt[0].z + 1.0;
    fillHistograms = false;
  }

  cout << endl << "printEventSummary  (simEvt)" << endl;
  cout << "---------------------------" << endl;
  cout << "event counter = " << eventcounter_ << "   " << message << endl;
  cout << "run " << run_ << "    event " << event_ << endl;
  cout << "---------------------------" << endl;

  map<unsigned int, unsigned int> rvmatch;  // reco vertex matched to sim vertex  (sim to rec)
  map<unsigned int, unsigned int> svmatch;  // sim vertex matched to rec vertex  (rec to sim)
  map<unsigned int, double> nmatch;         // highest number of truth-matched tracks of ev found in a recvtx
  map<unsigned int, double> wnmatch;        // highest sum of weights of truth-matched tracks of ev found in a recvtx
  map<unsigned int, double> purity;  // highest purity of a rec vtx (i.e. highest number of tracks from the same simvtx)
  map<unsigned int, double> wpurity;  // same for the sum of weights

  for (vector<pair<double, unsigned int>>::iterator itrec = zrecv.begin(); itrec != zrecv.end(); itrec++) {
    svmatch[itrec->second] = NOT_MATCHED;
    purity[itrec->second] = 0.;
    wpurity[itrec->second] = 0.;
  }

  if (zrecv.size() > 20) {
    cout << " full dump dropped because we have more than 20 vertices   nrec = " << zrecv.size() << endl;

    /* ?????

    // truthMatchedVertexTracks[irecvtx]=list of rec tracks that vertex for which we have truth matched simtracks
    // (not necessarily all from the same simvertex)
    map<unsigned int, int> truthMatchedVertexTracks;

    for (vector<pair<double, unsigned int>>::iterator itrec = zrecv.begin(); itrec != zrecv.end(); itrec++) {
      truthMatchedVertexTracks[itrec->second] =
          getTruthMatchedVertexTracks(recVtxs->at(itrec->second)).size();  // FIXME replace consistently
    }
    */
    
    
    for (vector<pair<double, unsigned int>>::iterator itsim = zsimv.begin(); itsim != zsimv.end(); itsim++) {
      // itsim->first = z of the simvx
      // itsim->second= index of the simvtx
      if ((itsim->first < zmin) || (itsim->first > zmax))
        continue;
      SimEvent* ev = &(simEvt[itsim->second]);
      rvmatch[itsim->second] = NOT_MATCHED;

      nmatch[itsim->second] = 0;   // highest number of truth-matched tracks of ev found in a recvtx
      wnmatch[itsim->second] = 0;  // highest sum of weights of truth-matched tracks of ev found in a recvtx
      //	double matchpurity=0;//,matchwpurity=0;//,matchpurity2=0;

      // compare this sim vertex to all recvertices:
      for (vector<pair<double, unsigned int>>::iterator itrec = zrecv.begin(); itrec != zrecv.end(); itrec++) {
        // itrec->first  = z coordinate of the recvtx
        // itrec->second = index of the recvtx
        unsigned int irecvtx = itrec->second;
        const reco::Vertex* v = &(recVtxs->at(irecvtx));

        // count tracks found in both, sim and rec
        double n = 0, wt = 0;
        for (trackit_t tv = v->tracks_begin(); tv != v->tracks_end(); tv++) {
          if (ev->hasRecoTrack(*tv)) {
            n++;
            wt += v->trackWeight(*tv);
          }
        }

        // consider for matching if reasonably close in z
        double deltaz = fabs(v->z() - itsim->first);
        if ((deltaz < (5 * v->zError())) && (deltaz < 0.5)) {
          // match by number of tracks
          if (n > nmatch[itsim->second]) {
            nmatch[itsim->second] = n;
            rvmatch[itsim->second] = itrec->second;
          }

          if (n > purity[itrec->second]) {
            purity[itrec->second] = n;
            svmatch[itrec->second] = itsim->second;
          }
        }

        // match by weight
        if (wt > wnmatch[itrec->second]) {
          wnmatch[itrec->second] = wt;
        }

        if (wt > wpurity[itrec->second]) {
          wpurity[itrec->second] = wt;
        }

      }  // end of reco vertex loop
    }
  } else {
    cout << " z[cm]       rec -->   ";
    cout.precision(4);
    for (vector<pair<double, unsigned int>>::iterator itrec = zrecv.begin(); itrec != zrecv.end(); itrec++) {
      if ((itrec->first < zmin) || (itrec->first > zmax))
        continue;
      cout << setw(8) << fixed << itrec->first;
    }
    cout << endl;

    cout << "                        ";
    for (vector<pair<double, unsigned int>>::iterator itrec = zrecv.begin(); itrec != zrecv.end(); itrec++) {
      if ((itrec->first < zmin) || (itrec->first > zmax))
        continue;
      cout << setw(7) << fixed << recVtxs->at(itrec->second).tracksSize();
      if (itrec->second == 0) {
        cout << "*";
      } else {
        cout << " ";
      }
    }
    cout << "   rec tracks" << endl;

    if(tracking_truth_available_){
    
    // truthMatchedVertexTracks[irecvtx]=list of rec tracks that vertex for which we have truth matched simtracks
    // (not necessarily all from the same simvertex)
    map<unsigned int, int> truthMatchedVertexTracks;
    
    for (vector<pair<double, unsigned int>>::iterator itrec = zrecv.begin(); itrec != zrecv.end(); itrec++) {
      truthMatchedVertexTracks[itrec->second] =
	getTruthMatchedVertexTracks(recVtxs->at(itrec->second)).size();  // FIXME replace consistently
    }
    
    cout << "                        ";
    for (vector<pair<double, unsigned int>>::iterator itrec = zrecv.begin(); itrec != zrecv.end(); itrec++) {
      if ((itrec->first < zmin) || (itrec->first > zmax))
	continue;
      cout << setw(7) << fixed << recvmatch[itrec->second].sumwnt
	   << " ";  // may contain tracks that have weight < 0.5 after the vertex fit
    }
    cout << "   truth matched " << endl;
    
    cout << "sim ------- trk  prim ----" << endl;
    
    for (vector<pair<double, unsigned int>>::iterator itsim = zsimv.begin(); itsim != zsimv.end(); itsim++) {
      // itsim->first = z of the simvx
      // itsim->second= index of the simvtx
      if ((itsim->first < zmin) || (itsim->first > zmax))
        continue;
      SimEvent* ev = &(simEvt[itsim->second]);
      rvmatch[itsim->second] = NOT_MATCHED;
      
      cout.precision(4);
      if (itsim->second == 0) {
        cout << setw(8) << fixed << ev->z << ")*" << setw(5) << ev->rtk.size() << setw(5) << ev->rtkprim.size() << "  | ";
      } else {
        cout << setw(8) << fixed << ev->z << ") " << setw(5) << ev->rtk.size() << setw(5) << ev->rtkprim.size() << "  | ";
      }

      nmatch[itsim->second] = 0;   // highest number of truth-matched tracks of ev found in a recvtx
      wnmatch[itsim->second] = 0;  // highest sum of weights of truth-matched tracks of ev found in a recvtx
      double matchpurity = 0;      //,matchwpurity=0;//,matchpurity2=0;

      // compare this sim vertex to all recvertices:
      for (vector<pair<double, unsigned int>>::iterator itrec = zrecv.begin(); itrec != zrecv.end(); itrec++) {
        // itrec->first  = z coordinate of the recvtx
        // itrec->second = index of the recvtx
        unsigned int irecvtx = itrec->second;
        const reco::Vertex* v = &(recVtxs->at(irecvtx));

        // count tracks found in both, sim and rec
        double n = 0, wt = 0;
        for (trackit_t tv = v->tracks_begin(); tv != v->tracks_end(); tv++) {
          if (ev->hasRecoTrack(*tv)) {
            n++;
            wt += v->trackWeight(*tv);
          }
        }

        if ((itrec->first >= zmin) && (itrec->first <= zmax)) {
          if (n > 0) {
            cout << setw(7) << int(n) << " ";
          } else {
            cout << "        ";
          }
        }

        // consider for matching if reasonably close in z
        double deltaz = fabs(v->z() - itsim->first);
        if ((deltaz < (5 * v->zError())) && (deltaz < 0.5)) {
          // match by number of tracks
          if (n > nmatch[itsim->second]) {
            nmatch[itsim->second] = n;
            rvmatch[itsim->second] = itrec->second;
            //matchpurity2=matchpurity;
            matchpurity = n / truthMatchedVertexTracks[itrec->second];
            //matchwpurity=wt/truthMatchedVertexTracks[itrec->second];
          }

          if (n > purity[itrec->second]) {
            purity[itrec->second] = n;
            svmatch[itrec->second] = itsim->second;
          }
        }

        // match by weight
        if (wt > wnmatch[itrec->second]) {
          wnmatch[itrec->second] = wt;
        }

        if (wt > wpurity[itrec->second]) {
          wpurity[itrec->second] = wt;
        }

      }  // end of reco vertex loop

      // summary of this sim vertex:
      cout << "  | " << setw(1) << ev->nwosmatch << "," << setw(1) << ev->nwntmatch;
      cout << " | ";
      if (nmatch[itsim->second] > 0) {
        if (matchpurity >= 0.5) {
          cout << "found  ";
        } else if (matchpurity >= 0.3) {
          cout << "ugly   ";
        } else {
          cout << "merged ";
        }
        cout << endl;
      } else {
        // sim vertex not matched to any rec vertex
        if (ev->rtk.size() == 0) {
          cout << "invisible" << endl;
        } else if (ev->rtk.size() == 1) {
          cout << "single track " << endl;
        } else if (ev->rec != NOT_MATCHED) {
          cout << "poor,  quality " << ev->matchQuality << "  for zrec " << recVtxs->at(ev->rec).z() << endl;
        } else {
          cout << "lost " << endl;
        }
      }
    }
    cout << "---------------------------" << endl;

    //  the purity of the reconstructed vertex
    cout << "               purity   ";
    for (vector<pair<double, unsigned int>>::iterator itrec = zrecv.begin(); itrec != zrecv.end(); itrec++) {
      if ((itrec->first < zmin) || (itrec->first > zmax))
        continue;
      cout << setw(7) << fixed << purity[itrec->second] / truthMatchedVertexTracks[itrec->second];
      if (itrec->second == 0) {
        cout << "*";
      } else {
        cout << " ";
      }  // flag the tagged vertex
    }
    cout << endl;

    //  test
    cout << "               purity*  ";
    for (vector<pair<double, unsigned int>>::iterator itrec = zrecv.begin(); itrec != zrecv.end(); itrec++) {
      if ((itrec->first < zmin) || (itrec->first > zmax))
        continue;
      if (recvmatch[itrec->second].sumwnt > 0) {
        cout << setw(7) << fixed << recvmatch[itrec->second].maxwnt / float(recvmatch[itrec->second].sumwnt);
      } else {
        cout << "   -   ";
      }
      if (itrec->second == 0) {
        cout << "*";
      } else {
        cout << " ";
      }  // flag the tagged vertex
    }
    cout << endl;

    // rec vertex classification: the good, the bad, and the ugly
    cout << "                     |   ";
    for (vector<pair<double, unsigned int>>::iterator itrec = zrecv.begin(); itrec != zrecv.end(); itrec++) {
      if ((itrec->first < zmin) || (itrec->first > zmax))
        continue;

      if ((svmatch[itrec->second] != NOT_MATCHED) && (rvmatch[svmatch[itrec->second]] == itrec->second)) {
        if ((purity[itrec->second] / truthMatchedVertexTracks[itrec->second]) >= 0.5) {
          cout << "   ok   ";
        } else {
          cout << "  ugly  ";
        }
      } else {
        if (recvmatch[itrec->second].split_from() >= 0) {
          cout << " split  ";
        } else {
          cout << "  junk  ";
        }
      }
    }
    cout << endl;

  } // !tracking_truth_available_

  
    // wos matching // FIXME this is duplicated, can we re-use tpmatch here?
    cout << "        matched vtx z| ";
    for (vector<pair<double, unsigned int>>::iterator itrec = zrecv.begin(); itrec != zrecv.end(); itrec++) {
      if ((itrec->first < zmin) || (itrec->first > zmax))
        continue;
      if (recvmatch[itrec->second].maxwos > 0) {
        cout << setw(8) << fixed << simEvt[recvmatch[itrec->second].wosmatch].z;
      } else if (recvmatch[itrec->second].maxwnt > 0) {
        cout << setw(8) << fixed << simEvt[recvmatch[itrec->second].wntmatch].z;
      } else {
        cout << "    -    ";
      }
    }
    cout << endl;
    cout << "---------------------------" << endl;

  }  // end of the sim vs rec table dump

  // print recvertices that were not successfully matched
  cout << "list of junk vertices" << endl;
  for (vector<pair<double, unsigned int>>::iterator itrec = zrecv.begin(); itrec != zrecv.end(); itrec++) {
    if ((itrec->first < zmin) || (itrec->first > zmax))
      continue;

    unsigned int iv = itrec->second;

    if ((svmatch[iv] == NOT_MATCHED) || (rvmatch[svmatch[iv]] != iv)) {
      cout << "zrec=" << setw(9) << fixed << setprecision(4) << itrec->first << "  ";
      if (recvmatch[iv].maxwos > 0) {
        cout << "maxwos = " << setw(10) << fixed << setprecision(1) << recvmatch[iv].maxwos << " / " << setw(10)
             << fixed << setprecision(1) << recvmatch[iv].sumwos << "   zwosmatch =" << setw(9) << fixed
             << setprecision(4) << simEvt[recvmatch[iv].wosmatch].z << "  ";
      }
      if (recvmatch[iv].maxwnt > 0) {
        cout << "maxnt = " << setw(3) << setprecision(0) << recvmatch[iv].maxwnt << " /" << setw(3) << setprecision(0)
             << recvmatch[iv].sumwnt << "  z=" << setw(8) << fixed << setprecision(4)
             << simEvt[recvmatch[iv].wntmatch].z << "  ";
      }
      cout << "  fake=" << recvmatch[iv].is_fake();
      cout << "  qual=" << recvmatch[iv].matchQuality;
      cout << "  split=" << recvmatch[iv].split_from();
      cout << endl;
    }
  }

  if (!fillHistograms)
    return;

  // FIXME this code needs an overhaul
  // list problematic tracks
  for (vector<pair<double, unsigned int>>::iterator itsim = zsimv.begin(); itsim != zsimv.end(); itsim++) {
    SimEvent* ev = &(simEvt[itsim->second]);

    //for (vector<TransientTrack>::iterator te = ev->tk.begin(); te != ev->tk.end(); te++) {
    for (auto tk : ev->rtk){
      const reco::Track& RTe = *tk.trk;

      unsigned int ivassign = NOT_ASSIGNED;  // will become the index of the vertex to which a track was assigned

      for (vector<pair<double, unsigned int>>::iterator itrec = zrecv.begin(); itrec != zrecv.end(); itrec++) {
        const reco::Vertex* v = &(recVtxs->at(itrec->second));

        for (trackit_t tv = v->tracks_begin(); tv != v->tracks_end(); tv++) {
          const reco::Track& RTv = *(tv->get());
          if (RTe.vz() == RTv.vz()) {
            ivassign = itrec->second;
          }
        }
      }

      double tantheta = tan(tk.theta);//(te->stateAtBeamLine().trackStateAtPCA()).momentum().theta());
      //reco::BeamSpot beamspot = (te->stateAtBeamLine()).beamSpot();
      reco::BeamSpot beamspot = (tk.tt->stateAtBeamLine()).beamSpot();
      double dz2 = pow(RTe.dzError(), 2) + pow(beamspot.BeamWidthX() / tantheta, 2);

      if ((ivassign == ev->rec) && (ev->matchQuality > 0)) {
        Fill(h, "correctlyassigned", RTe.eta(), RTe.pt());
        Fill(h, "correctlyassignedlogpt", RTe.eta(), log(RTe.pt()) / log(10.));
        Fill(h, "ptcat", RTe.pt());
        Fill(h, "etacat", RTe.eta());
        if (RTe.pt() > 2)
          Fill(h, "etacatpt2", RTe.eta());
        Fill(h, "phicat", RTe.phi());
        Fill(h, "dzcat", sqrt(dz2));
      } else {
        Fill(h, "misassigned", RTe.eta(), RTe.pt());
        Fill(h, "misassignedlogpt", RTe.eta(), log(RTe.pt()) / log(10.));
        Fill(h, "ptmis", RTe.pt());
        Fill(h, "etamis", RTe.eta());
        if (RTe.pt() > 2)
          Fill(h, "etamispt2", RTe.eta());
        Fill(h, "phimis", RTe.phi());
        Fill(h, "dzmis", sqrt(dz2));
      }
    }  // next simvertex-track

  }  //next simvertex
}
/***************************************************************************************/

/***************************************************************************************/
void PrimaryVertexAnalyzer4PU::ihmatch(
    const reco::VertexCollection* recVtxs,
    std::vector<SimEvent>& simEvt,
    std::vector<RSmatch>& recvmatch  // vector of RSmatch objects, parallel to recvertices, filled in tpmatch
    )
/***************************************************************************************/
{
  //implement the invisible hand algorithm
  // J.J. Kosowsky and A.L.Yuille,  Neural Networks, Vol 7,  No.3 pp 477-490, 1994

  unsigned int Nr = recVtxs->size();  // objects, a
  unsigned int Ns = simEvt.size();    // persons, i
  unsigned int N = max(Nr, Ns);
  // association benefits A[i][a] = recvmatch[a].wnt[i]  if it exists
  //  auto A = new double[N][N];        // A[i][a] association benefit ( = -association Energy)
  auto A = [N, Ns, Nr](auto& r2s, unsigned int a, unsigned int i) {
    if (a >= Nr)
      return 0.;
    if (i >= Ns)
      return 0.;
    if (r2s[a].wnt.find(i) == r2s[a].wnt.end()) {
      return 0.;
    } else {
      return r2s[a].wnt[i];
    }
  };

  //double Tmin = 0.01;
  double T = 1.;
  //  double coolingfactor = 0.1;

  std::vector<double> P(N);
  for (unsigned int a = 0; a < Nr; a++) {
    P[a] = 1. / Nr;
  }
  std::vector<double> s(N);
  std::vector<double> rho(N);  // = exp(beta*P[a])
  for (unsigned int a = 0; a < Nr; a++) {
    rho[a] = exp(1. / (T * Nr));
  }

  //  while( )
  //    {
  double beta = 1. / T;
  for (unsigned int m = 0; m < 10; m++) {
    for (unsigned int a = 0; a < Nr; a++) {
      s[a] = 0;
      for (unsigned int i = 0; i < Ns; i++) {
        double Zi = 0;
        for (unsigned int b = 0; b < Nr; b++) {
          Zi += rho[b] * exp(beta * A(recvmatch, i, b));
        }
        s[a] += exp(beta * A(recvmatch, i, a)) / Zi;
      }
    }

    for (unsigned int a = 0; a < Nr; a++) {
      rho[a] = rho[a] * s[a];
    }
  }

  //T = T * coolingfactor;
  for (unsigned int a = 0; a < Nr; a++) {
    double Sa = 0;
    unsigned int ia = 10000;
    for (unsigned int i = 0; i < Ns; i++) {
      double Zi = 0;
      for (unsigned int b = 0; b < Nr; b++) {
        Zi += rho[b] * exp(beta * A(recvmatch, i, b));
      }
      double Sia = rho[a] * exp(beta * A(recvmatch, i, a)) / Zi;
      if (Sia > Sa) {
        Sa = Sia;
        ia = i;
      }
    }

    cout << "a = " << a << "  rho_a = " << rho[a] << "  i[a] = " << ia << "  Sa=" << Sa << endl;
  }

  //    }
}

/***************************************************************************************/
std::vector<PrimaryVertexAnalyzer4PU::RSmatch> PrimaryVertexAnalyzer4PU::tpmatch(const reco::VertexCollection* recVtxs,
                                                                                 std::vector<SimEvent>& simEvt,
                                                                                 Tracks& tracks)
/***************************************************************************************/
{
  // collects truth-matched track information for matching rec and sim vertices
  // two different matches are used to determine the dominant sim vertex in a recvertex
  // wos match  = sum of "weight-over-sigma**2"
  // ntrk match = number of tracks (implicit with weight > 0.5), weighted with pt (cut-off at 1Gev)
  // this information is filled into both, recvmatch objects(|| recvertices) and simevts
  // a rec-to-sim match is made by identifying the the dominant sim vertex
  //   (wosmatch and wntmatch)
  // note: this is not a one-to-one match, a simvertex can dominate multiple recvtxs
  // simEvt.nwosmatch  and simEvt.nwtmatch keep count the number of recvertices dominated by a simvertex

  bool DEBUG = false;
  unsigned int ivdebug = 0;  // rec
  unsigned int ievdebug = 0; // sim
  if(DEBUG){ cout << "tpmatch  run:event = " << run_ << ":" << event_ << "  ivdebug=" << ivdebug << " ievdebug=" << ievdebug << endl;}

  std::vector<RSmatch> recvmatch;  // vector of RSmatch objects, parallel to recvertices
  recvmatch.clear();

  // clear old sim->rec pointers, only valid for one vertex collection
  for (vector<SimEvent>::iterator ev = simEvt.begin(); ev != simEvt.end(); ev++) {
    ev->clear_matching_info();
  }

  for (unsigned int iv = 0; iv < recVtxs->size(); iv++) {
    const reco::Vertex* v = &(recVtxs->at(iv));
    RSmatch M;  // holds information about match-candidate sim vertices
    if (v->isFake()) {
      recvmatch.push_back(M); //just a filler to keep the lists parallel
      continue;
    }
    unsigned int iev = 0;  // simevent index
    for (auto ev = simEvt.begin(); ev != simEvt.end(); ev++) {
      double evwos = 0;       // weight over sigma**2  of simEvt ev in the current recvtx
      double evwnt = 0;       // weighted number of tracks
      unsigned int evnt = 0;  // number of tracks from simEvt ev in the current recvtx

      assert(iev == ev->index);

      for (trackit_t tv = v->tracks_begin(); tv != v->tracks_end(); tv++) {
	if(v->trackWeight(*tv) < 0.5) continue;
	
        const RecoTrack t = tracks.from_key(tv->key());

        if (t.matched && (*t.simEvt == *ev)) {
          double dz2_beam = pow(vertexBeamSpot_.BeamWidthX() * cos(t.phi) / tan(t.theta), 2) +
                            pow(vertexBeamSpot_.BeamWidthY() * sin(t.phi) / tan(t.theta), 2);
          double dz2 = pow(t.dz, 2) + dz2_beam + pow(0.0020, 2); // added 20 um, some tracks have crazy small resolutions
          double wos = v->trackWeight(*tv) / dz2;
          if (f4D_ && t.has_timing && (t.dt > 0)) {
            wos = wos / erf(t.dt / sigmaT_);
          }
          double wnt = v->trackWeight(*tv) * min(t.pt, 1.0);
          M.addTrack(iev, wos, wnt);   // fill track(wos)  rec vtx <- sim vtx, increments sumnt and sumwos
          ev->addTrack(iv, wos, wnt);  // fill track(wos)  sim vtx -> rec vtx
          evwos += wos;
          evwnt += wnt;
          evnt++;
	  if (DEBUG &&((iev==0) && ((iv==0) || (iv==2)))){
	    cout << " iev=" << iev << "  iv=" << iv << "    wos=" << wos << "   evwos= " << evwos << "   dz2= "  << dz2 << "  w=" <<  v->trackWeight(*tv) << endl;
	  }
	  
        }
      }
      if(DEBUG && (iv == ivdebug) && (iev == ievdebug)){
	cout << "tpmatch rec " << iv << "  sim " << iev << "  evwos = " << evwos<< " M.maxwos=" << M.maxwos << "   evnt=" << evnt << endl;
      }
      if(DEBUG && (iv == 2) && (iev == ievdebug)){
	cout << "tpmatch rec " << iv << "  sim " << iev << "  evwos = " << evwos<< " M.maxwos=" << M.maxwos << "   evnt=" << evnt << endl;
      }
      // require 2 tracks for a wos-match
      if ((evwos > 0) && (evwos > M.maxwos) && (evnt > 1)) {
        M.wosmatch = iev;
        M.maxwos = evwos;
        M.maxwosnt = evnt;
      }

      // weighted track counting match, require at least one track
      if ((evnt > 0) && (evwnt > M.maxwnt)) {
        M.wntmatch = iev;
        M.maxwnt = evwnt;
      }

      iev++;
      
    } // end of simvertex loop

    // now wosmatch  holds the index of the dominant simvertex for this recvtx and maxwos the wos value
    if(DEBUG && (iv==ivdebug)){
      cout << "tpmatch: recvtx " << ivdebug << ", M.maxwos=" << M.maxwos << "   M.wosmatch="  << M.wosmatch << endl;
    }
    if(DEBUG && (iv==2)){
      cout << "tpmatch: recvtx " << ivdebug << ", M.maxwos=" << M.maxwos << "   M.wosmatch="  << M.wosmatch << endl;
    }
    if (M.maxwos > 0) {
      simEvt.at(M.wosmatch).nwosmatch++;  // count the recvertices dominated by a simvertex
      simEvt.at(M.wosmatch).wos_dominated_recv.push_back(iv);
      assert(iv < recVtxs->size());
    }

    if (M.maxwnt > 0) {
      simEvt.at(M.wntmatch).nwntmatch++;  // count the recvertices dominated by a simvertex
      //simEvt.at(M.wntmatch).wnt_dominated_recv.push_back(iv);
    }

    recvmatch.push_back(M);

  }  // end of recvertex loop

  return recvmatch;
}




void PrimaryVertexAnalyzer4PU::wos_match(std::vector<RSmatch> & recvmatch,
                               const reco::VertexCollection* recVtxs,
                               std::vector<SimEvent>& simEvt,
                               Tracks& tracks)
{
  bool DEBUG = false;
  //make rec <-> sim associations using the information filled by tpmatch
  // in contrast to tpmatch, this tries to make a one-to-one match
  //
  // first make only associations when a recvertex is dominated by a simvertex that is not
  // already assigned to another rec vertex
  // dominated = either by wos or (weighted) nt
  // notes:
  //  * unless a rec vertex consists entirely of fake tracks, it must be dominated
  //    by a simvertex
  //  * a rec-vertex that does not get sim vertex in this way is either a split-vertex
  //    in the sense that the dominating sim vertex also dominates another rec vertex
  //    or it consists entirely of unmatched tracks
  // * a merged vertex would be a simvertex that shares (enough? the majority of its?) tracks with a recvertex
  //   but another sim vertex dominates that vertex
  //
  // finally match rec vertices that are not necessarily dominated by a sim vertex,
  // but the dominating sim-vertex has been matched already and there is a simvertex that
  // really wants to have that rec-vertex. This corresponds to a "small" vertex close to "big" one
  // where the small vertex is contaminated a lot by the big one, but still recognizable
  //

  // reset
  if(DEBUG) {cout << "wos_match  recVtxs->size() = " << recVtxs->size()  << endl;}
  for (unsigned int iv = 0; iv < recVtxs->size(); iv++) {
    recvmatch[iv].sim = NOT_MATCHED;
    recvmatch[iv].matchQuality = 0;
  }
  for (unsigned int iev = 0; iev < simEvt.size(); iev++){
    simEvt[iev].rec = NOT_MATCHED;
    simEvt[iev].matchQuality = 0;
  }


  if(DEBUG){
    cout << "DEBUG: wos_match nwosmatch[0] = " << simEvt.at(0).nwosmatch << endl;
  }
  // when two or more rec vertices are dominated by the same sim vertex,
  // assign the rec vertex that got more from that sim
  for (unsigned int rank = 1; rank < 8; rank++)
    {
      for (unsigned int iev = 0; iev < simEvt.size(); iev++)
	{
	  if (simEvt.at(iev).rec != NOT_MATCHED) continue;
	  if (simEvt.at(iev).nwosmatch == 0) continue;
	  if (simEvt.at(iev).nwosmatch > rank) continue;
	  
	  unsigned int iv = NOT_MATCHED;
	      
	  for (unsigned int k = 0; k < simEvt.at(iev).wos_dominated_recv.size(); k++)
	    {
	      unsigned int rec = simEvt.at(iev).wos_dominated_recv.at(k);
	      if (recvmatch.at(rec).sim != NOT_MATCHED) continue; // already matched
	      if (fabs(simEvt.at(iev).z - recVtxs->at(rec).z()) > zWosMatchMax_) continue;// insanely far away
	      if ( (iv == NOT_MATCHED) || (simEvt.at(iev).wos.at(rec) > simEvt.at(iev).wos.at(iv)) )
		{
		  iv = rec;
		}
	    }
	  
	  if (iv != NOT_MATCHED)
	    {
	      recvmatch.at(iv).sim = iev;
	      simEvt.at(iev).rec = iv;
	      recvmatch.at(iv).matchQuality = rank;
	      simEvt.at(iev).matchQuality = rank;
	      if(DEBUG) {cout << "wos_match : match made  [" << iev << "]  <-> (" << iv << ")   rank=" << rank << endl;}
	    }
	}
    }
      

  // by now we have exhausted the rec vertices that are dominated by an unmatched simvertex
  // have we?
  for(unsigned int iv = 0; iv < recVtxs->size(); iv++) {
    if ((recvmatch.at(iv).sim == NOT_MATCHED) && (recvmatch.at(iv).maxwos >0)){
      if (simEvt.at(recvmatch.at(iv).wosmatch).rec == NOT_MATCHED){
	cout << "wos_match :  unmatched [" << recvmatch[iv].wosmatch << "]   dominantes unmatched  ("<< iv <<")   ????"<<endl;
	reportEvent("wos_match :  unmatched dominates matched");
	}
      }
    }


  // give vertices a chance that have a lot of overlap, but are still recognizably
  // caused by a specific simvertex (without being classified as dominating)
  // like a small peak sitting on the flank of a larger nearby peak

  unsigned int ntry = 0;
  while (ntry++ < 10)
    {
      unsigned nmatch = 0;
      for (unsigned int sim = 0; sim < simEvt.size(); sim++)
	{
	  if ((simEvt.at(sim).rec != NOT_MATCHED) || (simEvt.at(sim).wos.size() == 0))
	    continue;
	
	  // ok, single simvertex sim, who is your your favorite rec vertex?
	  unsigned int rec = NOT_MATCHED;
	  for (auto rv : simEvt.at(sim).wos) {
	    if ((rec == NOT_MATCHED) || (rv.second > simEvt.at(sim).wos.at(rec))) {
	      rec = rv.first;
	    }
	  }

	  if (rec == NOT_MATCHED){
	    cout << "wos_match : last hope failed for (" << rec << ")" << endl;
	    // second chance if wos didn't work?
	    for (auto rv : simEvt.at(sim).wnt) {
	      if ((rec == NOT_MATCHED) || (rv.second > simEvt.at(sim).wnt.at(rec))) {
		rec = rv.first;
	      }
	    }
	  }
	
	if (rec == NOT_MATCHED)
	  continue;  // should not happen
	if (recvmatch.at(rec).sim != NOT_MATCHED)
	  continue;  // already gone
	
	// do you, recvertex rec, take this simvertex sim as your lawful wedded truthmatch?
	unsigned int rec2sim = NOT_MATCHED;
	for (auto sv : recvmatch.at(rec).wos) {
	  if (simEvt.at(sv.first).rec != NOT_MATCHED)
	    continue;  // already used
	  if ((rec2sim == NOT_MATCHED) || (sv.second > recvmatch.at(rec).wos.at(rec2sim))) {
	    rec2sim = sv.first;
	  }
	}
	
	if (sim == rec2sim) {
	  // I do
	  recvmatch.at(rec).sim = sim;
	  recvmatch.at(rec).matchQuality = 8;
	  simEvt.at(sim).rec = rec;
	  simEvt.at(sim).matchQuality = 8;
	  //cout << "wos_match : late match  [" << sim << "]   <--> (" << rec << ")" <<  "    try=" << ntry << endl;
	  nmatch++;
	}
      }  //sim loop
      // cout << "wos_match late matches " << nmatch <<  "  in try " << ntry << endl;
      if (nmatch == 0) {
	break;
      }
    }  // ntry
  
  if (DEBUG) {
    for (unsigned int iv = 0; iv < recVtxs->size(); iv++) {
      if (iv == 0) {
        dumpThisEvent_ = true;
        cout << "wos_match    recvtx = " << setw(4) << iv << "   (sim=" << recvmatch[iv].sim << ")" << endl;
        cout << " splitfrom= " << recvmatch[iv].split_from() << endl;
        cout << " maxwos = " << recvmatch[iv].maxwos << endl;
        cout << " sumwos = " << recvmatch[iv].sumwos << endl;

        for (auto w = recvmatch[iv].wos.begin(); w != recvmatch[iv].wos.end(); w++) {
          cout << "matching (wos)  simevent = " << setw(4) << w->first << "  zsim = " << setw(8) << fixed << setprecision(4)
               << simEvt.at(w->first).z << "  wos=" << setw(10) << w->second << endl;
        }

        for (unsigned int iev = 0; iev < simEvt.size(); iev++) {
          for (auto rv = simEvt.at(iev).wos.begin(); rv != simEvt.at(iev).wos.end(); rv++) {
            if (rv->first == iv) {
              cout << "matching (wos)  simevent " << iev << "    zsim= " << setw(8) << fixed << setprecision(4)
                   << simEvt.at(iev).z << "    zrec= " << setw(8) << fixed << setprecision(4)
                   << recVtxs->at((*rv).first).z() << "  nt=" << (*rv).second << endl;
            }
          }
        }

	for (unsigned int iev = 0; iev < simEvt.size(); iev++) {
          for (auto rv = simEvt.at(iev).wnt.begin(); rv != simEvt.at(iev).wnt.end(); rv++) {
            if (rv->first == iv) {
              cout << "matching (wnt)  simevent " << iev << "    zsim= " << setw(8) << fixed << setprecision(4)
                   << simEvt.at(iev).z << "    zrec= " << setw(8) << fixed << setprecision(4)
                   << recVtxs->at((*rv).first).z() << "  nt=" << (*rv).second << endl;
            }
          }
        }
      }
    }
  }  // <<<< debugging

  /*
  // last chance for vertices that are not matched by tracks
  // but are reasonably close and not splitters
  for(unsigned int iev=0; iev<simEvt.size(); iev++)
    {
      if( simEvt[iev].rec == NOT_MATCHED )
	{
	  for(unsigned int iv=0; iv<recVtxs->size(); iv++)
	    {
	      
	      if ((recvmatch[iv].sim == NOT_MATCHED) && (recvmatch[iv].split_from() < 0))
		{
		  double dz = simEvt.at(iev).z-recVtxs->at(iv).z();
		  if( ( fabs(dz) < 0.1 )
		      ||( fabs(dz) < (2*recVtxs->at(iv).zError()) )
		      )
		  {
		    recvmatch[iv].matchQuality = 10;
		    recvmatch[iv].sim = iev;
		    simEvt.at(iev).rec = iv;
		    simEvt.at(iev).matchQuality = 10;
		  }
		}
	    }
	}
    }
  */

  // for convenience, store some values from the matched rec vertex with the sim vertex
  for (unsigned int iev = 0; iev < simEvt.size(); iev++) {
    int iv = simEvt.at(iev).rec;
    if (iv != NOT_MATCHED) {
      simEvt.at(iev).ndof = recVtxs->at(iv).ndof();
      simEvt.at(iev).zrec = recVtxs->at(iv).z();
    } else {
      simEvt.at(iev).ndof = 0.;
      simEvt.at(iev).zrec = 1000.;
    }
  }

}
/***************************************************************************************/




void PrimaryVertexAnalyzer4PU::nwt_match(std::vector<RSmatch> & recvmatch,
                               const reco::VertexCollection* recVtxs,
                               std::vector<SimEvent>& simEvt,
                               Tracks& tracks)
{
  bool DEBUG = false;
  //make rec <-> sim associations using the information filled by tpmatch
  // in contrast to tpmatch, this tries to make a one-to-one match
  //
  // first make only associations when a recvertex is dominated by a simvertex that is not
  // already assigned to another rec vertex
  // dominated = either by wos or (weighted) nt
  // notes:
  //  * unless a rec vertex consists entirely of fake tracks, it must be dominated
  //    by a simvertex
  //  * a rec-vertex that does not get sim vertex in this way is either a split-vertex
  //    in the sense that the dominating sim vertex also dominates another rec vertex
  //    or it consists entirely of unmatched tracks
  // * a merged vertex would be a simvertex that shares (enough? the majority of its?) tracks with a recvertex
  //   but another sim vertex dominates that vertex
  //
  // finally match rec vertices that are not necessarily dominated by a sim vertex,
  // but the dominating sim-vertex has been matched already and there is a simvertex that
  // really wants to have that rec-vertex. This corresponds to a "small" vertex close to "big" one
  // where the small vertex is contaminated a lot by the big one, but still recognizable
  //
  
  //start with the unambigous (nwosmatch==1)
  for (unsigned int iv = 0; iv < recVtxs->size(); iv++) {
    if ((recvmatch[iv].sim == NOT_MATCHED) && (recvmatch[iv].maxwos > 0)) {
      unsigned int cand = recvmatch[iv].wosmatch;  // simEvt index
      if ((simEvt.at(cand).rec == NOT_MATCHED) &&
          (simEvt.at(cand).nwosmatch == 1) &&  // reject insanely misreconstructed vertices
          (fabs(simEvt.at(cand).z - recVtxs->at(iv).z()) < zWosMatchMax_)) {
        recvmatch[iv].matchQuality = 1;
        recvmatch[iv].sim = cand;
        simEvt.at(cand).rec = iv;
        simEvt.at(cand).matchQuality = 1;
      }
    }
  }

  // second chance, track counting match
  for (unsigned int iv = 0; iv < recVtxs->size(); iv++) {
    if ((recvmatch[iv].sim == NOT_MATCHED) && (recvmatch[iv].maxwnt > 0)) {
      unsigned int cand = recvmatch[iv].wntmatch;
      if ((simEvt.at(cand).rec == NOT_MATCHED) && (simEvt.at(cand).nwntmatch == 1) &&
          (fabs(simEvt.at(cand).z - recVtxs->at(iv).z()) < 1.0))  // reject insanely misreconstructed  vertices
      {
        recvmatch[iv].matchQuality = 2;
        recvmatch[iv].sim = cand;
        simEvt.at(cand).rec = iv;
        simEvt.at(cand).matchQuality = 2;
      }
    }
  }

  // when two or more rec vertices are dominated by one sim vertex,
  // assign the rec vertex that got more tracks (from that sim vertex) to the sim
  for (unsigned int iev = 0; iev < simEvt.size(); iev++) {
    if ((simEvt.at(iev).rec != NOT_MATCHED) && (simEvt.at(iev).nwosmatch > 1)) {
      double nt = 0;
      unsigned int iv = 0;

      for (auto rv : simEvt.at(iev).wnt) {
        if ((recvmatch[rv.first].sim == NOT_MATCHED) && (rv.second > nt) &&
            (fabs(simEvt.at(iev).z - recVtxs->at(rv.first).z()) < zWosMatchMax_)) {
          nt = rv.second;
          iv = rv.first;
        }
      }

      if (nt > 0) {
        if (DEBUG) {
          cout << "match made  rec " << iv << " ==  sim " << iev << endl;
        }
        recvmatch[iv].matchQuality = 5;
        recvmatch[iv].sim = iev;
        simEvt.at(iev).rec = iv;
        simEvt.at(iev).matchQuality = 5;
      }
    }
  }

  // when two or more rec vertices are dominated by one sim vertex, assign the rec vertex with more tracks (from that sim vertex) to the sim
  for (unsigned int iev = 0; iev < simEvt.size(); iev++) {
    if ((simEvt.at(iev).rec == NOT_MATCHED) && (simEvt.at(iev).nwntmatch > 1)) {
      double nt = 0;
      unsigned int iv = 0;
      for (auto rv : simEvt.at(iev).wnt) {
        if ((recvmatch[rv.first].sim == NOT_MATCHED) && (rv.second > nt) &&
            (fabs(simEvt.at(iev).z - recVtxs->at(rv.first).z()) < zWosMatchMax_)) {
          nt = rv.second;
          iv = rv.first;
        }
      }

      if (nt > 0) {
        recvmatch[iv].matchQuality = 6;
        recvmatch[iv].sim = iev;
        simEvt.at(iev).rec = iv;
        simEvt.at(iev).matchQuality = 6;
      } else {
        cout << "BUG" << endl;
      }
    }
  }

  // give vertices a chance that have a lot of overlap, but are still recognizably
  // caused by a specific simvertex
  // like a small peak sitting on the flank of a larger nearby peak

  unsigned int ntry = 0;
  while (ntry++ < 5) {
    unsigned nmatch = 0;
    for (unsigned int sim = 0; sim < simEvt.size(); sim++) {
      if ((simEvt.at(sim).rec != NOT_MATCHED) || (simEvt.at(sim).wos.size() == 0))
        continue;

      // ok, single simvertex sim, who is your your favorite rec vertex
      unsigned int rec = NOT_MATCHED;
      for (auto rv : simEvt.at(sim).wnt) {
        if ((rec == NOT_MATCHED) || (rv.second > simEvt.at(sim).wnt[rec])) {
          rec = rv.first;
        }
      }

      if (rec == NOT_MATCHED)
        continue;  // should not happen
      if (recvmatch[rec].sim != NOT_MATCHED)
        continue;  // already gone

      // do you, recvertex rec, take this simvertex sim as your lawful wedded truthmatch?
      unsigned int rec2sim = NOT_MATCHED;
      for (auto sv : recvmatch[rec].wos) {
        if (simEvt[sv.first].rec != NOT_MATCHED)
          continue;  // already used
        if ((rec2sim == NOT_MATCHED) || (sv.second > recvmatch[rec].wos[rec2sim])) {
          rec2sim = sv.first;
        }
      }

      if (sim == rec2sim) {
        // I do
        recvmatch[rec].sim = sim;
        recvmatch[rec].matchQuality = 8;
        simEvt.at(sim).rec = rec;
        simEvt.at(sim).matchQuality = 8;
        //cout << " late match  [" << sim << "]   <--> (" << rec << ")" << endl;
        nmatch++;
      }
    }  //sim loop
    //cout << "late matches " << nmatch <<  "  in try " << ntry << endl;
    if (nmatch == 0) {
      break;
    }
  }  // ntry

  if (DEBUG) {
    for (unsigned int iv = 0; iv < recVtxs->size(); iv++) {
      if (iv == 68) {
        dumpThisEvent_ = true;
        cout << "matching    recvtx = " << setw(4) << iv << "   (sim=" << recvmatch[iv].sim << ")" << endl;
        cout << " splitfrom= " << recvmatch[iv].split_from() << endl;
        cout << " maxwos = " << recvmatch[iv].maxwos << endl;
        cout << " sumwos = " << recvmatch[iv].sumwos << endl;

        for (auto w = recvmatch[iv].wos.begin(); w != recvmatch[iv].wos.end(); w++) {
          cout << "sim " << setw(4) << w->first << "  zsim = " << setw(8) << fixed << setprecision(4)
               << simEvt.at(w->first).z << "  wos=" << setw(10) << w->second << endl;
        }

        for (unsigned int iev = 0; iev < simEvt.size(); iev++) {
          for (auto rv = simEvt.at(iev).wnt.begin(); rv != simEvt.at(iev).wnt.end(); rv++) {
            if (rv->first == iv) {
              cout << "matching (4)  simevent " << iev << "    zsim= " << setw(8) << fixed << setprecision(4)
                   << simEvt.at(iev).z << "    zrec= " << setw(8) << fixed << setprecision(4)
                   << recVtxs->at((*rv).first).z() << "  nt=" << (*rv).second << endl;
            }
          }
        }
      }
    }
  }  // <<<< debugging

  /*
  // last chance for vertices that are not matched by tracks
  // but are reasonably close and not splitters
  for(unsigned int iev=0; iev<simEvt.size(); iev++)
    {
      if( simEvt[iev].rec == NOT_MATCHED )
	{
	  for(unsigned int iv=0; iv<recVtxs->size(); iv++)
	    {
	      
	      if ((recvmatch[iv].sim == NOT_MATCHED) && (recvmatch[iv].split_from() < 0))
		{
		  double dz = simEvt.at(iev).z-recVtxs->at(iv).z();
		  if( ( fabs(dz) < 0.1 )
		      ||( fabs(dz) < (2*recVtxs->at(iv).zError()) )
		      )
		  {
		    recvmatch[iv].matchQuality = 10;
		    recvmatch[iv].sim = iev;
		    simEvt.at(iev).rec = iv;
		    simEvt.at(iev).matchQuality = 10;
		  }
		}
	    }
	}
    }
  */

  // for convenience, store some values from the matched rec vertex with the sim vertex
  for (unsigned int iev = 0; iev < simEvt.size(); iev++) {
    int iv = simEvt.at(iev).rec;
    if (iv != NOT_MATCHED) {
      simEvt.at(iev).ndof = recVtxs->at(iv).ndof();
      simEvt.at(iev).zrec = recVtxs->at(iv).z();
    } else {
      simEvt.at(iev).ndof = 0.;
      simEvt.at(iev).zrec = 1000.;
    }
  }

}
/***************************************************************************************/




/***************************************************************************************/
void PrimaryVertexAnalyzer4PU::analyzeRecVertexComposition(std::map<std::string, TH1*>& h,
                                                           const reco::Vertex& v,
							   Tracks & tracks,
                                                           RSmatch& rs,
                                                           vector<SimEvent>& simEvt,
							   float npu)
// analyzes the composition of reco vertices in terms simEvts
// uses trkkey2simevt_ filled by matchTP
/***************************************************************************************/
{
  // track counters per simEvt
  VertexCounter nt;     // all
  VertexCounter ntMTD;  // in MTD acceptace

  // weighted sums (aka purities)
  double sum_pt(0), sum_pt_majority(0), sum_pt_minority(0), sum_pt_unmatched(0);
  double sum_nt(0), sum_nt_majority(0), sum_nt_minority(0), sum_nt_unmatched(0);
  double sum_wt(0), sum_wt_majority(0), sum_wt_minority(0), sum_wt_unmatched(0);

  // count tracks in the rec vertex
  for (trackit_t t = v.tracks_begin(); t != v.tracks_end(); t++) {
    Float_t wt = v.trackWeight(*t);
    double wpt = min(10., t->get()->pt());
    unsigned int tk_sim = tracks.simevent_index_from_key(t->key());
    bool MTD = (t->get()->pt() > 0.9) && (fabs(t->get()->eta()) < 2.7);

    if (wt > 0.5) {
      
      if(tk_sim != NOT_MATCHED){
	nt.count(tk_sim);
	if (MTD) {
	  ntMTD.count(tk_sim);
	}
      }
      
      if (tk_sim == rs.sim) { // so majority actually stands for "track matched to the same simevent as the vertex"
        sum_pt_majority += wpt;  
        sum_nt_majority += 1.;
        sum_wt_majority += wt;
      } else if (tk_sim == NOT_MATCHED) {
        sum_pt_unmatched += wpt;
        sum_nt_unmatched += 1.;
        sum_wt_unmatched += wt;
      } else {
        sum_pt_minority += wpt;
        sum_nt_minority += 1.;
        sum_wt_minority += wt;
      }
    }
  }

  
  // fill histograms
  if (rs.is_signal()) {
    Fill(h, "MTDTDR", sum_nt_minority);
  }

  Fill(h, "nsimevt", float(nt.nkey()), rs.is_signal());
  Fill(h, "nsimevt_nt2", float(nt.nkey(2)), rs.is_signal());
  Fill(h, "nsimevtMTD", float(ntMTD.nkey()), rs.is_signal());

  sum_pt = sum_pt_majority + sum_pt_minority + sum_pt_unmatched;
  sum_nt = sum_nt_majority + sum_nt_minority + sum_nt_unmatched;
  sum_wt = sum_wt_majority + sum_wt_minority + sum_wt_unmatched;

  if (false && (sum_pt_majority == 0) && (sum_pt > 0)) {
    cout << "analyzeRecVertexComposition, vertex without majority tracks ?" << endl;
    cout << "z = " << v.z() << endl;
    cout << "nt = " << sum_nt_majority << "," << sum_nt_minority << "," << sum_nt_unmatched << endl;
    cout << "match = " << rs.sim << endl;
    cout << "wosmatch = " << rs.wosmatch << " maxwos =" << rs.maxwos << " sumwos =" << rs.sumwos << endl;
    cout << "wntmatch = " << rs.wntmatch << " maxnt=" << rs.maxwnt << " sumnt=" << rs.sumwnt << endl;
    cout << "quality = " << rs.matchQuality << endl;
    unsigned int cand = rs.wosmatch;
    cout << "simEvt.at(cand).rec = " << simEvt.at(cand).rec << "    zrec=" << simEvt.at(cand).zrec << endl;
    cout << "simEvt.at(cand).nwosmatch = " << simEvt.at(cand).nwosmatch << endl;
    cout << "simEvt.at(cand).z = " << simEvt.at(cand).z << endl;
    cout << "dz = " << simEvt.at(cand).z - v.z() << " cutoff = " << zWosMatchMax_ << endl;
    cout << endl;
  }

  Fill(h, "pt_majority_frac", sum_pt_majority / sum_pt, rs.is_signal());
  Fill(h, "pt_minority_frac", sum_pt_minority / sum_pt, rs.is_signal());
  Fill(h, "pt_unmatched_frac", sum_pt_unmatched / sum_pt, rs.is_signal());

  Fill(h, "pt_majority_frac_vsz", v.z(), sum_pt_majority / sum_pt, rs.is_signal());
  Fill(h, "pt_minority_frac_vsz", v.z(), sum_pt_minority / sum_pt, rs.is_signal());
  Fill(h, "pt_unmatched_frac_vsz", v.z(), sum_pt_unmatched / sum_pt, rs.is_signal());

  Fill(h, "nt_majority_frac", sum_nt_majority / sum_nt, rs.is_signal());
  Fill(h, "nt_minority_frac", sum_nt_minority / sum_nt, rs.is_signal());
  Fill(h, "nt_unmatched_frac", sum_nt_unmatched / sum_nt, rs.is_signal());

  Fill(h, "nt_majority_frac_vsz", v.z(), sum_nt_majority / sum_nt, rs.is_signal());
  Fill(h, "nt_minority_frac_vsz", v.z(), sum_nt_minority / sum_nt, rs.is_signal());
  Fill(h, "nt_unmatched_frac_vsz", v.z(), sum_nt_unmatched / sum_nt, rs.is_signal());

  Fill(h, "wt_majority_frac", sum_wt_majority / sum_wt, rs.is_signal());
  Fill(h, "wt_minority_frac", sum_wt_minority / sum_wt, rs.is_signal());
  Fill(h, "wt_unmatched_frac", sum_wt_unmatched / sum_wt, rs.is_signal());

  if(sum_nt > 0){
    float purity = sum_nt_majority / sum_nt;
    Fill(h, "vtxtrkpurity", purity, rs.is_signal());             // same as nt_majority_frac
    Fill(h, "vtxtrkpurityvsz", v.z(), purity, rs.is_signal());    // average purity vs z
    Fill(h, "vtxtrkpurityvspu", npu, purity, rs.is_signal());
  }

  float sum_nt_m = sum_nt_majority + sum_nt_minority;
  if(sum_nt_m > 0){
    float purity_m = sum_nt_majority / sum_nt_m;
    Fill(h, "vtxtrkpuritym", purity_m, rs.is_signal());             // same as nt_majority_frac
    Fill(h, "vtxtrkpuritymvsz", v.z(), purity_m, rs.is_signal());    // average purity vs z
    Fill(h, "vtxtrkpuritymvspu", npu, purity_m, rs.is_signal());
  }
  
  // note that this efficiency is biased, because we fill it for vertices that have been reconstructed only
  double nmatch = sum_nt_majority;
  if(rs.sim != NOT_MATCHED){
    auto numtk = simEvt[rs.sim].rtk.size();           // denominator for assignment efficiency
    if (numtk > 0) {
      Fill(h, "trkAssignmentEfficiency", std::min(double(nmatch) / numtk, 1.001), rs.is_signal());  // also x-check matchVtxFraction
      Fill(h, "trkAssignmentEfficiencyvspu", npu, double(nmatch)/numtk, rs.is_signal());  // also x-check matchVtxFraction
      Fill(h, "trkAssignmentEfficiencyvsntrk", std::min(double(numtk), 50.5), nmatch / numtk, rs.is_signal());
    }
    auto numprimtk = simEvt[rs.sim].rtkprim.size();
    if (numprimtk > 0) {
      Fill(h, "primtrkAssignmentEfficiency", min(nmatch / numprimtk, 1.001), rs.is_signal());
    }
  }
}
/***************************************************************************************/




/***************************************************************************************/
void PrimaryVertexAnalyzer4PU::analyzeVertexMergeRateTP(std::map<std::string, TH1*>& h,
                                                        const reco::VertexCollection* recVtxs,
                                                        vector<SimEvent>& simEvt,
                                                        std::vector<RSmatch>& recvmatch,
                                                        const string message)
/***************************************************************************************/
// definitions ?
//  * a sim vertex that was not found but would have been reconstructed if there had not been other vertices (nearby)?
//  * a rec vertex that is matched to >= 2 sim vertices (dqm) but what is the merge rate then? <> merge-fraction ?
//
{
  auto nv = recVtxs->size();
  
  for (unsigned int iv = 0; iv < nv; iv++) {
    
    // only consider selected vertices
    if (!select(recVtxs->at(iv)))
      continue;


    unsigned int nsim = 0;      // count sim vertices that make their largest contribution in this recvtx
    unsigned int nsim_any = 0;  // count sim vertices that make a contribution to this vertex
    unsigned int nsim_any_not_matched = 0;  // count sim vertices that make a contribution but are not matched to this vertex

    //cout << setw(4) << iv << "  -> [" << recvmatch[iv].wos.size() << "]  ";
    for (auto s : recvmatch[iv].wos) {
      auto sim = s.first;
      nsim_any++;
      // only count if that sim vertex was unmatched or matched to this rec vertex
      if ((simEvt[sim].rec == NOT_MATCHED) || (simEvt[sim].rec == iv)) {
        nsim_any_not_matched++;
        // count if this rec vertex received the largest contribution of that simvtx
        if (simEvt[sim].max_nt_vtx() == iv) {
          nsim++;
        }
      }
    }

    //    cout << " ==>  nsim=" << nsim << " nsim_any=" << nsim_any << " nsim_any_not_matched=" << nsim_any_not_matched << endl;
    
    double z = recVtxs->at(iv).z();
    bool signal = recvmatch[iv].is_signal();
    Fill(h, "nsim_recvtx_sel", float(nsim), signal);
    if (nsim_any>10.) nsim_any=9.9;
    Fill(h, "nsim_any_recvtx_sel", float(nsim_any), signal);
    Fill(h, "nsim_any_not_matched_recvtx_sel", float(nsim_any_not_matched), signal);
    Fill(h, "nsimvsz_recvtx_sel", z, float(nsim), signal);

    // z-dependent merge rates (in some sense)
    Fill(h, "zrec_recvtx_sel", z, signal);  // denominator for "zrec_recvtx_selnsim*"
    if (nsim == 0) {
      Fill(h, "zrec_recvtx_selnsim0", z, signal);  // no sim vertex really wants this vertex
    } else if (nsim == 1) {
      Fill(h, "zrec_recvtx_selnsim1", z, signal);
    } else {
      Fill(h, "zrec_recvtx_selnsimgt1", z, signal);
    }
  }

  
  // fake rate and fake fraction vs distance to the signal vertex
  if ((simEvt.size()>0) && (simEvt[0].is_signal()) && (simEvt[0].matched()))
    {
      double zsimsignal = simEvt[0].z;
      double zrecsignal = simEvt[0].zrec;
      double recsignal =  simEvt[0].rec;

      Fill(h,"countzrecsignal",1); // denominator
      
      for (unsigned int iv = 0; iv < recVtxs->size(); iv++) {
	if (!select(recVtxs->at(iv)) || (iv == recsignal) )
	  continue;
	double dzrecsim = recVtxs->at(iv).z() - zsimsignal;
	double dzrecrec = recVtxs->at(iv).z() - zrecsignal;
	bool fake = (recvmatch[iv].sim == NOT_MATCHED);
	Fill(h,"zrecsimsignal_unmatchedPUfraction", dzrecsim, fake ? 1. : 0.); // 
	Fill(h,"zrecrecsignal_unmatchedPUfraction", dzrecrec, fake ? 1. : 0.); //
	Fill(h,"zrecsimsignal_unmatchedPUfraction_hr", dzrecsim, fake ? 1. : 0.); // 
	Fill(h,"zrecrecsignal_unmatchedPUfraction_hr", dzrecrec, fake ? 1. : 0.); //
	if ( fake )
	  {
	    // a unmatched/fake vertex
	    Fill(h,"zrecsimsignal_unmatchedPU", dzrecsim); // 
	    Fill(h,"zrecrecsignal_unmatchedPU", dzrecrec); //
	    Fill(h,"zrecsimsignal_unmatchedPU_hr", dzrecsim); // 
	    Fill(h,"zrecrecsignal_unmatchedPU_hr", dzrecrec); //
	    if (std::abs(dzrecsim) < 0.0010){//0.0050){
	      reportVertex(recVtxs->at(iv), Form("fake vertex on top of signal at z=%f in %s", recVtxs->at(iv).z(), message.c_str()), dump_fake_vertex_on_top_of_signal_);
	    }
	    
	  }
      }
    
  }

  // another approach, count simvertices that have no associated recvertex.
  // in order to separate merging and other inefficiencies, count in bins of
  // some measure of local density of other vertices, i.e. the distance to the nearest
  // recvertex or simvtx? distance to the signal vertex !





  for (unsigned int sim1 = 1; sim1 < simEvt.size(); sim1++)
    {
      double dzmin = 1000.;
      // note that we are not putting any restrictions on the first sim vertex here other than that it is NOT the signal vertex
      for (unsigned int sim2 = 1; sim2 < simEvt.size(); sim2++)
	{
	  if (sim1 == sim2) continue;
	  double dz = simEvt.at(sim2).z - simEvt.at(sim1).z;
	  if (std::abs(dz) < std::abs(dzmin)){dzmin = dz;}
	  Fill(h, "effvsdzsimPU_PU", dz, simEvt.at(sim2).matched() ? 1. : 0.);
	  Fill(h, "effvsdzsimPU_PU_hr", dz, simEvt.at(sim2).matched() ? 1. : 0.);
	  if(simEvt.at(sim1).matched())
	    {
	      Fill(h, "effvsdzsimmatchedPU_PU", dz, simEvt.at(sim2).matched() ? 1. : 0.);
	      Fill(h, "effvsdzsimmatchedPU_PU_hr", dz, simEvt.at(sim2).matched() ? 1. : 0.);
	      
	      bool sel = simEvt.at(sim2).matched() && select(recVtxs->at(simEvt.at(sim2).rec)) ;
	      Fill(h, "effselvsdzsimmatchedPU_PU", dz, sel ? 1. : 0.);
	      Fill(h, "effselvsdzsimmatchedPU_PU_hr", dz, sel ? 1. : 0.);
	    }
	  else
	    {
	      Fill(h, "effvsdzsimunmatchedPU_PU", dz, simEvt.at(sim2).matched() ? 1. : 0.);
	      Fill(h, "effvsdzsimunmatchedPU_PU_hr", dz, simEvt.at(sim2).matched() ? 1. : 0.);

	      bool sel = simEvt.at(sim2).matched() && select(recVtxs->at(simEvt.at(sim2).rec)) ;
	      Fill(h, "effselvsdzsimunmatchedPU_PU", dz, sel ? 1. : 0.);
	      Fill(h, "effselvsdzsimunmatchedPU_PU_hr", dz, sel ? 1. : 0.);
	    }
	}
      // efficiency vs distance to the nearest other PU sim vertex (excluding signal vertices)
      Fill(h, "effallvsdzsimminPU_PU", dzmin, simEvt.at(sim1).matched() ? 1. :0);
      bool sel1 = simEvt.at(sim1).matched() && select(recVtxs->at(simEvt.at(sim1).rec)) ;
      Fill(h, "effselvsdzsimminPU_PU", dzmin, sel1 ? 1. :0);
    }
  
  // efficiency of pu vtx finding vs distance to the signal vertex
  if ((simEvt.size()>0) && (simEvt.at(0).matched()))
    {
      assert( simEvt.at(0).is_signal());
      for (unsigned int sim = 1; sim < simEvt.size(); sim++)
	{
	  double dz = simEvt.at(sim).z - simEvt.at(0).z;
	  Fill(h, "effvsdzsimsignal_PU", dz, simEvt.at(sim).matched() ? 1. : 0.);
	  Fill(h, "effvsdzsimsignal_PU_hr", dz, simEvt.at(sim).matched() ? 1. : 0.);
	  
	  bool sel = simEvt.at(sim).matched() && select(recVtxs->at(simEvt.at(sim).rec)) ;
	  Fill(h, "effselvsdzsimsignal_PU", dz, sel ? 1. : 0.);
	  Fill(h, "effselvsdzsimsignal_PU_hr", dz, sel ? 1. : 0.);
	}
    }
  
}//analyzeVertexMergeRateTP


/***************************************************************************************/





/***************************************************************************************/
void PrimaryVertexAnalyzer4PU::analyzeVertexCollectionTP(std::map<std::string, TH1*>& h,
                                                         const reco::VertexCollection* recVtxs,
                                                         Tracks& tracks,
                                                         vector<SimEvent>& simEvt,
                                                         std::vector<RSmatch>& recvmatch,
                                                         const string message)
/***************************************************************************************/
{
  if (verbose_) {
    cout << "analyzeVertexCollectionTP, simEvts= " << simEvt.size() << endl;
  }

  if (simEvt.size() == 0)
    return;

  //EncodedEventId iSignal = simEvt[0].eventId;

  // count pile-up
  int npu0 = simEvt.size();
  Fill(h, "npu0", npu0);
  float npu = simEvt.size();

  if (eventSummaryCounter_++ < nEventSummary_) {
    printEventSummary(h, recVtxs, tracks, simEvt, recvmatch, message);
  }

  //
  unsigned int ntpfake = 0;
  unsigned int ntpfake4 = 0;
  unsigned int ntpfound = 0, ntpfound4 = 0;
  unsigned int ntpfakesel = 0, ntprealsel = 0, ntpsplitsel = 0, ntpotherfakesel = 0, ntpsel = 0;
  unsigned int ntpsplitselfromsignal = 0, ntpsplitselfrompu = 0;

  auto nv = recVtxs->size();
  for (unsigned int iv = 0; iv < nv; iv++) {

    if (select(recVtxs->at(iv))) {

      // (signed) distance to the nearest selected other vertex
      // used for signed quantities in fillVertexHistos
      double deltaz = 1000.;
      for(unsigned int jv = 0; jv < recVtxs->size(); jv++)
	{
	  if ( (jv != iv) && (fabs(recVtxs->at(jv).z() - recVtxs->at(iv).z()) < fabs(deltaz)) && select(recVtxs->at(jv)) )
	    {
	      deltaz = recVtxs->at(jv).z() - recVtxs->at(iv).z();
	    }
	}
      
      ntpsel++;
      if (recvmatch[iv].is_fake()) {
        ntpfakesel++;
      }
      if (recvmatch[iv].is_real()) {
        ntprealsel++;
	fillVertexHistosMatched(h, "matchedvtxsel",  &(recVtxs->at(iv)), tracks, recvmatch[iv], simEvt, iv, deltaz);
	if(recvmatch[iv].is_signal()){
	  fillVertexHistosMatched(h, "matchedsignalvtxsel",  &(recVtxs->at(iv)), tracks, recvmatch[iv], simEvt, iv, deltaz);
	  for (trackit_t t = recVtxs->at(iv).tracks_begin(); t != recVtxs->at(iv).tracks_end(); t++) {
	    auto tk = tracks.from_ref(*t);
	    if (tk.dz < 0.01) {
	      fillTrackHistos(h, "matchedsignalvtxseldriver", tk, &(recVtxs->at(iv)));
	    }
	  }
	  if (iv == 0){
	    reportEvent(Form("signal vertex selected at index [%d]", iv), false);
	  }else{
	    std::string comment = "?";
	    if (recvmatch[0].is_real()){ comment = "real vertex at [0]";}
	    else if(recvmatch[0].other_fake()){ comment = "other fake vertex at [0]";}
	    else if(recvmatch[0].split_from() > 0){ comment = "split from " + std::to_string(recvmatch[0].split_from());}
	    else if(recvmatch[0].is_fake()){ comment = "fake vertex at [0]";}
	    else if(recvmatch[0].is_signal()){ comment = "signal vertex at [0]?????";}
	    reportEvent(Form("signal vertex selected at index [%d]    %s", iv,comment.c_str()), false);
	  }
	}else{
	  fillVertexHistosMatched(h, "matchedpuvtxsel",  &(recVtxs->at(iv)), tracks, recvmatch[iv], simEvt, iv, deltaz);
	  if(iv == 0){
	    reportEvent(Form("PU vertex selected at index [0]"), false);
	  }
	}
      }

      if (recvmatch[iv].split_from() >= 0) {
        ntpsplitsel++;
	fillVertexHistos(h, "splitvtxsel",  &(recVtxs->at(iv)), tracks, iv, deltaz);
      }

      if (recvmatch[iv].split_from() == 0) {
        ntpsplitselfromsignal++;
	fillVertexHistos(h, "splitvtxselfromsignal",  &(recVtxs->at(iv)), tracks, iv, deltaz);
	reportEvent(Form("split signal vertex selected at index [%d]", iv), false);

	for (trackit_t t = recVtxs->at(iv).tracks_begin(); t != recVtxs->at(iv).tracks_end(); t++) {
	  auto tk = tracks.from_ref(*t);
	  if (tk.dz < 0.0100) {
	    fillTrackHistos(h, "splitvtxselfromsignaldriver", tk, &(recVtxs->at(iv)));
	  }
        }
	
      }
      if (recvmatch[iv].split_from() > 0) {
        ntpsplitselfrompu++;
	fillVertexHistos(h, "splitvtxselfrompu",  &(recVtxs->at(iv)), tracks, iv, deltaz);
      }
      if (recvmatch[iv].other_fake()) {
        ntpotherfakesel++;
	fillVertexHistos(h, "otherfakevtxsel",  &(recVtxs->at(iv)), tracks, iv, deltaz);
      }
    }

    double z = recVtxs->at(iv).z();
    double vxx = recVtxs->at(iv).covariance(iX, iX) + pow(vertexBeamSpot_.BeamWidthX(), 2);
    double vyy = recVtxs->at(iv).covariance(iY, iY) + pow(vertexBeamSpot_.BeamWidthY(), 2);
    ;
    double vxy = recVtxs->at(iv).covariance(iX, iY);
    double dx = recVtxs->at(iv).x() - vertexBeamSpot_.x(z);
    double dy = recVtxs->at(iv).y() - vertexBeamSpot_.y(z);
    double D = vxx * vyy - vxy * vxy;
    double c2xy = pow(dx, 2) * vyy / D + pow(dy, 2) * vxx / D - 2 * dx * dy * vxy / D;
    Fill(h, "vtxcxy2", c2xy);

    double vtx_trec = recVtxs->at(iv).t();
    double vtx_terr = recVtxs->at(iv).tError();

    // this is obsolete, use what is filled by   fillVertexHistosMatched(h, "matchedvtx"...) instead
    if (recvmatch[iv].is_real() && select(recVtxs->at(iv))) {
      double vtx_tsim = simEvt[recvmatch[iv].sim].t;
      Fill(h, "trecvtx_selmatched", vtx_trec);
      Fill(h, "tsimvtx_selmatched", vtx_tsim);
      Fill(h, "tresvtx_selmatched", vtx_trec - vtx_tsim);
      Fill(h, "tpullvtx_selmatched", (vtx_trec - vtx_tsim) / vtx_terr);
      auto xres = recVtxs->at(iv).x() - simEvt[recvmatch[iv].sim].x;
      auto yres = recVtxs->at(iv).y() - simEvt[recvmatch[iv].sim].y;
      auto zres = recVtxs->at(iv).z() - simEvt[recvmatch[iv].sim].z;
      bool is_signal = recvmatch[iv].is_signal();
      Fill(h, "xresvtx_selmatched", xres, is_signal);
      Fill(h, "yresvtx_selmatched", yres, is_signal);
      Fill(h, "zresvtx_selmatched", zres, is_signal);
    }

    if (recvmatch[iv].is_real()) {
      fillVertexHistosMatched(h, "matchedvtx", &recVtxs->at(iv), tracks, recvmatch[iv], simEvt, iv);
    }

    if (recvmatch[iv].is_fake()) {
      ntpfake++;
      if (recVtxs->at(iv).ndof() > 4) {
        ntpfake4++;
      }

      Fill(h, "vtxcxy2Fake", c2xy);
      Fill(h, "vtxprobxyFake", TMath::Prob(c2xy, 2));
      Fill(h, "ndofFake", recVtxs->at(iv).ndof());
      Fill(h, "logndofFake", recVtxs->at(iv).ndof() > 0 ? log(recVtxs->at(iv).ndof())/log(10.) : -10.);

      if (recVtxs->at(iv).ndof() > 50) {
        reportVertex(recVtxs->at(iv), Form("big fake vertex ndof=%f , split_from=%d", 
					   recVtxs->at(iv).ndof(), recvmatch[iv].split_from()), dump_big_fakes_);
      }
    } else {
      Fill(h, "vtxcxy2Matched", c2xy);
      Fill(h, "vtxprobxyMatched", TMath::Prob(c2xy, 2));
      Fill(h, "ndofMatched", recVtxs->at(iv).ndof());
      ntpfound++;

      if (select(recVtxs->at(iv))) {
        ntpfound4++;
      }
    }

  }  // end of recvtx loop


  Fill(h, "ntpsplitselfromsignal", float(ntpsplitselfromsignal));
  Fill(h, "ntpsplitselfrompu", float(ntpsplitselfrompu));
  Fill(h, "ntpfake0", float(ntpfake));
  Fill(h, "ntpfake4", float(ntpfake4));
  Fill(h, "ntpfound0", float(ntpfound));
  Fill(h, "ntpfound4", float(ntpfound4));
  Fill(h, "ntpfake4prof", npu, float(ntpfake4));
  Fill(h, "ntpfound4prof", npu, float(ntpfound4));

  Fill(h, "ntpselvssimPU", simPU_, float(ntpsel));
  Fill(h, "ntpfakeselvssimPU", simPU_, float(ntpfakesel));
  Fill(h, "ntprealselvssimPU", simPU_, float(ntprealsel));
  Fill(h, "ntpsplitselvssimPU", simPU_, float(ntpsplitsel));
  Fill(h, "ntpotherfakeselvssimPU", simPU_, float(ntpotherfakesel));




  for(unsigned int iv = 0; iv < recVtxs->size(); iv++){
    if( recvmatch[iv].is_signal() ){
      Fill(h, "index_signal", float(iv));
    }else if( recvmatch[iv].is_real() ){
      Fill(h, "index_pu", float(iv));
    }else if( recvmatch[iv].split_from() == 0) {
      Fill(h, "index_splitfromsignal", float(iv));
    }else if( recvmatch[iv].split_from() > 0) {
      Fill(h, "index_splitfrompu", float(iv));
    }else{
      Fill(h, "index_otherfake", float(iv));
    }
  }


  // ranking by track sumpt2:
  std::vector< std::pair<double, unsigned int> > ranking;
  for (unsigned int iv = 0; iv < recVtxs->size(); iv++) { // why did it start at 1 ?
    if (select(recVtxs->at(iv))){
      ranking.push_back( std::make_pair(-vertex_sumpt2(recVtxs->at(iv)),iv) );
    }
  }

  stable_sort(ranking.begin(), ranking.end());
  for(unsigned int n = 0; n < ranking.size(); n++){
    unsigned int iv = ranking[n].second;
    if( recvmatch[iv].is_signal() ){
      Fill(h, "trksumpt2rank_signal", float(n));
    }else if( recvmatch[iv].is_real() ){
      Fill(h, "trksumpt2rank_pu", float(n));
    }else if( recvmatch[iv].split_from() == 0) {
      Fill(h, "trksumpt2rank_splitfromsignal", float(n));
    }else if( recvmatch[iv].split_from() > 0) {
      Fill(h, "trksumpt2rank_splitfrompu", float(n));
    }else{
      Fill(h, "trksumpt2rank_otherfake", float(n));
    }
  }




  // for the classification of unmatched rec- and sim-vertices
  // rec first: how many sim vertices do I need to get 50% of my wos
  for (unsigned int iv = 0; iv < recVtxs->size(); iv++)
    {
      if (select(recVtxs->at(iv))) {
	
	  // sort the wos values
	  std::vector<double> wos;
	  double sumwos = 0.;
	  for (auto it : recvmatch[iv].wos) {
	    sumwos += it.second;
	    wos.push_back(it.second);
	  }
	  stable_sort(wos.rbegin(), wos.rend());  // reverse

	  // find the number of sim vertices needed to get 70 %
	  // interpretation for fakes:
	  // a number like 1 means splitter
	  // a number like 2 means mergers (unlikely to be classified as fakes)
	  // a large number means garbage collection
	  double cumwos = 0.;
	  unsigned n = 0;
	  for(; n < wos.size(); n++){
	    cumwos += wos[n];
	    if (cumwos > (0.7 * sumwos)) break;
	  }
	  // for n==0: cumwos==recvmatch[iv].sumwomaxwos, sumwos==recvmatch[iv].sumwos
	  
	  n = n > 4 ? 4 : n;// fill overflow in the last bin
	  double v = log(float(sumwos))/log(10.);
	  if (recvmatch[iv].is_fake()){
	    Fill(h, "wornk_unmatchedsel", v, float(n));
	    if(recvmatch[iv].split_from() >= 0){
	      Fill(h, "wornk_splitsel", v, float(n));
	    }else{
	      Fill(h, "wornk_othersel", v, float(n));
	    }
	  }else{
	    Fill(h, "wornk_matchedsel", v, float(n));
	  }
	  
      }
    }
	
  // now something similar for sim vertices
  // here we want to distinguish clear mergers from 'just below threshold'
  // 
  for (auto ev = simEvt.begin(); ev != simEvt.end(); ev++)
    {
    
      // sort the wos values
      std::vector<double> wos;
      double sumwos = 0.;
      for (auto it : ev->wos) {
	sumwos += it.second;
	wos.push_back(it.second);
      }
      stable_sort(wos.rbegin(), wos.rend());  // reverse

      // find the number of rec vertices that receive 70 %
      // a number like 1 means swallowed/merged by a single vertex
      // a large number means it was eaten by the pack
      double cumwos = 0.;
      unsigned n = 0;
      for(; n < wos.size(); n++){
	cumwos += wos[n];
	if (cumwos > (0.7 * sumwos)) break;
      }
      n = n > 4 ? 4 : n;// fill overflow in the last bin
      double v = log(float(sumwos))/log(10.);
      if (ev->matched()){
	Fill(h, "wornk_matchedsim", v, float(n));
      }else{
	Fill(h, "wornk_unmatchedsim", v, float(n));
      }
    }

  // determine the recvertex that is the signalvertex (if there is one)
  unsigned int signalv = 0;
  bool has_signal = false;
  for (unsigned int iv = 0; iv < recVtxs->size(); iv++) {  // was size-1, why?
    if (select(recVtxs->at(iv)) && recvmatch[iv].is_real() && recvmatch[iv].is_signal()) {
      signalv = iv;
      has_signal = true;
      break;
    }
  }

  // properties of vertices that are split from the signal vertex
  if (has_signal) {
    for (unsigned int jv = 0; jv < recVtxs->size(); jv++) {
      if (recvmatch[jv].split_from() == 0) {
        Fill(h, "ndofsignalsplit", recVtxs->at(jv).ndof());

	if (select(recVtxs->at(signalv)) && select(recVtxs->at(jv))){
          Fill(h, "zdiffrecsignalsel", fabs(recVtxs->at(signalv).z() - recVtxs->at(jv).z()));
	}

        if (select(recVtxs->at(jv))) {
          Fill(h, "zdiffrecselsignalsplit", fabs(recVtxs->at(signalv).z() - recVtxs->at(jv).z()));
          Fill(h, "zdiffrecselsignalsplit-dzbin", fabs(recVtxs->at(signalv).z() - recVtxs->at(jv).z()));
          // add more, e.g. sumpt, fraction of pt
        }
      }
    }
  }

  // distance histogram for pairs of selected reco vertices
  for (unsigned int iv = 0; iv < recVtxs->size(); iv++) {
    if (select(recVtxs->at(iv))) {
      for (unsigned int jv = 0; jv < recVtxs->size(); jv++) {
        if ((!(jv == iv)) && select(recVtxs->at(jv))) {
          double dz = recVtxs->at(iv).z() - recVtxs->at(jv).z();
          if (recvmatch[jv].is_real() && recvmatch[iv].is_real()) {
            Fill(h, "zdiffrecselrealreal-dzbin", std::abs(dz));
          } else {
            Fill(h, "zdiffrecselanyfake-dzbin", std::abs(dz));
          }
        }
      }
    }
  }

  // histogram the neighborhood of the signal vertex
  if (has_signal) {
    for (unsigned int jv = 0; jv < recVtxs->size(); jv++) {

      if ((!(jv == signalv)) && select(recVtxs->at(jv)) && (!recvmatch[jv].is_signal())) {
        double dz = recVtxs->at(signalv).z() - recVtxs->at(jv).z();
        if (recvmatch[jv].is_real()) {
          Fill(h, "zdiffrecselsignalrealpu", std::abs(dz));
          if(jv == 1) Fill(h, "zdiffrecselsignalrealpuidx1", std::abs(dz));
        } else {
          Fill(h, "zdiffrecselsignalfake", std::abs(dz));
          if (fabs(dz) < 0.0) {  // only useful at low pu
	    reportEvent(Form(" suspected signal splitting :  z(signal)= %10.4f  z(fake)= %10.4f  distance=%10.4f",recVtxs->at(signalv).z(),recVtxs->at(jv).z() , fabs(dz)));
          }
	  if (recvmatch[jv].split_from() == 0){
	    Fill(h, "zdiffrecselsignalsplitfromsignal", std::abs(dz));
	  }else if(recvmatch[jv].split_from() > 0){
	    Fill(h, "zdiffrecselsignalsplitfrompu", std::abs(dz));
	  }else{
	    Fill(h, "zdiffrecselsignalotherfake", std::abs(dz));
	  }
        }
      }

      // more or less the same, but older
      if ((!(jv == signalv)) && select(recVtxs->at(jv)) && (!recvmatch[jv].is_signal())) {
        double dz = recVtxs->at(signalv).z() - recVtxs->at(jv).z();
        if (recvmatch[jv].is_real()) {
          Fill(h, "zdiffrecselsignalreal", std::abs(dz));
          Fill(h, "zdiffrecselsignalreal-dzbin", std::abs(dz));
        } else {
          Fill(h, "zdiffrecselsignalfake", std::abs(dz));
          Fill(h, "zdiffrecselsignalfake-dzbin", std::abs(dz));
        }
      }
    }
  }

  // histogram the neighborhood of a PU vertex
  for (unsigned int iv = 0; iv < recVtxs->size() - 1; iv++) {
    if ( select(recVtxs->at(iv)) && recvmatch[iv].is_real() && (!recvmatch[iv].is_signal())) {
      for (unsigned int jv = 0; jv < recVtxs->size(); jv++) {
        if ((!(jv == iv)) && select(recVtxs->at(jv)) && (!recvmatch[jv].is_signal())) {
          double dz = recVtxs->at(iv).z() - recVtxs->at(jv).z();
          if (recvmatch[jv].is_real()) {
            Fill(h, "zdiffrecselPUrealPU", std::abs(dz));
          } else {
            Fill(h, "zdiffrecselPUfake", std::abs(dz));
          }
        }
      }
    }
  }

  // all pairs of reconstructed vertices
  timer_start("TP_all_pairs_loop");
  for (unsigned int iv = 0; iv < recVtxs->size() - 1; iv++) {
    if (recVtxs->at(iv).ndof() > 4.) {
      double mindistance_realreal = 1e10;

      for (unsigned int jv = iv; jv < recVtxs->size(); jv++) {
        if ((!(jv == iv)) && select(recVtxs->at(jv))) {
          double dz = recVtxs->at(iv).z() - recVtxs->at(jv).z();
          if (recvmatch[iv].is_real() && recvmatch[jv].is_real()) {
            Fill(h, "zdiffrecselfound", dz);
            if (fabs(dz) < fabs(mindistance_realreal)) {
              mindistance_realreal = dz;
            }
          } else if (recvmatch[iv].is_fake() && recvmatch[jv].is_fake()) {
            Fill(h, "zdiffrecselfakefake", dz);
          }
        }
      }

      double mindistance_fakereal = 1e10;
      double mindistance_realfake = 1e10;
      double zf = 0, zr = 0;
      for (unsigned int jv = 0; jv < recVtxs->size(); jv++) {
        if ((!(jv == iv)) && select(recVtxs->at(jv))) {
          double dz = recVtxs->at(iv).z() - recVtxs->at(jv).z();

          if (recvmatch[iv].is_fake() && recvmatch[jv].is_real()) {
            Fill(h, "zdiffrecselfakereal", dz);
            if (fabs(dz) < fabs(mindistance_fakereal)) {
              mindistance_fakereal = dz;
              zr = recvmatch[iv].is_real() ? recVtxs->at(iv).z() : recVtxs->at(jv).z();
              zf = recvmatch[iv].is_fake() ? recVtxs->at(iv).z() : recVtxs->at(jv).z();
            }
          }

          if (recvmatch[iv].is_real() && recvmatch[jv].is_fake() && (fabs(dz) < fabs(mindistance_realfake))) {
            mindistance_realfake = dz;
          }
        }
      }

      if (recvmatch[iv].is_real()) {
        Fill(h, "zdiffmin4realreal", fabs(mindistance_realreal));
        Fill(h, "zdiffmin4realfake", fabs(mindistance_realfake));
        if (recvmatch[iv].is_signal()) {
          Fill(h, "zdiffmin4signalreal", fabs(mindistance_realreal));
          Fill(h, "zdiffmin4signalfake", fabs(mindistance_realfake));
        }
      } else if (recvmatch[iv].is_fake()) {
        Fill(h, "zdiffmin4fakereal", fabs(mindistance_fakereal));
        if (fabs(mindistance_fakereal) < 0.00) {
          cout << "run : event = " << run_ << ":" << event_
               << " suspected vertex splitting :  z(real)=" << setprecision(4) << zr << "   z(fake)=" << zf << endl;
        }
      }
    }
  }
  timer_stop("TP_all_pairs_loop");


  // properties of the first/second vertex (in the sorted list)
  if ( (nv > 1) && select(recVtxs->at(0)) && select(recVtxs->at(1)) ){
    double zdiff = recVtxs->at(1).z() - recVtxs->at(0).z();
    Fill(h, "zdiffrec10-selsel", zdiff);
    if (recvmatch[0].is_signal()){
      Fill(h, "zdiffrec10-signalsel", zdiff);
      if (recvmatch[1].is_real()){ Fill(h, "zdiffrec10-signalreal", zdiff);} // redundant and obsolete
      if (recvmatch[1].is_real()){ Fill(h, "zdiffrec10-signalrealpu", zdiff);}
      if (recvmatch[1].is_fake()){ Fill(h, "zdiffrec10-signalfake", zdiff);}
      if (recvmatch[1].split_from() == 0){ 
	Fill(h, "zdiffrec10-signalsplitfromsignal", zdiff);
      }else if (recvmatch[1].split_from() > 0){ 
	Fill(h, "zdiffrec10-signalsplitfrompu", zdiff);
      }else{
	Fill(h, "zdiffrec10-signalotherfake", zdiff);
      }

      if (recvmatch[1].is_real() &&  (fabs(zdiff) < 0.05)) reportEvent("real pu close to signal promoted to nr 2", false);
    }
  }


  // efficiency histograms for simEvts
  for (unsigned int iev = 0; iev < simEvt.size(); iev++) {
    bool is_signal = (iev == 0);

    //simEvt matching done in tpmatch
    float wall = 0., wsel=0., wsel3sigma=0.;
    if (simEvt.at(iev).rec != NOT_MATCHED) {

      wall = 1.;

      // selected 
      if( select(recVtxs->at(simEvt.at(iev).rec)) ){
	wsel = 1.;
	if (fabs(recVtxs->at(simEvt.at(iev).rec).z() - simEvt.at(iev).z) < 3.* recVtxs->at(simEvt.at(iev).rec).zError()){
	  wsel3sigma = 1.;
	}
      }
    }

    Fill(h, "effallTP", wall, is_signal);
    Fill(h, "effselTP", wsel, is_signal);
    Fill(h, "effsel3sigma", wsel3sigma, is_signal);
    Fill(h, "effallvsnsimevt", npu, wall, is_signal);
    Fill(h, "effselvsnsimevt", npu, wsel, is_signal);
    Fill(h, "effsel3sigmavsnsimevt", npu, wsel3sigma, is_signal);
    Fill(h, "effallvsnrectp", float(simEvt.at(iev).rtk.size()), wall, is_signal);
    Fill(h, "effselvsnrectp", float(simEvt.at(iev).rtk.size()), wsel, is_signal);
    Fill(h, "effallvsngentp", float(simEvt.at(iev).nChTP), wall, is_signal);
    Fill(h, "effselvsngentp", float(simEvt.at(iev).nChTP), wsel, is_signal);
    Fill(h, "effallvssumpt", float(simEvt.at(iev).sumpt), wall, is_signal);
    Fill(h, "effselvssumpt", float(simEvt.at(iev).sumpt), wsel, is_signal);
    Fill(h, "effallvssqrtsumpt2", sqrt(float(simEvt.at(iev).sumpt2)), wall, is_signal);
    Fill(h, "effselvssqrtsumpt2", sqrt(float(simEvt.at(iev).sumpt2)), wsel, is_signal);
    //Fill(h, "effallvstrksqrtsumpt2", sqrt(), wall, is_signal); FIXME, todo
    //Fill(h, "effselvstrksqrtsumpt2", sqrt(), wsel, is_signal);

    if ((wall == 0) && (is_signal)) {
      reportEvent(
		  Form("signal vertex not TP matched (z=%8.4f    charged trks=%d)",
		       simEvt.at(iev).z, (int)simEvt.at(iev).rtk.size()),
		  dump_signal_vertex_not_tpmatched_);
      
      cout << " signal vertex not found :" << message << "   run : event = " << run_ << " : " << event_ << endl;
      cout << "sim " << iev << "  z=" << simEvt.at(iev).z << " ntk=" << simEvt.at(iev).rtk.size()
           << "  nChTP=" << simEvt.at(iev).rtk.size() << endl;
    } else if (verbose_ && (wall == 0) && (simEvt.at(iev).rtk.size() > 20)) {
       reportEvent(
		   Form("big PU vertex not found (TP) (z=%8.4f    charged trks=%d)",
			simEvt.at(iev).z, (int)simEvt.at(iev).rtk.size()),
		   false);
      cout << " big PU vertex not found :" << message << "  ";
      cout << "sim " << iev << "  z=" << simEvt.at(iev).z << " ntk=" << simEvt.at(iev).rtk.size()
           << "  nChTP=" << simEvt.at(iev).rtk.size() << endl;
    }


  }

  // some more vertex quantities
  for (vector<SimEvent>::iterator ev = simEvt.begin(); ev != simEvt.end(); ev++) {
    Fill(h, "nwosmatch", float(ev->nwosmatch));
    Fill(h, "Tc", ev->Tc, ev == simEvt.begin());
    Fill(h, "Chisq", ev->chisq, ev == simEvt.begin());
    if (ev->chisq > 0)
      Fill(h, "logChisq", log(ev->chisq), ev == simEvt.begin());
    Fill(h, "dzmax", ev->dzmax, ev == simEvt.begin());
    Fill(h, "dztrim", ev->dztrim, ev == simEvt.begin());
    Fill(h, "m4m2", ev->m4m2, ev == simEvt.begin());
    if (ev->Tc > 0) {
      Fill(h, "logTc", log(ev->Tc) / log(10.), ev == simEvt.begin());
    }

    // and the same for pairs of vertices
    for (vector<SimEvent>::iterator ev2 = ev + 1; ev2 != simEvt.end(); ev2++) {
      vector<RecoTrack> xt;
      if ((ev->rtkprimsel.size() > 0) && (ev2->rtkprimsel.size() > 0) &&
          (ev->rtkprimsel.size() + ev2->rtkprimsel.size()) > 1) {
        xt.insert(xt.end(), ev->rtkprimsel.begin(), ev->rtkprimsel.end());
        xt.insert(xt.end(), ev2->rtkprimsel.begin(), ev2->rtkprimsel.end());
        double xTc, xChsq, xDzmax, xDztrim, xm4m2;
        getTc(xt, xTc, xChsq, xDzmax, xDztrim, xm4m2);
        if (xTc > 0) {
          Fill(h, "xTc", xTc, ev == simEvt.begin());
          Fill(h, "logxTc", log(xTc) / log(10), ev == simEvt.begin());
          Fill(h, "xChisq", xChsq, ev == simEvt.begin());
          if (xChsq > 0) {
            Fill(h, "logxChisq", log(xChsq), ev == simEvt.begin());
          };
          Fill(h, "xdzmax", xDzmax, ev == simEvt.begin());
          Fill(h, "xdztrim", xDztrim, ev == simEvt.begin());
          Fill(h, "xm4m2", xm4m2, ev == simEvt.begin());
        }
      }
    }
  }

  // vertex pairs with >=4 tracks
  for (vector<SimEvent>::iterator ev1 = simEvt.begin(); ev1 != simEvt.end(); ev1++) {
    if (ev1->tk.size() < 4)
      continue;
    for (vector<SimEvent>::iterator ev2 = ev1 + 1; ev2 != simEvt.end(); ev2++) {
      if (ev2->tk.size() < 4)
        continue;
      double deltazsim = ev2->z - ev1->z;
      Fill(h, "zdiffsimallTP", deltazsim);

      if ((ev1->rec > 0) && (ev2->rec > 0)) {
        // both sim vertices of this pair were found
        Fill(h, "zdiffsimfoundTP", deltazsim);
        if ((ev2->ndof > 4.) && (ev1->ndof > 4.)) {
          Fill(h, "zdiffsimfound4TP", deltazsim);
          if (ev1 == simEvt.begin()) {
            Fill(h, "zdiffsimfound4SignalTP", deltazsim);
          }
        }
        double deltazrec = ev2->zrec - ev1->zrec;
        Fill(h, "zdiffrecvssimTP", deltazsim, deltazrec);
        Fill(h, "zdiffrecvsdsimTP", deltazsim, deltazrec - deltazsim);
        Fill(h, "zdiffrecvsdsimTPprof", deltazsim, deltazrec - deltazsim);
      }

      /* nothing about mergers, yet */
    }
  }

  timer_start("TP_track_vertex_loop");
  // track properties, fake and real
  for (unsigned int iv = 0; iv < recVtxs->size(); iv++) {
    for (trackit_t t = recVtxs->at(iv).tracks_begin(); t != recVtxs->at(iv).tracks_end(); t++) {
      auto tk = tracks.from_ref(*t);
      if (recvmatch[iv].is_real()) {
        if (t->get()->dzError() < 0.01) {
          fillTrackHistos(h, "realvtxdriver", tk, &(recVtxs->at(iv)));
        }
      } else {
        if (t->get()->dzError() < 0.01) {
          fillTrackHistos(h, "fakevtxdriver", tk, &(recVtxs->at(iv)));
        }
      }
    }
  }
  timer_stop("TP_track_vertex_loop");

  
  for (unsigned int iv = 0; iv < recVtxs->size(); iv++) {
    if (select(recVtxs->at(iv))) {
      analyzeRecVertexComposition(h, recVtxs->at(iv), tracks, recvmatch[iv], simEvt, npu);
    }
  }


  // unbiased track assignment efficiency
  // (unbiased in contrast to the conditional assignment efficiency where the existence of the rec vertex is required)
  for(unsigned int i = 0; i < tracks.size(); i++){
    auto tk = tracks(i);
    if (tk.matched){
      unsigned int iv = tk.simEvt->rec;
      auto numtk = tk.simEvt->rtk.size();
      
      double w = 0.;
      if ((iv != NOT_MATCHED ) && (select(recVtxs->at(iv))) ){
	w = 1.;
      }
      Fill(h, "utrkAssignmentEfficiency", w, tk.simEvt->is_signal());
      Fill(h, "utrkAssignmentEfficiencyvspu", npu, w, tk.simEvt->is_signal()); 
      Fill(h, "utrkAssignmentEfficiencyvsntrk", std::min(double(numtk), 50.5), w, tk.simEvt->is_signal());
      if (tk.is_primary){
	Fill(h, "uprimtrkAssignmentEfficiency", w, tk.simEvt->is_signal());
      }
    }
  }
}

/***************************************************************************************/





/***************************************************************************************/
void PrimaryVertexAnalyzer4PU::analyzeVertexCollectionPtvis(std::map<std::string, TH1*>& h,
                                                         const reco::VertexCollection* recVtxs,
                                                         Tracks& tracks,
                                                         vector<SimEvent>& simEvt,
                                                         std::vector<RSmatch>& recvmatch,
                                                         const string message)
/***************************************************************************************/
{
  unsigned int k = 0;
  for(; k < recvmatch.size(); k++){
    if(recvmatch[k].sim == 0) break;
  }
  if(k == recvmatch.size()){
    report_counted("analyzeVertexCollectionPtvis : Signal vertex not found", 1);
    return;
  }
  
  auto v = recVtxs->at(k);
  double zvtx = v.z();
  double density = simPU_ * exp(-0.5*pow((zvtx - vertexBeamSpot_.z0()) / vertexBeamSpot_.sigmaZ(),2))/sqrt(2*3.14159)/vertexBeamSpot_.sigmaZ(); 

  vector<double> ptvis = {0,0,0,0,0};
  vector<double> px =  {0,0,0,0,0};
  vector<double> py =  {0,0,0,0,0};
  

  for (unsigned int i = 0; i < tracks.size(); i++)
    {
    
      RecoTrack tk = tracks(i);
      double tkpx  = tk.pt * cos(tk.phi);
      double tkpy  = tk.pt * sin(tk.phi);
      
      if (fabs(tk.z - zvtx) < 0.1){
	ptvis[0] += tk.pt;
	px[0] += tkpx;
	py[0] += tkpy;
      }
      if (fabs(tk.z - zvtx) < 0.2){
	ptvis[1] += tk.pt;
	px[1] += tkpx;
	py[1] += tkpy;
      }
      double zpull = (tk.z - zvtx) / tk.dz;
      if (fabs(zpull) < 3.){
	ptvis[2] += tk.pt;
	px[2] += tkpx;
	py[2] += tkpy;
      }

      double w3 = 1;
      double Z4 = 0;
      for(unsigned int iv=0; iv < recVtxs->size(); iv++)
	{
	  if(select(recVtxs->at(iv))){
	    if (fabs(recVtxs->at(iv).z() -tk.z) <  fabs(zvtx - tk.z)){
	      w3 = 0;
	    }
	      
	    double p = (recVtxs->at(iv).z() - tk.z) / tk.dz;
	    if (fabs(p) < 5.)
	      {
		Z4 +=  exp(-0.5* p * p);
	      }
	  }
	}

      if (w3 > 0){// nearest
	ptvis[3] += tk.pt;
	px[3] += tkpx;
	py[3] += tkpy;
      }
      
      if ((Z4 > 0) && (zpull < 5)){
	ptvis[4] += exp(-0.5*zpull*zpull) / Z4 * tk.pt;
	px[4] += exp(-0.5*zpull*zpull) / Z4 * tkpx;
	py[4] += exp(-0.5*zpull*zpull) / Z4 * tkpy;
      }
      
    }// tracks

  //now histogram ratio to "true" ptvs vs pu and z
  double ptvis0 = simEvt[0].ptvis;

  for(unsigned int i = 0; i < 5; i++){
    Fill(h, Form("ptvis/ptvis%d", i), ptvis[i] / ptvis0);
    Fill(h, Form("ptvis/ptvis%dvsz", i), zvtx,ptvis[i] / ptvis0 );
    Fill(h, Form("ptvis/ptvis%dvspu", i), simPU_, ptvis[i] / ptvis0);
    Fill(h, Form("ptvis/ptvis%dvsdensity", i), density,  ptvis[i] / ptvis0); 


    double ptmiss = sqrt(px[i] * px[i] + py[i] * py[i]);
    double ptmiss_sim = sqrt(pow(simEvt[0].pxvis, 2) + pow(simEvt[0].pyvis, 2));
    Fill(h, Form("ptvis/ptmiss_sim%d", i), ptmiss_sim);
    Fill(h, Form("ptvis/ptmiss%d", i), ptmiss);
    Fill(h, Form("ptvis/ptmiss%dvsz", i), zvtx, ptmiss );
    Fill(h, Form("ptvis/ptmiss%dvspu", i), simPU_, ptmiss);
    Fill(h, Form("ptvis/ptmiss%dvsdensity", i), density,  ptmiss);
  }

}
/***************************************************************************************/



/***************************************************************************************/

void PrimaryVertexAnalyzer4PU::analyzeVertexCollectionDQMMC(std::map<std::string, TH1*>& h,
                                                            const reco::VertexCollection* recVtxs,
                                                            std::vector<simPrimaryVertex>& simpv,
                                                            std::vector<SimPart>& tsim,
                                                            const std::string message)
/*
The logic for the merge rate goes as follows.

Reco and sim vertices are currently matched if their z distance is < 3*sigma_z and < 1 mm (sigma_z = reco uncertainty).

Reco vertex is classified as merged if it is matched to >= 2 sim vertices.

             N(merged reco vertices)
Merge rate = -----------------------
               N(all reco vertices)

 */
{
  // dqm style z-matching
  std::vector<std::vector<std::pair<double, unsigned int>>> r2s;
  std::vector<std::vector<std::pair<double, unsigned int>>> s2r;
  for (unsigned int j = 0; j < simpv.size(); j++) {
    std::vector<std::pair<double, unsigned int>> v;
    s2r.push_back(v);
  }

  for (unsigned int k = 0; k < recVtxs->size(); k++) {
    std::vector<std::pair<double, unsigned int>> v;
    r2s.push_back(v);
  }

  for (unsigned int k = 0; k < recVtxs->size(); k++) {
    if (recVtxs->at(k).ndof() < 4.)
      continue;
    double zrec = recVtxs->at(k).z();
    double sigmaz = recVtxs->at(k).zError();
    double dzmatch = std::min(0.1, 3 * sigmaz);
    for (unsigned int j = 0; j < simpv.size(); j++) {
      double dz = std::abs(simpv[j].z - zrec);
      if (dz < dzmatch) {
        s2r[j].push_back(make_pair(dz, k));
        r2s[k].push_back(make_pair(dz, j));
      }
    }
  }

  for (unsigned int j = 0; j < simpv.size(); j++) {
    std::stable_sort(s2r[j].begin(), s2r[j].end(), std::less<std::pair<double, unsigned int>>());
  }

  for (unsigned int k = 0; k < recVtxs->size(); k++) {
    std::stable_sort(r2s[k].begin(), r2s[k].end(), std::less<std::pair<double, unsigned int>>());
  }

  // now fill
  for (unsigned int k = 0; k < recVtxs->size(); k++) {
    if (recVtxs->at(k).ndof() < 4.)
      continue;
    double zfill = 0;
    if (r2s[k].size() > 0) {
      zfill = simpv[r2s[k][0].second].closest_vertex_distance_z;
    }

    if (zfill > 0) {
      Fill(h, "mergerate_denominator_dqm", zfill);
      Fill(h, "mergerate_denominator_lin", zfill);
    }

    if ((zfill > 0) && (r2s[k].size() > 1)) {
      Fill(h, "mergerate_numerator_dqm", zfill);
      Fill(h, "mergerate_numerator_lin", zfill);
    }
  }

  for (unsigned int j1 = 0; j1 < simpv.size() - 1; j1++) {
    if (s2r[j1].size() == 0)
      continue;
    for (unsigned int j2 = j1 + 1; j2 < simpv.size(); j2++) {
      if (s2r[j2].size() == 0)
        continue;
      double dz = std::abs(simpv[j1].z - simpv[j2].z);
      Fill(h, "mergerate_denominator", dz);  // two sim vertices compatible with a rec vertex
      if (s2r[j1][0].second == s2r[j2][0].second) {
        Fill(h, "mergerate_numerator", dz);  // both matched to the same vertex
      }
    }
  }
}
/***************************************************************************************/

/***************************************************************************************
*  for samples with rec tracks, but without sim tracks (for PU)                        *
***************************************************************************************/
void PrimaryVertexAnalyzer4PU::analyzeVertexCollectionSimPvNoSimTracks(std::map<std::string, TH1*>& h,
                                                                       const reco::VertexCollection* recVtxs,
                                                                       Tracks& tracks,
                                                                       std::vector<simPrimaryVertex>& simpv,
                                                                       const std::string message) {
  if(verbose_){
    cout << "PrimaryVertexAnalyzer4PU::analyzeVertexCollectionSimPvNoTracks " << message << endl;
    cout << "analyzeVertexCollectionSimPv simpv.size = " << simpv.size() << endl;
  }
  Fill(h, "simPU", float(simPU_));

  //
  for (unsigned int iv = 0; iv < recVtxs->size(); iv++) {
    reco::Vertex vrec = recVtxs->at(iv);
    if (!select(vrec, 0))
      continue;
    double ptmax2 = vertex_ptmax2(vrec);

    // for this reconstructed primary find the distance to the nearest sim vertex
    double dzmin = 1e100;
    for (unsigned int idx = 1; idx < simpv.size(); idx++) {
      double dz = simpv[idx].z - vrec.z();
      if (std::abs(dz) < std::abs(dzmin)) {
        dzmin = dz;
      }

      // also fill dz histos
      Fill(h, "dzrecsim", dz);
      if (ptmax2 < 0.4) {
        Fill(h, "dzrecsimptmax2lt04", dz);
      } else {
        Fill(h, "dzrecsimptmax2gt04", dz);
      }
    }

    Fill(h, "dzrecsimmin", dzmin);
    if (ptmax2 < 0.4) {
      Fill(h, "dzrecsimminptmax2lt04", dzmin);
    } else {
      Fill(h, "dzrecsimminptmax2gt04", dzmin);
    }

    if (std::abs(dzmin) < 0.1) {
      Fill(h, "zmatched01", vrec.z());
      Fill(h, "ptmax2matched01", ptmax2);
    }
    if (std::abs(dzmin) < 0.2) {
      Fill(h, "zmatched02", vrec.z());
      Fill(h, "ptmax2matched02", ptmax2);
    }

    // fake?, no real vertex within 1 (or 2 ) mm
    if (std::abs(dzmin) > 0.1) {
      Fill(h, "zunmatched01", vrec.z());
      Fill(h, "ptmax2unmatched01", ptmax2);
      if (ptmax2 < 0.4) {
        Fill(h, "zunmatched01ptmax2lt04", vrec.z());
      }
    }
    if (std::abs(dzmin) > 0.2) {
      Fill(h, "zunmatched02", vrec.z());
      Fill(h, "ptmax2unmatched02", ptmax2);
      if (ptmax2 < 0.4) {
        Fill(h, "zunmatched02ptmax2lt04", vrec.z());
      }
    }
  }

  // fill track property histograms for selected rec tracks that are far away from any simulated vertex
  //  for(auto t = recTrks.begin(); t !=recTrks.end(); t++)
  for (unsigned int i = 0; i < tracks.size(); i++) {
    RecoTrack& tk = tracks(i);
    if (theTrackFilter(*(tk.tt))) {
      double z0 = tk.z;
      double dzmin = 1e10;
      for (unsigned int idx = 1; idx < simpv.size(); idx++) {
        double dz = simpv[idx].z - z0;
        if (std::abs(dz) < std::abs(dzmin)) {
          dzmin = dz;
        }
      }
      if (std::abs(dzmin) > 0.3) {
        fillTrackHistos(h, "faraway", tk);
      }
    }
  }
}
/***************************************************************************************/


/***************************************************************************************/

void PrimaryVertexAnalyzer4PU::analyzeVertexCollectionSimPv(std::map<std::string, TH1*>& h,
                                                            const reco::VertexCollection* recVtxs,
                                                            Tracks& tracks,
                                                            std::vector<simPrimaryVertex>& simpv,
                                                            std::vector<SimPart>& tsim,
                                                            const std::string message) {
  /* 
     not using tracking particles 
     kind of obsolete, revive for samples without TP
   */
  if(verbose_){
    cout << "PrimaryVertexAnalyzer4PU::analyzeVertexCollectionSimPv " << message << endl;
    cout << " eventSummaryCounter_= " << eventSummaryCounter_ << "  nEventSummary_ =" << nEventSummary_ << endl;
    cout << "analyzeVertexCollectionSimPv simpv.size = " << simpv.size() << endl;
  }


  int nrectrks = tracks.size();
  int nrecvtxs = recVtxs->size();

  if (simpv.size() == 0) {
    return;
  }

  int npu0 = simpv.size();
  Fill(h, "npu", float(npu0));
  Fill(h, "nsimtrkSimpv", float(tsim.size()));

  bool debug = false;

  // find the reco vertex matching the signal via truth matched tracks (don't have those for pu)
  int* rectosim = NULL;
  if (MC_)
    rectosim = supf(tsim, *tracks.trackCollectionH.product());
  double zsimsignal = simpv[0].z;

  vector<double> wos0(nrecvtxs);  // truth matched signal tracks contributing to a vertex (to identify splitting)
  vector<double> wos1(nrecvtxs);  // only tracks from primary particles -> for identifying the signal vertex
  for (unsigned int it = 0; it < tracks.size(); it++) {
    auto t = tracks.ref(it);
    if (rectosim[it] >= 0) {
      unsigned int iv = trkkey2recvtx_[t.key()];  // FIXME, store this in the RecoTrack instead
      double dz2 = pow(t.get()->dzError(), 2);
      double wos = recVtxs->at(iv).trackWeight(t) / dz2;
      wos0[iv] += wos;
      if (tsim[rectosim[it]].type == 0) {
        wos1[iv] += wos;
      }
    }
  }

  vector<double> sumwos(nrecvtxs);  // all tracks contributing to a vertex
  for (unsigned int iv = 0; iv < recVtxs->size(); iv++) {
    reco::Vertex vrec = recVtxs->at(iv);
    if (vrec.isFake())
      continue;
    for (trackit_t tv = vrec.tracks_begin(); tv != vrec.tracks_end(); tv++) {
      sumwos[iv] += vrec.trackWeight(*tv) / pow(tv->get()->dzError(), 2);
    }
    if (sumwos[iv] > 0) {
      wos0[iv] /= sumwos[iv];
    } else {
      wos0[iv] = 0;  // should be anyway
    }
  }

  // find the signal recvertex as the one with the highest contribution from primary signal tracks
  // require that those tracks make a significant contribution to the vertex
  // when many non-primary or non-matched tracks are present, this can fail,
  // in that case accept the rec vertex when it is closer than 100 um to the sim vertex
  // or should this be the highest relative contribution ? or some combination?
  unsigned int vsignal = 0;
  bool found_signal = false;
  double wosmax = 0;
  for (unsigned int iv = 0; iv < recVtxs->size(); iv++) {
    if ((wos1[iv] > wosmax)) {
      reco::Vertex vrec = recVtxs->at(iv);
      if ((wos1[iv] > 0.5 * sumwos[iv]) ||
          ((wos1[iv] > 0.1 * sumwos[iv]) && (std::abs(vrec.z() - zsimsignal) < 4 * vrec.zError()))) {
        vsignal = iv;
        wosmax = wos1[iv];
        found_signal = true;
      }
    }
  }

  //
  for (unsigned int iv = 0; iv < recVtxs->size(); iv++) {
    reco::Vertex vrec = recVtxs->at(iv);
    if (vrec.isFake())
      continue;
    if (debug) {
      cout << "VVVVV " << setw(10) << setprecision(4) << fixed << vrec.z() << " " << setw(10) << setprecision(1)
           << wos0[iv] * sumwos[iv] << " " << setw(10) << sumwos[iv] << "  " << setprecision(3) << wos0[iv] << "  "
           << setprecision(3) << wos1[iv] / sumwos[iv] << "  ndof=" << vrec.ndof();
      if ((found_signal) && (iv == vsignal)) {
        cout << "  Signal" << endl;
      } else {
        cout << endl;
      }
    }
  }

  double zrecsignal = 999.;
  double ntrkacc = simpv[0].nTrk;
  
  if (found_signal) {
    zrecsignal = recVtxs->at(vsignal).z();
    Fill(h, "effSignal", 1.);
    Fill(h, "effSignalvsNpu0", npu0, 1.);
    Fill(h, "effSignalvsnsimtrkacc", ntrkacc, 1.);
    if (recVtxs->at(vsignal).ndof() > 4) {
      //fillVertexHistos(h, "signalvtx4", &(recVtxs->at(vsignal)), tracks);
      Fill(h, "effSignal4", 1.);
      Fill(h, "effSignal4vsNpu0", npu0, 1.);
      Fill(h, "effSignal4vsnsimtrkacc", ntrkacc, 1.);
    } else {
      reportEvent(
          Form("signal vertex not selected (z=%8.4f    ndof=%f)", recVtxs->at(vsignal).z(), recVtxs->at(vsignal).ndof()),
          true);
      Fill(h, "effSignal4", 0.);
      Fill(h, "effSignal4vsNpu0", npu0, 0.);
      Fill(h, "effSignalvsnsimtrkacc", ntrkacc, 0.);
      Fill(h, "effSignal4vsnsimtrkacc", ntrkacc, 0.);
    }
  } else {
    found_signal = false;
    if ((simpv[0].type == 1) && (simpv[0].nGenTrk == 0)) {
      // neutrino gun
    } else {
      reportEvent(Form("signal vertex not matched at all, type=%d, nGenTrk=%d", simpv[0].type, simpv[0].nGenTrk), true);
      Fill(h, "effSignal", 0.);
      Fill(h, "effSignalvsNpu0", npu0, 0.);
      Fill(h, "effSignal4", 0.);
      Fill(h, "effSignal4vsNpu0", npu0, 0.);
    }
  }

  for (unsigned int iv = 0; iv < recVtxs->size(); iv++) {
    if ((!(iv == vsignal) && (wos0[iv] > 0.5) && (recVtxs->at(iv).ndof() > 4))) {
      cout << "split or split/merge vertex z=" << setw(10) << setprecision(4) << recVtxs->at(iv).z()
           << "  wos0 =" << setprecision(3) << wos0[iv] << endl;  // todo, make sure there is no real vertex right there
    }
  }

  // match in z, like printEventSummary does, but pull out the signal
  vector<pair<double, unsigned int>> zrecv;
  for (unsigned int idx = 0; idx < recVtxs->size(); idx++) {
    if (found_signal && (idx == vsignal))
      continue;
    zrecv.push_back(make_pair(recVtxs->at(idx).z(), idx));
  }
  stable_sort(zrecv.begin(), zrecv.end(), lt);

  // same for simulated vertices, don't include the signal
  vector<pair<double, unsigned int>> zsimv;
  for (unsigned int idx = 1; idx < simpv.size(); idx++) {
    zsimv.push_back(make_pair(simpv[idx].z, idx));
  }
  stable_sort(zsimv.begin(), zsimv.end(), lt);

  unsigned int idxrec = 0;
  unsigned int idxsim = 0;
  double zmatch = 0.05;
  int nsplit0 = 0;
  int nsplit4 = 0;
  int nfake4 = 0;
  int nfake0 = 0;
  int nlostpu = 0;
  int nfoundpu4 = 0;
  int nfoundpu0 = 0;

  while ((idxrec < zrecv.size()) || (idxsim < zsimv.size())) {
    double ndof = 0;
    double pxy = 1.;
    bool sel = false;
    if (idxrec < zrecv.size()) {
      sel = select(recVtxs->at(zrecv[idxrec].second));
      ndof = recVtxs->at(zrecv[idxrec].second).ndof();
      pxy = vertex_pxy(recVtxs->at(zrecv[idxrec].second));
    }

    double dznearestrec = 0.;
    if (idxrec < zrecv.size()) {
      dznearestrec = 1000.;
      if ((idxrec > 0) && (std::abs(zrecv[idxrec].first - zrecv[idxrec - 1].first) < std::abs(dznearestrec))) {
        dznearestrec = zrecv[idxrec].first - zrecv[idxrec - 1].first;
      }
      if (((idxrec + 1) < zrecv.size()) &&
          (std::abs(zrecv[idxrec].first - zrecv[idxrec + 1].first) < std::abs(dznearestrec))) {
        dznearestrec = zrecv[idxrec].first - zrecv[idxrec + 1].first;
      }
    }

    if ((idxrec < zrecv.size()) && (idxsim < zsimv.size()) &&
        (abs(zrecv[idxrec].first - zsimv[idxsim].first) < (zmatch + recVtxs->at(zrecv[idxrec].second).zError())) &&
        (((idxsim + 1) == simpv.size()) ||
         (fabs(zrecv[idxrec].first - zsimv[idxsim].first) < fabs(zrecv[idxrec].first - zsimv[idxsim + 1].first))) &&
        (((idxrec + 1) == recVtxs->size()) ||
         (fabs(zrecv[idxrec].first - zsimv[idxsim].first) < fabs(zrecv[idxrec + 1].first - zsimv[idxsim].first)))) {
      nfoundpu0++;
      //fillVertexHistos(h, "puvtx0", &(recVtxs->at(zrecv[idxrec].second)), tracks, dznearestrec);

      if (sel) {
        nfoundpu4++;
        Fill(h, "zdiffrecselsignalrealSimpv", std::abs(zrecv[idxrec].first - zrecsignal));
        Fill(h, "zdiffrecselsignalrealSimpv-dzbin", std::abs(zrecv[idxrec].first - zrecsignal));
        if (pxy > 0.01)
          Fill(h, "zdiffrecselpxysignalrealSimpv", std::abs(zrecv[idxrec].first - zrecsignal));
        if (pxy > 0.01)
          Fill(h, "zdiffrecselpxysignalrealSimpv-dzbin", std::abs(zrecv[idxrec].first - zrecsignal));
        Fill(h, "zdiffsimfoundselsignalSimpva", std::abs(zsimv[idxsim].first - zsimsignal));
        Fill(h, "zdiffsimfoundselsignalSimpv", zsimv[idxsim].first - zsimsignal);
        //fillVertexHistos(h, "puvtx4", &(recVtxs->at(zrecv[idxrec].second)), tracks, dznearestrec);
      }
      if (debug) {
        cout << "matched       " << setw(10) << setprecision(4) << zrecv[idxrec].first << endl;
      }
      idxrec++;
      idxsim++;

    } else if (((idxrec < zrecv.size()) && (idxsim < zsimv.size()) && (zrecv[idxrec].first < zsimv[idxsim].first)) ||
               ((idxrec < zrecv.size()) && (idxsim == zsimv.size()))) {
      // recvertex does not match any: fake or split
      if (found_signal && (sel)) {
        Fill(h, "zdiffrecselsignalfakeSimpv", std::abs(zrecv[idxrec].first - zrecsignal));
        Fill(h, "zdiffrecselsignalfakeSimpv-dzbin", std::abs(zrecv[idxrec].first - zrecsignal));
        if (pxy > 0.01)
          Fill(h, "zdiffrecselpxysignalfakeSimpv", std::abs(zrecv[idxrec].first - zrecsignal));
        if (pxy > 0.01)
          Fill(h, "zdiffrecselpxysignalfakeSimpv-dzbin", std::abs(zrecv[idxrec].first - zrecsignal));
        if (wos0[zrecv[idxrec].second] > 0.2) {
          if (debug) {
            reportEvent(Form("split from signal  z=%10.4f", zrecv[idxrec].first), verbose_);
          }
          //double deltaz = zrecv[idxrec].first - zrecsignal;
          //fillVertexHistos(h, "splitvtx4", &(recVtxs->at(zrecv[idxrec].second)), tracks, deltaz, verbose_);
          Fill(h, "zdiffrecselsignalsplitSimpv", std::abs(zrecv[idxrec].first - zrecsignal));
          Fill(h, "zdiffrecselsignalsplitSimpv-dzbin", std::abs(zrecv[idxrec].first - zrecsignal));
          Fill(h, "ndofsignalsplitSimpv", ndof);
          if (pxy > 0.01) {
            Fill(h, "zdiffrecselpxysignalsplitSimpv", std::abs(zrecv[idxrec].first - zrecsignal));
            Fill(h, "zdiffrecselpxysignalsplitSimpv-dzbin", std::abs(zrecv[idxrec].first - zrecsignal));
            Fill(h, "ndofpxysignalsplitSimpv", ndof);
          }

        } else {
          nfake4++;
	  //          fillVertexHistos(h, "otherfakevtx", &(recVtxs->at(zrecv[idxrec].second)), tracks);
          if (debug) {
            cout << "fake          " << setw(10) << setprecision(4) << zrecv[idxrec].first << endl;
          }
        }
      }

      if (found_signal && (ndof > 0.)) {
        nsplit0++;
        //double deltaz = zrecv[idxrec].first - zrecsignal;
        if (wos0[zrecv[idxrec].second] > 0.2) {
          //fillVertexHistos(h, "splitvtx0", &(recVtxs->at(zrecv[idxrec].second)), tracks, deltaz);
        } else {
          nfake0++;
        }
      }

      if (idxrec < zrecv.size()) {
        idxrec++;
      } else {
        reportEvent("FIXME : you should never be here", false);
      }

    } else if (((idxrec < zrecv.size()) && (idxsim < zsimv.size()) && (zrecv[idxrec].first > zsimv[idxsim].first)) ||
               ((idxrec == zrecv.size()) && (idxsim < zsimv.size()))) {
      // lost
      nlostpu++;
      idxsim++;
    } else {
      cout << "oooops" << endl;
      break;
    }
  }

  if (MC_)
    delete[] rectosim;

  Fill(h, "nsplitvtx0vsnpuprof", npu0, float(nsplit0));
  Fill(h, "nsplitvtx4vsnpuprof", npu0, float(nsplit4));
  Fill(h, "nfakevtx0vsnpuprof", npu0, float(nfake0));
  Fill(h, "nfakevtx4vsnpuprof", npu0, float(nfake4));
  Fill(h, "nlostpuvtxvsnpuprof", npu0, float(nlostpu));
  Fill(h, "nfoundpuvtx0vsnpuprof", npu0, float(nfoundpu0));
  Fill(h, "nfoundpuvtx4vsnpuprof", npu0, float(nfoundpu4));

  // vertex matching and efficiency bookkeeping
  for (std::vector<simPrimaryVertex>::iterator vsim = simpv.begin(); vsim != simpv.end(); vsim++) {
    bool is_signal = vsim == simpv.begin();

    Fill(h, "nsimtrk", static_cast<double>(vsim->nGenTrk), is_signal);          // nsimtrkSignal or nsimtrkPU
    Fill(h, "nsimtrkacc", static_cast<double>(vsim->nTrk), is_signal);          // nsimtrkSignal or nsimtrkPU
    Fill(h, "nsimtrkaccprim", static_cast<double>(vsim->nTrkPrim), is_signal);  // nsimtrkSignal or nsimtrkPU
    Fill(h, "nsimtrkaccsec", static_cast<double>(vsim->nTrkSec), is_signal);    // nsimtrkSignal or nsimtrkPU

    // look for a matching reconstructed vertex
    vsim->recVtx = NULL;

    // find the nearest recvertex  (multiple sims may be mapped to the same rec)
    for (auto vrec = recVtxs->begin(); vrec != recVtxs->end(); ++vrec) {
      if (!(vrec->isFake())) {
        if ((vsim->recVtx == NULL) || ((fabs(vsim->recVtx->position().z() - vsim->z) > fabs(vrec->z() - vsim->z)))) {
          vsim->recVtx = &(*vrec);
        }
      }
    }

    if (vsim->recVtx == NULL) {
      cout << "PrimaryVertexAnalyzer4PU:analyzerVertexCollectionSimPV : no recvtx ? " << endl;
      return;
    }

    // histogram properties of matched vertices
    if (vsim->recVtx && (fabs(vsim->recVtx->z() - vsim->z * simUnit_) < zmatch_)) {
      double w4 = vsim->recVtx->ndof() > 4 ? 1. : 0;

      Fill(h, "matchedVtxNdof", vsim->recVtx->ndof());
      Fill(h, "resxvsNdofprof", vsim->recVtx->ndof(), vsim->recVtx->x() - vsim->x * simUnit_);
      Fill(h, "resyvsNdofprof", vsim->recVtx->ndof(), vsim->recVtx->y() - vsim->y * simUnit_);
      Fill(h, "resxvsNdofSpread", vsim->recVtx->ndof(), vsim->recVtx->x() - vsim->x * simUnit_);
      Fill(h, "resyvsNdofSpread", vsim->recVtx->ndof(), vsim->recVtx->y() - vsim->y * simUnit_);

      // residuals an pulls with respect to simulated vertex
      if (vsim->recVtx->ndof() > 4) {
        Fill(h, "resx", vsim->recVtx->x() - vsim->x * simUnit_, is_signal);
        Fill(h, "resy", vsim->recVtx->y() - vsim->y * simUnit_, is_signal);
        Fill(h, "resz", vsim->recVtx->z() - vsim->z * simUnit_, is_signal);
        Fill(h, "resz10", vsim->recVtx->z() - vsim->z * simUnit_, is_signal);
        Fill(h, "pullx", (vsim->recVtx->x() - vsim->x * simUnit_) / vsim->recVtx->xError(), is_signal);
        Fill(h, "pully", (vsim->recVtx->y() - vsim->y * simUnit_) / vsim->recVtx->yError(), is_signal);
        Fill(h, "pullz", (vsim->recVtx->z() - vsim->z * simUnit_) / vsim->recVtx->zError(), is_signal);
      }

      if (vsim->recVtx->ndof() > 50) {
        Fill(h, "resx50", vsim->recVtx->x() - vsim->x * simUnit_);
        Fill(h, "resy50", vsim->recVtx->y() - vsim->y * simUnit_);
        Fill(h, "resz50", vsim->recVtx->z() - vsim->z * simUnit_);
      }

      // efficiency with zmatch within 500 um (or whatever zmatch is)
      Fill(h, "eff0", 1., is_signal);
      Fill(h, "eff4", w4, is_signal);

      Fill(h, "effvsnrectrk", nrectrks, 1.);
      Fill(h, "effvsz", vsim->z * simUnit_, 1.);
      Fill(h, "effvsz2", vsim->z * simUnit_, 1.);
      Fill(h, "effvspt_hat", vsim->pt_hat, 1.);
      if (is_signal) {
        Fill(h, "effvsnsimvtx", simpv.size(), 1.);
        Fill(h, "effvsnrecvtx", nrecvtxs, 1.);
      }

      if (vsim->type == 1) {  // full (i.e. not just PUInfo)
        Fill(h, "effvsnsimtrk", vsim->nGenTrk, 1.);
        Fill(h, "effvsptsq", vsim->ptsq, 1.);
        Fill(h, "effvsr", sqrt(vsim->x * vsim->x + vsim->y * vsim->y) * simUnit_, 1.);
      }

    } else {  // no matching rec vertex found for this simvertex

      bool plapper = veryverbose_ && vsim->nGenTrk;
      if (plapper) {
        // be quiet about invisble vertices
        std::cout << "primary not found " << message << " " << eventcounter_ << "  x=" << vsim->x * simUnit_
                  << "  y=" << vsim->y * simUnit_ << " z=" << vsim->z * simUnit_ << " nGenTrk=" << vsim->nGenTrk
                  << std::endl;
      }

      //int mistype=0;
      if (vsim->recVtx) {
        if (plapper) {
          std::cout << "nearest recvertex at " << vsim->recVtx->z()
                    << "   dz=" << vsim->recVtx->z() - vsim->z * simUnit_ << std::endl;
        }

        if (fabs(vsim->recVtx->z() - vsim->z * simUnit_) < 0.2) {
          Fill(h, "effvsz2", vsim->z * simUnit_, 1.);
        }

        if (fabs(vsim->recVtx->z() - vsim->z * simUnit_) < 0.5) {
          if (plapper) {
            std::cout << "type 1, lousy z vertex" << std::endl;
          }
          Fill(h, "zlost1", vsim->z * simUnit_, 1.);
          //mistype=1;
        } else {
          if (plapper) {
            std::cout << "type 2a no vertex anywhere near" << std::endl;
          }
          //mistype=2;
        }
      } else {  // no recVtx at all
        //mistype=2;
        if (plapper) {
          std::cout << "type 2b, no vertex at all" << std::endl;
        }
      }

      Fill(h, "eff0", 0., is_signal);
      Fill(h, "eff4", 0., is_signal);

      Fill(h, "effvsnsimtrk", float(vsim->nGenTrk), 0.);
      Fill(h, "effvsnsimtrkacc", float(vsim->nGenTrk), 0.);
      Fill(h, "effvsnrectrk", nrectrks, 0.);
      Fill(h, "effvspt_hat", vsim->pt_hat, 0.);
      Fill(h, "effvsz", vsim->z * simUnit_, 0.);
      if (is_signal) {
        Fill(h, "effvsnsimvtx", simpv.size(), 0.);
        Fill(h, "effvsnrecvtx", nrecvtxs, 0.);
      }
      if (vsim->type == 1) {
        Fill(h, "effvsptsq", vsim->ptsq, 0.);
        Fill(h, "effvsr", sqrt(vsim->x * vsim->x + vsim->y * vsim->y) * simUnit_, 0.);
      }

    }  // no recvertex for this simvertex

  }  // vsim loop

  int nrecvtxs4 = 0;
  for (reco::VertexCollection::const_iterator vrec = recVtxs->begin(); vrec != recVtxs->end(); ++vrec) {
    if (!(vrec->isFake()) && (vrec->ndof() > 4)) {
      nrecvtxs4++;
    }
  }
  Fill(h, "nrecvsnpu", float(npu0), float(nrecvtxs));
  Fill(h, "nrec4vsnpu", float(npu0), float(nrecvtxs4));
  Fill(h, "nrec4vsnpuprof", float(npu0), float(nrecvtxs4));

  // look for rec vertices with no matching sim vertex
  for (reco::VertexCollection::const_iterator vrec = recVtxs->begin(); vrec != recVtxs->end(); ++vrec) {
    simPrimaryVertex* match = NULL;
    double zmax = zmatch_;
    if ((3 * vrec->zError()) > zmatch_)
      zmax = 3 * vrec->zError();

    if (!(vrec->isFake())) {
      for (std::vector<simPrimaryVertex>::iterator vsim = simpv.begin(); vsim != simpv.end(); vsim++) {
        if ((vsim->recVtx == &(*vrec)) &&
            ((match == NULL) || (fabs(vrec->position().z() - vsim->z) < fabs(vrec->position().z() - match->z)))) {
          match = &(*vsim);
        }
      }

      if ((match == NULL) || (fabs(vrec->position().z() - match->z) > zmax)) {
        Fill(h, "fakeVtxZ", vrec->z());
        if (vrec->ndof() >= 0.5)
          Fill(h, "fakeVtxZNdofgt05", vrec->z());
        if (vrec->ndof() >= 2.0)
          Fill(h, "fakeVtxZNdofgt2", vrec->z());
        if (vrec->ndof() >= 4.0)
          Fill(h, "fakeVtxZNdofgt4", vrec->z());
        if (vrec->ndof() >= 8.0)
          Fill(h, "fakeVtxZNdofgt8", vrec->z());
        Fill(h, "fakeVtxNdof", vrec->ndof());
        Fill(h, "fakeVtxNtrk", vrec->tracksSize());
      }
    }
  }

  // compare the signal vertex with the nearest rec vertex
  double deltaznearest = 9999.;
  int indexnearest = -1, idx = 0;
  for (reco::VertexCollection::const_iterator vrec = recVtxs->begin(); vrec != recVtxs->end(); ++vrec) {
    //if( !(vrec->isFake()) && (vrec->ndof()>4) ) {
    if (!(vrec->isFake()) && (vrec->ndof() > 0)) {
      Double_t dz = vrec->z() - (*simpv.begin()).z * simUnit_;
      if (fabs(dz) < fabs(deltaznearest)) {
        deltaznearest = dz;
        indexnearest = idx;
      }
    }
    idx++;
  }

  Fill(h, "zdistancenearest", deltaznearest);
  Fill(h, "abszdistancenearest", fabs(deltaznearest));
  Fill(h, "indexnearest", float(indexnearest));

  if (eventSummaryCounter_++ < nEventSummary_) {
    printEventSummary(h, recVtxs, tracks, simpv, message);
  }

  // isolated simulated vertices and vertex pairs

  // sort sim vertices in z, now include the signal
  //  vector< pair<double,unsigned int> >  zsimv;
  zsimv.clear();
  for (unsigned int idx = 0; idx < simpv.size(); idx++) {
    zsimv.push_back(make_pair(simpv[idx].z, idx));
  }
  stable_sort(zsimv.begin(), zsimv.end(), lt);

  // pairs of simulated vertices vs pairs of reconstructed vertices
  if (simpv.size() > 1) {
    for (std::vector<simPrimaryVertex>::iterator vsim1 = simpv.begin(); vsim1 != (simpv.end() - 1); vsim1++) {
      if (vsim1->nGenTrk < 4)
        continue;
      for (std::vector<simPrimaryVertex>::iterator vsim2 = vsim1 + 1; vsim2 != simpv.end(); vsim2++) {
        if ((vsim2->nGenTrk < 4))
          continue;
        double deltazsim = vsim2->z - vsim1->z;
        Fill(h, "zdiffsimall", deltazsim);
        if ((vsim1->recVtx == NULL) || (vsim2->recVtx == NULL))
          continue;
        double deltazrec = vsim2->recVtx->position().z() - vsim1->recVtx->position().z();
        if (vsim2->recVtx == vsim1->recVtx) {
          // merged or lost for some other reason
          Fill(h, "zdiffsimmerge", deltazsim);
        } else {
          // separated
          Fill(h, "zdiffrecvssim", deltazsim, deltazrec);
          if (select(*(vsim1->recVtx)) && select(*(vsim2->recVtx)))
            Fill(h, "zdiffrecselvssim", deltazsim, deltazrec);
          if ((vsim1->recVtx->ndof() > 12) && (vsim2->recVtx->ndof() > 12))
            Fill(h, "zdiffrec12vssim", deltazsim, deltazrec);
          Fill(h, "zdiffsimfound", deltazsim);
          Fill(h, "zdiffrecvsdsim", deltazsim, deltazrec - deltazsim);
          Fill(h, "zdiffrecvsdsimprof", deltazsim, deltazrec - deltazsim);
        }
      }
    }

    // look for isolated pairs of simvertices, then count rec vertices inside the isolated interval

    double ziso = 0.5;
    for (unsigned int idxsim = 0; idxsim < simpv.size() - 1; idxsim++) {
      if (((idxsim == 0) || ((idxsim > 0) && (zsimv[idxsim].first - zsimv[idxsim - 1].first > ziso))) &&
          ((idxsim + 1 >= zsimv.size() - 1) ||
           ((idxsim + 1 < zsimv.size() - 1) && (zsimv[idxsim + 2].first - zsimv[idxsim + 1].first > ziso)))) {
        if ((simpv[zsimv[idxsim].second].nGenTrk > 4) && (simpv[zsimv[idxsim + 1].second].nGenTrk > 4)) {
          double dzsim = zsimv[idxsim + 1].first - zsimv[idxsim].first;
          double zmin = zsimv[idxsim].first - 0.5 * ziso;
          double zmax = zsimv[idxsim + 1].first + 0.5 * ziso;
          int nreciso = 0;
          std::vector<double> zrec;
          for (reco::VertexCollection::const_iterator vrec = recVtxs->begin(); vrec != recVtxs->end(); ++vrec) {
            if (!(vrec->isFake()) && (vrec->ndof() > 4) && (vrec->z() > zmin) && (vrec->z() < zmax)) {
              nreciso++;
              zrec.push_back(vrec->z());
            }
          }
          Fill(h, "zdiffsimisoall", dzsim);
          if (nreciso == 0)
            Fill(h, "zdiffsimiso0", dzsim);
          if (nreciso == 1)
            Fill(h, "zdiffsimiso1", dzsim);
          if (nreciso == 2)
            Fill(h, "zdiffsimiso2", dzsim);
          if (nreciso == 3)
            Fill(h, "zdiffsimiso3", dzsim);
          if (nreciso == 2) {
            double dzrec = fabs(zrec[1] - zrec[2]);
            Fill(h, "zdiffreciso2", dzrec);
            Fill(h, "dzrecvssimiso2", fabs(dzsim), dzrec);
          }
        }
      }
    }

    // single isolated
    for (unsigned int idxsim = 0; idxsim < simpv.size(); idxsim++) {
      if ((simpv[zsimv[idxsim].second].nGenTrk > 4) &&
          ((idxsim == 0) || ((idxsim > 0) && (zsimv[idxsim].first - zsimv[idxsim - 1].first > ziso))) &&
          ((idxsim >= zsimv.size() - 1) ||
           ((idxsim < zsimv.size() - 1) && (zsimv[idxsim + 1].first - zsimv[idxsim].first > ziso)))) {
        double zmin = zsimv[idxsim].first - 0.5 * ziso;
        double zmax = zsimv[idxsim].first + 0.5 * ziso;
        int nreciso = 0;
        for (reco::VertexCollection::const_iterator vrec = recVtxs->begin(); vrec != recVtxs->end(); ++vrec) {
          if (!(vrec->isFake()) && (vrec->ndof() > 4) && (vrec->z() > zmin) && (vrec->z() < zmax)) {
            nreciso++;
            Fill(h, "zreciso", vrec->z() - zsimv[idxsim].first);
          }
          Fill(h, "nreciso", float(nreciso));
        }
      }
    }

  }  //simpv>1
}
/***************************************************************************************/

/***************************************************************************************/
void PrimaryVertexAnalyzer4PU::analyzeVertexRecoCPUTime(std::map<std::string, TH1*>& h,
                                                        const reco::VertexCollection* recVtxs,
                                                        const std::string message)
/***************************************************************************************/
{
  double tclu = 0;
  double tfit = 0;
  int nsel = 0;
  bool found_timing_info = false;


  for (reco::VertexCollection::const_iterator v = recVtxs->begin(); v != recVtxs->end(); ++v) {
    if (v->isFake()) {
      tclu = v->covariance(iX, iX);
      tfit = v->covariance(iY, iY);
      found_timing_info = true;
      //std::cout << "extracted timing info " << message << "  tclu = "  << tclu << " ms   tfit =" << tfit << " ms" << endl;
    } else {
      if (select(*v)) {
        nsel++;
      }
    }
  }

  if (found_timing_info) {
    Fill(h, "cputime/tcluvsSimPU", simPU_, tclu);
    Fill(h, "cputime/tfitvsSimPU", simPU_, tfit);
    Fill(h, "cputime/tcluvsLPU", lumiPU_, tclu);
    Fill(h, "cputime/tfitvsLPU", lumiPU_, tfit);
    Fill(h, "cputime/tcluvsnsel", nsel, tclu);
    Fill(h, "cputime/tfitvsnsel", nsel, tfit);
  }
}
/***************************************************************************************/

//******* the following code does not require MC and will/should work for data, requires tracks **********
/***************************************************************************************/
void PrimaryVertexAnalyzer4PU::analyzeVertexCollectionReco(std::map<std::string, TH1*>& h,
                                                           const reco::VertexCollection* recVtxs,
                                                           Tracks& tracks,
                                                           const std::string message)
/***************************************************************************************/
{
  //  cout << "PrimaryVertexAnalyzer4PU::analyzeVertexCollectionReco " << message << endl;
  int nrectrks = tracks.size();

  // repeat some counting here from analyzeVertexCollectionNotrack, needed for eventclassification and correlation
  int nrecvtx = 0;
  int nselvtx = 0;
  int nvtxselgt1sigmaz = 0;
  for (unsigned int iv = 0; iv < recVtxs->size(); iv++) {
    const reco::Vertex* v = &(recVtxs->at(iv));
    if (!(v->isFake()) && (v->ndof() > 0)) {
      nrecvtx++;
    }

    fillVertexHistos(h, "rec", &(*v), tracks, iv);
    if(iv == 0){fillVertexHistos(h, "tagged", &(*v), tracks, iv);}
    
    if (select(*v)) {
      nselvtx++;
      fillVertexHistos(h, "selected", &(*v), tracks, iv);
      if (std::abs(v->position().z() - vertexBeamSpot_.z0() - dzb_) > sigmaZ_) {
        nvtxselgt1sigmaz++;
      }
    }
  }


  // distance to first vertex and ndof for different positions in the vertex list
  if( (recVtxs->size() > 1) && (select(recVtxs->at(0))) ){
    double z0 = recVtxs->at(0).position().z();
    Fill(h, "ndof-vtx0", recVtxs->at(0).ndof());
    Fill(h, "logsumpt2-vtx0", log(vertex_sumpt2(recVtxs->at(0))));
    Fill(h, "logndof-vtx0", log(recVtxs->at(0).ndof())/log(10.));

    for (unsigned int iv = 1; iv < recVtxs->size(); iv++) {
      const reco::Vertex* v = &(recVtxs->at(iv));
      if(select(*v)){
	if(iv < 10){
	    Fill(h, Form("zdiffrec-vtx%d", iv), v->position().z() - z0);
	    Fill(h, Form("ndof-vtx%d", iv), v->ndof());
	    Fill(h, Form("logndof-vtx%d", iv), log(v->ndof())/log(10.));
	    Fill(h, Form("logsumpt2-vtx%d", iv), log(vertex_sumpt2(*v)));
	}else{
	  Fill(h, "zdiffrec-vtxgt9", v->position().z() - z0);
	  Fill(h, "ndof-vtxgt9", v->ndof());
	  Fill(h, "logndof-vtxgt9", log(v->ndof())/log(10.));
	  Fill(h, "logsumpt2-vtxgt9", log(vertex_sumpt2(*v)));
	}
      }
    }
  }



  bool is_hiPU = lumiPU_ > 50;
  //bool is_tail =  (lumiPU_ > 40) && (nselvtx > (0.8*lumiPU_ + 3*sqrt(0.8*lumiPU_)));
  bool is_tail = (lumiPU_ > 40) && (nvtxselgt1sigmaz > (0.24 * lumiPU_ + 3 * sqrt(0.24 * lumiPU_)));

  // -----------------  reconstructed tracks  ------------------------
  // the list of selected tracks can only be correct if the selection parameterset  and track collection
  // is the same that was used for the reconstruction

  for (unsigned int i = 0; i < tracks.size(); i++) {
    RecoTrack tk = tracks(i);
    RefToBase<reco::Track> t = tracks.ref(i);

    if ((recVtxs->size() > 0) && (recVtxs->begin()->isValid())) {
      fillTrackHistos(h, "trkall", tk, &(*recVtxs->begin()));
      //fillTrackClusterHistos(h, "trkall", *t, &(*recVtxs->begin()));
    } else {
      fillTrackHistos(h, "trkall", tk);
      //fillTrackClusterHistos(h, "trkall", *t);
    }

    if (MC_)
    // histogram the fraction of truth-matched tracks vs z, requires TP and should not be here!!!
    {
      if (tracks(i).matched) {
        Fill(h, "matchedallfractionvsz", t->vz(), 1.);
        Fill(h, "unmatchedallfractionvsz", t->vz(), 0.);
      } else {
        Fill(h, "unmatchedallfractionvsz", t->vz(), 1.);
        Fill(h, "matchedallfractionvsz", t->vz(), 0.);
      }
    }

    if (MC_ && f4D_ && tk.matched && tk.has_timing)
      {
	if ((tk.dt >0) && (tk.dt < 0.1) && (fabs(tk.t-tk.tsim) > 5.*tk.dt))
	  {
	    fillTrackHistos(h, "MTDtail", tk);
	  }
      }
  }

  // ---- selected tracks

  int nseltrks = 0;
  int nseltrksptlt04 = 0;
  int nseltrksptlt03 = 0;
  int nseltrksptlt02 = 0;

  for (unsigned int i = 0; i < tracks.size(); i++) {
    
    RefToBase<reco::Track> t = tracks.ref(i);
    TransientTrack* tt = tracks(i).tt;
    RecoTrack tk = tracks(i);

    if (theTrackFilter(*tt)) {
      fillTrackHistos(h, "trksel", tk, &(*recVtxs->begin()));
      nseltrks++;
      if (t->pt() < 0.4)
        nseltrksptlt04++;
      if (t->pt() < 0.3)
        nseltrksptlt03++;
      if (t->pt() < 0.2)
        nseltrksptlt02++;

      auto nmiss = t->hitPattern().numberOfLostHits(HitPattern::MISSING_INNER_HITS); 
      if(nmiss == 0){
	fillTrackHistos(h, "missing_inner0", tk, &(*recVtxs->begin()));
      }else if(nmiss == 1){
	fillTrackHistos(h, "missing_inner1", tk, &(*recVtxs->begin()));
      }else if (nmiss == 2){
	fillTrackHistos(h, "missing_inner2", tk, &(*recVtxs->begin()));
      }

      auto nbarrel = t->hitPattern().pixelBarrelLayersWithMeasurement();
      if(nbarrel <2){
	fillTrackHistos(h, "nbarrel_lt2", tk, &(*recVtxs->begin()));
      }else if(nbarrel ==2){
	fillTrackHistos(h, "nbarrel_eq2", tk, &(*recVtxs->begin()));
      }      

      if((abs(tk.z) < 4) && (tk.dz < 0.01) && (abs(tk.eta) > 1.4) && (abs(tk.eta) < 1.8)){
	fillTrackHistos(h, "highetadriver", tk, &(*recVtxs->begin()));
      }
      if((abs(tk.z) < 4) && (tk.dz < 0.01)){
	fillTrackHistos(h, "alletadriver", tk, &(*recVtxs->begin()));
      }

      if (is_hiPU) {
        fillTrackHistos(h, "thipu", tk, &(*recVtxs->begin()));
      }  // all tracks here, not only those that are in a vertex (->"hipu" without a "t")

      if (is_tail) {
        fillTrackHistos(h, "ttail", tk, &(*recVtxs->begin()));
        double z = (tt->stateAtBeamLine().trackStateAtPCA()).position().z();
        if (fabs(z - vertexBeamSpot_.z0() - dzb_) > sigmaZ_)
          fillTrackHistos(h, "ttailzgt1", tk, &(*recVtxs->begin()));
      }

      if (MC_) {
        if (tk.matched) {
          fillTrackHistos(h, "seltpmatched", tk);
	  if (abs(tk.pt > 10.) && (abs(tk.eta)<0.5))
	    fillTrackHistos(h, "hiptcentral", tk);
          if (abs(tk.t) > 1.0)
            fillTrackHistos(h, "seltpmatched_tgt1", tk);
	  if (tk.simEvt == 0){
            fillTrackHistos(h, "seltpmatchedSignal", tk);
	  }else{
            fillTrackHistos(h, "seltpmatchedPU", tk);
	  }
          Fill(h, "matchedselfractionvsz", tk.z, 1.);
          Fill(h, "unmatchedselfractionvsz", tk.z, 0.);
        } else {
          Fill(h, "matchedselfractionvsz", tk.z, 0.);
          //Fill(h, "unmatchedselfractionvsz", t->vz(), 1.);
          Fill(h, "unmatchedselfractionvsz", tk.z, 1.);
          fillTrackHistos(h, "seltpunmatched", tk);
        }
      }

      int foundinvtx = 0;
      int nvtemp = -1;
      for (reco::VertexCollection::const_iterator v = recVtxs->begin(); v != recVtxs->end(); ++v) {
        nvtemp++;
        if ((v->isFake()) || (v->ndof() < -2))
          break;
        for (trackit_t tv = v->tracks_begin(); tv != v->tracks_end(); tv++) {
          if (((**tv).vz() == t->vz() && ((**tv).phi() == t->phi()))) {
            foundinvtx++;
          }
        }
      }

      if (foundinvtx == 0) {
        fillTrackHistos(h, "sellost", tk);
      } else if (foundinvtx > 1) {
        //	cout << "track found in " << foundinvtx << " vertices" endl;
      }
    }
  }

  if (nseltrks < minNumberOfSelTrks_) {
    cout << "Event without tracks skipped   run= " << run_ << " ls = " << luminosityBlock_ << " BX=" << bunchCrossing_
         << "  instBXLumi_=" << 1.e3 * instBXLumi_ << " expected PU = " << lumiPU_ << "    selected tracks " << nseltrks
         << " rec tracks = " << nrectrks << std::endl;
    return;
  }

  // event
  Fill(h, "nseltrk", nseltrks);
  Fill(h, "nrectrk", nrectrks);
  Fill(h, "nseltrkvsLPUprof", lumiPU_, float(nseltrks));
  Fill(h, "nrectrkvsLPUprof", lumiPU_, float(nrectrks));
  Fill(h, "nseltrkptlt04vsLPUprof", lumiPU_, float(nseltrksptlt04));
  Fill(h, "nseltrkptlt03vsLPUprof", lumiPU_, float(nseltrksptlt03));
  Fill(h, "nseltrkptlt02vsLPUprof", lumiPU_, float(nseltrksptlt02));

  Fill(h, "nrecvtx4vsnrectrk", float(nrectrks), float(nselvtx));
  Fill(h, "nrecvtx4vsnrectrkprof", float(nrectrks), float(nselvtx));
  Fill(h, "nrecvtx4vsnseltrk", float(nseltrks), float(nselvtx));
  Fill(h, "nrecvtx4vsnseltrkprof", float(nseltrks), float(nselvtx));

  if (is_tail) {
    Fill(h, "nrecvtx4vsnseltrk_tail", float(nseltrks), float(nselvtx));
    Fill(h, "nrecvtx4vsnseltrkprof_tail", float(nseltrks), float(nselvtx));
  }
  if (is_hiPU) {
    Fill(h, "nrecvtx4vsnseltrk_hipu", float(nseltrks), float(nselvtx));
    Fill(h, "nrecvtx4vsnseltrkprof_hipu", float(nseltrks), float(nselvtx));
  }

  if (nrecvtx > 0) {
    Fill(h, "eff0vsntrec", nrectrks, 1.);
    Fill(h, "eff0vsntsel", nseltrks, 1.);
  } else {
    Fill(h, "eff0vsntrec", nrectrks, 0.);
    Fill(h, "eff0vsntsel", nseltrks, 0.);
  }

  // properties of events without a vertex
  if ((nrecvtx == 0) || (recVtxs->begin()->isFake())) {
    Fill(h, "nrectrk0vtx", nrectrks);
    Fill(h, "nseltrk0vtx", nseltrks);
  }

  for (reco::VertexCollection::const_iterator v = recVtxs->begin(); v != recVtxs->end(); ++v) {
    // count vertex tracks
    double npt1 = 0, ntrkwgt05 = 0;
    for (trackit_t t = v->tracks_begin(); t != v->tracks_end(); t++) {
      if (v->trackWeight(*t) > 0.5) {
        ntrkwgt05++;
        Fill(h, "trackwtgt05vsz", v->position().z(), v->trackWeight(*t));
      }
      if ((**t).pt() > 1.0)
        npt1++;
      Fill(h, "trackwtvsz", v->position().z(), v->trackWeight(*t));
    }
    if (v->ndof() > 4) {
      Fill(h, "ntrkpt1vsz", v->position().z(), npt1);
      Fill(h, "ntrkwgt05vsz", v->position().z(), ntrkwgt05);
      Fill(h, "ftrkwgt05vsz", v->position().z(), ntrkwgt05 / v->tracksSize());
    }

    Fill(h, "nbtksinvtx", v->tracksSize());
    if (instBXLumi_ > 0) {
      Fill(h, "nbtksinvtxvsL", 1e3 * instBXLumi_, float(v->tracksSize()));
      Fill(h, "nbtksinvtxvsLPU", lumiPU_, float(v->tracksSize()));
      if (v->ndof() > 4)
        Fill(h, "nbtksinvtx4vsL", 1e3 * instBXLumi_, float(v->tracksSize()));
      if (v->ndof() > 4)
        Fill(h, "nbtksinvtx4vsLPU", lumiPU_, float(v->tracksSize()));
    }
    Fill(h, "nbtksinvtx2", v->tracksSize());

    double sumw = 0.5 * (v->ndof() + 2);
    if (v->ndof() > 4) {
      Fill(h, "sumwvsz", v->position().z(), sumw);
      Fill(h, "sumntkvsz", v->position().z(), (double)v->tracksSize());
      Fill(h, "sumwoverntk", sumw / v->tracksSize());
      Fill(h, "vtxndfoverntk", v->ndof() / v->tracksSize());
      Fill(h, "vtxndf2overntk", (v->ndof() + 2) / v->tracksSize());
      Fill(h, "sumwoverntkvsz", v->position().z(), sumw / v->tracksSize());
      Fill(h, "sumwoverntkvsz4", v->position().z(), sumw / v->tracksSize());
      if (ntrkwgt05 > 0) {
        Fill(h, "sumwoverntkwgt05vsz", v->position().z(), sumw / ntrkwgt05);
      }
    }
    if ((v->ndof() > 4) && (v->ndof() < 10)) {
      Fill(h, "sumwoverntkvszlo", v->position().z(), sumw / v->tracksSize());
    }
    if ((v->ndof() > 4) && (v->ndof() < 10)) {
      Fill(h, "sumwoverntkvszlo", v->position().z(), sumw / v->tracksSize());
    }
    if (v->ndof() > 20) {
      Fill(h, "sumwoverntkvszhi", v->position().z(), sumw / v->tracksSize());
    }
    Fill(h, "ntrkvsz", v->position().z(), float(v->tracksSize()));
  }

}  // end of analyzeVertexCollectionReco
/***************************************************************************************/


 
/***************************************************************************************/
void PrimaryVertexAnalyzer4PU::analyzeVertexCollectionRecoNoTracks(std::map<std::string, TH1*>& h,
                                                                   const reco::VertexCollection* recVtxs,
                                                                   const std::string message)
//******* for data/MC with vertices but not tracks, e.g. MINIAOD
/***************************************************************************************/
{
  // make a z-sorted list for gap histogramming and neighbour searches
  vector<pair<double, unsigned int>> zrecv;
  for (unsigned int idx = 0; idx < recVtxs->size(); idx++) {
    if (select(recVtxs->at(idx))) {
      zrecv.push_back(make_pair(recVtxs->at(idx).z(), idx));
    }
  }
  stable_sort(zrecv.begin(), zrecv.end(), lt);

  // count vertices outside +/1 sigma z
  int nvtxselgt1sigmaz = 0;
  for (unsigned int idx = 0; idx < zrecv.size(); idx++) {
    if (std::abs(zrecv[idx].first - vertexBeamSpot_.z0() - dzb_) > sigmaZ_) {
      nvtxselgt1sigmaz++;
    }
  }

  //bool is_tail = (lumiPU_ > 40) && (zrecv.size() > (0.8*lumiPU_ + 3*sqrt(0.8*lumiPU_)));
  bool is_tail = (lumiPU_ > 40) && (nvtxselgt1sigmaz > (0.24 * lumiPU_ + 3 * sqrt(0.24 * lumiPU_)));
  bool is_hiPU = (lumiPU_ > 50);

  Fill(h, "bunchCrossing", float(bunchCrossing_));
  Fill(h, "bunchCrossing_PU", float(bunchCrossing_), 1. / lumiPU_);

  for (unsigned int idx = 0; idx < zrecv.size(); idx++) {
  /*
  //unused if we don't fill NoTracks here , FIXME do it, but only when no tracks are available?
    float dznearest = 0; 
    if (zrecv.size() > 1) {
      if (idx == 0) {
        dznearest = zrecv[0].first - zrecv[1].first;
      } else if (idx == (zrecv.size() - 1)) {
        dznearest = zrecv[idx].first - zrecv[idx - 1].first;
      } else {
        if (fabs(zrecv[idx].first - zrecv[idx - 1].first) < fabs(zrecv[idx].first - zrecv[idx + 1].first)) {
          dznearest = zrecv[idx].first - zrecv[idx - 1].first;
        } else {
          dznearest = zrecv[idx].first - zrecv[idx + 1].first;
        }
      }
    }

    const reco::Vertex* v = &(recVtxs->at(zrecv[idx].second));
    fillVertexHistosNoTracks(h, "sel", v, dznearest);  FIXME avoid calling it twice 
    */

    // this is the test for the "empty" vertices, remove it
  for (unsigned int idx = 0; idx < zrecv.size(); idx++) {
    if ( (recVtxs->at(idx).ndof() > 4.09) && (recVtxs->at(idx).ndof() < 4.11) 
	 && (recVtxs->at(idx).tracksSize() == 0)
	 && (recVtxs->at(idx).chi2() > 0.9)
	 && (recVtxs->at(idx).chi2() < 1.1)){
      Fill(h, "indexempty", idx);
    }
  }

    /*
    if (is_tail)
      fillVertexHistosNoTracks(h, "tail", v, dznearest);
    if (is_hiPU)
      fillVertexHistosNoTracks(h, "hipu", v, dznearest);
    double z = zrecv[idx].first;
    if (std::abs(z - vertexBeamSpot_.z0() - dzb_) > sigmaZ_) {
      fillVertexHistosNoTracks(h, "selzgt1", v, dznearest);
      if (is_tail)
        fillVertexHistosNoTracks(h, "tailzgt1", v, dznearest);
    }
    */
  }

  // count and z-map gaps of various sizes
  int ng1 = 0, ng05 = 0, ng03 = 0, ng02 = 0;

  if (zrecv.size() > 0) {
    Fill(h, "zofselvtx", zrecv[0].first);
    for (unsigned int idx = 1; idx < zrecv.size(); idx++) {
      double dz = fabs(zrecv[idx].first - zrecv[idx - 1].first);
      double zg = 0.5 * (zrecv[idx].first + zrecv[idx - 1].first);
      if (fabs(dz) < 0.02)
        Fill(h, "zofzgap02mm", zg);
      if (fabs(dz) < 0.03)
        Fill(h, "zofzgap03mm", zg);
      if (fabs(dz) < 0.05)
        Fill(h, "zofzgap05mm", zg);
      if (fabs(dz) < 0.1)
        Fill(h, "zofzgap1mm", zg);
      Fill(h, "zofselvtx", zrecv[idx].first);
      if (is_tail) {
        if (fabs(dz) < 0.02)
          Fill(h, "zofzgap02mm_tail", zg);
        if (fabs(dz) < 0.03)
          Fill(h, "zofzgap03mm_tail", zg);
        if (fabs(dz) < 0.05)
          Fill(h, "zofzgap05mm_tail", zg);
        if (fabs(dz) < 0.1)
          Fill(h, "zofzgap1mm_tail", zg);
        Fill(h, "zofselvtx_tail", zrecv[idx].first);
      }
    }

    ng1 = 1;
    ng05 = 1;
    ng03 = 1;
    ng02 = 1;
    double z0_1 = zrecv[0].first;
    double z0_05 = zrecv[0].first;
    double z0_03 = zrecv[0].first;
    double z0_02 = zrecv[0].first;
    for (unsigned int idx = 1; idx < zrecv.size(); idx++) {
      double z1 = zrecv[idx].first;
      if ((z1 - z0_1) > 0.1) {
        z0_1 = z1;
        ng1++;
      }

      if ((z1 - z0_05) > 0.05) {
        z0_05 = z1;
        ng05++;
      }

      if ((z1 - z0_03) > 0.03) {
        z0_03 = z1;
        ng03++;
      }

      if ((z1 - z0_02) > 0.02) {
        z0_02 = z1;
        ng02++;
      }
    }
  }
  // done with gap counting

  // estimate the event sigmaz, use the beamspot position
  double sigmaz_event = 0;
  if (zrecv.size() > 0) {
    double sumz2 = 0;
    for (unsigned int idx = 0; idx < zrecv.size(); idx++) {
      sumz2 += pow(zrecv[idx].first - vertexBeamSpot_.z0(), 2);
    }
    sigmaz_event = sqrt(sumz2 / zrecv.size());  // not robustified yet
    Fill(h, "sigmaz_event", sigmaz_event);
    Fill(h, "sigmaz_event_over_sigmaZ", sigmaz_event / sigmaZ_);
    Fill(h, "sigmaz_event_vs_beam", sigmaZ_, sigmaz_event);
    Fill(h, "sigmaz_event_vs_bx", bunchCrossing_, sigmaz_event / sigmaZ_);
    if (is_tail) {
      Fill(h, "sigmaz_event_vs_beam_tail", sigmaZ_, sigmaz_event);
      Fill(h, "sigmaz_event_vs_bx_tail", bunchCrossing_, sigmaz_event / sigmaZ_);
    }
  }

  // z coordinate of the tagged vertex
  double ztag = 1000.;
  if (!(recVtxs->begin()->isFake())) {
    ztag = recVtxs->begin()->position().z();
    Fill(h, "vtxz_tag", ztag);
  }

  // count vertices above some ndof thresholds
  int nrec = 0, nrec0 = 0, nrec2 = 0, nrec3 = 0, nrec4 = 0, nrec5 = 0, nrec6 = 0, nrec7 = 0, nrec8 = 0;
  for (reco::VertexCollection::const_iterator v = recVtxs->begin(); v != recVtxs->end(); ++v) {
    if (!(v->isFake()) && v->ndof() > 0) {
      nrec++;
      if (v->ndof() > 0)
        nrec0++;
      if (v->ndof() > 2)
        nrec2++;
      if (v->ndof() > 3)
        nrec3++;
      if (v->ndof() > 4)
        nrec4++;
      if (v->ndof() > 5)
        nrec5++;
      if (v->ndof() > 6)
        nrec6++;
      if (v->ndof() > 7)
        nrec7++;
      if (v->ndof() > 8)
        nrec8++;
    }
  }
  Fill(h, "nrecvtx", nrec);
  Fill(h, "nrecvtx2", nrec2);
  Fill(h, "nrecvtx3", nrec3);
  Fill(h, "nrecvtx4", nrec4);
  Fill(h, "nrecvtx5", nrec5);
  Fill(h, "nrecvtx6", nrec6);
  Fill(h, "nrecvtx7", nrec7);
  Fill(h, "nrecvtx8", nrec8);

  if (instBXLumi_ > 0) {
    Fill(h, "nrecvtxvsLPU", lumiPU_, float(nrec));
    Fill(h, "nrecvtxvsLPUprof", lumiPU_, float(nrec));
    Fill(h, "nrecvtx4vsLPU", lumiPU_, float(nrec4));
    Fill(h, "nrecvtx4vsLPUprof", lumiPU_, float(nrec4));
    Fill(h, "nrecvtx4vsavgLPU", avglumiPU_, float(nrec4));
    Fill(h, "nrecvtx4vsavgLPUprof", avglumiPU_, float(nrec4));
    Fill(h, "LPUnrecvtx4vsLPUprof", lumiPU_, lumiPU_);  // keep track of the proper bin centers
  }

  Fill(h, "nevtLPU", lumiPU_);
  Fill(h, "sigmaZ", sigmaZ_);

  //  properties of (valid) vertices
  double ndof2 = -10, ndof1 = -10, zndof1 = 0, zndof2 = 0;

  int nv = 0;
  for (reco::VertexCollection::const_iterator v = recVtxs->begin(); v != recVtxs->end(); ++v) {
    nv++;

    if (v->isFake()) {
      Fill(h, "isFake", 1.);
    } else {
      Fill(h, "isFake", 0.);
    }
    if (v->isFake() || ((v->ndof() < -1) && (v->ndof() > -3))) {
      Fill(h, "isFake1", 1.);
    } else {
      Fill(h, "isFake1", 0.);
    }

    if ((v->isFake()) || (v->ndof() < -1))
      continue;

    //fillVertexHistosNoTracks(h, "rec", &(*v), nv-1); // FIXME

    if (v->ndof() > 4) {
      Fill(h, "ndofvsz4", v->position().z(), v->ndof());

      Fill(h, "vtxndofvsLS", lsglb_, v->ndof());
      if (fabs(v->position().z() - vertexBeamSpot_.z0() - dzb_) > sigmaZ_) {
        Fill(h, "vtxndofzgt1vsLS", lsglb_, v->ndof());
      }
      if (v == recVtxs->begin()) {
        //fillVertexHistosNoTracks(h, "tagged", &(*v), nv-1); FIXME avoid calling twice
        Fill(h, "vtxndf_tagged", v->ndof());
      }
    }

    // low multiplicity tagged vertex in high pile-up
    if ((v->ndof() > 4) && (v == recVtxs->begin()) && (recVtxs->size() > 20)) {
      Fill(h, "zrec_tagged_hpu", v->position().z());
      if (v->ndof() < 20) {
        Fill(h, "zrec_tagged_lowndof_hpu", v->position().z());
        // quiet for now reportEvent(Form("low multiplicity tagged vertex in high pile-up  ndof=%4.1f",v->ndof()));
      }
    }

    if (v->ndof() > ndof1) {
      ndof2 = ndof1;
      zndof2 = zndof1;
      ndof1 = v->ndof();
      zndof1 = v->position().z();
    } else if (v->ndof() > ndof2) {
      ndof2 = v->ndof();
      zndof2 = v->position().z();
    }

    if (v->ndof() > 4) {
      Fill(h, "trkchi2overndof", v->chi2() / v->ndof());
    }

    Fill(h, "vtxprob", ChiSquaredProbability(v->chi2(), v->ndof()));
    Fill(h, "ndofvsz0", v->position().z(), v->ndof());

    if (v->ndof() < 4.) {
      Fill(h,
           "zrecBeamPullndoflt4",
           (v->position().z() - vertexBeamSpot_.z0()) / sqrt(pow(v->zError(), 2) + pow(sigmaZ_, 2)));
    }
    if (v->ndof() > 4.0) {  // enter only vertices that really contain tracks
      Fill(h, "xrec", v->position().x());
      Fill(h, "yrec", v->position().y());
      Fill(h, "zrec", v->position().z());
      Fill(h, "xrec1", v->position().x());
      Fill(h, "yrec1", v->position().y());
      Fill(h, "zrec1", v->position().z());
      Fill(h, "xrec2", v->position().x());
      Fill(h, "yrec2", v->position().y());
      Fill(h, "zrec2", v->position().z());
      Fill(h, "xrec3", v->position().x());
      Fill(h, "yrec3", v->position().y());
      Fill(h, "zrec3", v->position().z());
      Fill(h, "zrec3a", v->position().z());
      Fill(h, "xrecb", v->position().x() - vertexBeamSpot_.x0());
      Fill(h, "yrecb", v->position().y() - vertexBeamSpot_.y0());
      Fill(h, "zrecb", v->position().z() - vertexBeamSpot_.z0());

      // dx,dy taking the beam slope into account
      double dx = v->position().x() - vertexBeamSpot_.x(v->position().z()) - dxb_;
      double dy = v->position().y() - vertexBeamSpot_.y(v->position().z()) - dyb_;

      Fill(h, "xrecBeam", dx);
      Fill(h, "yrecBeam", dy);
      Fill(h, "zrecBeam", v->position().z() - vertexBeamSpot_.z0() - dzb_);

      Fill(h, "xrecBeamPull", dx / sqrt(pow(v->xError(), 2) + pow(vertexBeamSpot_.BeamWidthX(), 2)));
      Fill(h, "yrecBeamPull", dy / sqrt(pow(v->yError(), 2) + pow(vertexBeamSpot_.BeamWidthY(), 2)));
      Fill(h, "zrecBeamPull", (v->position().z() - vertexBeamSpot_.z0()) / sqrt(pow(v->zError(), 2) + pow(sigmaZ_, 2)));
      Fill(h, "zrecBeamPull0", (v->position().z() - vertexBeamSpot_.z0()) / sigmaZ_);
      if (v->ndof() < 10.) {
        Fill(h, "zrecBeamPullndof4-10", (v->position().z() - vertexBeamSpot_.z0()) / sigmaZ_);
      }
      if (v->ndof() > 12.) {
        Fill(h, "zrecBeamPull12", (v->position().z() - vertexBeamSpot_.z0()) / sigmaZ_);
      }

      if (DO_BEAMSPOT_ANALYSIS) {
        Fill(h, "dxbinavg", v->xError(), v->xError());
        Fill(h, "dybinavg", v->yError(), v->yError());
        Fill(h, "dx2binavg", pow(v->xError(), 2), pow(v->xError(), 2));
        Fill(h, "dy2binavg", pow(v->yError(), 2), pow(v->yError(), 2));
        Fill(h, "xrecBeamvsdx", v->xError(), dx);
        Fill(h, "yrecBeamvsdy", v->yError(), dy);
        Fill(h, "xrecBeamvsdxprof", v->xError(), dx);
        Fill(h, "yrecBeamvsdyprof", v->yError(), dy);
        Fill(h, "xrecBeam2vsdx2prof", pow(v->xError(), 2), pow(dx, 2));
        Fill(h, "yrecBeam2vsdy2prof", pow(v->yError(), 2), pow(dy, 2));

        Fill(h, "xrecBeamvszprof", v->position().z(), dx);
        Fill(h, "yrecBeamvszprof", v->position().z(), dy);

        Fill(h, "xrecBeamvsNdofprof", v->ndof(), dx);
        Fill(h, "yrecBeamvsNdofprof", v->ndof(), dy);
      }

      if (bunchCrossing_ > 0) {
        Fill(h, "zvsls", float(luminosityBlock_), v->position().z());
        Fill(h, "zbeamvsls", float(luminosityBlock_), vertexBeamSpot_.z0());
      }

    }  // ndof>4

    double z0 = v->position().z() - vertexBeamSpot_.z0() - dzb_;

    // reference shape for zdiff4 and other z-histograms
    if ( select(*v) ) {
      Fill(h, "zrecsel", z0);
      Fill(h, "zrecselvsLPU", lumiPU_, z0);
      Fill(h, "zrecselvsordprof", float(nv), z0);
    }


    //  properties of (valid) neighbour vertices
    reco::VertexCollection::const_iterator v1 = v;
    v1++;
    for (; v1 != recVtxs->end(); v1++) {
      if ((v1->isFake()) || (v1->ndof() < -1))
        continue;
      double z1 = v1->position().z() - vertexBeamSpot_.z0() - dzb_;

      npair_++;
      Fill(h, "zdiffrec", z1 - z0);
      Fill(h, "zdiffrechr", z1 - z0);
      if (v == recVtxs->begin()) {
        Fill(h, "zdiffrec_tagged", z1 - z0);
        Fill(h, "zdiffrechr_tagged", z1 - z0);
      }

      // lower ndof of the pair
      double ndoflow = v1->ndof();
      double ndofhi = v->ndof();
      if (v1->ndof() > v->ndof()) {
        ndofhi = v1->ndof();
        ndoflow = v->ndof();
      } else {
        ndofhi = v->ndof();
        ndoflow = v1->ndof();
      }

      // central vertices, avoid acceptance issues
      if ((ndoflow > 4) && (fabs(z0) < 1.)) {
        Fill(h, "dzreccentral", fabs(z0 - z1));
        Fill(h, "ndofcentral", fabs(z0 - z1), v->ndof());
        Fill(h, "ndoflocentral", fabs(z0 - z1), ndoflow);
        Fill(h, "ndofhicentral", fabs(z0 - z1), ndofhi);
        Fill(h, "ndofsumcentral", fabs(z0 - z1), ndofhi + ndoflow);
        if (v == recVtxs->begin()) {
          Fill(h, "dzreccentral_tagged", fabs(z0 - z1));
          Fill(h, "ndofcentral_tagged", fabs(z0 - z1), v->ndof());
          Fill(h, "ndoflocentral_tagged", fabs(z0 - z1), ndoflow);
          Fill(h, "ndofhicentral_tagged", fabs(z0 - z1), ndofhi);
          Fill(h, "ndofsumcentral_tagged", fabs(z0 - z1), ndofhi + ndoflow);
          Fill(h, "ndofothercentral_tagged", fabs(z0 - z1), v1->ndof());
        }
      }

      Fill(h, "zndof2", zndof2);

      double zbar = 0.5 * (z1 + z0);
      double zbarp = zbar / sigmaZ_;
      if (DO_ZDIFFVSZ_ANALYSIS) {
        Fill(h, "zdiffvsz", z1 - z0, zbar);
        Fill(h, "zdiffvszp", z1 - z0, zbarp);

        if (fabs(z1 - z0) < 0.2) {
          Fill(h, "zbarFakeEnriched", zbar);
          if (ndoflow > 5)
            Fill(h, "zbarFakeEnriched5", zbar);
          if (ndoflow > 2)
            Fill(h, "zbarFakeEnriched2", zbar);
        }
        if ((fabs(z1 - z0) > 2.0) && (v->ndof() > 10) && (v1->ndof() > 10)) {
          Fill(h, "zbarFakeDepleted", zbar);
        }  // just for comparison , pure real
      }

      if (select(*v) && select(*v1)) {
        Fill(h, "zdiffrecsel", z1 - z0);
        Fill(h, "zdiffrecselhr", z1 - z0);

        if (v == recVtxs->begin()) {
          Fill(h, "zdiffrecsel_tagged", z1 - z0);
          Fill(h, "zdiffrecselhr_tagged", z1 - z0);
          Fill(h, "ndofother_tagged", fabs(z0 - z1), v1->ndof());
          if (v->ndof() < 20) {
            Fill(h, "zdiffrecselhr_tagged_lowndof", z1 - z0);
            Fill(h, "ndofother_tagged_lowndof", fabs(z0 - z1), v1->ndof());
          }
        }
        if (abs(0.5 * (z0 + z1) - vertexBeamSpot_.z0() - dzb_) > sigmaZ_) {
          Fill(h, "zdiffrecselhr_gt1sigma", z1 - z0);
        } else {
          Fill(h, "zdiffrecselhr_lt1sigma", z1 - z0);
        }
      }

      if (select(*v) && (!select(*v1))) {
        Fill(h, "zdiffrecseldeselhr", z1 - z0);
      } else if (!(select(*v)) && select(*v1)) {
        Fill(h, "zdiffrecseldeselhr", z0 - z1);
      } else if (!(select(*v)) && (!select(*v1))){
        Fill(h, "zdiffrecdeselhr", z0 - z1);
      }

      if (select(*v) && select(*v1)) {

        if (DO_ZDIFFVSZ_ANALYSIS) {
          Fill(h, "zdiffvszsel", z1 - z0, zbar);
          Fill(h, "zdiffvszselhr", z1 - z0, zbar);
          if (is_hiPU)
            Fill(h, "zdiffvszselhr_hipu", z1 - z0, zbar);
          Fill(h, "zdiffvszpsel", z1 - z0, zbarp);
          if (nrec4 == 2)
            Fill(h, "zdiffvszselNv2", z1 - z0, zbar);
        }

        if (fabs(zbar - vertexBeamSpot_.z0() - dzb_) > sigmaZ_) {
          Fill(h, "zdiffrecselhr_zgt1", z1 - z0);
        }
        if (is_tail) {
          Fill(h, "zdiffrecselhr_tail", z1 - z0);
          if (fabs(zbar - vertexBeamSpot_.z0() - dzb_) > sigmaZ_) {
            Fill(h, "zdiffrecselhr_tailzgt1", z1 - z0);
          }
        }
        if (is_hiPU) {
          Fill(h, "zdiffrecselhr_hipu", z1 - z0);
        }
        if ((fabs(z1) < 5.) && (fabs(z0) < 5)) {
          Fill(h, "zdiffrecselhrzlt5", z1 - z0);
        }
        Fill(h, "zdiffrec4f", z1 - z0);
        if (instBXLumi_ > 0) {
          Fill(h, "zdiffrecselhrvsL", 1e3 * instBXLumi_, z1 - z0);
        }

        if (DO_ZVSZ_ANALYSIS) {
          Fill(h, "zvszrecsel", v->position().z(), v1->position().z());
          Fill(h, "zvszrecselref", v->position().z());
          Fill(h, "zvszrecselref", v1->position().z());
          Fill(h, "pzvspzsel", TMath::Freq(z0 / sigmaZ_), TMath::Freq(z1 / sigmaZ_));
        }
      }

      if (DO_ZDIFFVSZ_ANALYSIS) {
        if ((v->ndof() > 6) && (v1->ndof() > 6)) {
          Fill(h, "zdiffvsz6", z1 - z0, zbar);
          Fill(h, "zdiffrec6", z1 - z0);
          Fill(h, "zdiffvszp6", z1 - z0, zbarp);
          if (nrec6 == 2)
            Fill(h, "zdiffvsz6Nv2", z1 - z0, zbar);
        }

        if ((v->ndof() > 7) && (v1->ndof() > 7)) {
          Fill(h, "zdiffvsz7", z1 - z0, zbar);
          Fill(h, "zdiffrec7", z1 - z0);
          Fill(h, "zdiffvszp7", z1 - z0, zbarp);
          if (nrec7 == 2)
            Fill(h, "zdiffvsz7Nv2", z1 - z0, zbar);
        }

        if ((v->ndof() > 8) && (v1->ndof() > 8)) {
          Fill(h, "zdiffvsz8", z1 - z0, zbar);
          Fill(h, "zdiffrec8", z1 - z0);
          Fill(h, "zdiffvszp8", z1 - z0, zbarp);
          if (nrec8 == 2)
            Fill(h, "zdiffvsz8Nv2", z1 - z0, zbar);
        }
        if ((v->ndof() > 12) && (v1->ndof() > 12)) {
          Fill(h, "zdiffvsz12", z1 - z0, zbar);
          Fill(h, "zdiffrec12", z1 - z0);
          Fill(h, "zdiffvszp12", z1 - z0, zbarp);
        }
        if ((v->ndof() > 20) && (v1->ndof() > 20)) {
          Fill(h, "zdiffrec20", z1 - z0);
        }
      }
    }

    // is it isolated?
    double deltaz = 1e10;
    for (reco::VertexCollection::const_iterator v1 = recVtxs->begin(); v1 != recVtxs->end(); ++v1) {
      if (v->position().z() == v1->position().z())
        continue;
      if (fabs(v->position().z() - v1->position().z()) < fabs(deltaz))
        deltaz = v->position().z() - v1->position().z();
    }
    if (fabs(deltaz) > 2.0)
      Fill(h, "vtxndfIso", v->ndof());

  }  // vertex loop (v)

  if (true) {
    double zrange = 30.;

    float nvtxsel = 0;
    float nvtxsel_2 = 0;  // select level 2
    float nvtxsel_3 = 0;  // select level 2
    float nvtxsel1sigmaz = 0;
    float nvtxsel2sigmaz = 0;
    float nvtxsel3sigmaz = 0;
    float nvtxselgt1sigmaz = 0;
    float nvtxselgt2sigmaz = 0;
    float nvtxselrgt200 = 0;
    float nvtxselrlt200 = 0;
    float nvtxselndof10 = 0;
    float nvtxselndof20 = 0;
    float nvtxselpxytight = 0;
    float nvtxselaptsumgt2 = 0;
    float nvtxselaptsumlt2 = 0;
    float nvtxselptmax2gt04 = 0;
    float nvtxselptmax2lt04 = 0;

    float nvtx[nzbins_];
    std::fill_n(nvtx, nzbins_, 0);

    //  clear the dzzbin_event histograms
    if (DO_DENSITY_ANALYSIS) {
      for (unsigned int i = 0; i < dzbins_.size() - 1; i++) {
        h[Form("dzzbin_event_recsel%d", i)]->Reset();
        h[Form("dzzbin_event_recsel_2_%d", i)]->Reset();
        h[Form("dzzbin_event_recreject%d", i)]->Reset();
      }
      h["inclusive_event"]->Reset();
      h["inclusive_event_2_"]->Reset();
      h["inclusive_event_ndoflt10"]->Reset();
      h["inclusive_event_ndofgt10"]->Reset();
      h["inclusive_event_ndofgt20"]->Reset();
    }

    //--------------------- vertex pair loop
    for (auto v1 = recVtxs->begin(); v1 != recVtxs->end(); ++v1) {
      if (v1->isFake())
        continue;
      double z1 = v1->position().z() - vertexBeamSpot_.z0() - dzb_;
      if (fabs(z1) > zrange)
        continue;
      if (select(*v1, 2)) {
        Fill(h, "inclusive_event_2_", fabs(z1));
      }
      if (select(*v1)) {
        Fill(h, "inclusive_event", fabs(z1));  // remember the factor of two!! v17 and higher
        if (v1->ndof() < 10)
          Fill(h, "inclusive_event_ndoflt10", fabs(z1));
        if (v1->ndof() > 10)
          Fill(h, "inclusive_event_ndofgt10", fabs(z1));
        if (v1->ndof() > 20)
          Fill(h, "inclusive_event_ndofgt20", fabs(z1));
        nvtxsel++;
        if (fabs(z1) < 1. * sigmaZ_) {
          nvtxsel1sigmaz++;
        } else {
          nvtxselgt1sigmaz++;
        }
        if (fabs(z1) < 2. * sigmaZ_) {
          nvtxsel2sigmaz++;
        } else {
          nvtxselgt2sigmaz++;
        }
        if (fabs(z1) < 3. * sigmaZ_)
          nvtxsel3sigmaz++;
        if (v1->ndof() > 10.)
          nvtxselndof10++;
        if (v1->ndof() > 20.)
          nvtxselndof20++;
        if (vertex_pxy(*v1) > 0.1)
          nvtxselpxytight++;
        if (vertex_r(*v1) > 0.02) {
          nvtxselrgt200++;
        } else {
          nvtxselrlt200++;
        }
        if (vertex_aptsum(*v1) > 2.0) {
          nvtxselaptsumgt2++;
        } else {
          nvtxselaptsumlt2++;
        }
        if (vertex_ptmax2(*v1) > 0.4) {
          nvtxselptmax2gt04++;
        } else {
          nvtxselptmax2lt04++;
        }

        Fill(h, "vtxz_recsel", z1);
        Fill(h, "vtxzhr_recsel", z1);
        if (is_hiPU) {
          Fill(h, "vtxz_recsel_hipu", z1);
          Fill(h, "vtxzhr_recsel_hipu", z1);
        }
        Fill(h, "vtxzbybx_recsel", z1, (float)bunchCrossing_);
        int zbin1 = nzbins_ * (z1 + zbinmax_) / (2 * zbinmax_);
        if ((zbin1 >= 0) && (zbin1 < nzbins_)) {
          nvtx[zbin1]++;
        }
      }

      if (select(*v1, 2)) {
        nvtxsel_2++;
        Fill(h, "vtxz_recsel_2", z1);
        Fill(h, "vtxzhr_recsel_2", z1);
        if (is_hiPU)
          Fill(h, "vtxz_recsel_2_hipu", z1);
        if (is_hiPU)
          Fill(h, "vtxzhr_recsel_2_hipu", z1);
      }

      if (select(*v1, 3)) {
        nvtxsel_3++;
        Fill(h, "vtxz_recsel_3", z1);
        Fill(h, "vtxzhr_recsel_3", z1);
        if (is_hiPU)
          Fill(h, "vtxz_recsel_3_hipu", z1);
        if (is_hiPU)
          Fill(h, "vtxzhr_recsel_3_hipu", z1);
      }

      if (v1 == recVtxs->end())
        continue;
      for (auto v2 = v1 + 1; v2 != recVtxs->end(); ++v2) {
        if (v2->isFake())
          continue;
        double z2 = v2->position().z() - vertexBeamSpot_.z0() - dzb_;
        if (fabs(z2) > zrange)
          continue;
        double zbar = std::abs(0.5 * (z1 + z2));  // note that by filling abs(z) we get a factor two in density
        double dz = std::abs(z1 - z2);

        /*
	if( select( *v1 ) && select( *v2 )){
	  Fill(h, "dzz_recsel", dz,  zbar );
	}else{
	  Fill(h, "dzz_recreject", dz, zbar );
	}
	*/

        if (DO_DENSITY_ANALYSIS) {
          for (unsigned int dzbin = 0; dzbin < dzbins_.size() - 1; dzbin++) {
            if ((dz >= dzbins_[dzbin]) && (dz < dzbins_[dzbin + 1])) {
              if (select(*v1) && select(*v2)) {
                Fill(h, Form("dzzbin_event_recsel%d", dzbin), zbar);
              } else {
                Fill(h, Form("dzzbin_event_recreject%d", dzbin), zbar);
              }

              if (select(*v1, 2) && select(*v2, 2)) {
                Fill(h, Form("dzzbin_event_recsel_2_%d", dzbin), zbar);
              }
            }
          }
        }
      }
    }
    //--------------------- end of vertex pair loop

    if (DO_DENSITY_ANALYSIS) {
      // now fill the dzzbin projections from per-event "histograms"
      for (unsigned int dzbin = 0; dzbin < dzbins_.size() - 1; dzbin++) {
        TH1* hrecsel = h[Form("dzzbin_event_recsel%d", dzbin)];
        TH1* hrecsel_2 = h[Form("dzzbin_event_recsel_2_%d", dzbin)];
        TH1* hrecreject = h[Form("dzzbin_event_recreject%d", dzbin)];
        double dz = 0.5 * (dzbins_[dzbin] + dzbins_[dzbin + 1]);
        double frho = TMath::Gaus(dz, 0., 2 * sigmaZ_, false);
        for (unsigned int b = 1; b < 201; b++) {
          double zbar = hrecsel->GetBinCenter(b);
          double wzbar = hrecsel->GetBinWidth(b);
          double rho = lumiPU_ * TMath::Gaus(zbar, 0., sigmaZ_, true);  // beam center already subtracted
          Fill(h, Form("dzzbin_recsel%d", dzbin), rho, hrecsel->GetBinContent(b) / wzbar);
          Fill(h, Form("dzzbin_recsel_2_%d", dzbin), rho, hrecsel_2->GetBinContent(b) / wzbar);
          Fill(h, Form("dzzbin2_recsel%d", dzbin), rho * frho, hrecsel->GetBinContent(b) / wzbar);
          Fill(h, Form("dzzbin_recreject%d", dzbin), rho, hrecreject->GetBinContent(b) / wzbar);
        }
      }

      // same for the single inclusive test histo
      for (unsigned int b = 1; b < 201; b++) {
        double zbar = h["inclusive_event"]->GetBinCenter(b);
        double wzbar = h["inclusive_event"]->GetBinWidth(b);
        double rho = lumiPU_ * TMath::Gaus(zbar, 0., sigmaZ_, true);  // beam center already subtracted
        Fill(h,
             "inclusive",
             rho,
             0.5 * h["inclusive_event"]->GetBinContent(b) /
                 wzbar);  // the factor 0.5 compensate histogramming |z| instead of z
        Fill(h,
             "inclusive_2_",
             rho,
             0.5 * h["inclusive_event_2_"]->GetBinContent(b) /
                 wzbar);  // the factor 0.5 compensate histogramming |z| instead of z
        Fill(h, "inclusive_ndoflt10", rho, 0.5 * h["inclusive_event_ndoflt10"]->GetBinContent(b) / wzbar);
        Fill(h, "inclusive_ndofgt10", rho, 0.5 * h["inclusive_event_ndofgt10"]->GetBinContent(b) / wzbar);
        Fill(h, "inclusive_ndofgt20", rho, 0.5 * h["inclusive_event_ndofgt20"]->GetBinContent(b) / wzbar);
      }
    }

    Fill(h, "lPUbybx", (float)bunchCrossing_, lumiPU_);
    Fill(h, "nselvtxbybx", (float)bunchCrossing_, nvtxsel);
    Fill(h, "nselvtxbyLS", (float)luminosityBlock_, nvtxsel);
    Fill(h, "nselvtx", nvtxsel);
    Fill(h, "LPU", lumiPU_);
    Fill(h, "lBX", 1.e3 * instBXLumi_);
    Fill(h, "nlPU", lumiPU_);
    Fill(h, "lPULS", lsglb_, lumiPU_);
    Fill(h, "sigmazLS", lsglb_, sigmaZ_);
    Fill(h, "zbeamLS", lsglb_, vertexBeamSpot_.z0());
    if (lumiPU_ > 0) {
      Fill(h, "nselvtxoverPUvsLS", lsglb_, nvtxsel / lumiPU_);
      if (nvtxsel > 30) {
        Fill(h, "f10vsLS", lsglb_, 1. - nvtxselndof10 / nvtxsel);
        Fill(h, "f20vsLS", lsglb_, 1. - nvtxselndof20 / nvtxsel);
      }

      if (lumiPU_ > 40) {
        if (is_tail) {
          Fill(h, "nvtxtailvsLS", lsglb_, 1.);
        } else {
          Fill(h, "nvtxtailvsLS", lsglb_, 0.);
        }
      }
      if (lumiPU_ > 60.) {
        Fill(h, "nselvtxoverPU60-80vsLS", lsglb_, nvtxsel / lumiPU_);
      } else if (lumiPU_ > 40.) {
        Fill(h, "nselvtxoverPU40-60vsLS", lsglb_, nvtxsel / lumiPU_);
      } else if (lumiPU_ > 20.) {
        Fill(h, "nselvtxoverPU20-40vsLS", lsglb_, nvtxsel / lumiPU_);
      } else {
        Fill(h, "nselvtxoverPU00-20vsLS", lsglb_, nvtxsel / lumiPU_);
      }
    }

    Fill(h, "LPUvsLPU", lumiPU_, lumiPU_);
    Fill(h, "LPUvsavgLPU", avglumiPU_, lumiPU_);
    Fill(h, "avgLPUvsLPU", lumiPU_, avglumiPU_);
    Fill(h, "nselvtxvsLPU", lumiPU_, nvtxsel);
    Fill(h, "nselvtxvssimPU", simPU_, nvtxsel);
    Fill(h, "nselvtx_2vsLPU", lumiPU_, nvtxsel_2);
    Fill(h, "nselvtxvsLPUpxytight", lumiPU_, nvtxselpxytight);
    Fill(h, "nselvtxvsLPUndof10", lumiPU_, nvtxselndof10);
    Fill(h, "nselvtxvsLPUndoflt10", lumiPU_, nvtxsel - nvtxselndof10);
    Fill(h, "nselvtxvsLPUndof20", lumiPU_, nvtxselndof20);
    Fill(h, "nselvtxvsLPU1sigmaz", lumiPU_, nvtxsel1sigmaz);
    Fill(h, "nselvtxvsLPU2sigmaz", lumiPU_, nvtxsel2sigmaz);
    Fill(h, "nselvtxvsLPU3sigmaz", lumiPU_, nvtxsel3sigmaz);
    Fill(h, "nselvtxvsLPUgt1sigmaz", lumiPU_, nvtxselgt1sigmaz);
    Fill(h, "nselvtxvsLPUgt2sigmaz", lumiPU_, nvtxselgt2sigmaz);
    Fill(h, "nselvtxvsLPUaptsumgt2", lumiPU_, nvtxselaptsumgt2);
    Fill(h, "nselvtxvsLPUaptsumlt2", lumiPU_, nvtxselaptsumlt2);
    Fill(h, "nselvtxvsLPUptmax2gt04", lumiPU_, nvtxselptmax2gt04);
    Fill(h, "nselvtxvsLPUptmax2lt04", lumiPU_, nvtxselptmax2lt04);

    //Fill(h, "nselvtxgt1vslt1sigmaz", nvtxsel1sigmaz, nvtxselgt1sigmaz);  already existed as "nvtxgt1vslt1sigmaz"
    Fill(h, "nselvtxvsavgLPU", avglumiPU_, nvtxsel);
    Fill(h, "nselvtx1sigmaz", nvtxsel1sigmaz);
    Fill(h, "nselvtx2sigmaz", nvtxsel2sigmaz);
    /*
    Fill(h, "ng1vsLPU", lumiPU_, float(ng1));
    Fill(h, "ng05vsLPU", lumiPU_, float(ng05));
    Fill(h, "ng03vsLPU", lumiPU_, float(ng03));
    Fill(h, "ng02vsLPU", lumiPU_, float(ng02));
    Fill(h, "ng1vsLPUprof", lumiPU_, float(ng1));
    Fill(h, "ng05vsLPUprof", lumiPU_, float(ng05));
    Fill(h, "ng03vsLPUprof", lumiPU_, float(ng03));
    Fill(h, "ng02vsLPUprof", lumiPU_, float(ng02));
    */
    if (firstBXinTrain_) {
      Fill(h, "nselvtxvsLPUfirstBX", lumiPU_, nvtxsel);
    } else {
      Fill(h, "nselvtxvsLPUnotfirstBX", lumiPU_, nvtxsel);
    }
    if (bunchCrossing_ < 2) {
      Fill(h, "nselvtxvsLPUBX0", lumiPU_, nvtxsel);
    }

    Fill(h, "nselvtxvsLPUprof", lumiPU_, nvtxsel);
    Fill(h, "nselvtxvssimPUprof", simPU_, nvtxsel);
    Fill(h, "LPUnselvtxvsLPUprof", lumiPU_, lumiPU_);
    Fill(h, "nselvtxvsLPU1sigmazprof", lumiPU_, nvtxsel1sigmaz);
    Fill(h, "nselvtxvsLPU2sigmazprof", lumiPU_, nvtxsel2sigmaz);
    Fill(h, "nselvtxvsLPU3sigmazprof", lumiPU_, nvtxsel3sigmaz);
    Fill(h, "nselvtxvsavgLPUprof", avglumiPU_, nvtxsel);
    Fill(h, "nselvtxvsLPUrlt200", lumiPU_, nvtxselrlt200);
    Fill(h, "nselvtxvsLPUrgt200", lumiPU_, nvtxselrgt200);

    //Fill(h, "nselvtxvsLPU1sigmazprof", lumiPU_, nvtxsel1sigmaz); // duplicate

    //if (nvtxsel>70){
    //  reportEvent(Form("high-nvtx-event  nvtxsel = %5.1f",nvtxsel),false);
    //}

    Fill(h, "nvtxgt1vslt1sigmaz", float(nvtxsel1sigmaz), float(nvtxselgt1sigmaz));

    /*
    if ((nvtxsel>70) && (lumiPU_<50)){
      reportEvent(Form("nvtx-too-high-for-pu  nvtxsel = %5.1f",nvtxsel), false);
    }

    if ((lumiPU_ >40 ) && (nvtxsel < 7)) {
      reportEvent(Form("low-nvtx-tail  nvtxsel = %5.1f  expected = %5.1f",nvtxsel, lumiPU_ ), false);
    }
    */

    Fill(h, "nvtxgt1vslt1sigmaz", float(nvtxsel1sigmaz), float(nvtxselgt1sigmaz));
    if ((nvtxsel1sigmaz < 2) && (nvtxselgt1sigmaz > 12)) {
      reportEvent(Form("no lt1sigma vertices  n(<1sigma)= %3d ,  n(>1sigma) = %3d, sigmaZ=%5.2f",
                       int(nvtxsel1sigmaz),
                       int(nvtxselgt1sigmaz),
                       sigmaZ_));
    }
    if ((nvtxsel1sigmaz > 2) && (nvtxsel1sigmaz < 20) && (nvtxselgt1sigmaz > 25 + 0.5 * nvtxsel1sigmaz)) {
      reportEvent(Form("anomalous lt1sigma vs gt1sigma  n(< 1sigma)= %3d ,  n(>1sigma) = %3d, sigmaZ=%5.2f",
                       int(nvtxsel1sigmaz),
                       int(nvtxselgt1sigmaz),
                       sigmaZ_));
    }

    if (is_tail) {
      //reportEvent(Form("high-nvtx-tail  event nvtxsel = %5.1f,   expected = %5.1f   sigmaz_event=%5.2f  napt=%5.1f",nvtxsel,  lumiPU_ * 0.8, sigmaz_event, nvtxselaptsumlt2), false);
      Fill(h, "vtxz_tag_tail", ztag);
      Fill(h, "nvtxgt1vslt1sigmaz_tail", float(nvtxsel1sigmaz), float(nvtxselgt1sigmaz));
      //if(( nvtxselaptsumlt2 > 5) &&  (nvtxselgt1sigmaz > (0.24*lumiPU_ + 5*sqrt(0.24*lumiPU_)))){
      //if(( nvtxselaptsumlt2 > 4)){
      //	reportEvent(Form("very-high-nvtx-tail  event nvtxsel = %5.1f,   expected = %5.1f   sigmaz_event=%5.2f  napt=%5.1f",nvtxsel,  lumiPU_ * 0.8, sigmaz_event, nvtxselaptsumlt2), true);      }
    }

    if (is_hiPU) {
      Fill(h, "nvtxgt1vslt1sigmaz_hipu", float(nvtxsel1sigmaz), float(nvtxselgt1sigmaz));
    }
  }
}  // end of analyzeVertexCollectionRecoNoTracks

//define this as a plug-in
DEFINE_FWK_MODULE(PrimaryVertexAnalyzer4PU);
