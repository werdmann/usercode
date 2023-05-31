#include "usercode/PrimaryVertexAnalyzer/interface/PrimaryVertexAnalyzer4PU.h"
#include "usercode/PrimaryVertexAnalyzer/interface/FFA.h"
//  std::cout << "XDBG " << __func__ << " : " << __LINE__ << std::endl;

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "MagneticField/Engine/interface/MagneticField.h"

// reco track and vertex
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/deltaPhi.h"
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
#include <TRandom.h>
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

// some of this is probably redundant or unused
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
#include <algorithm>


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
      vecPileupSummaryInfoToken_(consumes<std::vector<PileupSummaryInfo>>(edm::InputTag(std::string("slimmedAddPileupInfo")))),
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
      lumiInfoToken_(consumes<LumiInfo>(iConfig.getUntrackedParameter<edm::InputTag>("lumiInfoTag"))),
      pdtToken_(esConsumes()),
      transientTrackBuilderToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
      trackerTopologyToken_(esConsumes())
{

  // v34 and higher, avoid biases
  minNumberOfRecTrks_ = 0.;  // FIXME make this configurable or maybe even better obsolete (are these empty BXs?)
  minNumberOfSelTrks_ = 0.;  // FIXME make this configurable

  // definition of visible vertices
  etaMaxVisible_ = 4.0;
  ptMinVisible_ = 0.2;
  numTrkHitsVisible_ = 4;

  fill_track_histos_ = iConfig.getUntrackedParameter<bool>("fill_track_histos", false);
  selNdofNoBS_ = iConfig.getUntrackedParameter<double>("selNdof", 4.);
  std::cout << "PrimaryVertexAnalyzer4PU: selNDof_(noBS) = " << selNdofNoBS_ << std::endl;
  selNdofWithBS_ = iConfig.getUntrackedParameter<double>("selNdofWithBS", 7.);
  std::cout << "PrimaryVertexAnalyzer4PU: selNDofWithBS_ = " << selNdofWithBS_ << std::endl;
  selNdof_ = selNdofNoBS_; // to be changed later according to the collection name
  ndof0trk_ = 0.;

  min_trk_in_vtx_weight_ = 0.2;

  f4D_ = (iConfig.getUntrackedParameter<bool>("f4D", true));
  frefit_ = (iConfig.getUntrackedParameter<bool>("frefit", false));
  if (f4D_) {
    std::cout << "PrimaryVertexAnalyzer4PU: TrackTimeResosLabel "
              << iConfig.getUntrackedParameter<edm::InputTag>("TrackTimeResosLabel") << std::endl;
    trkTimesToken_ = consumes<edm::ValueMap<float>>(iConfig.getUntrackedParameter<edm::InputTag>("TrackTimesLabel"));
    trkTimeResosToken_ =
      consumes<edm::ValueMap<float>>(iConfig.getUntrackedParameter<edm::InputTag>("TrackTimeResosLabel"));  
    trkTimeQualityThreshold_ = iConfig.getUntrackedParameter<double>("TrackTimeQualityThreshold", 0.8);
    trkTimeQualityToken_ =  consumes<edm::ValueMap<float> >(iConfig.getUntrackedParameter<edm::InputTag>("TrackTimeQualityMapLabel"));//mtdTrackQualityMVA:mtdQualMVA"

    MTD_pathlength_Token_ = consumes<edm::ValueMap<float> >(edm::InputTag("trackExtenderWithMTD:generalTrackPathLength")); 
    MTD_time_Token_ = consumes<edm::ValueMap<float> >(edm::InputTag("trackExtenderWithMTD:generalTracktmtd"));
    MTD_timeerror_Token_ = consumes<edm::ValueMap<float> >(edm::InputTag("trackExtenderWithMTD:generalTracksigmatmtd")); 
    MTD_momentum_Token_ = consumes<edm::ValueMap<float> >(edm::InputTag("trackExtenderWithMTD:generalTrackp")); 
  } else {
    std::cout << "PrimaryVertexAnalyzer4PU: no timing" << std::endl;
  }



  MINIAOD_ = iConfig.getUntrackedParameter<bool>("MINIAOD", false);
 
  if(MINIAOD_){
    report_counted("MINIAOD ",1);
    theTracksToken_= consumes<edm::View<pat::PackedCandidate> >(edm::InputTag("packedPFCandidates"));
    theLostTracksToken_= consumes<edm::View<pat::PackedCandidate> >(edm::InputTag("lostTracks"));
    theGenParticlesToken_= consumes<edm::View<pat::PackedGenParticle> >(edm::InputTag("packedGenParticles"));
    theGenParticlesXyz0Token_= consumes<GenEventVertex>(edm::InputTag("genParticles","xyz0","HLT"));// 
    theGenParticlesT0Token_= consumes<float>(edm::InputTag("genParticles","t0","HLT"));//
    thePrunedGenParticlesToken_ = consumes<edm::View<reco::GenParticle> >(edm::InputTag("prunedGenParticles"));
  }

  
  reco_vertex_collections_ = iConfig.getParameter<std::vector<edm::InputTag>>("vertexRecoCollections");

  for (auto const& l : reco_vertex_collections_) {
    auto l_encode = l.encode();
    std::replace(l_encode.begin(), l_encode.end(), ':', '_');
    std::cout << "PrimaryVertexAnalyzer4PU: vertex collection [" << l_encode << "]" << std::endl;
    vertexCollectionLabels_.push_back(l_encode);

    auto token = edm::EDGetTokenT<reco::VertexCollection>(consumes<reco::VertexCollection>(edm::InputTag(l)));
    vertexCollectionTokens_[l_encode] = token;
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
    cout << "PrimaryVertexAnalyzer4PU: veryverbose = " << veryverbose_ << endl;
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
    cout << "PrimaryVertexAnalyzer4PU: zmatch =" << zmatch_ << endl;
    cout << "     cut-off for matching sim to reco by z" << endl;
    cout << endl;
  }

  zWosMatchMax_ = iConfig.getUntrackedParameter<double>("zMatchMax", 1.);
  if (verbose_) {
    cout << "PrimaryVertexAnalyzer4PU: zMatchMax = " << zWosMatchMax_ << endl;
    cout << "     cut-off for insane recvertex <-> simvertex  matches" << endl;
    cout << "     (TrackingParticles, matching by weight/sigma^2) " << endl;
    cout << "     default is 1 cm, configurable for exotic events where" << endl;
    cout << "     all tracks appear far away from the vertex   " << endl;
    // such as LongLivedChi0ToNuLL_MSquark-1000_MChi-148_TuneZ2Star_8TeV-pythia6
    cout << endl;
  }
  
  zwindow_sigmas_ = iConfig.getUntrackedParameter<double>("zMatchWindow_sigmas", 3.0);
  if (verbose_) {
    cout << "PrimaryVertexAnalyzer4PU: zwindow_sigmas_ = " << zwindow_sigmas_ << endl;
    cout << "     for z-coordinate based recvertex <-> simvertex  matchin" << endl;
    cout << "     window size in multiples of  the reconstructed sigma(z)" << endl;
    cout << endl;
  }

  RECO_ = iConfig.getUntrackedParameter<bool>("RECO", false);
  if (verbose_) {
    cout << "PrimaryVertexAnalyzer4PU: RECO = " << RECO_ << endl;
    cout << "      use RECO information (pixel hits and some trackextra)" << endl;
    cout << endl;
  }

  autoDumpCounter_ = iConfig.getUntrackedParameter<int>("autodump", 0);
  if (verbose_) {
    cout << "PrimaryVertexAnalyzer4PU: autodump = " << autoDumpCounter_ << endl;
    cout << "      dump detailed information about the first <autodump> events " << endl;
    cout << endl;
  }

  trkhiptmin_ = 3.0;
  trkloptmax_ = 1.0;
  trkcentraletamax_ = 1.5;
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
  addn(h, new TProfile("tfitvsnsel", "vertex fitting cpu-time vs nvertex", 300, 0., 300., 0., 1.e8));
  addn(h, new TProfile("ttimevsnsel", "vertex timing cpu-ime vs nvertex", 300, 0., 300., 0., 1.e8));

  addn(h, new TProfile("tcluvsLPU", "clustering time vs pu", 300, 0., 300., 0., 1.e8));
  addn(h, new TProfile("tfitvsLPU", "vertex fitting cpu-time vs pu", 300, 0., 300., 0., 1.e8));
  addn(h, new TProfile("ttimevsLPU", "vertex timing cpu-time vs pu", 300, 0., 300., 0., 1.e8));

  addn(h, new TProfile("tcluvsSimPU", "clustering time vs #simvtx", 300, 0., 300., 0., 1.e8));
  addn(h, new TProfile("tfitvsSimPU", "vertex fitting cpu-time vs #simvtx", 300, 0., 300., 0., 1.e8));
  addn(h, new TProfile("ttimevsSimPU", "vertex timing cpu-time vs #simvtx", 300, 0., 300., 0., 1.e8));
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
  const float etarange = 4.;
  const int netabin = 40.;

  //const int nzbin_normal = 30; until used
  const int nzbin_normal_fine = 60;
  const float zrange_normal = 15.;

  const int nrecmax = 300;
  const int ntrkmax = 6000.;

  int nvtxbin = 150;
  float nvtxrange = 300.;

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
  add(h, new TH1F("zdiffrecselsignalfakev2", "z-distance between reconstructed Signal and fake selected vertices", nbinzdiffrec, 0., 2.)); // eventually get rid of this
  add(h, new TH1F("zdiffsimselsignalrealpu", "simulated z-distance between reconstructed Signal and real selected vertices", nbinzdiffrec, 0., 2.));
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

  add(h, new TH1F("ntpreal", "found vertices", 300, 0., 300.));
  add(h, new TH1F("ntprealsel", "found selected vertices", 300, 0., 300.));
  add(h, new TH1F("ntpfake", "fake vertices", 300, 0., 300.));
  add(h, new TH1F("ntpfakesel", "fake selected vertices", 300, 0., 300.));
  add(h, new TH1F("ntpotherfake", "other fake", 10, 0., 10.));
  add(h, new TH1F("ntpotherfakesel", "selected other fake", 10, 0., 10.));
  add(h, new TH1F("ntpsplit", "split", 10, 0., 10.));
  add(h, new TH1F("ntpsplitsel", "selected split", 10, 0., 10.));
  add(h, new TH1F("ntpsplitfromsignal", "reconstruced split from signal", 10, 0., 10.));
  add(h, new TH1F("ntpsplitfrompu", "reconstructed split from PU", 100, 0., 100.));
  add(h, new TH1F("ntpsplitselfromsignal", "selected split from signal", 10, 0., 10.));
  add(h, new TH1F("ntpsplitselfrompu", "selected split from PU", 100, 0., 100.));

  // proviles vs pu for selected vertices
  add(h, new TProfile("ntprecselvssimPU", "selected vertices", nvtxbin, 0., nvtxrange, 0., nvtxrange * 2.));
  add(h, new TProfile("ntpfakeselvssimPU", "fake selected vertices", nvtxbin, 0., nvtxrange, 0., nvtxrange * 2.));
  add(h, new TProfile("ntprealselvssimPU", "matched selected vertices", nvtxbin, 0., nvtxrange, 0., nvtxrange * 2.));
  add(h, new TProfile("ntpsplitselvssimPU", "split selected vertices", nvtxbin, 0., nvtxrange, 0., nvtxrange * 2.));
  add(h, new TProfile("ntpotherfakeselvssimPU", "other fake selected vertices", nvtxbin, 0., nvtxrange, 0., nvtxrange * 2.));
  add(h, new TProfile("ntpfakeselratevssimPU", "fake selected vertex rate", nvtxbin, 0., nvtxrange, 0., 2.));

  // same as the previous block, but for all reconstructed vertices
  add(h, new TProfile("ntprecvssimPU", "reconstructed vertices", nvtxbin, 0., nvtxrange, 0., nvtxrange * 2.));
  add(h, new TProfile("ntpfakevssimPU", "fake reconstructed vertices", nvtxbin, 0., nvtxrange, 0., nvtxrange * 2.));
  add(h, new TProfile("ntprealvssimPU", "matched reconstructed vertices", nvtxbin, 0., nvtxrange, 0., nvtxrange * 2.));
  add(h, new TProfile("ntpsplitvssimPU", "split reconstructed vertices", nvtxbin, 0., nvtxrange, 0., nvtxrange * 2.));
  add(h, new TProfile("ntpotherfakevssimPU", "other fake reconstructed vertices", nvtxbin, 0., nvtxrange, 0., nvtxrange * 2.));
  add(h, new TProfile("ntpfakeratevssimPU", "fake reconstructed vertex rate", nvtxbin, 0., nvtxrange, 0., 2.));


  add(h, new TH2F("wornk_matchedsel","n70% vs sumwos", 100, -1, 7., 5, 0., 5.));
  add(h, new TH2F("wornk_unmatchedsel","n70% vs sumwos", 100, -1., 7., 5, 0., 5.));
  add(h, new TH2F("wornk_splitsel","n70% vs sumwos", 100, -1., 7., 5, 0., 5.));
  add(h, new TH2F("wornk_othersel","n70% vs sumwos", 100, -1., 7., 5, 0., 5.));
  add(h, new TH2F("wornk_matchedsim","n70% vs sumwos", 100, -1., 7., 5, 0., 5.));
  add(h, new TH2F("wornk_unmatchedsim","n70% vs sumwos", 100, -1., 7., 5, 0., 5.));

  addSP(h, new TProfile("unmatchedFractionVsWeight", "fraction of unmatched tracks vs weight", 101, 0., 1.01, 0., 2.));
  addSP(h, new TProfile("correctVtxFractionVsWeight", "fraction of correctly assigned tracks vs weight", 101, 0., 1.01, 0., 2.));
  addSP(h, new TProfile("correctVtxFractionAllVsWeight", "fraction of matched and correctly assigned tracks vs weight", 101, 0., 1.01, 0., 2.));
  addSP(h, new TH1F("trkWeightCorrectVtx", "weight of correctly assigned tracks", 101, 0., 1.01));
  addSP(h, new TH1F("trkWeightIncorrectVtx", "weight of incorrectly assigned tracks", 101, 0., 1.01));
  addSP(h, new TH1F("trkWeightUnmatched", "weight of unmatched", 101, 0., 1.01));

  addSP(h, new TH1F("vtxtrkpurity", "track purity in vertex", 101, 0., 1.01));
  addSP(h, new TProfile("vtxtrkpurityvspu", "track purity in vertex vs PU", npubin2, 0., npumax, 0., 2.));
  addSP(h, new TProfile("vtxtrkpurityvsz", "track purity", nzbin_wide_fine, -zrange_wide, zrange_wide, 0., 2.));
  addSP(h, new TProfile("vtxtrkpurityvsdz", "track purity", 100, -1., 1., 0., 2.));
  addSP(h, new TProfile("vtxtrkpurityvsdzpull", "track purity", 100, -10., 10., 0., 2.));
  addSP(h, new TProfile("vtxtrkpurityvsdzerror", "track purity vs zError", 100, 0, 0.1, 0., 2.));
  
  addSP(h, new TH1F("vtxtrkallpurity", "track purity in vertex (all weights)", 101, 0., 1.01));
  addSP(h, new TProfile("vtxtrkallpurityvspu", "track purity in vertex vs PU (all weights)", npubin2, 0., npumax, 0., 2.));
  addSP(h, new TProfile("vtxtrkallpurityvsputest", "track purity in vertex vs PU (all weights, test)", npubin2, 0., npumax, 0., 2.));
  addSP(h, new TProfile("vtxtrkallpurityvsz", "track purity (all weights)", nzbin_wide_fine, -zrange_wide, zrange_wide, 0., 2.));

  addSP(h, new TH1F("vtxtrktimingpurity", "timing track purity in vertex", 101, 0., 1.01));
  addSP(h, new TProfile("vtxtrktimingpurityvspu", "timing track purity in vertex vs PU", npubin2, 0., npumax, 0., 2.));
  addSP(h, new TProfile("vtxtrktimingpurityvsz", "timing track purity in vertex", nzbin_wide_fine, -zrange_wide, zrange_wide, 0., 2.));
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
  addSP(h, new TH1F("utrkAssignmentEfficiency", "track to any vertex assignment efficiency", 101, 0., 1.01));
  addSP(h, new TH1F("uprimtrkAssignmentEfficiency", "track to any vertex assignment efficiency", 101, 0., 1.01));
  addSP(h, new TProfile("utrkAssignmentEfficiencyvspu", "track to any vertex assignment efficiency vs PU", npubin2, 0., npumax, 0., 2.));
  addSP(h, new TProfile("utrkAssignmentEfficiencyvszerror", "track to any vertex assignment efficiency vs zerror", 100, 0, 0.1, 0., 2.));
  //
  addSP(h, new TH1F("utrkCorrectAssignmentEfficiency", "track to correct vertex assignment efficiency", 101, 0., 1.01));
  addSP(h, new TH1F("uprimtrkCorrectAssignmentEfficiency", "track to correct vertex assignment efficiency", 101, 0., 1.01));
  addSP(h, new TProfile("utrkCorrectAssignmentEfficiencyvspu", "track to correct vertex assignment efficiency vs PU", npubin2, 0., npumax, 0., 2.));
  addSP(h, new TProfile("utrkCorrectAssignmentEfficiencyvsntrk", "track to correct vertex assignment efficiency vs number of tracks",
			101, 0., 101, 0., 2.));  // last bin is overflow
  addSP(h, new TProfile("utrkCorrectAssignmentEfficiencyvszerror", "track to correct vertex assignment efficiency vs zerror", 200, 0, 0.2, 0., 2.));
  addSP(h, new TProfile("utrkCorrectAssignmentEfficiencyvsdz", "track to correct vertex assignment efficiency vs dz", 100, -1., 1., 0., 2.));
  addSP(h, new TProfile("utrkCorrectAssignmentEfficiencyvsdzpull", "track to correct vertex assignment efficiency vs dz pull", 100, -10., 10., 0., 2.));

  {// analyzeVertexTrackAssociation
    addSP(h, new TH1F("trkVtxAssocEffic_TPMatchedTracks_tpPt", "trkVtxAssocEffic_TPMatchedTracks_tpPt", 100, 0., 5.));
    addSP(h, new TH1F("trkVtxAssocEffic_TPMatchedTracks_tpEta", "trkVtxAssocEffic_TPMatchedTracks_tpEta", 160, -4., 4.));
    addSP(h, new TH1F("trkVtxAssocEffic_TPMatchedTracks_tpEta_tpPt000to001", "trkVtxAssocEffic_TPMatchedTracks_tpEta_tpPt000to001", 160, -4., 4.));
    addSP(h, new TH1F("trkVtxAssocEffic_TPMatchedTracks_tpEta_tpPt001to003", "trkVtxAssocEffic_TPMatchedTracks_tpEta_tpPt001to003", 160, -4., 4.));
    addSP(h, new TH1F("trkVtxAssocEffic_TPMatchedTracks_tpEta_tpPt003to010", "trkVtxAssocEffic_TPMatchedTracks_tpEta_tpPt003to010", 160, -4., 4.));
    addSP(h, new TH1F("trkVtxAssocEffic_TPMatchedTracks_tpEta_tpPt010toINF", "trkVtxAssocEffic_TPMatchedTracks_tpEta_tpPt010toINF", 160, -4., 4.));

    addSP(h, new TH1F("trkVtxAssocEffic_TPMatchedTracksWithRecoVtx_tpPt", "trkVtxAssocEffic_TPMatchedTracksWithRecoVtx_tpPt", 100, 0., 5.));
    addSP(h, new TH1F("trkVtxAssocEffic_TPMatchedTracksWithRecoVtx_tpEta", "trkVtxAssocEffic_TPMatchedTracksWithRecoVtx_tpEta", 160, -4., 4.));
    addSP(h, new TH1F("trkVtxAssocEffic_TPMatchedTracksWithRecoVtx_tpEta_tpPt000to001", "trkVtxAssocEffic_TPMatchedTracksWithRecoVtx_tpEta_tpPt000to001", 160, -4., 4.));
    addSP(h, new TH1F("trkVtxAssocEffic_TPMatchedTracksWithRecoVtx_tpEta_tpPt001to003", "trkVtxAssocEffic_TPMatchedTracksWithRecoVtx_tpEta_tpPt001to003", 160, -4., 4.));
    addSP(h, new TH1F("trkVtxAssocEffic_TPMatchedTracksWithRecoVtx_tpEta_tpPt003to010", "trkVtxAssocEffic_TPMatchedTracksWithRecoVtx_tpEta_tpPt003to010", 160, -4., 4.));
    addSP(h, new TH1F("trkVtxAssocEffic_TPMatchedTracksWithRecoVtx_tpEta_tpPt010toINF", "trkVtxAssocEffic_TPMatchedTracksWithRecoVtx_tpEta_tpPt010toINF", 160, -4., 4.));

    addSP(h, new TH1F("trkVtxAssocEffic_TPMatchedTracksWithCorrectRecoVtx_tpPt", "trkVtxAssocEffic_TPMatchedTracksWithCorrectRecoVtx_tpPt", 100, 0., 5.));
    addSP(h, new TH1F("trkVtxAssocEffic_TPMatchedTracksWithCorrectRecoVtx_tpEta", "trkVtxAssocEffic_TPMatchedTracksWithCorrectRecoVtx_tpEta", 160, -4., 4.));
    addSP(h, new TH1F("trkVtxAssocEffic_TPMatchedTracksWithCorrectRecoVtx_tpEta_tpPt000to001", "trkVtxAssocEffic_TPMatchedTracksWithCorrectRecoVtx_tpEta_tpPt000to001", 160, -4., 4.));
    addSP(h, new TH1F("trkVtxAssocEffic_TPMatchedTracksWithCorrectRecoVtx_tpEta_tpPt001to003", "trkVtxAssocEffic_TPMatchedTracksWithCorrectRecoVtx_tpEta_tpPt001to003", 160, -4., 4.));
    addSP(h, new TH1F("trkVtxAssocEffic_TPMatchedTracksWithCorrectRecoVtx_tpEta_tpPt003to010", "trkVtxAssocEffic_TPMatchedTracksWithCorrectRecoVtx_tpEta_tpPt003to010", 160, -4., 4.));
    addSP(h, new TH1F("trkVtxAssocEffic_TPMatchedTracksWithCorrectRecoVtx_tpEta_tpPt010toINF", "trkVtxAssocEffic_TPMatchedTracksWithCorrectRecoVtx_tpEta_tpPt010toINF", 160, -4., 4.));

    addSP(h, new TH1F("trkVtxAssocPurity_TracksOfRecoVtx_trkPt", "trkVtxAssocPurity_TracksOfRecoVtx_trkPt", 100, 0., 5.));
    addSP(h, new TH1F("trkVtxAssocPurity_TracksOfRecoVtx_trkEta", "trkVtxAssocPurity_TracksOfRecoVtx_trkEta", 160, -4., 4.));
    addSP(h, new TH1F("trkVtxAssocPurity_TracksOfRecoVtx_trkEta_trkPt000to001", "trkVtxAssocPurity_TracksOfRecoVtx_trkEta_trkPt000to001", 160, -4., 4.));
    addSP(h, new TH1F("trkVtxAssocPurity_TracksOfRecoVtx_trkEta_trkPt001to003", "trkVtxAssocPurity_TracksOfRecoVtx_trkEta_trkPt001to003", 160, -4., 4.));
    addSP(h, new TH1F("trkVtxAssocPurity_TracksOfRecoVtx_trkEta_trkPt003to010", "trkVtxAssocPurity_TracksOfRecoVtx_trkEta_trkPt003to010", 160, -4., 4.));
    addSP(h, new TH1F("trkVtxAssocPurity_TracksOfRecoVtx_trkEta_trkPt010toINF", "trkVtxAssocPurity_TracksOfRecoVtx_trkEta_trkPt010toINF", 160, -4., 4.));

    addSP(h, new TH1F("trkVtxAssocPurity_TracksOfRecoVtxMatchedToSimVtx_trkPt", "trkVtxAssocPurity_TracksOfRecoVtxMatchedToSimVtx_trkPt", 100, 0., 5.));
    addSP(h, new TH1F("trkVtxAssocPurity_TracksOfRecoVtxMatchedToSimVtx_trkEta", "trkVtxAssocPurity_TracksOfRecoVtxMatchedToSimVtx_trkEta", 160, -4., 4.));
    addSP(h, new TH1F("trkVtxAssocPurity_TracksOfRecoVtxMatchedToSimVtx_trkEta_trkPt000to001", "trkVtxAssocPurity_TracksOfRecoVtxMatchedToSimVtx_trkEta_trkPt000to001", 160, -4., 4.));
    addSP(h, new TH1F("trkVtxAssocPurity_TracksOfRecoVtxMatchedToSimVtx_trkEta_trkPt001to003", "trkVtxAssocPurity_TracksOfRecoVtxMatchedToSimVtx_trkEta_trkPt001to003", 160, -4., 4.));
    addSP(h, new TH1F("trkVtxAssocPurity_TracksOfRecoVtxMatchedToSimVtx_trkEta_trkPt003to010", "trkVtxAssocPurity_TracksOfRecoVtxMatchedToSimVtx_trkEta_trkPt003to010", 160, -4., 4.));
    addSP(h, new TH1F("trkVtxAssocPurity_TracksOfRecoVtxMatchedToSimVtx_trkEta_trkPt010toINF", "trkVtxAssocPurity_TracksOfRecoVtxMatchedToSimVtx_trkEta_trkPt010toINF", 160, -4., 4.));

    addSP(h, new TH1F("trkVtxAssocPurity_TracksOfRecoVtxWithCorrectSimVtxMatch_trkPt", "trkVtxAssocPurity_TracksOfRecoVtxWithCorrectSimVtxMatch_trkPt", 100, 0., 5.));
    addSP(h, new TH1F("trkVtxAssocPurity_TracksOfRecoVtxWithCorrectSimVtxMatch_trkEta", "trkVtxAssocPurity_TracksOfRecoVtxWithCorrectSimVtxMatch_trkEta", 160, -4., 4.));
    addSP(h, new TH1F("trkVtxAssocPurity_TracksOfRecoVtxWithCorrectSimVtxMatch_trkEta_trkPt000to001", "trkVtxAssocPurity_TracksOfRecoVtxWithCorrectSimVtxMatch_trkEta_trkPt000to001", 160, -4., 4.));
    addSP(h, new TH1F("trkVtxAssocPurity_TracksOfRecoVtxWithCorrectSimVtxMatch_trkEta_trkPt001to003", "trkVtxAssocPurity_TracksOfRecoVtxWithCorrectSimVtxMatch_trkEta_trkPt001to003", 160, -4., 4.));
    addSP(h, new TH1F("trkVtxAssocPurity_TracksOfRecoVtxWithCorrectSimVtxMatch_trkEta_trkPt003to010", "trkVtxAssocPurity_TracksOfRecoVtxWithCorrectSimVtxMatch_trkEta_trkPt003to010", 160, -4., 4.));
    addSP(h, new TH1F("trkVtxAssocPurity_TracksOfRecoVtxWithCorrectSimVtxMatch_trkEta_trkPt010toINF", "trkVtxAssocPurity_TracksOfRecoVtxWithCorrectSimVtxMatch_trkEta_trkPt010toINF", 160, -4., 4.));

    addSP(h, new TProfile("trkVtxAssocPurity_vs_pu", "track-vertex purity vs PU", npubin2, 0., npumax, 0., 2.));
    addSP(h, new TProfile("trkVtxAssocPurity_vs_vz", "track-vertex purity vs vz", nzbin_wide_fine, -zrange_wide, zrange_wide, 0., 1.));

    addSP(h, new TProfile("trkVtxAssocPurityWithoutFakeRecoVtxs_vs_pu", "track-vertex purity vs PU", npubin2, 0., npumax, 0., 2.));
    addSP(h, new TProfile("trkVtxAssocPurityWithoutFakeRecoVtxs_vs_vz", "track-vertex purity vs vz", nzbin_wide_fine, -zrange_wide, zrange_wide, 0., 1.));
  }

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
  addSP(h, new TProfile("effallvssimpu", "efficiency (all) vs #PU", 30, 0., 300., 0., 1.));
  addSP(h, new TProfile("effselvssimpu", "efficiency (selected) vs #PU", 30, 0., 300., 0., 1.));
  addSP(h, new TProfile("effsel3sigmavssimpu", "efficiency (selected) vs #PU", 30, 0., 300., 0., 1.));

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
  add(h, new TH1F("zdiffsimfoundselTP", "delta-zsim of found simulated distance with at least 4 tracks, selected", 4000, -10., 10));
  add(h, new TH1F("zdiffsimfoundselSignalTP", "delta-zsim of found simulated distance Signal-PU, ntrk>3, selected", 4000, -10., 10));
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

  add(h, new TH1F("nselvtx", "# selected vertices", 300, 0., 300.));
  add(h, new TH1F("nselvtx1sigmaz", "# selected vertices |z|<1 sigma", 300, 0., 300.));
  add(h, new TH1F("nselvtx2sigmaz", "# selected vertices |z|<2 sigma", 300, 0., 300.));
  add(h, new TH1F("LPU", "expected PU", 300, 0., 300.));

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

  add(h, new TH1F("zvtx_rec", "z", 400, -20., 20.));
  add(h, new TH1F("zvtx_sel0", "z", 400, -20., 20.));
  add(h, new TH1F("zvtx_sel1", "z", 400, -20., 20.));
  add(h, new TH1F("zvtx_sel2", "z", 400, -20., 20.));
  add(h, new TH1F("zvtx_sel3", "z", 400, -20., 20.));

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



  double* dzbins = &dzbins_[0];  // I love C++
  unsigned int ndzbin = dzbins_.size() - 1;


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


  //  for analyzeVertexCollectionSimTracks , i.e. we have sim tracks for the signal vertex and only mediocre truth-matching
  add(h, new TH1F("sigwosfracvsdzsignal", "signal-wos fraction vs distance to signal", 400, -1.0, 1.0));
  add(h, new TH1F("sigwntfracvsdzsignal", "signal-wnt fraction vs distance to signal", 400, -1.0, 1.0));
  add(h, new TH1F("sigwosfracmax", "highest signal-wos fraction of any rec vertex", 20, 0.0, 1.0));
  add(h, new TH1F("sigwntfracmax", "highest signal-wnt fraction of any rec vertex ", 20, 0.0, 1.0));
  add(h, new TH1F("sigwosfractagged", "signal-wos fraction of the tagged rec vertex", 20, 0.0, 1.0));
  add(h, new TH1F("sigwntfractagged", "signal-wnt fraction of the tagged rec vertex ", 20, 0.0, 1.0));
  add(h, new TProfile("nrecwithsigwos", "rec vertices with signal-wos above threshold vs threshold", 20, 0., 1., 0, 1000));
  add(h, new TProfile("nrecwithsigwnt", "rec vertices with signal-wnt above threshold vs threshold", 20, 0., 1., 0, 1000));
  

  // for analyzeVertexCollectionZmatching
  add(h, new TProfile("zmatcheffvspu", "eff vs sim pu", 300, 0., 300., 0., 2.));
  add(h, new TProfile("zmatchambigvspu", "ambiguity vs sim pu", 300, 0., 300., 0., 2.));
  add(h, new TProfile("zmatchfakevspu", "fake fraction vs sim pu", 300, 0., 300., 0., 2.));
  add(h, new TProfile("zmatchfakeallpu", "fake fraction vs sim pu", 300, 0., 300., 0., 2.));
  add(h, new TH1F("zmatchnsimmatch", "sim matches per rec vtx", 20, 0., 20.));
  add(h, new TH1F("zmatchnrecmatch", "rec matches per sim vtx", 20, 0., 20.));

  add(h, new TProfile("zcmatcheffvspu", "eff vs sim pu", 300, 0., 300., 0., 2.));
  add(h, new TProfile("zcmatchambigvspu", "ambiguity vs sim pu", 300, 0., 300., 0., 2.));
  add(h, new TProfile("zcmatchfakevspu", "fake fraction vs sim pu", 300, 0., 300., 0., 2.));
  add(h, new TH1F("zcmatchnsimmatch", "sim matches per rec vtx", 20, 0., 20.));
  add(h, new TH1F("zcmatchnrecmatch", "rec matches per sim vtx", 20, 0., 20.));

  add(h, new TProfile("ztpmatcheffvspu", "efficiency vs sim pu", 300, 0., 300., 0., 2.));
  add(h, new TProfile("ztpmatchambigvspu", "ambiguity vs sim pu", 300, 0., 300., 0., 2.));
  add(h, new TProfile("ztpmatchfakevspu", "fake fraction vs sim pu", 300, 0., 300., 0., 2.));
  add(h, new TH1F("ztpmatchnsimmatch", "sim matches per rec vtx", 20, 0., 20.));
  add(h, new TH1F("ztpmatchnrecmatch", "rec matches per sim vtx", 20, 0., 20.));

  add(h, new TProfile("zrandomeffvspu", "random point eff vs sim pu", 300, 0., 300., 0., 2.));
  add(h, new TProfile("zrandomfakevspu", "random point eff vs sim pu", 300, 0., 300., 0., 2.));

  // augmented z-matching
  add(h, new TProfile("FFAzmatcheffvspu", "eff vs sim pu", 300, 0., 300., 0., 2.));
  add(h, new TProfile("FFAzmatchfakevspu", "fake fraction vs sim pu", 300, 0., 300., 0., 2.));
  add(h, new TProfile("FFAztpmatcheffvspu", "eff vs sim pu", 300, 0., 300., 0., 2.));
  add(h, new TProfile("FFAztpmatchfakevspu", "fake fraction vs sim pu", 300, 0., 300., 0., 2.));
  add(h, new TProfile("FFAzcmatcheffvspu", "eff vs sim pu", 300, 0., 300., 0., 2.));
  add(h, new TProfile("FFAzcmatchfakevspu", "fake fraction vs sim pu", 300, 0., 300., 0., 2.));
  // including selection
  add(h, new TProfile("FFAzmatchseleffvspu", "eff vs sim pu", 300, 0., 300., 0., 2.));
  add(h, new TProfile("FFAzmatchselfakevspu", "fake fraction vs sim pu", 300, 0., 300., 0., 2.));
  add(h, new TProfile("FFAztpmatchseleffvspu", "eff vs sim pu", 300, 0., 300., 0., 2.));
  add(h, new TProfile("FFAztpmatchselfakevspu", "fake fraction vs sim pu", 300, 0., 300., 0., 2.));
  add(h, new TProfile("FFAzcmatchseleffvspu", "eff vs sim pu", 300, 0., 300., 0., 2.));
  add(h, new TProfile("FFAzcmatchselfakevspu", "fake fraction vs sim pu", 300, 0., 300., 0., 2.));

  // histograms of vertex properties, for fillVertexHistos(*)
  vector<string> vtypes = {"recvtx",
                           "selectedvtx",
                           "taggedvtx",
			   "matchedvtx",
			   "matchedvtxsel",
			   "fakevtxsel",
			   "fakevtx",
			   "splitvtxsel",
			   "otherfakevtxsel",
			   "splitvtxselfromsignal",
			   "splitvtxselfrompu"
  };
  unsigned int nvtype = vtypes.size();
  for (unsigned int t = 0; t < nvtype; t++) {
    string st = vtypes[t];
    dir->mkdir(vtypes[t].c_str())->cd();
    //addn :  add a histogram in a subdirectory and use the subdirectory name in the map key
    addn(h, new TH1F("c2xy", "c2xy", 100, 0., 10.));
    addn(h, new TProfile("c2xyvsntrk", "c2xy vs ntrk", 100, 0., 100., 0, 10.));
    addn(h, new TH1F("probxy", "probxy", 100, 0., 1.));
    addn(h, new TH1F("chi2", "chi**2", 100, 0., 100.));
    addn(h, new TH1F("chi2overntk", "chi**2", 100, 0., 4.));
    addn(h, new TH1F("xrecbeam", "xvtx-xbeam", 100, -0.025, 0.025));
    addn(h, new TH1F("yrecbeam", "yvtx-ybeam", 100, -0.025, 0.025));
    addn(h, new TH2F("xyrecbeam", "xyvtx-ybeam", 25, -0.025, 0.025, 25, -0.025, 0.025));
    addn(h, new TH1F("r", "r", 100, 0., 0.1));
    addn(h, new TH1F("index", "index", 200, 0., 200.));
    addn(h, new TH2F("logndofvsindex", "log ndof vs index", 200, 0., 200., 100, 0., 3.));
    addn(h, new TProfile("ndofvsindex", "mean ndof vs index", 200, 0., 200.,  0., 1e9));
    addn(h, new TH1F("zpullbeam", "z / sigmaz(Beam)", 200, -5., 5.));
    addn(h, new TH1F("ndof", "ndof", 500, 0., 100));
    addn(h, new TProfile("ndofvspu", "ndof vs pu", 300, 0., 300, 0., 1000));
    addn(h, new TH1F("numtrk", "number of tracks", 1000, 0., 1000.));  // log_10 : 1..1000.
    addn(h, new TH1F("logndof", "log ndof", 100, 0., 3.));  // log_10 : 1..1000.
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
    addn(h, new TH1F("logerrz", "log(z-error/um)", 100, 0, 5));

    addnSP(h, new TH1F("sumpt2", "sum pt**2", 200, 0., 1000));
    addnSP(h, new TH1F("logsumpt2", "log sum pt**2", 200, -1., 9.));
    addnSP(h, new TH1F("sumpt", "sum |pt|", 200, 0., 200));
    addnSP(h, new TH1F("logsumpt", "log sum |pt|", 60, -1., 5.));
    addnSP(h, new TH1F("ptmax2", "second highest pt", 100, 0., 10.));
    addn(h, new TH1F("zvtx", "z", 400, -20., 20.));
    addn(h, new TH1F("xvtx", "x", 400, -1., 1.));
    addn(h, new TH1F("yvtx", "y", 400, -1., 1.));

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
    addnSP(h, new TH1F("zrecsimHR","zrec - zsim", 400, -0.02, 0.02));
    addnSP(h, new TH1F("xrecsim","xrec - xsim", 100, -0.01, 0.01));
    addnSP(h, new TH1F("yrecsim","yrec - ysim", 100, -0.01, 0.01));
    addnSP(h, new TH1F("urecsim","urec - usim", 100, -0.01, 0.01));
    addnSP(h, new TH1F("vrecsim","vrec - vsim", 100, -0.01, 0.01));
    addnSP(h, new TH1F("xrecsimHR","xrec - xsim", 400, -0.01, 0.01));
    addnSP(h, new TH1F("yrecsimHR","yrec - ysim", 400, -0.01, 0.01));
    addnSP(h, new TH1F("urecsimHR","urec - usim", 400, -0.01, 0.01));
    addnSP(h, new TH1F("vrecsimHR","vrec - vsim", 400, -0.01, 0.01));
    addnSP(h, new TH1F("xrecerr","xrec uncertainty", 100, 0.0, 0.01));
    addnSP(h, new TH1F("yrecerr","yrec uncertainty", 100, 0.0, 0.01));
    addnSP(h, new TH1F("zrecsimpull","(zrec - zsim)/error", 100, -10, 10));
    addnSP(h, new TH1F("xrecsimpull","(xrec - xsim)/error", 100, -10, 10));
    addnSP(h, new TH1F("yrecsimpull","(yrec - ysim)/error", 100, -10, 10));
    addnSP(h, new TH1F("urecsimpull","(urec - usim)/error", 100, -10, 10));
    addnSP(h, new TH1F("vrecsimpull","(vrec - vsim)/error", 100, -10, 10));
    addnSP(h, new TH2F("xrecsimpullvsntrk","(xrec - xsim)/error", 20, 0, 100, 40, -10., 10));
    addnSP(h, new TH2F("yrecsimpullvsntrk","(yrec - ysim)/error", 20, 0, 100, 40, -10., 10));
    addnSP(h, new TH2F("zrecsimpullvsntrk","(zrec - zsim)/error", 20, 0, 100, 40, -10., 10));
    addnSP(h, new TH2F("urecsimpullvsntrk","(urec - usim)/error", 20, 0, 100, 40, -10., 10));
    addnSP(h, new TH2F("vrecsimpullvsntrk","(vrec - vsim)/error", 20, 0, 100, 40, -10., 10));
    addnSP(h, new TProfile("xrecsimvsxsimbeam","(xrec-xsim) vs (xsim-xbeam)", 50, -0.0050, 0.0050, -0.01, 0.01));
    addnSP(h, new TProfile("yrecsimvsysimbeam","(yrec-ysim) vs (ysim-ybeam)", 50, -0.0050, 0.0050, -0.01, 0.01));
    addnSP(h, new TProfile("xrecsimvsxsimbeamntklt5","(xrec-xsim) vs (xsim-xbeam)", 50, -0.0050, 0.0050, -0.01, 0.01));
    addnSP(h, new TProfile("yrecsimvsysimbeamntklt5","(yrec-ysim) vs (ysim-ybeam)", 50, -0.0050, 0.0050, -0.01, 0.01));
    addnSP(h, new TProfile("xrecsimvsxsimbeamntkge5","(xrec-xsim) vs (xsim-xbeam)", 50, -0.0050, 0.0050, -0.01, 0.01));
    addnSP(h, new TProfile("yrecsimvsysimbeamntkge5","(yrec-ysim) vs (ysim-ybeam)", 50, -0.0050, 0.0050, -0.01, 0.01));
    addnSP(h, new TH1F("nbsimtksinvtx", "simulated tracks in vertex", 200, -0.5, 199.5));
    addnSP(h, new TH1F("xsim", "simulated x", 100, -0.1, 0.1));
    addnSP(h, new TH1F("ysim", "simulated y", 100, -0.1, 0.1));
    addnSP(h, new TH1F("zsim", "simulated z", 120, -30., 30.));
    addnSP(h, new TH1F("dzminsim", "simulated min(dz)", 120, 0., 12.));
    addnSP(h, new TH1F("logpthatsim", "simulated log(pt-hat)", 120, 0., 10.));
    addnSP(h, new TH1F("logsumptsim", "simulated log(sum-pt)", 120, 0., 10.));
    addnSP(h, new TH1F("logsumpt2sim", "simulated log(sum-pt2)", 120, 0., 12.));

    // resolution vs sumpt2
    double log_pt2_bins[16] = {
      0.0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0, 1000.0, 2000.0, 5000.0, 10000.0};
    addnSP(h, new TProfile("xrecsimvssumptprof","x resolution vs sumpt2", 15, log_pt2_bins, -1., 1.));
    addnSP(h, new TH2F("xrecsimvssumpt","x resolution vs sumpt2", 15, log_pt2_bins, 400, -0.05, 0.05));
    addnSP(h, new TProfile("xrecsimvssumptsimprof","x resolution vs sumpt2", 15, log_pt2_bins, -1., 1.));
    addnSP(h, new TH2F("xrecsimvssumptsim","x resolution vs sumpt2", 15, log_pt2_bins, 400, -0.05, 0.05));

    addnSP(h, new TProfile("yrecsimvssumptprof","y resolution vs sumpt2", 15, log_pt2_bins, -1., 1.));
    addnSP(h, new TH2F("yrecsimvssumpt","y resolution vs sumpt2", 15, log_pt2_bins, 400, -0.05, 0.05));
    addnSP(h, new TProfile("yrecsimvssumptsimprof","y resolution vs sumpt2", 15, log_pt2_bins, -1., 1.));
    addnSP(h, new TH2F("yrecsimvssumptsim","y resolution vs sumpt2", 15, log_pt2_bins, 400, -0.05, 0.05));

    addnSP(h, new TProfile("zrecsimvssumptprof","z resolution vs sumpt2", 15, log_pt2_bins, -1., 1.));
    addnSP(h, new TH2F("zrecsimvssumpt","z resolution vs sumpt2", 15, log_pt2_bins, 400, -0.05, 0.05));
    addnSP(h, new TProfile("zrecsimvssumptsimprof","z resolution vs sumpt2", 15, log_pt2_bins, -1., 1.));
    addnSP(h, new TH2F("zrecsimvssumptsim","z resolution vs sumpt2", 15, log_pt2_bins, 400, -0.05, 0.05));
    
    double log_ntrk_bins[25] = {0.,   2.0,  4.0,  6.0,  8.0,  10.,  12.0, 14.0, 16.0, 18.0,  22.0,  26.0, 30.0,
                             35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 70.0, 80.0, 90.0, 100.0, 150.0, 200.0};
    addnSP(h, new TProfile("xrecsimvsntrkprof","x resolution vs number of vertex tracks", 24, log_ntrk_bins, -1., 1.));
    addnSP(h, new TH2F("xrecsimvsntrk","x resolution vs number of vertex tracks", 24, log_ntrk_bins, 400, -0.05, 0.05));
    addnSP(h, new TProfile("xrecsimvsnsimtrkprof","x resolution vs number of vertex tracks", 24, log_ntrk_bins, -1., 1.));
    addnSP(h, new TH2F("xrecsimvsnsimtrk","z resolution vs number of vertex tracks", 24, log_ntrk_bins, 400, -0.05, 0.05));

    addnSP(h, new TProfile("yrecsimvsntrkprof","y resolution vs number of vertex tracks", 24, log_ntrk_bins, -1., 1.));
    addnSP(h, new TH2F("yrecsimvsntrk","y resolution vs number of vertex tracks", 24, log_ntrk_bins, 400, -0.05, 0.05));
    addnSP(h, new TProfile("yrecsimvsnsimtrkprof","y resolution vs number of vertex tracks", 24, log_ntrk_bins, -1., 1.));
    addnSP(h, new TH2F("yrecsimvsnsimtrk","z resolution vs number of vertex tracks", 24, log_ntrk_bins, 400, -0.05, 0.05));

    addnSP(h, new TProfile("zrecsimvsntrkprof","z resolution vs number of vertex tracks", 24, log_ntrk_bins, -1., 1.));
    addnSP(h, new TH2F("zrecsimvsntrk","z resolution vs number of vertex tracks", 24, log_ntrk_bins, 400, -0.05, 0.05));
    addnSP(h, new TProfile("zrecsimvsnsimtrkprof","z resolution vs number of vertex tracks", 24, log_ntrk_bins, -1., 1.));
    addnSP(h, new TH2F("zrecsimvsnsimtrk","z resolution vs number of vertex tracks", 24, log_ntrk_bins, 400, -0.05, 0.05));
 
    if (f4D_) {
      // for matched vertices with timing
      addnSP(h, new TH1F("tsim", "simulated t", 400, -2., 2.));
      addnSP(h, new TH1F("ntimingvtx", "#timing tracks", 200, 0., 200.));
      addnSP(h, new TH1F("ntimingqual05vtx", "#timing tracks above quality threshold 0.5", 200, 0., 200.));
      addnSP(h, new TH1F("ntimingqual08vtx", "#timing tracks above quality threshold 0.8", 200, 0., 200.));

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

      // for the various vertex_time_fromtracks versions
      vector<string> suffixes = {"","_qual","_pid","_qual_pid", "_qual_pid_new"};
      for(auto & suffix : suffixes){
	addnSP(h, new TH1F(("trecsim_fromtracks" + suffix).c_str(), ("vertex time residual from tracks "+ suffix).c_str(), 200, -0.1, 0.1));
	addnSP(h, new TH1F(("trecerr_fromtracks" + suffix).c_str(), ("vertex time error from tracks"+ suffix).c_str(), 500, 0., 0.1));
	addnSP(h, new TH1F(("trecsimpull_fromtracks" + suffix).c_str(), ("tpull from tracks" + suffix).c_str(), 200, -10., 10.));
	//references, vertex time for the same set of vertices
	addnSP(h, new TH1F(("trecsim_withtracks" + suffix).c_str(), ("vertex time residual from vertex " + suffix).c_str(), 200, -0.1, 0.1)); 
	addnSP(h, new TH1F(("trecerr_withtracks" + suffix).c_str(), ("vertex time error from vertex" + suffix).c_str(), 500, 0., 0.1));
	addnSP(h, new TH1F(("trecsimpull_withtracks" + suffix).c_str(), ("tpull from vertex" + suffix).c_str(), 200, -10., 10.));
      }
    }
    dir->cd();
  }// vtypes

  //  histograms of track quality for fillTrackHistos and fillVertexHistos
  //  note that track histograms are automatically also created and filled for all vtypes !
  vector<string> ttypes = {"trkall",
                           "trksel",
                           "trkseltiming",
                           "sellost",
                           "selused",
                           "trkselhipt",
                           "trkselhiptfwd",
                           "trkselhiptcentral",
			   "trksellopt",
			   "trkselloptfwd",
			   "trkselloptcentral",
                           "unmatchedVtx",
                           "seltpmatched",
			   "seltpmatchedhipt",
                           "seltpmatchedSignal",
                           "seltpmatchedPU",
                           "seltpunmatched",
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
    addn(h, new TH1F("phi", "phi", 320, -3.14159, 3.14159));
    addn(h, new TH1F("eta", "eta ", netabin, -etarange, etarange));
    addn(h, new TH1F("pt", "pt ", 100, 0., 5.));
    addn(h, new TH1F("logpt", "logpt ", 60, -1., 5.));
    addn(h, new TH2F("logpteta", "logpt vs eta ", netabin, -etarange, etarange, 60, -1., 5.));
    addn(h, new TH2F("z-eta", "z vs eta ", netabin, -etarange, etarange, 60, -15., 15.));
    addn(h, new TH2F("z-pt", "pt vs z", 60, -15., 15., 100, 0., 5.));
    addn(h, new TH2F("z-logpt", "log pt vs z", 60, -15., 15., 60, -1., 5.));
    addn(h, new TH1F("ptfwd", "pt (forward)", 100, 0., 5.));
    addn(h, new TH1F("ptcentral", "pt (central)", 100, 0., 5.));
    addn(h, new TH1F("found", "found hits", 40, 0., 40.));
    addn(h, new TH1F("lost", "lost hits", 20, 0., 20.));
    addn(h, new TH1F("validfraction", "fraction of valid hits", 50, 0., 1.));
    addn(h, new TH1F("ndoftrk", "track fit ndof", 100, 0., 100.));
    addn(h, new TH1F("chi2trk", "track chi2", 200, 0., 100.));
    addn(h, new TH1F("nchi2trk", "normalized track chi2", 100, 0., 20.));
    addn(h, new TProfile("nchi2trkvsz", "normalized track chi2", 120, -30., 30., 0., 20.));
    addn(h, new TH1F("rstart", "start radius", 100, 0., 20.));
    addn(h, new TH1F("expectedInner", "expected inner hits ", 10, 0., 10.));
    addn(h, new TH1F("expectedOuter", "expected outer hits ", 10, 0., 10.));
    addn(h, new TH1F("logtresxy", "log10(track r-phi resolution/um)", 100, 0., 5.));
    addn(h, new TH1F("logtresz", "log10(track z resolution/um)", 100, 0., 5.));
    addn(h, new TH1F("tpullxy", "track r-phi pull", 100, -10., 10.));
    addn(h, new TProfile("tpullxyvsz", "track r-phi pull", 120, -30., 30., -10., 10.));
    addn(h, new TProfile("tpullxyvseta", "track r-phi pull", netabin, -etarange, etarange, -10., 10.));
    addn(h, new TProfile("tpullzvsz", "track z pull", 120, -30., 30., -10., 10.));
    addn(h, new TProfile("tpullzvseta", "track z pull", netabin, -etarange, etarange, -10., 10.));
    addn(h, new TH2F("lvseta", "cluster length vs eta", netabin, -etarange, etarange, 20, 0., 20));
    addn(h, new TH2F("lvstanlambda", "cluster length vs tan lambda", 60, -6., 6., 20, 0., 20));
    addn(h, new TH1F("longestbarrelhit", "longest barrel cluster", 40, 0., 40.));
    // note : histograms involving the relation wrt a rec vertex are usually filled with the tagged vertex
    // for pile-up this is usually the wrong one, expect a peak on top of a more or less flat backgound
    // useful mostly for MC without PU
    addn(h, new TH1D("zrestrk", "z-residuals (track vs vertex)", 200, -2., 2.));
    addn(h,
        new TH2F("zrestrkvsphi", "z-residuals (track - vertex)", 12, -3.14159, 3.14159, 100, -1., 1.));
    addn(h, new TH2F("zrestrkvseta", "z-residuals (track - vertex)", 16, -etarange, etarange, 100, -1., 1.));
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
            "zpulltrkvseta", "normalized z-residuals (track - vertex)", 12, -etarange, etarange, 100, -5., 5.));
    addn(h,
        new TH2F(
            "zpulltrkvsz", "normalized z-residuals (track - vertex) vs z", 100, -20., 20., 100, -5., 5.));
    addn(h, new TH1D("zpulltrk", "normalized z-residuals (track vs vertex)", 100, -5., 5.));
    addn(h, new TH1D("zerrtrk", "z-resolution (excluding beam)", 500, 0., 0.5));
    addn(h, new TH1D("logzerrtrk", "log(z-resolution)", 100, 0., 5.));
    addn(h, new TH1D("zerrtrk_withbeam", "z-resolution (including beam)", 500, 0., 0.5));
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

    // fillTrackHistosMatched
    //  - tracks with MC truth (matched to TrackingParticles)
    const double pull_cut_off = 10.; // for profile histograms
    const double zrecsim_cut_off = 2.; 
    addnSP(h, new TH1D("tkzrecsim", "zrec-zsim", 200, -0.1, 0.1));
    addnSP(h, new TH1D("tkzrecsimpull", "(zrec-zsim) / #sigma_{z}", 300, -15., 15.));
    addnSP(h, new TProfile("tkzrecsimpullsqvseta", "zrec-zsim pull**2 vs eta", netabin, -etarange, etarange, 0., pull_cut_off * pull_cut_off));
    addnSP(h, new TProfile("tkzrecsimsqvseta", "(zrec-zsim)**2 vs eta", netabin, -etarange, etarange, 0., zrecsim_cut_off*zrecsim_cut_off));
    addnSP(h, new TProfile("tkzrecsimsqvsetahipt", "(zrec-zsim)**2 vs eta (high pt)", netabin, -etarange, etarange, 0., zrecsim_cut_off*zrecsim_cut_off));
    addnSP(h, new TH2F("tkzrecsimpullvseta", "zrec-zsim pull vs eta", netabin, -etarange, etarange, 200, -20., 20.));
    addnSP(h, new TProfile("tkzrecsimvsz",   "zrec-zsim vs z", 100, -15., 15., -zrecsim_cut_off, zrecsim_cut_off));
    addnSP(h, new TProfile("tkzrecsimvseta", "zrec-zsim vs eta", netabin, -etarange, etarange, -zrecsim_cut_off, zrecsim_cut_off));
    addnSP(h, new TProfile("tkzrecsimvsetaz", "zrec-zsim vs eta with z-flip", netabin, -etarange, etarange, -zrecsim_cut_off, zrecsim_cut_off));
    addnSP(h, new TH2F("tkzrecsimvseta2d", "zrec-zsim vs eta", netabin, -etarange, etarange, 400, -2., 2.));
    addnSP(h, new TH2F("tkzrecsimvseta2dhipt", "zrec-zsim vs eta (high pt)", netabin, -etarange, etarange, 400, -0.4, 0.4));

    addnSP(h, new TH2F("tkptrecsimrelvssimeta2d", "(pTrec-pTsim)/pTrec vs eta_{TP}", netabin, -etarange, etarange, 300, -0.3, 0.3));
    addnSP(h, new TH2F("tkptrecsimpullvssimeta2d", "pTrec-pTsim pull vs eta_{TP}", netabin, -etarange, etarange, 100, -15., 15.));
    addnSP(h, new TH2F("tkdptrelrecvssimeta2d", "#sigma(pTrec)/pTrec vs eta_{TP}", netabin, -etarange, etarange, 100, -3., 3.));
    addnSP(h, new TH2F("tketarecsimvssimeta2d", "etarec-etasim vs eta_{TP}", netabin, -etarange, etarange, 100, -0.05, 0.05));
    addnSP(h, new TH2F("tkphirecsimvssimeta2d", "phirec-phisim vs eta_{TP}", netabin, -etarange, etarange, 100, -0.05, 0.05));
    addnSP(h, new TH2F("tkzrecsimvssimeta2d", "zrec-zsim vs eta_{TP}", netabin, -etarange, etarange, 200, -0.2, 0.2));
    addnSP(h, new TH2F("tkzrecsimpullvssimeta2d" , "zrec-zsim pull vs eta_{TP}", netabin, -etarange, etarange, 300, -15., 15.));
    addnSP(h, new TH2F("tkdzrecvssimeta2d", "zError vs eta_{TP}", netabin, -etarange, etarange, 200, -0.2, 0.2));


    addnSP(h, new TH1D("tktrecsim", "trec-tsim", 200, -0.3, 0.3));
    addnSP(h, new TH1D("tktrecsimpull", "(trec-tsim)/terr", 200, -10., 10.));
    addnSP(h, new TH1D("tktrecsim_pid", "trec-tsim (pid)", 200, -0.3, 0.3));
    addnSP(h, new TH1D("tktrecsimpull_pid", "(trec-tsim)/terr (pid)", 200, -10., 10.));
    addnSP(h, new TH1D("tktrecsimpullwide", "trec-tsim", 200, -20., 20.));
    addnSP(h, new TH2F("tktrecsimvseta2d", "trec-tsim vs eta", netabin, -etarange, etarange, 200, -1.0, 1.0));
    addnSP(h, new TProfile("tktrecsimpullsqvserr", "((trec-tsim)/terr)**2 vs terr", 50, 0., 0.5, 0., 100.0));
    addnSP(h, new TH2F("tktrecsimpullvserr", "(trec-tsim)/terr vs terr", 50, 0., 0.5, 50, 0., 100.0));
    addnSP(h, new TProfile("tkzrecsimvslogpt",   "(zrec-zsim) vs log pt", 100, -1., 5., -pull_cut_off, pull_cut_off));

    addnSP(h, new TH2F("tkpidvsetalogpt", "all particle types", netabin/2, 0, etarange, 30, -1., 2.));
    addnSP(h, new TH2F("tkpidpionvsetalogpt", "pions", netabin/2, 0, etarange, 30, -1., 2.));
    addnSP(h, new TH2F("tkpidkaonvsetalogpt", "kaons", netabin/2, 0, etarange, 30, -1., 2.));
    addnSP(h, new TH2F("tkpidprotonvsetalogpt", "protons", netabin/2, 0, etarange, 30, -1., 2.));
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
  add(h, new TH1F("dzrecsimtaggedsignalnarrowovfl", "z(tagged,rec) - z(signal,sim)", 200, 0, 0.1));
  add(h, new TH1F("dzrecsimtaggedsignalnarrow", "z(tagged,rec) - z(signal,sim)", 200, -0.1, 0.1));
  add(h, new TH1F("dzrecsimtaggedsignalwide", "z(tagged,rec) - z(signal,sim)", 200, -5., 5.));
  add(h, new TProfile("tag01mmvspu","tagged vertex within 0.1 mm of signal vs pu",vmax, 0., vmax, 0., 1.));
  add(h, new TProfile("tag02mmvspu","tagged vertex within 0.2 mm of signal vs pu",vmax, 0., vmax, 0., 1.));
  add(h, new TProfile("tag05mmvspu","tagged vertex within 0.5 mm of signal vs pu",vmax, 0., vmax, 0., 1.));
  add(h, new TProfile("tag1mmvspu","tagged vertex within 1 mm of signal vs pu",vmax, 0., vmax, 0., 1.));
  add(h, new TProfile("tag2mmvspu","tagged vertex within 2 mm of signal vs pu",vmax, 0., vmax, 0., 1.));
  add(h, new TProfile("tag5mmvspu","tagged vertex within 5 mm of signal vs pu",vmax, 0., vmax, 0., 1.));
  add(h, new TH2F("indexvsdistancerank","rec index vs signal distance index", 10, 0., 10., 10, 0., 10.));
  const std::vector<std::string>t = {"_tagged","_nottagged"};
  for(auto & tag : t){
    add(h, new TH2F(("nrecvsnsimwithin5mm"+tag).c_str(), "#rec vs #sim within 5mm of the signal vertex", 10, 0, 10., 10, 0., 10.));
    add(h, new TH2F(("nrecvsnsimwithin2mm"+tag).c_str(), "#rec vs #sim within 2mm of the signal vertex", 10, 0, 10., 10, 0., 10.));
    add(h, new TH2F(("nrecvsnsimwithin1mm"+tag).c_str(), "#rec vs #sim within 1mm of the signal vertex", 10, 0, 10., 10, 0., 10.));
    add(h, new TH2F(("nrecvsnsimwithin05mm"+tag).c_str(), "#rec vs #sim within 0.5mm of the signal vertex", 10, 0, 10., 10, 0., 10.));
    add(h, new TH2F(("nrecvsnsimwithin02mm"+tag).c_str(), "#rec vs #sim within 0.2mm of the signal vertex", 10, 0, 10., 10, 0., 10.));
    add(h, new TH2F(("nrecvsnsimwithin01mm"+tag).c_str(), "#rec vs #sim within 0.1mm of the signal vertex", 10, 0, 10., 10, 0., 10.));
    add(h, new TH1F(("nrecwithin5mm"+tag).c_str(), "#rec within 5 mm of the signal vertex", 50, 0, 50));
    add(h, new TH1F(("nsimwithin5mm"+tag).c_str(), "#sim within 5 mm of the signal vertex", 50, 0, 50));
    add(h, new TH1F(("nrecwithin1mm"+tag).c_str(), "#rec within 1 mm of the signal vertex", 50, 0, 50));
    add(h, new TH1F(("nsimwithin1mm"+tag).c_str(), "#sim within 1 mm of the signal vertex", 50, 0, 50));
    add(h, new TH2F(("indexvsdistancerank"+tag).c_str(),"rec index vs signal distance index", 10, 0., 10., 10, 0., 10.));
  }
  add(h, new TH1F("dzrecsim", "distance of sim vertex", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("dzrecsimHR", "distance of sim vertex", nbinzdiffrec, -0.2, 0.2));
  add(h, new TH1F("dzrecsimSignal", "distance between signal and rec vertices", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("dzrecsimSignalHR", "distance between signal and rec vertices", nbinzdiffrec, -0.2, 0.2));
  add(h, new TH1F("dzrecsimSignalVHR", "distance between signal and rec vertices", nbinzdiffrec, -0.05, 0.05));
  /* do we care?
  add(h, new TH1F("dzrecsimPU", "distance between PU and rec vertices", nbinzdiffrec, -2., 2.));
  add(h, new TH1F("dzrecsimPUHR", "distance between PU and rec vertices", nbinzdiffrec, -0.2, 0.2));
  add(h, new TH1F("dzrecsimPUVHR", "distance between PU and rec vertices", nbinzdiffrec, -0.05, 0.05));
  */
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
  add(h, new TH1F("sumpt2", "sumpt2 of simulated tracks", 100, 0., 100.));
  add(h, new TH1F("sumpt2Signal", "sumpt2 of simulated tracks in Signal events", 100, 0., 200.));
  add(h, new TH1F("sumpt2PU", "sumpt2 of simulated tracks in PU events", 100, 0., 200.));
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
  add(h, new TH1F("matchVtxFractionCumSignal", "fraction of sim vertexes track found in a recvertex", 101, 0, 1.01));
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
  add(h, new TH2F("correctlyassigned", "pt and eta of correctly assigned tracks", netabin, -etarange, etarange, 100, 0, 10.));
  add(h, new TH2F("misassigned", "pt and eta of mis assigned tracks", netabin, -etarange, etarange, 100, 0, 10.));
  add(h, new TH2F("correctlyassignedlogpt", "log10 pt and eta of correctly assigned tracks", netabin, -etarange, etarange, 40, -1., 3.));
  add(h, new TH2F("misassignedlogpt", "log10 pt and eta of mis assigned tracks", netabin, -etarange, etarange, 40, -1., 3.));

  add(h, new TH1F("ptcat", "pt of correctly assigned tracks", 100, 0, 10.));
  add(h, new TH1F("etacat", "eta of correctly assigned tracks", netabin, -etarange, etarange));
  add(h, new TH1F("etacatpt2", "eta of correctly assigned tracks pt>2GeV", netabin, -etarange, etarange));
  add(h, new TH1F("phicat", "phi of correctly assigned tracks", 100, -3.14159, 3.14159));
  add(h, new TH1F("dzcat", "dz of correctly assigned tracks", 100, 0., 1.));

  add(h, new TH1F("ptmis", "pt of mis-assigned tracks", 100, 0, 10.));
  add(h, new TH1F("etamis", "eta of mis-assigned tracks", netabin, -etarange, etarange));
  add(h, new TH1F("etamispt2", "eta mis-correctly assigned tracks pt>2GeV", netabin, -etarange, etarange));
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
  return h;
}


void PrimaryVertexAnalyzer4PU::bookTrackHistograms(const char * directory_name)
// book histograms for truth-matched tracks,
// filled in analyzeTracksTP
{
  const float etarange = 4.;
  const int netabin = 40.;

  rootFile_->cd();
  TDirectory* dir = rootFile_->mkdir(directory_name);
  dir->cd();

  // this prim-only group may not be needed
  add(hTrk, new TH1F("zpulltrksec", "reconstructed z- generated z / error for non-primary tracks", 200, -10., 10.));

  // these are filled
  add(hTrk, new TH1F("zpulltrk_primselmatched", "reconstructed z- generated z for primary tracks", 200, -10., 10.));
  add(hTrk,
      new TH1F(
          "zpulltrkt_primselmatched", "reconstructed z- generated z for primary tracks with timing", 200, -10., 10.));
  add(hTrk, new TH1F("zpulltrkprimsel", "reconstructed z- generated z for primary tracks with timing", 200, -10., 10.));
  add(hTrk, new TH1F("zrestrk_primselmatched", "reconstructed z- generated z for primary tracks", 200, -0.2, 0.2));
  add(hTrk, new TH2F("zpulltrkprimselvseta", "", netabin, -etarange, etarange, 200, -10., 10.));


  //>>>>>>>>>>>>>>>>>
  for(auto bin = 0u; bin < trkdzbin_.size(); bin++){
    add(hTrk, new TH2F(Form("zpulltrkprimselvseta_%s", trkdzbin_[bin].c_str()), "", netabin, -etarange, etarange, 200, -10., 10.));
    add(hTrk, new TH2F(Form("zpulltrkprimselbpxlt2vseta_%s", trkdzbin_[bin].c_str()), "", netabin, -etarange, etarange, 200, -10., 10.));
    add(hTrk, new TH2F(Form("zpulltrkprimselbpxgt2vseta_%s", trkdzbin_[bin].c_str()), "", netabin, -etarange, etarange, 200, -10., 10.));
    add(hTrk, new TH2F(Form("zpulltrkprimselvslogpt_%s", trkdzbin_[bin].c_str()), "", 40, -1., 3., 200, -10., 10.));
    add(hTrk, new TProfile(Form("ztailtprimselvslogpt_%s",trkdzbin_[bin].c_str()) , "", 40, -1., 3., 0., 2.));
    add(hTrk, new TProfile(Form("ztailtprimselvslogpt_etahi_%s",trkdzbin_[bin].c_str()) , "", 40, -1., 3., 0., 2.));
    add(hTrk, new TProfile(Form("ztailtprimselvslogpt_etalo_%s",trkdzbin_[bin].c_str()) , "", 40, -1., 3., 0., 2.));
    add(hTrk, new TH2F(Form("ztailtprimselvslogpteta_%s", trkdzbin_[bin].c_str()), "", netabin, -etarange, etarange, 40, -1., 3.));
    add(hTrk, new TH2F(Form("tprimselvslogpteta_%s", trkdzbin_[bin].c_str()), "", netabin, -etarange, etarange, 40, -1., 3.));
  }

  add(hTrk, new TH2F("zpulltrkprimselvslogpt", "", 40, -1., 3., 200, -10., 10.));

  add(hTrk, new TH2F("ztailtprimselvslogpteta", "", netabin, -etarange, etarange, 40, -1., 3.));
  add(hTrk, new TH2F("tprimselvslogpteta", "", netabin, -etarange, etarange, 40, -1., 3.));
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

  // FIXME, add pt bins
  add(hTrk, new TH1F("etaprim", "eta of primary tracks", 100, -4., 4.));
  add(hTrk, new TH1F("etaprim_Pt000to001", "eta of primary tracks pt < 1", 100, -4., 4.));
  add(hTrk, new TH1F("etaprim_Pt001to003", "eta of primary tracks 1 < pt < 3", 100, -4., 4.));
  add(hTrk, new TH1F("etaprim_Pt003to010", "eta of primary tracks 3 < pt < 10", 100, -4., 4.));
  add(hTrk, new TH1F("etaprimsel", "eta of primary selected tracks", 100, -4., 4.));
  add(hTrk, new TH1F("etaprimsel_Pt000to001", "eta of primary selected tracks pt < 1", 100, -4., 4.));
  add(hTrk, new TH1F("etaprimsel_Pt001to003", "eta of primary selected tracks 1 < pt < 3", 100, -4., 4.));
  add(hTrk, new TH1F("etaprimsel_Pt003to010", "eta of primary selected tracks 3 < pt < 10", 100, -4., 4.));

  add(hTrk, new TH1F("d0pullprim", "d0/error for primary tracks", 200, -10., 10.));
  add(hTrk, new TH1F("d0pullprimsel", "d0/error for primary selected tracks", 200, -10., 10.));
  add(hTrk, new TH1F("d0pullsec", "d0/error for non-primary tracks", 200, -10., 10.));
  add(hTrk, new TH2F("zpullvsd0pullprim", "z pull vs d0-pull for primary tracks", 100, 0., 10., 100, 0., 10.));
  add(hTrk, new TH2F("zpullvsd0pullprimsel", "z pull vs d0-pull for primary selected tracks", 100, 0., 10., 100, 0., 10.));
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
  addSP(hsimPV, new TH1F("xsim", "simulated x", 100, -0.1, 0.1));
  addSP(hsimPV, new TH1F("ysim", "simulated y", 100, -0.1, 0.1));
  addSP(hsimPV, new TH1F("zsim", "simulated z", 120, -30., 30.));
  addSP(hsimPV, new TH1F("zvtx_sim", "simulated z", 400, -20., 20.)); // binning compatible with zvtx_rec at al
  addSP(hsimPV, new TH1F("dzminsim", "simulated min(dz)", 120, 0., 12.));
  addSP(hsimPV, new TH1F("xsimb", "simulated x", 100, -0.01, 0.01));  // 0.01cm = 100 um
  addSP(hsimPV, new TH1F("ysimb", "simulated y", 100, -0.01, 0.01));
  addSP(hsimPV, new TH1F("zsimb", "simulated z", 120, -30., 30.));
  if (f4D_) addSP(hsimPV, new TH1F("tsim", "simulated t", 400, -2., 2.));

  addSP(hsimPV, new TH1F("xbeamspot", "beamspot x", 100, -1., 1.));
  addSP(hsimPV, new TH1F("ybeamspot", "beamspot y", 100, -1., 1.));
  addSP(hsimPV, new TH1F("zbeamspot", "beamspot z", 100, -5., 5));
  addSP(hsimPV, new TH1F("wxbeamspot", "beamspot sigma x", 100, 0., 0.02));
  addSP(hsimPV, new TH1F("wybeamspot", "beamspot sigma y", 100, 0., 0.02));
  addSP(hsimPV, new TH1F("sigmaZbeamspot", "beamspot sigma z", 100, 0., 10.));
  addSP(hsimPV, new TH1F("nsimvtx", "# of simulated vertices", 250, 0., 250.));
  addSP(hsimPV, new TH1F("nbsimtksinvtx", "simulated tracks in vertex", 200, -0.5, 199.5));

  addSP(hsimPV, new TH1F("logpthatsim", "simulated log(pt-hat)", 120, 0., 10.));
  addSP(hsimPV, new TH1F("logsumptsim", "simulated log(sum-pt)", 120, 0., 10.));
  addSP(hsimPV, new TH1F("logsumpt2sim", "simulated log(sum-pt2)", 120, 0., 12.));
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
    //iSetup.getData(pdt_);
    pdt_ = iSetup.getHandle(pdtToken_);
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

    Fill(hsimPV, "xbeamspot", vertexBeamSpot_.x0());
    Fill(hsimPV, "wxbeamspot", vertexBeamSpot_.BeamWidthX());
    Fill(hsimPV, "ybeamspot", vertexBeamSpot_.y0());
    Fill(hsimPV, "wybeamspot", vertexBeamSpot_.BeamWidthY());
    Fill(hsimPV, "zbeamspot", vertexBeamSpot_.z0());
    Fill(hsimPV, "sigmaZbeamspot", vertexBeamSpot_.sigmaZ());
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

/********************************************************************************************************/
bool PrimaryVertexAnalyzer4PU::get_miniaod_tracks(const edm::EventSetup& iSetup,
                                                  const edm::Event& iEvent,
						  const std::string &miniaod_vertexcollection_label,
						  Tracks& tracks) 
/********************************************************************************************************/
{   // 
   // more info https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2017
   //
   // code from Kirill
   Handle<edm::View<pat::PackedCandidate> > tracksPackedHandle;
   iEvent.getByToken(theTracksToken_, tracksPackedHandle);

   Handle<edm::View<pat::PackedCandidate> > lostTracksPackedHandle;
   iEvent.getByToken(theLostTracksToken_, lostTracksPackedHandle);

   
   // Create pseudo-track collection
   edm::View<pat::PackedCandidate> tracksPacked = (*tracksPackedHandle.product());
   edm::View<pat::PackedCandidate> lostTracksPacked = (*lostTracksPackedHandle.product());

   if(verbose_){
     cout << "get_miniaod_tracks  found " << tracksPacked.size() << " packed tracks and  " << lostTracksPacked.size() << " lost tracks"<<endl;
   }
   
   // loop over both types in one go

   tracks.reserve(tracksPacked.size() + lostTracksPacked.size());  // without this the track pointers become invalid

   for(size_t it = 0; it < tracksPacked.size() + lostTracksPacked.size(); it++){
     const pat::PackedCandidate &trkPacked = it < tracksPacked.size() ? tracksPacked[it] : lostTracksPacked[it - tracksPacked.size()];
     if (trkPacked.charge() == 0.) continue; // this is a neutral, not a track
     auto tk = MTrack(tracks.size(), trkPacked, it, theB_, vertexBeamSpot_, miniaod_vertexcollection_label, false);  // false = no
     if(tk.has_transienttrack()){
       tk._selected = theTrackFilter(tk.transientTrack());
     }
     tracks.push_back(tk);
   }

   return true;
}
/********************************************************************************************************/



   
/********************************************************************************************************/
bool PrimaryVertexAnalyzer4PU::get_reco_and_transient_tracks(const edm::EventSetup& iSetup,
                                                             const edm::Event& iEvent,
                                                             Tracks& tracks) {
/********************************************************************************************************/
  // requires beamspot
  if (!iEvent.getByToken(edmView_recoTrack_Token_, tracks.trackCollectionH)) {
    if (verbose_) {
      report_counted( " no reco tracks found, bailing out", 10);
    }
    return false;
  }

  const View<reco::Track>* recTrks = tracks.trackCollectionH.product();

  edm::Handle<edm::ValueMap<float>> trackTimesH;
  edm::Handle<edm::ValueMap<float>> trackTimeResosH;

  bool trktime = f4D_
    && iEvent.getByToken(trkTimesToken_, trackTimesH)
    && iEvent.getByToken(trkTimeResosToken_, trackTimeResosH);
  

  std::vector<reco::TransientTrack> t_tks;
  if (trktime) {
    t_tks = (*theB_).build(tracks.trackCollectionH, vertexBeamSpot_, *(trackTimesH.product()), *(trackTimeResosH.product()));
  }else { 
    t_tks = (*theB_).build(tracks.trackCollectionH, vertexBeamSpot_);
  }


  if (trktime){
    report_counted("PrimaryVertexAnalyzer4PU::get_reco_and_transient_tracks timing found",1);
    //cout << "pathlength " << have_MTD_pathlength << "  mtdtime " << have_MTD_time << "   mtdtimerror" << have_MTD_timeerror << "   mtdmomentum " << have_MTD_momentum << endl;
  }else{
    report_counted("PrimaryVertexAnalyzer4PU::get_reco_and_transient_tracks timing requested but not found",1);
  }

  bool have_timequality = false;
  bool have_MTD_pathlength = false;
  bool have_MTD_time = false;
  bool have_MTD_timeerror = false;
  bool have_MTD_momentum = false;
  edm::Handle<edm::ValueMap<float>> MTD_pathlength_H, MTD_time_H, MTD_timeerror_H, MTD_momentum_H;
  edm::Handle<edm::ValueMap<float>> trackTimeQualityH;

  if (f4D_){
    have_timequality = iEvent.getByToken(trkTimeQualityToken_, trackTimeQualityH);
    have_MTD_pathlength = iEvent.getByToken(MTD_pathlength_Token_, MTD_pathlength_H);
    have_MTD_time = iEvent.getByToken(MTD_time_Token_, MTD_time_H);
    have_MTD_timeerror = iEvent.getByToken(MTD_timeerror_Token_, MTD_timeerror_H);
    have_MTD_momentum = iEvent.getByToken(MTD_momentum_Token_, MTD_momentum_H);
  }

  /* fill private track container and make it globally available */
  for (View<reco::Track>::size_type i = 0; i < recTrks->size(); ++i) {

    const reco::TransientTrack* tt = &(t_tks[i]);
    unsigned int key1 = t_tks[i].trackBaseRef().key();
    RefToBase<reco::Track> trb(tracks.trackCollectionH, i);
    unsigned int key2 = trb.key();
    // paranoia is in bloom
    if (key1 != key2) {
      cout << "get_reco_and_transient_tracks : key confusion" << endl;
    }

    auto tk = MTrack(i, &(recTrks->at(i)), t_tks[i], key1, trktime);

    if(trktime && have_timequality){
      tk._timeQuality = (*trackTimeQualityH)[tt->trackBaseRef()];
      if( tk._timeQuality < 0.001){
	tk._has_timing = false;
      }
    }
    
    if(have_MTD_time){
      tk._MTD_time = (*MTD_time_H)[tt->trackBaseRef()];
      //cout << "MTD time = " << tk.MTD_time << ", has_timing=" << tk.has_timing << endl;
    }
    if(have_MTD_timeerror){
      tk._MTD_timeerror = (*MTD_timeerror_H)[tt->trackBaseRef()];
      //cout << "MTD time error= " << tk.MTD_timeerror << "    track time error = " << tk.dt <<  endl;
    }else{
      tk._MTD_timeerror = 1.e10;
    }

    if(have_MTD_momentum){
      tk._MTD_momentum = (*MTD_momentum_H)[tt->trackBaseRef()];
      //cout << "MTD momentum = " << tk.MTD_momentum << "  vs " << tk.pt / std::abs(tan(tk.theta)) << endl;
    }else{
      tk._MTD_momentum = tk.pt() / std::abs(tan(tk.theta()));
    }
    if(have_MTD_pathlength){
      tk._MTD_pathlength = (*MTD_pathlength_H)[tt->trackBaseRef()];
      //cout << "MTD pathlength = " << tk.MTD_pathlength << endl;
    }else{
      tk._MTD_pathlength = -1.;
    }

    if (have_MTD_time && have_MTD_momentum && have_MTD_pathlength && (tk.MTD_pathlength() > 0) ){
      tk.th[0] = tk.get_t_pid( 0.139570 );
      tk.th[1] = tk.get_t_pid( 0.493677 );
      tk.th[2] = tk.get_t_pid( 0.938272 );
      //cout << "t(track)=" <<  tk.t << " t(pi) = " << tk.th[0] << " t(K) = " << tk.th[1] << "  t(p) = " << tk.th[2] << endl;
    }else{
      for(unsigned int j = 0; j < 3; j++){ tk.th[j] = 0; }
    }
    tk._selected = theTrackFilter(*tt);

    tracks.push_back(tk);
  }

  /* moved to a separate method 
  // set-up the reco  track --> vertex pointers

  for (auto label : vertexCollectionLabels_) {

    if (recVtxs_[label] == NULL) continue;

    unsigned int iv = 0;
    for (auto recv : *recVtxs_[label]) {
      if (recv.tracksSize() > 0){
	for (trackit_t t = recv.tracks_begin(); t != recv.tracks_end(); t++) {
          auto & tk = tracks.from_key(t->key());
          const auto weight = recv.trackWeight(*t);
          assert(tk.key() == t->key());
          if( tk.get_recv(label) == NO_RECVTX ||  tk.get_weight(label) < weight){
            tk._recv[label] = iv;
            tk._weight[label] =  weight;
            assert(tracks.from_key(t->key()).get_weight(label) == weight);
            assert(tracks.from_key(t->key()).get_recv(label) == iv);
          }
	}
      }
      iv++;
    }
  }
  */
  
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
/********************************************************************************************************/


   
/********************************************************************************************************/
void  PrimaryVertexAnalyzer4PU::fill_track_to_vertex_pointers(Tracks& tracks) {
/********************************************************************************************************/
  for (auto label : vertexCollectionLabels_) {

    if (recVtxs_[label] == NULL) continue;
    
    unsigned int iv = 0;
    for (auto recv : *recVtxs_[label]) {
      if (recv.tracksSize() > 0){
	for (trackit_t t = recv.tracks_begin(); t != recv.tracks_end(); t++) {
          auto & tk = tracks.from_key(t->key());
          const auto weight = recv.trackWeight(*t);
          assert(tk.key() == t->key());
          if( tk.get_recv(label) == NO_RECVTX ||  tk.get_weight(label) < weight){
            tk._recv[label] = iv;
            tk._weight[label] =  weight;
            assert(tracks.from_key(t->key()).get_weight(label) == weight);
            assert(tracks.from_key(t->key()).get_recv(label) == iv);
          }
	}
      }
      iv++;
    }
  }
}
/********************************************************************************************************/



  
/********************************************************************************************************/
void PrimaryVertexAnalyzer4PU::add_timing_to_vertex_collection(const string & label, Tracks& tracks) {
/********************************************************************************************************
 * does the pid timing fit for each vertex one-by-one, 3d part is not touched
*********************************************************************************************************/

    if (recVtxs_[label] == NULL){
      return;
    }
    
    auto new_collection = new reco::VertexCollection;
    
    for (auto recv : *recVtxs_[label]) {
      
      auto result = vertex_time_from_tracks_pid(recv, tracks, 0.8, false);
      
      if (result.successful()){
	auto err = recv.covariance4D();
	err(3, 3) = result.tError() * result.tError();
	reco::Vertex new_recv = Vertex(recv.position(),  err, result.t(), recv.chi2(), recv.ndof(), recv.tracksSize());
	// copy the tracks and weights
	for(auto tk : recv.tracks()){   // recv.tracks() ->  std::vector<TrackBaseRef>;
	  new_recv.add(tk, recv.trackWeight(tk));
	}
	new_collection->push_back(new_recv);
      }else{
	new_collection->push_back(recv);
      }
    } // vertices
    
    delete recVtxs_[label];
    recVtxs_[label] = new_collection;
    
}
/********************************************************************************************************/
  




/********************************************************************************************************/
void PrimaryVertexAnalyzer4PU::refit_recvertices_after_timing(Tracks& tracks, double min_tk_vtx_weight) {
/********************************************************************************************************
 * does the pid timing fit for each vertex one-by-one
 * removes tracks that are rejected by the timing fit
 * first observation:
 *   small improvement of the purity but also significant reduction of assignment efficiency
 *   in cases where a "small" vertex is swamped by tracks from a "big" neighbour this doesn't help
 *    -> need a multi-vertex approach here
*********************************************************************************************************/

  for( auto label : vertexCollectionLabels_){

    if (recVtxs_[label] == NULL){
      continue;
    }
    
    auto new_collection = new reco::VertexCollection;
    
    for (auto recv : *recVtxs_[label]) {
      
      auto result = vertex_time_from_tracks_pid(recv, tracks, 0.8, false);
      
      // select the tracks
      vector<TransientTrack> vtxtracks;
      unsigned int tkidx = 0;
      for(auto tv = recv.tracks_begin(); tv != recv.tracks_end(); tv++){
	const reco::TrackRef trackRef = tv->castTo<reco::TrackRef>();
	
	TransientTrack  transientTrack = theB_->build(trackRef); 
	transientTrack.setBeamSpot(vertexBeamSpot_);
	if(result.successful()){
	  // use only accepted tracks
	  if ( result.tk_tweight[tkidx] > min_tk_vtx_weight) {
	    vtxtracks.push_back(transientTrack);
	  }
	} else {
	  // take all (may reject the vertex later)
	  vtxtracks.push_back(transientTrack);
	}
	tkidx++;
      }
      
      // refit
      AdaptiveVertexFitter  theFitter(GeometricAnnealing(2.5)); //;algorithm.fitter = new AdaptiveVertexFitter(GeometricAnnealing(algoconf->getParameter<double>("chi2cutoff")));
      TransientVertex v = theFitter.vertex(vtxtracks);

      if (v.isValid()){

	if (result.successful()){
	  // make it a 4d vertex
	  auto err = v.positionError().matrix4D();
	  err(3, 3) = result.tError() * result.tError();
	  auto trkWeightMap3d = v.weightMap();  // copy the 3d-fit weights
	  v = TransientVertex(v.position(), result.t(), err, v.originalTracks(), v.totalChiSquared(), v.degreesOfFreedom());
	  v.weightMap(trkWeightMap3d);
	}

	reco::Vertex new_recv = v;
	std::cout << " refitted   (" <<  result.status << ")" 
		  <<  "     z : " << fixed << setw(8) << setprecision(4) << recv.z() << " -> " << fixed << setw(8) << setprecision(4)  << new_recv.z() 
		  <<  "     tracks : " << fixed << setw(4) << recv.tracksSize() << " -> " <<  fixed << setw(4) << new_recv.tracksSize()
		  <<  "     ndof : "  << fixed << setw(6) << setprecision(2) << recv.ndof()  << " -> "  << fixed << setw(6) << setprecision(2) << new_recv.ndof();
	if (result.not_converged()) {
	  cout << "   rejected " << endl;
	}else{
	  cout << "   accepted " << endl;
	  new_collection->push_back(new_recv);
	}
      }else{
	// refit failed, discard or copy old version
	std::cout << " refit failed (" << result.status << ")"
		  <<  "    z : " << fixed << setw(8) << setprecision(4) << recv.z() 
		  <<  "    tracks : "  << fixed << setw(4) << recv.tracksSize() << " -> "  << fixed << setw(4) << vtxtracks.size() 
		  <<  "    ndof : "  << fixed << setw(6) << setprecision(2) << recv.ndof() 
		  << "   rejected " << std::endl;
	//new_collection->push_back(recv);
      } // valid fit ?
      
    } // vertices
    
    std::cout << "about to replace collection " << label << "  " << recVtxs_[label]->size() << " ->" << new_collection->size() << std::endl;
    delete recVtxs_[label];
    recVtxs_[label] = new_collection;
    std::cout << "done " << std::endl;
    
  }// collections
}
/********************************************************************************************************/



/************************************************************************************************************/
void PrimaryVertexAnalyzer4PU::getSimEvents_pu(PileupSummaryInfo& puInfo, std::vector<SimEvent>& simEvents)
/********************************************************************************************************
 * append the pile-up vertices to the simEvent list 
 * hard scatter assumed to be already filled in the first entry
 ************************************************************************************************************/
{
    report_counted("adding simEvts from puInfo", 1);
    if(simEvents.size() != 1 ){
      cout << "getSimEvents_pu : Warning !!!!  the size of simEvents is " <<  simEvents.size()  << " instead of 1" << endl;
    }

    for (unsigned int i = 0; i < puInfo.getPU_zpositions().size(); i++) {
      double t = 0;
      if (puInfo.has_times()) {
	t = puInfo.getPU_times()[i];  // appears to be in ns already
      }
      if (veryverbose_) {
	cout << setw(4) << i << "  z=" << setw(8) << fixed << setprecision(4) << puInfo.getPU_zpositions()[i]
	     << "  t=" << setw(7) << fixed << setprecision(4 ) << t 
	     << "  pthat= " << scientific << puInfo.getPU_pT_hats()[i] 
	     << " sumpt hi" << puInfo.getPU_sumpT_highpT()[i]                     // 0
	     << " sumpt lo" << puInfo.getPU_sumpT_lowpT()[i]                      // 0
	     << " ntrk hi" << fixed << puInfo.getPU_ntrks_highpT()[i]             // nonsense
	     << " ntrk lo" << fixed << puInfo.getPU_ntrks_lowpT()[i]              // nonsense
	     << endl;
      }

      SimEvent e(i+1);
      e.z = puInfo.getPU_zpositions()[i];
      e.x = vertexBeamSpot_.x(e.z);
      e.y = vertexBeamSpot_.y(e.z);
      e.t = t;

      e.type = FROM_PU_SUMMARY;  // partial
      e.pt_hat = puInfo.getPU_pT_hats()[i];

      simEvents.push_back(e);
    }
}
/********************************************************************************************************/






/********************************************************************************************************/
bool PrimaryVertexAnalyzer4PU::get_MC_truth(const edm::Event& iEvent,
                                            Tracks& tracks,
                                            bool bPuInfo,
                                            PileupSummaryInfo& puInfo,
                                            std::vector<SimEvent>& simEvt)
/********************************************************************************************************/
{
  //  std::string mcproduct = "generator";       // starting with 3_1_0 pre something
  tracking_truth_available_ = false;         // meaning trackingparticles

  // genParticles for AOD et al:  https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookGenParticleCandidate
  Handle<reco::GenParticleCollection> genParticlesH;
  bool bgenParticles = iEvent.getByToken(genParticleCollection_Token_, genParticlesH);
  if (bgenParticles) report_counted("found genParticles", 1);

  Handle<HepMCProduct> evtMC;              // generator level
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

  

  if (gotTP)
    {
      edm::Handle<reco::RecoToSimCollection> recoToSimH;
      iEvent.getByToken(recoToSimAssociationToken_, recoToSimH);
      tp_r2s_ = recoToSimH.product();
      tracking_truth_available_ = true;
      simEvt = getSimEvents_tp(TPCollectionH, tracks);
      analyzeTracksTP(tracks, simEvt);
    }
  else
    {
      // look for mc info for the signal vertex

      if( bSimTrks && bSimVtxs ){
	simEvt = getSimEvents_simtrks(simTrks, simVtxs, tracks);
	//tsim = PrimaryVertexAnalyzer4PU::getSimTrkParameters(simTrks, simVtxs, simUnit_);
	MC_ = true;
      }

      else if (iEvent.getByToken(edmHepMCProductToken_, evtMC)) {
	if(verbose_){
	  cout << "all we have is  hepmc" << endl;
	}
	//simEvt = getSimEvents_hepmc(evtMC); FIXME  convert obsolete getsimpv
	MC_ = true;
      }

      else if (bgenParticles){
	if(verbose_){
	  cout << "making sim events from genparticles" << endl;
	  //simEvt = getSimEvents_hepmc(evtMC); FIXME  convert obsolete getsimpv
	  MC_= false;
	}

      }

      else if (MINIAOD_){
	simEvt = getSimEvents_miniaod(iEvent, tracks);
	MC_ = (simEvt.size() >0);
      }

      //  throw in some pile-up
      if (bPuInfo) {
	report_counted("filling simEvents from signal MC  + puInfo", 1);
	getSimEvents_pu(puInfo, simEvt);
	if (verbose_) {
	  cout << "PileupSummaryInfo  nPU=" << puInfo.getPU_NumInteractions() << endl;
	}
      }
    }
  return MC_;
}//get_MC_truth
/***********************************************************************************/


/***********************************************************************************/
void PrimaryVertexAnalyzer4PU::fill_simvtx_histos(std::vector<SimEvent>& simEvts)
// fill some basic sim vertex quantities not using any reco stuff
// booked in bookSimPVHistograms, get their own subdirectory "simpvs"
// called by analyze, information filled by various getSimEvt
{
  Fill(hsimPV, "nsimvtx", simEvts.size());
  for (auto const& e : simEvts) {

    if (e.type == FROM_TRACKING_TRUTH or e.type == FROM_PU_SUMMARY) {

      Fill(hsimPV, "zvtx_sim", e.z * simUnit_, e.is_signal());
      Fill(hsimPV, "zsim", e.z * simUnit_, e.is_signal());
      Fill(hsimPV, "zsimb", e.z * simUnit_ - vertexBeamSpot_.z0(), e.is_signal());
      Fill(hsimPV, "dzminsim", e.dzmin * simUnit_, e.is_signal());
      if (f4D_) Fill(hsimPV, "tsim", e.t, e.is_signal());

      if (e.type == FROM_TRACKING_TRUTH) {
        Fill(hsimPV, "nbsimtksinvtx", float(e.nGenTrk), e.is_signal());
        Fill(hsimPV, "xsim", e.x * simUnit_, e.is_signal());
        Fill(hsimPV, "ysim", e.y * simUnit_, e.is_signal());
        Fill(hsimPV, "xsimb", e.x * simUnit_ - vertexBeamSpot_.x(e.z), e.is_signal());
        Fill(hsimPV, "ysimb", e.y * simUnit_ - vertexBeamSpot_.y(e.z), e.is_signal());

        Fill(hsimPV, "logpthatsim", std::log10(e.pt_hat), e.is_signal());
        Fill(hsimPV, "logsumptsim", std::log10(e.sumpt), e.is_signal());
        Fill(hsimPV, "logsumpt2sim", std::log10(e.sumpt2), e.is_signal());
      }
    }
  }
}

void PrimaryVertexAnalyzer4PU::beginJob() {
  matchsummaries_ = 4;   // number of match summaries to be dumped

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


  // do some normalization / cumulation etc before saving
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
  for (auto it : pvtimers_ ){
    std::cout << setw(50) << it.first   << std::setw(15) << std::fixed << std::setprecision(3) << (it.second).duration * 1e-6<< " s" << "  counter = " << (it.second).counter <<std::endl;
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

  simPU_ = 0.;

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
            (pdgCode == -321) || (pdgCode == -3222) || (pdgCode == 3312)) {
          Q = -1;
        } else if ((pdgCode == -11) || (pdgCode == -13) || (pdgCode == -15) || (pdgCode == 211) || (pdgCode == 2212) ||
                   (pdgCode == 321) || (pdgCode == 3222) ||(pdgCode == -3312)) {
          Q = 1;
        } else {
          //std::cout << pdgCode << " " <<std::endl;
        }
        math::XYZTLorentzVectorD p(t->momentum().x(), t->momentum().y(), t->momentum().z(), t->momentum().e());
        if ((Q != 0) && (p.pt() > 0.1) && (fabs(t->momentum().eta()) < 5.0) && fabs(v.z() * simUnit < 20) &&
            (sqrt(v.x() * v.x() + v.y() * v.y()) < 10.)) {

	  // FIXME, switch to the constructors of SimPart here
	  
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
	  sp.pt = 0;
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
	  sp.tvtx = 0;
	  sp.charge = Q;
	  sp.ldec =0 ;

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
void PrimaryVertexAnalyzer4PU::analyzeVertexCollectionSimTracks(
							    std::map<std::string, TH1*>& h,
							    MVertexCollection& vertexes,
							    Tracks & tracks,
							    std::vector<SimEvent>& simEvt,
							    const std::string message)
/*************************************************************************************************************/
{
  if (verbose_){
    cout << "analyzeVertexCollectionSimTracks " << message << "  simEvt.size()=" << simEvt.size() << " vertexes.size=" << vertexes.size() << endl;
  }
  //
  if (simEvt.size() ==0){
    return;
  }
  if (vertexes.size() == 0){
    return;
  }
  //
  //
  auto & e = simEvt[0];
  e.clear_matching_info();

  unsigned int nvtxtrks = 0;
  for(auto v : vertexes){
    v.clear_matching_info();
    nvtxtrks += v.tracksSize();
  }
  if (nvtxtrks == 0){
    report_counted("analyzeVertexCollectionSimTracks  " + vertexes.label() + " : no vertices have tracks attached", 10);
    return;
  }

  vector<unsigned int> r2s = supfT(e.parts, tracks); // parallel to tracks

  // analog to tpmatch, but we only have truth-matching for the signal tracks
  for(auto & v : vertexes){
    if( ! (select(v) && (fabs(v.z() - e.z) < 1.0) )) continue;
    for(auto & tv : v.tracks){
      assert( tv->index() < tracks.size());
      double wt = v.trackWeight(tv);
      if(wt < min_trk_in_vtx_weight_) continue;
      double wos = wt / pow(v.zError(),2); // FIXME (align with tpmatch, include time)
      double wnt = wt * min(tv->pt(), 1.0);  // (truncated-)pt-weighted track count
      if (r2s[tv->index()] == NOT_MATCHED_TK_SIM){
	// track not matched (at least no to event[0], don't know about others)
	// we do nothing ?
      }else{
	// track matched to event[0]
	e.addTrack(v.index(), wos, wnt);          // fill vertex wos/wnt list, increment sumwmt and sumwos
	v.add_truthmatched_track(0, wos, wnt);   // fill track(wos)  rec vtx <- sim vtx, increments sumnt and sumwos
      }
    }
  }

  // normalize and select a 'best' match
  unsigned int bestmatch =  NOT_MATCHED_VTX_SIM;
  for(unsigned int iv = 0; iv < vertexes.size(); iv++){
    MVertex & v = vertexes[iv];
    v.sigwosfrac = v.wos[0] / e.sumwos;
    v.sigwntfrac = v.wnt[0] / e.sumwnt;
    if ((v.sigwosfrac > 0.5) && (v.sigwntfrac > 0.5)){
      bestmatch = iv;
    }
  }

  if (bestmatch != NOT_MATCHED_VTX_SIM){
    vertexes[bestmatch].sim = 0;
    e.rec = bestmatch;
  }
    
  
  if(verbose_){
    cout << endl << "miniaod signal event matching, z(signal) = " << e.z << endl;
    cout << "best match " << bestmatch << endl;
    for(auto & vtxwos : e.wos){
      auto iv = vtxwos.first;
      cout << iv << "  wos = " << vertexes(iv).sigwosfrac  << "  wnt = " << vertexes(iv).sigwntfrac 
	   << " zrec =  " << vertexes(iv).z() 
	   << " zrec-zsignal   " << vertexes(iv).z() - e.z
	   << endl;
    }
    cout << endl;
  }


  Fill(h, "sigwosfractagged", vertexes(0).sigwosfrac);
  Fill(h, "sigwntfractagged", vertexes(0).sigwntfrac);
  
  double sigwosfracmax = 0., sigwntfracmax = 0.;
  for(const auto & v : vertexes){
    if(v.sigwosfrac > sigwosfracmax){ sigwosfracmax = v.sigwosfrac;}
    if(v.sigwntfrac > sigwntfracmax){ sigwntfracmax = v.sigwntfrac;}
  }
  Fill(h, "sigwosfracmax", sigwosfracmax);
  Fill(h, "sigwntfracmax", sigwntfracmax);

  // how much is the signal vertex spread out in z
  // i.e. histogram the wos/wnt fractions of vertices as a function of their distance to the true signal
  for(auto & vtxwos : e.wos){
    auto iv = vtxwos.first;
    Fill(h, "sigwosfracvsdzsignal", vertexes(iv).z() - e.z, vertexes(iv).sigwosfrac );
  }
  for(auto & vtxwnt : e.wnt){
    auto iv = vtxwnt.first;
    Fill(h, "sigwntfracvsdzsignal", vertexes(iv).z() - e.z, vertexes(iv).sigwntfrac );
  }
  
  // how many rec vertics have a wos value above some threshold
  const unsigned int nbin = 20;// number of bins in nrecwithsigwos
  vector<unsigned int> nv(nbin, 0.);
  for(auto & vtxwos : e.wos){
    auto iv = vtxwos.first;
    if (vertexes(iv).sigwosfrac == 0) continue;
    // count and cumulate
    unsigned int bin = vertexes(iv).sigwosfrac * nbin;
    if (bin > nbin){ bin = nbin;} // if sigwosfrac is == 1.0
    for(unsigned int b = 0; b < bin; b++){ nv[b]++; }
  }
  // then fill
  for(unsigned int b=0; b<nbin; b++){
    Fill(h, "nrecwithsigwos", float(b + 0.5)/nbin, float(nv.at(b)));
  }
  
  
}
/*************************************************************************************************************/



/*************************************************************************************************************/
vector<unsigned int> PrimaryVertexAnalyzer4PU::supfT(std::vector<SimPart>& simtrks, const Tracks& trks)
// track rec to sim matching for simParts  (e.g. from miniaod candidates)
/*************************************************************************************************************/
{
  bool debug_supf = false;
  unsigned int ndumpcov = 0;

  unsigned int nrec = trks.size();
  vector<unsigned int> rectosim(nrec, NOT_MATCHED_TK_SIM);

  unsigned int nsim = simtrks.size();
  if (nsim == 0)
    return rectosim;
  
  for (unsigned int j = 0; j < nsim; j++){ simtrks[j].rec = NOT_MATCHED_TK_REC; }

  // get tracks ordered by volume in z-eta-phi space // other choices det(V), # hits
  std::vector<std::pair<double, unsigned int>> vol;
  unsigned int nsel = 0;
  for (auto tk : trks){
    if(tk.has_track()){
      vol.push_back(std::make_pair(tk.track().etaError() * tk.track().phiError() * tk.track().dszError(), tk.index()));
      nsel++;
    }
  }
  std::stable_sort(vol.begin(), vol.end(), std::less<std::pair<double, unsigned int>>());


  unsigned int nmatch = 0;
  const vector<double>chi2max = {10., 100., 1000.};  // relax for subsequent passes
  const vector<double> covscale = {0.5, 1.0, 1.0, 0.01, 0.01};  // downweight pt(0), dxy(3), and z(4)
  for (unsigned int pass = 0; pass < chi2max.size(); pass++) {

    for (const auto & vi : vol){

      unsigned int irec = vi.second;
      const auto & t = trks(irec).track();
      if (rectosim[irec] != NOT_MATCHED_TK_SIM) continue; // already matched

      // get track parameters and covariance matrix
      ParameterVector par = t.parameters();
      reco::TrackBase::CovarianceMatrix S = t.covariance();

      // pseudotracks sometimes come with > 1 correlations, especially at high pt
      //double t34 = t.covariance(3,4) / sqrt(t.covariance(3,3)*t.covariance(4,4));
      double t14 = t.covariance(1,4) / sqrt(t.covariance(1,1)*t.covariance(4,4));
      double t23 =  t.covariance(2,3) / sqrt(t.covariance(2,2)*t.covariance(3,3));
      string msg2 = Form("track pt= %8.2f ", t.pt());
      if( std::abs(t23) > 0.95){
	report_counted("fudged dphi    dxy correlation ", msg2, 10);
	S(2,3) = t23 / std::abs(t23) * 0.9 * sqrt(t.covariance(2,2)*t.covariance(3,3));
      }
      if( std::abs(t14) > 0.95){
	report_counted("fudged dlambda dz  correlation ",msg2, 10);
	S(1,4) = t14 / std::abs(t14) * 0.9 * sqrt(t.covariance(1,1)*t.covariance(4,4));
      }
      
      if (!S.Invert()) {
	report_counted("covariance matrix inversion failed for track ", 10);
        continue;
      }

      double mindiag = S(0,0);
      for(unsigned int a=1; a<5; a++){
	if (S(a,a) < mindiag){ mindiag = S(a,a);}
      }
      if(mindiag < 0){
	report_counted("inverse covariance hase negative diagonal values ", 10);
	if(debug_supf){
	  cout << "track.pt = " << t.pt() << endl;
	  cout << S << endl;
	  cout << "cov(3,4) = " << t.covariance(3,4) / sqrt(t.covariance(3,3)*t.covariance(4,4)) << endl;
	  cout << "cov(1,4) = " <<  t.covariance(1,4) / sqrt(t.covariance(1,1)*t.covariance(4,4)) << endl;
	  cout << "cov(2,3) = " <<  t.covariance(2,3) / sqrt(t.covariance(2,2)*t.covariance(3,3)) << endl;
	}
	continue;
      }

      // try the remaining sim tracks for this rec track
      unsigned int jmatch = NOT_MATCHED_TK_SIM;
      unsigned int jnextbest = NOT_MATCHED_TK_SIM;
      double cmatch = 1e10;
      double cnextbest = 1e10;

      for (unsigned int j = 0; j < nsim; j++) {
        if (simtrks[j].rec != NOT_MATCHED_TK_REC) continue;
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

	// are we even close in eta-phi?
	if( (std::abs(d[2]) > 0.2) || (std::abs(d[3])>0.2) ) continue;

        // matching chi**2
        double c = 0, cd = 0;
        for (int k = 0; k < 5; k++) {
	  cd += d[k] * d[k] / t.covariance()(k, k); // back-up for meaningless inverse
          for (int l = 0; l < 5; l++) {
            c += d[k] * d[l] * S(k, l) * covscale[k] * covscale[l];
          }
        }

        if (c <= 0) {
          report_counted("supfT: non-positive-definite covariance matrix encountered ", 10);
	  if( ndumpcov++ < 10){
	    cout << scientific << setw(8) << setprecision(2);
	    cout << "t.pt() = " << t.pt() << endl;
	    cout << "c = " << c << endl;
	    cout << "cd = " << cd << endl;
	    cout << t.covariance() << endl;
	    cout << "inverse:" <<endl;
	    cout << S << endl;
	  }
	  c=cd;
        }

        if ((debug_supf) && (pass==3)) {
          double c0 =   pow((par[0] - s.par[0]) / t.qoverpError(), 2) * 0.1  //!
                      + pow((par[1] - s.par[1]) / t.lambdaError(), 2)
	              + pow((par[2] - s.par[2]) / t.phiError(), 2) 
                      + pow((par[3] - s.par[3]) / t.dxyError(), 2) * 0.1;  //!
	  if ( std::abs((par[1] - s.par[1]) < 0.5)  && std::abs((par[2] - s.par[2]) < 0.5) ){
            cout << "pass " << pass << endl;
            cout << setw(4) << irec << " rec " << scientific << setw(10) << setprecision(2) << par << "    " 
                 << fixed << setprecision(3) << setw(8) << t.charge() << " " << t.pt() << " " 
		 << t.phi() << " " << t.eta() << endl;
            cout << setw(4) << j << " sim " << scientific << setw(10)  << setprecision(2) << s.par << "    " 
		 << fixed << setprecision(3) << setw(8) << simtrks[j].charge << " " << simtrks[j].pt 
		 << " "  << simtrks[j].phi << " " << simtrks[j].eta  << "  pdg=" << simtrks[j].pdg 
		 <<  " ---> C=" << scientific << c << endl;
            cout << "       "  << fixed << setw(10) << setprecision(4)
		 << setw(10) << d[0] << " " << setw(10) << d[1]<< " "  << setw(10) << d[2]<< " "  << setw(10) << d[3] << " " <<  setw(10) << d[4]
                 << "   match=" << match(simtrks[j].par, trks(irec).track().parameters()) 
		 << endl;
            cout << "       " 
		 << fixed << setw(10) << setprecision(1)
		 << setw(10) << (par[0] - s.par[0]) / t.qoverpError() << " "
                 << setw(10) << (par[1] - s.par[1]) / t.lambdaError() << " "
		 << setw(10) << (par[2] - s.par[2]) / t.phiError() << " "
                 << setw(10) << (par[3] - s.par[3]) / t.dxyError() << " "
		 << setw(10) << (par[4] - s.par[4]) / t.dszError()  << " "
		 << "   c0=" << scientific << c0
                 << endl
                 << endl;
          }
        }

        if ((jmatch == NOT_MATCHED_TK_SIM) || (c < cmatch)) {
          jnextbest = jmatch;
          cnextbest = cmatch;
          jmatch = j;
          cmatch = c;
        } else if ((jnextbest  == NOT_MATCHED_TK_SIM) || (c < cnextbest)) {
          jnextbest = j;
          cnextbest = c;
        }
      }

      if (debug_supf) {
      cout << setw(4) << irec << " rec --> " << setw(4) << jmatch;
        if (jmatch != NOT_MATCHED_TK_SIM) {
	  if(cmatch < chi2max[pass]) { cout <<"**";}else{cout << "  ";}
          cout << " cmatch = " << fixed << setw(8) << setprecision(1) << cmatch << "  next best= " << cnextbest;
        }
        cout << endl;// << endl << endl;
      }

      if ((jmatch != NOT_MATCHED_TK_SIM) && (cmatch < chi2max[pass])) {
	// for the last pass with insanely poor chi**2 require at least that nothing else is close
	if ((pass == chi2max.size()+1 ) && (cmatch > 0.1 * cnextbest)) continue;
        nmatch++;
        rectosim[irec] = jmatch;
        simtrks[jmatch].rec = irec;
      }

    }  // rec trk
  }    // pass

  if (verbose_) {
    std::cout << "simtracks (pt>1) without a matching rec track: " << std::endl;
    int nunmatched = 0;
    for (unsigned int j = 0; j < nsim; j++) {
      if (simtrks[j].rec == NOT_MATCHED_TK_REC) {
        nunmatched++;
        double pt = 1. / simtrks[j].par[reco::TrackBase::i_qoverp] / tan(simtrks[j].par[1]);
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
  return rectosim; 
}  //supfT

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

double PrimaryVertexAnalyzer4PU::vertex_sumw(const reco::Vertex& v) {
  double wsum = 0.;
  for (trackit_t t = v.tracks_begin(); t != v.tracks_end(); t++) {
    wsum += v.trackWeight(*t);
  }
  return wsum;
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


bool  PrimaryVertexAnalyzer4PU::uv_analysis(const MVertex & vtx, const SimEvent & simevt, double & du, double & dv, double & uError, double & vError, double &pol){
  du = 0;
  dv = 0;
  uError = 1.;
  vError = 1.;
  double pt_max = 0, px_max = 0, py_max = 0;
  double sumpxy = 0, sumpxx=0, sumpyy = 0;
  for (auto tk : vtx.tracks){
    double px = tk->px();
    double py = tk->py();
    sumpxy += px * py;
    sumpxx += px * px;
    sumpyy += py * py;
    if(tk->pt() > pt_max){
      pt_max = tk->pt();
      px_max = tk->px();
      py_max = tk->py();
    }
  }
  double nx = 0, ny = 0; 
  if(sumpxx != sumpyy){
    double tan2phi = 2 * sumpxy / (sumpxx - sumpyy);
    double phi0 = 0.5 * atan(tan2phi);
    double T0 = sumpyy*cos(phi0)*cos(phi0) - 2*sumpxy * cos(phi0) * sin(phi0) + sumpxx * sin(phi0) * sin(phi0);
    double phi1 = phi0 + 0.5 * 3.14159;
    double T1 = sumpyy*cos(phi1)*cos(phi1) - 2*sumpxy * cos(phi1) * sin(phi1) + sumpxx * sin(phi1) * sin(phi1);
    // make (nx, ny) point into the thrust direction (minimizes perpendicular momentum**2 sum)
    // a small T0/T1 should correspond to a pencil-like track list, T0 ~ T1 more spherical
    if(T0 < T1){
      nx = cos(phi0);
      ny = sin(phi0);
      pol = T0/T1;
    }else{
      nx = cos(phi1);
      ny = sin(phi1);
      pol = T1/T0;
    }
    /*
    std::cout << "T0= " << T0 << "   T1=" << T1 <<  "   phi0=" << phi0 << "   phi1= " << phi1 << std::endl;
    for(unsigned int j =0; j< 12; j++){
      double phi = 2 * 3.14159 * j / 12;
      double pp2 = 0;
      double pl2 = 1e-10;
      for (auto tk : vtx.tracks){
	pp2 += pow(tk->pt() * sin(tk->phi() - phi), 2);
	pl2 += pow(tk->pt() * cos(tk->phi() - phi), 2);
      }
      double T = sumpyy*cos(phi)*cos(phi) - 2*sumpxy * cos(phi) * sin(phi) + sumpxx * sin(phi) * sin(phi);
      std::cout << " T  ("<< phi << ")   =  "  << T <<  " pp2 = " << pp2 <<   "   pl2 = " << pl2 <<   " p/l = " << pp2/pl2 << std:: endl;
    }
    */
  }else{
    if (pt_max == 0) {
      return false;
    }
    nx = px_max / pt_max;
    ny = py_max / pt_max;
  }

  // (nx,ny) = thrust direction, minimizes the sum of perpendicular momentum**2 sum
  // u = (x,y) component along the thrust axis
  // v = (x,y) component perpendicular to the thrust axis
  double vxx = vtx.covariance(iX, iX);
  double vyy = vtx.covariance(iY, iY);
  double vxy = vtx.covariance(iX, iY);
  double dx = vtx.x() - simevt.x;
  double dy = vtx.y() - simevt.y;
  du = dx * nx + dy * ny; 
  dv = dx * ny - dy * nx;
  uError = sqrt(nx * nx * vxx + ny * ny * vyy + 2 * nx * ny * vxy);
  vError = sqrt(ny * ny * vxx + nx * nx * vyy - 2 * nx * ny * vxy);
  return true;
}


PrimaryVertexAnalyzer4PU::Vertex_time_result PrimaryVertexAnalyzer4PU::vertex_time_from_tracks(const reco::Vertex& v,
                                                       Tracks& tracks,
						       double minquality,
						       bool verbose) {

  PrimaryVertexAnalyzer4PU::Vertex_time_result result;

  const double min_trk_in_vtx_weight = 0.2;
  double tsum = 0;
  double wsum = 0;
  double w2sum = 0;
  double t = 0;
  for (auto tk = v.tracks_begin(); tk != v.tracks_end(); tk++) {
    if (v.trackWeight(*tk) > min_trk_in_vtx_weight) {
      auto trk = tracks.from_ref(*tk);
      if (trk.has_timing() && (trk.timeQuality() >= minquality)) {
        double w = v.trackWeight(*tk) / (trk.dt() * trk.dt());
        wsum += w;
        tsum += w * trk.t();
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
	    if (v.trackWeight(*tk) > min_trk_in_vtx_weight)
	      {
		auto trk = tracks.from_ref(*tk);
		if (trk.has_timing() && (trk.timeQuality() >= minquality))
		  {
		    double tpull = (trk.t() - t0) / trk.dt();
		    if (fabs(tpull) < 5.)
		      {
			double wt = 1./(1.+ exp(0.5 * tpull * tpull - 0.5 * 9.));
			double w = wt * v.trackWeight(*tk) / (trk.dt() * trk.dt());
			wsum += w;
			tsum += w * trk.t();
			w2sum += w * w * trk.dt() * trk.dt();
		      }
		  }
	      }
	  }

    	if (wsum > 0)
	  {
	    t = tsum / wsum;
	    if (fabs(t - t0) < 1e-3)
	      {
		//tError = sqrt(w2sum) / wsum;
		result.success(tsum / wsum, sqrt(w2sum)/ wsum, nit);
		return result;
	      }
	  }
	t0 = t;
      }
  }
  return result;
}


PrimaryVertexAnalyzer4PU::Vertex_time_result PrimaryVertexAnalyzer4PU::vertex_time_from_tracks_pid(const reco::Vertex& v,
							   Tracks& tracks,
							   double minquality,
							   bool verbose) {
  PrimaryVertexAnalyzer4PU::Vertex_time_result result;
  result.tk_tweight.assign(v.tracksSize(), 1.);
  double tsum = 0;
  double wsum = 0;
  double w0sum = 0;
  double w2sum = 0;
  double t;
  const double min_trk_in_vtx_weight = 0.2;
  double a[3] = {0.7,0.2,0.1};
  constexpr double cooling_factor = 0.5;


  if(verbose) {
    cout << "vertex_time_from_tracks_pid vtx x=" << v.x()
	 << " y=" << v.y()
	 << " z=" << v.z()
	 << " t=" << v.t()
         << endl;
  }

  // initial guess
  for (auto tk = v.tracks_begin(); tk != v.tracks_end(); tk++) {
    if (v.trackWeight(*tk) > min_trk_in_vtx_weight) {
      auto trk = tracks.from_ref(*tk);
      if (trk.has_timing() && (trk.timeQuality() >= minquality)) {
	double w = v.trackWeight(*tk) / (trk.MTD_timeerror()*trk.MTD_timeerror());
        wsum += w;
	for(unsigned int j=0; j < 3; j++){
	  tsum += w * trk.th[j] * a[j];
	}

	if(verbose){
          cout << "vertex_time_from_tracks_pid:     track"
               << " pt=" << trk.pt() << " eta=" << trk.eta() << " phi=" << trk.phi() << " t=" << trk.t()
               << " vtxWeight=" << v.trackWeight(*tk) << " time=" << trk.MTD_time() << " timeError=" << trk.MTD_timeerror()
               << " timeQuality=" << trk.timeQuality() << " pathLength=" << trk.MTD_pathlength() << " momentum=" << trk.MTD_momentum()
               << " timeHyp[pion]=" << trk.th[0] << " timeHyp[kaon]=" << trk.th[1]
               << " timeHyp[proton]=" << trk.th[2] << endl;
        }
      }
    }
  }

  if (wsum > 0) {

    if(verbose) {
      cout << "vertex_time_from_tracks_pid  minquality=" << minquality 
	   << " wsum = " << wsum 
	   << " tsum = " << tsum 
	   << " t0 = " << (wsum > 0 ? tsum/wsum : 0)
	   << " trec = " << v.t();
      cout << endl;
    }

    double t0 = tsum / wsum;
    int nit = 0;
    double beta = 1./256.;
    while ( (nit++) < 100)
      {
	tsum = 0;
	wsum = 0;
	w2sum = 0;
	w0sum = 0;
	double wtcsum = 0;
	unsigned int tkidx = 0;
	for (auto tk = v.tracks_begin(); tk != v.tracks_end(); tk++)
	  {
	    if (v.trackWeight(*tk) > min_trk_in_vtx_weight)
	      {
		auto trk = tracks.from_ref(*tk);
		if (trk.has_timing() && (trk.timeQuality() >= minquality))
		  {
		    double dt = trk.MTD_timeerror();
		    double e[3] = {0,0,0};
		    double Z = exp(-beta * 0.5* 3.* 3.);
		    for(unsigned int j = 0; j < 3; j++){
		      double tpull =  (trk.th[j] - t0) / dt;
		      e[j] = exp(- 0.5 * beta * tpull * tpull);
		      Z += a[j] * e[j];
		    }

		    double wtsum_trk = 0;
		    double wsum_trk = 0; // aka s10
		    double s11 = 0, s21 = 0;
		    for(unsigned int j = 0; j < 3; j++){
		      double wt = a[j] * e[j] / Z;
		      wtsum_trk += wt;
		      double w = wt / (dt * dt);
		      // for the vertex time
		      wsum += w * v.trackWeight(*tk);
		      tsum += w * v.trackWeight(*tk) * trk.th[j]; 
		      // for error propagation
		      wsum_trk += w;
		      s11 +=  (trk.th[j] - t0)  * w;
		      s21 +=  pow(trk.th[j] - t0, 2)  * w;
		    }
		    result.tk_tweight[tkidx] = wtsum_trk * v.trackWeight(*tk);
		    double dfodt = (wsum_trk + beta * (s11*s11 - s21/(dt*dt)) )* v.trackWeight(*tk); // aka f_i
		    w2sum += pow(dfodt * dt, 2);
		    w0sum += dfodt;
		    wtcsum += s21/(dt*dt);
		  }// track has timing
	      } // track weight > min
	    tkidx ++;
	  }
	if (wsum < 1e-10)
	  {
	    //report_counted("PrimaryVertexAnalyzer4PU::vertex_time_from_tracks_pid failed while iterating", 10);
	    return result;
	  }

	t = tsum / wsum;
	if(verbose){
	  cout << " XX " << nit << " T= " << 1/beta << "  t=" << t <<  "    t-t0=" <<  t-t0 << "  Tc=" << wtcsum/wsum << endl;
	  // dump the full monty
	  unsigned int tkidx = 0;
	  for (auto tk = v.tracks_begin(); tk != v.tracks_end(); tk++)
	    {
	      if (v.trackWeight(*tk) > min_trk_in_vtx_weight)
		{
		  auto trk = tracks.from_ref(*tk);
		  if (trk.has_timing() && (trk.timeQuality() >= minquality))
		  {
		    double dt = trk.MTD_timeerror();
		    double e[3] = {0,0,0};
		    double Z = exp(-beta * 0.5* 3.* 3.);
		    for(unsigned int j = 0; j < 3; j++){
		      double tpull =  (trk.th[j] - t0) / dt;
		      e[j] = exp(- 0.5 * beta * tpull * tpull);
		      Z += a[j] * e[j];
		    }
		    cout << "    " << fixed << setw(7) << setprecision(7) << trk.t() 
			 << "+/-" << fixed << setw(6) << setprecision(7) << dt;
		    for(unsigned int j = 0; j < 3; j++){
		      double wt = a[j] * e[j] / Z;
		      cout << "  [" <<  fixed << setw(1) << j << "]  " 
			   << fixed << setw(7) << setprecision(7) << trk.th[j] << "ns  " 
			   << fixed << setw(7) << setprecision(5) << wt;  
		    }
		    cout << "  wtot = " << fixed << setw(7) << setprecision(5) << result.tk_tweight[tkidx];  
		    if(trk.is_matched()) {
		      cout << "   m " << fixed << setw(8) << setprecision(7) << trk.get_t_pid()
			   << " gen " << fixed << setw(8) << setprecision(7) << trk.tsim();
		    }
		    cout << endl;
		  }
	      }
	      tkidx ++;
	    }// end of track dump loop
	}
	if ((fabs(t - t0) < 1e-4 / sqrt(beta)) && (beta >= 1.))
	  {
	    double tError = sqrt(w2sum) / w0sum;
	    if(verbose) {
	      cout << "vertex_time_from_tracks_pid         minquality=" << minquality 
		   << " tfit = " << t << " +/- " << tError
		   << " trec = " << v.t()
		   << " iteration " <<  nit
		   << " Tc = " <<  wtcsum/wsum ;
	      cout << endl;
	    }
	    result.success(t, tError, nit,  wtcsum/wsum );
	    return result;
	  }
	
	if ((fabs(t - t0) < 1e-3) && ((beta < 1.)))
	  { 
	    beta = std::min(1., beta / cooling_factor);
	  }

	t0 = t;

      }
    //report_counted("PrimaryVertexAnalyzer4PU::vertex_time_from_tracks_pid failed to converge", 10);
    result.status = 3;
  }else{
    //report_counted("PrimaryVertexAnalyzer4PU::vertex_time_from_tracks_pid has no track timing info", 10);
    result.status = 2;
  }
  return result;
}




PrimaryVertexAnalyzer4PU::Vertex_time_result PrimaryVertexAnalyzer4PU::vertex_time_from_tracks_pid_newton(const reco::Vertex& v,
							   Tracks& tracks,
							   double minquality, 
							   bool verbose) {
  /* determine the vertex time and best mass assignment using annealing 
     this version applies a Newton method to minimize F
  */
  PrimaryVertexAnalyzer4PU::Vertex_time_result result;
  double t;

  bool lverbose = verbose;
  double a_hyp[3] = {0.7,0.2,0.1};  // prior probabilities of mass hypotheses, FIXME : to be determined
  constexpr double cooling_factor = 0.5;
  constexpr double cut_off = 3.;
  constexpr double min_trk_in_vtx_weight = 0.2;


  // initial guess
  double S0 = 0, S1 = 0;
  double S0p = 0, S1p = 0, S2p = 0;
  for (auto tk = v.tracks_begin(); tk != v.tracks_end(); tk++) {
    if (v.trackWeight(*tk) > min_trk_in_vtx_weight) {
      auto trk = tracks.from_ref(*tk);
      if (trk.has_timing() && (trk.timeQuality() >= minquality)) {
	double q = 1. / (trk.MTD_timeerror()*trk.MTD_timeerror());
	double w = v.trackWeight(*tk) *q;
	double wp = w * q;
        S0 += w;
	S0p += wp;
	for(unsigned int j=0; j < 3; j++){
	  S1  += w  * a_hyp[j] * trk.th[j];
	  S1p += wp * a_hyp[j] * trk.th[j];
	  S2p += wp * a_hyp[j] * trk.th[j] * trk.th[j];
	}
      }
    }
  }
  double t0 = 0;
  double Tmax = 512.;
  double beta = 1 / Tmax;
  if (S0 > 0) {
    t0 = S1 / S0;
    Tmax = (S2p - 2*t0 * S1p + S0p * t0 * t0) / S0;
    int coolingsteps = 1 - int(std::log(Tmax) / std::log(cooling_factor));
    beta = std::pow(cooling_factor, coolingsteps);


    if(lverbose) {
      cout << "vertex_time_from_tracks_pid_newton  minquality=" << minquality 
	 << " sum(w) = " << S0 
	 << " sum(t) = " << S1
	 << " t0 = " << (S0 > 0 ? S1 / S0 : 0)
	 << " trec = " << v.t()
	 << " Tmax = " << Tmax
	 << " Tstart = " << 1/beta;
      cout << endl;
    }



    int nit = 0;
    int nit_newt = 0;
    double beta0 = beta / cooling_factor;  // previous temperature
    double Z0 = exp(-beta * 0.5* cut_off* cut_off); // update when T changes
    double F0 = 1, F1 = 0, F2 = 0;
    double Sw2 = 0;
    double Sdtwprime = 0;

    while ( (nit++) < 100)
      {

	unsigned int nT = 0;
	while(nT++ < 10){ // may have to repeat after T-changes if we did not end up in a minimum
	  S1 = 0;
	  S0 = 0; 
	  F0 = 1.; // really exp(-beta F0), temporary, not necessarily needed
	  Sw2 = 0;
	  Sdtwprime = 0;
	  for (auto tk = v.tracks_begin(); tk != v.tracks_end(); tk++)
	    {
	      if (v.trackWeight(*tk) > min_trk_in_vtx_weight)
		{
		  auto trk = tracks.from_ref(*tk);
		  if (trk.has_timing() && (trk.timeQuality() >= minquality))
		    {
		      double sigmat_trk = trk.MTD_timeerror();
		      double sigmat_trk_sq = sigmat_trk * sigmat_trk;
		      
		      double ae_hyp[3] = {0,0,0};
		      double deltat_hyp[3];
		      double Z_trk = Z0;
		      for(unsigned int j = 0; j < 3; j++){
			deltat_hyp[j] = trk.th[j] - t0;
			ae_hyp[j] = a_hyp[j] * exp(- 0.5 * beta * deltat_hyp[j] * deltat_hyp[j] / sigmat_trk_sq);
			Z_trk += ae_hyp[j];
		      }
		      
		      double Sdt2w_trk = 0, S0_trk = 0, S1_trk = 0, S0w_trk = 0;
		      double one_over_Z_sigma_sq = 1./(Z_trk * sigmat_trk_sq);
		      for(unsigned int j = 0; j < 3; j++){
			double w_hyp = ae_hyp[j] * one_over_Z_sigma_sq ;
			S0_trk += w_hyp;
			S1_trk += w_hyp * trk.th[j];
			S0w_trk += w_hyp * w_hyp;            // needed for error calculation only
			Sdt2w_trk += deltat_hyp[j] * deltat_hyp[j] * w_hyp;
		      }
		      S0  += v.trackWeight(*tk) * S0_trk;
		      S1  += v.trackWeight(*tk) * S1_trk;
		      Sw2 += v.trackWeight(*tk) * S0w_trk * sigmat_trk_sq; // error calculation only
		      auto Sdtw_trk = S1_trk - t0 * S0_trk;
		      Sdtwprime += v.trackWeight(*tk) * Sdt2w_trk / sigmat_trk_sq - Sdtw_trk * Sdtw_trk;
		      
		      F0 *= Z_trk;
		    }
		}
	    }
	  
	  Sdtwprime *= beta; 
	  F1 = t0 * S0 - S1;
	  F2 = S0 - Sdtwprime;

	  if (( F2 < 0 ) && (beta > beta0)){
	    if(lverbose){std::cout << "T-step rejected  " << nT << "   T(try) = " << 1/beta << "   F2="   << F2 << "  Tnew =  " << 1/sqrt(beta * beta0) <<  "  T0=" << 1/beta0<< std::endl;}
	    // T-step too big ?
	    beta = sqrt(beta * beta0);
	    Z0 = exp(-beta * 0.5 * cut_off * cut_off);
	  }else{
	    break; // T-step accepted
	  }

	} // one more time

	// T-step accepted (or gave up)

	double dF_rose = -F1*F1 / S0;
	double dF_newt = -F1*F1 / F2; // same as t0 - F1/F2 
	double t_rose = S1 / S0;
	double t_newt = (S1 - t0 * Sdtwprime) / (S0 - Sdtwprime);// same as t0 - F1/F2 
	if(lverbose){
	  std::cout << " t_rose = " << t_rose << " t_newton = " << t_newt  << "  dF_rose = " << dF_rose << "  dF_newton = " << dF_newt << std::endl;
	}

	if (dF_rose < dF_newt){
	  t = t_rose;
	  if(lverbose){std::cout << " newton step rejected     F0=" << -log(F0)/beta <<  "  F1=" << F1 << "  F2=" << F2 << "  Sw2= " <<  Sw2 <<  " tError=" <<sqrt(Sw2) / S0<< std::endl;}
	  // means F2 > F0, i.e. Sdtwprime < 0 
	}else{
	  if(lverbose){ std::cout << " newton step accepted     F0=" << -log(F0)/beta <<  "  F1=" << F1 << "  F2=" << F2 << "  Sw2= " <<  Sw2 <<  " tError=" <<sqrt(Sw2) / S0<< std::endl;}
	  t = t_newt;
	  nit_newt ++;
	}

	
	if(lverbose){
	  cout << " YY " << nit << " T= " << 1/beta << "  t=" << t <<  " +/- " << sqrt(Sw2)/S0 <<  "    t-t0=" <<  t-t0 << "  F0=" << -log(F0)/beta<< endl;
	  // dump the full monty
	  for (auto tk = v.tracks_begin(); tk != v.tracks_end(); tk++)
	    {
	      if (v.trackWeight(*tk) > min_trk_in_vtx_weight)
		{
		  auto trk = tracks.from_ref(*tk);
		  if (trk.has_timing() && (trk.timeQuality() >= minquality))
		  {
		    double sigmat_trk = trk.MTD_timeerror();
		    double ae[3] = {0,0,0};
		    double Z = Z0;
		    for(unsigned int j = 0; j < 3; j++){
		      double tpull =  (trk.th[j] - t0) / sigmat_trk;
		      ae[j] = a_hyp[j] * exp(- 0.5 * beta * tpull * tpull);
		      Z += ae[j];
		    }
		    cout << "    " << fixed << setw(7) << setprecision(4) << trk.t() 
			 << "+/-" << fixed << setw(6) << setprecision(4) << sigmat_trk;
		    for(unsigned int j = 0; j < 3; j++){
		      double wt = ae[j] / Z;
		      cout << "  [" <<  fixed << setw(1) << j << "]  " 
			   << fixed << setw(7) << setprecision(4) << trk.th[j] << "ns  " 
			   << fixed << setw(7) << setprecision(4) << wt;  
		    }
		    if(trk.is_matched()) {
		      cout << "   m " << fixed << setw(8) << setprecision(4) << trk.get_t_pid()
			   << " gen " << fixed << setw(8) << setprecision(4) << trk.tsim();
		    }
		    cout << endl;
		  }
	      }
	    }// end of track dump loop
	} // verbose

	if ((fabs(t - t0) < 1e-4 / sqrt(beta)) && (beta >= 1.))
	  {
	    double tError = sqrt(Sw2) / S0;
	    if(lverbose) {
	      cout << "vertex_time_from_tracks_pid_newton  minquality=" << minquality 
		   << " tfit = " << t << " +/- " << tError
		   << " trec = " << v.t()
		   << " iteration " <<  nit  << "   newton " << nit_newt;
	      cout << " Tfinal = " << 1./beta;
	      cout << endl;
	    }
	    result.success(t, tError, nit);
	    return result;
	  }
	
	if ((fabs(t - t0) < 1e-3) && ((beta < 1.)))
	  { 
	    beta0 = beta;
	    beta = std::min(1., beta / cooling_factor);
	    Z0 = exp(-beta * 0.5 * cut_off * cut_off);
	  }

	t0 = t;

      }
    result.status = 3;
    result.niteration = nit;
    //report_counted("PrimaryVertexAnalyzer4PU::vertex_time_from_tracks_pid_newton failed to converge", 10);
  }else{
    result.status = 2;
    //report_counted("PrimaryVertexAnalyzer4PU::vertex_time_from_tracks_pid_newton has no track timing info", 10);
  }
  return result;
}




PrimaryVertexAnalyzer4PU::Vertex_time_result PrimaryVertexAnalyzer4PU::vertex_time_from_tracks_analysis(
								PrimaryVertexAnalyzer4PU::Vertex_time_result  time_method(const reco::Vertex&, Tracks& , double , bool ), 
								const MVertex & v,
								Tracks& tracks,
								double minquality,
								const SimEvent & simevt, 
								const string label,
								std::map<std::string, TH1*>& h,
								const string vtype, 
								bool verbose){
  string suffix = (label=="") ? label : "_"+label;
  
  string timer_label = "vertex_time_fromtracks" + suffix;
  if (!verbose) timer_start(timer_label);
  auto result = time_method(*(v.recvtx), tracks, minquality, verbose);
  if (!verbose) timer_stop(timer_label);


  if (!verbose){
    if (result.successful())
      {
	// timing from tracks
	Fill(h, vtype + "/trecsim_fromtracks" + suffix, result.t() - simevt.t, simevt.is_signal());
	Fill(h, vtype + "/trecerr_fromtracks" + suffix, result.tError(), simevt.is_signal());
	Fill(h, vtype + "/trecsimpull_fromtracks" +suffix, (result.t() - simevt.t) / result.tError(), simevt.is_signal());
      
	// also fill histos with the default values for the same list of vertices for comparison
	if(v.has_timing()){
	  Fill(h, vtype + "/trecsim_withtracks" + suffix, v.t() - simevt.t, simevt.is_signal());
	  Fill(h, vtype + "/trecerr_withtracks" + suffix, v.tError(), simevt.is_signal());
	  Fill(h, vtype + "/trecsimpull_withtracks" + suffix, (v.t() - simevt.t) / v.tError(), simevt.is_signal());
	}
      }
    else
      {
	report_counted("fillVertexHistosMatched: "+label+ " timing from tracks failed",1);
      }
  }
  return result;
}



int multi_vertex_dump_counter_ = -1;
void PrimaryVertexAnalyzer4PU::mass_constrained_multi_vertex_time_from_tracks_pid(const std::string label,
							   Tracks& tracks,
							   double minquality,
							   bool verbose) {
  multi_vertex_dump_counter_++;
  // work on this collection
  auto & input_vertices = recVtxs_[label];


  // cluster properties and clusters cluster-sums
  unsigned int nvtx = input_vertices->size();
  std::vector<double> tsum(nvtx);
  std::vector<double> psum(nvtx);
  std::vector<double> wsum(nvtx);
  std::vector<double> w2sum(nvtx);
  std::vector<double> w0sum(nvtx);
  std::vector<double> vtx_rhot(nvtx);
  std::vector<double> vtx_t(nvtx);
  std::vector<double> vtx_tError(nvtx);
  std::vector<double> vtx_z(nvtx);     // input only

  

  // build a track table for the relevant tracks with timing
  std::vector<unsigned int> timing_tracklist;
  std::vector<unsigned int> timing_track_vtx_input;
  for(unsigned int k =0 ; k< nvtx; k++){
    auto & v = input_vertices->at(k);
    vtx_z[k] = v.z();
    for (auto tk = v.tracks_begin(); tk != v.tracks_end(); tk++) {
      if(v.trackWeight(*tk) < min_trk_in_vtx_weight_) continue; // just a test
      unsigned int itrk = tracks.index_from_ref(*tk);
      if (tracks[itrk].has_timing() && (tracks[itrk].timeQuality() >= minquality)){
	if (find(timing_tracklist.begin(), timing_tracklist.end(), itrk) ==  timing_tracklist.end()){
	  timing_tracklist.push_back(itrk);
	  timing_track_vtx_input.push_back(k);
	  vtx_rhot[k] += 1.;
	}else{
	  // should not happen currently
	  cout << "track assigned to multiple vertices ?" << endl;
	}
      }
    }
  }

  // FIXME use previously unassigned tracks ?

  // prepare arrays, pre-fill track z and z-error
  unsigned int ntrk = timing_tracklist.size();
  std::vector<double> Z(ntrk);
  std::vector<double> track_z_error(ntrk);
  std::vector<double> track_z(ntrk);
  std::vector<double> track_vtx_assignment_prob(ntrk);
  std::vector<unsigned int> track_vtx_assignment(ntrk);

  for(unsigned int i = 0; i < ntrk; i++){
    auto trk = tracks[timing_tracklist[i]];
    auto k = timing_track_vtx_input[i];
    track_z[i] = trk.z();
    track_z_error[i] = sqrt(pow(trk.dzError(), 2) + pow(0.0020 / tan(trk.theta()),2) + pow(0.0020,2));
    if(multi_vertex_dump_counter_ < 1){
      cout << fixed << setw(4) << setprecision(0) << i <<  ")"
	   << " k= " << fixed << setw(3) << setprecision(0) << k
	   << " zvtx = "  << fixed << setw(8) << setprecision(4) << vtx_z[k]
	   << " ztrk = "  << fixed << setw(8) << setprecision(4) << track_z[i] 
	   << " +/-"  << fixed << setw(7) << setprecision(4) << track_z_error[i] 
	   << " zpull = " << fixed << setw(5) << setprecision(2) << (track_z[i]-vtx_z[k])/track_z_error[i]
	   << "     ttrk[pi/K/p] = "   << fixed << setw(7) << setprecision(3) << trk.th[0] << " " << setw(7) << setprecision(3)<< trk.th[1] << " " << setw(7) << setprecision(3) << trk.th[2]
	   << "  +/-"  << fixed << setw(7) << setprecision(3) << trk.MTD_timeerror();
      if(trk.is_matched()){
	cout << "  zsim = " <<  fixed << setw(8) << setprecision(4) << trk.zsim()
	     << "  tsim = " <<  fixed << setw(7) << setprecision(3) << trk.tsim()
	     << "  m=" << trk.get_particle_mass();
      }
      cout << endl;
    }
  }


  // constants
  double a[3] = {0.7,0.2,0.1};
  constexpr double cooling_factor = 0.5;
  const double Ezmax = 50;                        // ignore everything above
  const double Etmax = 50;                        // ignore everything above
  constexpr double Ezcutoff = 0.5 * 3 * 3;        // const. term in partition function, z-only
  constexpr double Eztcutoff = 0.5 * (3*3 + 3*3); // zt

  // initial guess for cluster timing at high T

  // reset sums
  for(unsigned int k =0 ; k< nvtx; k++){
    psum[k] = 0;
    tsum[k] = 0;
    wsum[k] = 0;
    vtx_rhot[k] = vtx_rhot[k] / ntrk;
  }


  for(unsigned int i = 0; i < ntrk; i++){
    auto  trk = tracks[timing_tracklist[i]];
    double tErr = trk.MTD_timeerror();
    double tErrSq = tErr * tErr;
    double tbar = 0;
    for(unsigned int j=0; j < 3; j++){ tbar += a[j] * trk.th[j];  }

    // track assignment probability denominator (partition function)
    // at beta -> 0 this is z-only!
    double Ztrk = exp(-Ezcutoff); // cut-off 
    for(unsigned int k = 0; k < nvtx; k++){
      double Ez = 0.5 * pow((track_z[i] - vtx_z[k]) / track_z_error[i], 2); 
      if( Ez < Ezmax) { Ztrk += exp(-Ez) * vtx_rhot[k]; } 
    }

    // use Ztrk for the the assignment probabilities now
    if(Ztrk > exp(-Ezmax)){ // always true if Ezmax > Ezcutoff
      for(unsigned int k = 0; k < nvtx; k++){
	double Ez = 0.5 * pow((track_z[i] - vtx_z[k]) / track_z_error[i], 2); 
	if( Ez < Ezmax){
	  double p = vtx_rhot[k] * exp(-Ez) / Ztrk; 
	  double w = p / tErrSq; 
	  psum[k] += p;
	  wsum[k] += w;
	  tsum[k] += w * tbar;
	}
      }
    }
  } // track loop

  // evaluate cluster time estimates
  for(unsigned int k=0; k< nvtx; k++){
    vtx_rhot[k] =  psum[k] / ntrk;

    if(wsum[k] > 0){
      vtx_t[k] = tsum[k] / wsum[k];
    }else{
      vtx_t[k] = 0;
    }
      
    if(multi_vertex_dump_counter_ < 1){
      if(wsum[k] > 0 ){
	std::cout << "initial guess for vertex " << fixed << setw(4) << k << " at "  
		  << "  z=" << fixed << setprecision(4) << setw(8) << vtx_z[k] 
		  << "  t=" << fixed << setprecision(3) << setw(8) << vtx_t[k]
		  << "  rhot=" << fixed << setprecision(6) << setw(8) << vtx_rhot[k]
		  << "  t[org]=" << fixed << setprecision(3) << setw(7) << input_vertices->at(k).t()
		  << endl;
      }else{
	std::cout << "mass_constrained_multi_vertex_time_from_tracks_pid : no initial guess for vertex " 
		  << k << " at "  << vtx_z[k]   <<  "  rho_t="  << setprecision(6) << vtx_rhot[k] 
		  << "  t[org]=" << fixed << setprecision(3) << setw(7) << input_vertices->at(k).t()
		  << endl;
      }
    }

  }


  // annealing loop

  int nit = 0;
  double beta = 1. / 256.;
  double tres_sq = 0;
  while ( (nit++) < 200)  // sometimes converge veeery slowly
    {

      double Z0 = exp(-Eztcutoff);

      // first pass, partition function, Z[i]
      for(unsigned int i = 0; i < ntrk; i++){
	auto  trk = tracks[timing_tracklist[i]];
	double tErr = trk.MTD_timeerror();

	// probability denominator
	Z[i] =  Z0;
	for(unsigned int k = 0; k < nvtx; k++){
	  double Ez = 0.5 * pow((vtx_z[k] - track_z[i]) / track_z_error[i], 2); 
	  if( Ez > Ezmax) continue;
	  for(unsigned int j = 0; j < 3; j++){
	    double thpull =  (trk.th[j] - vtx_t[k]) / tErr;
	    double Eth = 0.5 * thpull * thpull;
	    if (fabs(beta * Eth) < Etmax){
	      Z[i] += a[j] * exp(- beta * Eth - Ez) * vtx_rhot[k];
	    }
	  }
	}
      }// 1st track loop


      // reset cluster sums
      for(unsigned int k =0 ; k< nvtx; k++){
	tsum[k] = 0;
	wsum[k] = 0;
	psum[k] = 0;
	w2sum[k] = 0;
	w0sum[k] = 0;
      }

      // second track loop, assignment and cluster sums
      for(unsigned int i=0; i< ntrk; i++){
	track_vtx_assignment_prob[i] = 0.;
	track_vtx_assignment[i] = 0;

	auto  trk = tracks[timing_tracklist[i]];
	double tErr = trk.MTD_timeerror();
	double tErrSq = tErr * tErr;

	for(unsigned int k = 0; k < nvtx; k++){
	  double Ez = 0.5 * pow((vtx_z[k] - track_z[i]) / track_z_error[i], 2); 
	  if( Ez > Ezmax) continue;

	  double tk_vtx_prob = 0; 
	  double wsum_trk = 0;
	  double s11=0, s21=0;
	  for(unsigned int j = 0; j < 3; j++){
	    double thpull =  (trk.th[j] - vtx_t[k]) / tErr;
	    double Eth = 0.5 * thpull * thpull;
	    if(fabs(beta * Eth) < Etmax){
	      double probh = a[j] * exp(- beta * Eth - Ez) / Z[i] * vtx_rhot[k];
	      psum[k] += probh;
	      tk_vtx_prob = max(tk_vtx_prob,  probh);
	      double w = probh / tErrSq;
	      // for the vertex time
	      tsum[k] += w * trk.th[j];
	      wsum[k] += w;
	      // error propagation
	      wsum_trk += w;
	      s11 +=  (trk.th[j] - vtx_t[k])  * w;// new
	      s21 +=  pow(trk.th[j] - vtx_t[k], 2)  * w;
	    }
	  }
	  //w2sum[k] += wsum_trk * wsum_trk * (dt * dt);  // naive error propagation 
	  double dfodt = wsum_trk + beta * (s11*s11 - s21/tErrSq);
	  w2sum[k] += pow(dfodt * tErr, 2);
	  w0sum[k] += dfodt;

	  if(tk_vtx_prob > track_vtx_assignment_prob[i]){
	    track_vtx_assignment_prob[i] = tk_vtx_prob;
	    track_vtx_assignment[i] = k;
	  }
	}// cluster k
      }//track i


      // update cluster times
      tres_sq = 0;
      double sum_rhot = 0.;
      for(unsigned int k = 0; k < nvtx; k++){
	vtx_rhot[k] = psum[k] / ntrk;
	if (wsum[k] < 1e-10) {
	  // essentially no timing for this cluster
	  //cout << "no timing " << k << "  wsum="  << scientific << wsum[k] << "  rho=" << scientific << vtx_rhot[k] << " t[org]=" << input_vertices->at(k).t() << endl;
	  vtx_t[k] = input_vertices->at(k).t();
	  vtx_tError[k] = 9.999;
	}else{
	  double t_new = tsum[k] / wsum[k];
	  tres_sq += pow(t_new - vtx_t[k], 2);
	  vtx_t[k] = t_new;
	  vtx_tError[k] = sqrt(w2sum[k]) / w0sum[k];
	  sum_rhot += vtx_rhot[k];
	}
      }
      //cout << "DEBUG sum_rhot=" << sum_rhot << endl;
      
      // at final temperature and converged?
      if ((tres_sq < 1e-8 / beta) && (beta >= 1.)){
	if(multi_vertex_dump_counter_ < 20){
	  cout << "multi converged" 
	       << " iteration " <<  nit
	       << " sum_rho = " << scientific << sum_rhot << fixed
	       << " tres = " << scientific << tres_sq << fixed
	       << endl;
	}
	

	if(multi_vertex_dump_counter_ < 1){
	  for(unsigned int k = 0; k < nvtx; k++){
	    cout << fixed << setw(3) << k << " ) " 
		 << "  z = " << fixed << setw(8) << setprecision(4) << vtx_z[k]
		 << "  t = " << fixed << setw(7) << setprecision(3) << vtx_t[k]
		 << "  +/- " << fixed << setw(7) << setprecision(3) << vtx_tError[k]
		 << "  rhot" << fixed << setw(8) << setprecision(5) << vtx_rhot[k]
		 << "       t[org]=" << fixed << setw(7) << setprecision(3) << input_vertices->at(k).t() 
		 << "  +/- " << fixed << setw(7) << setprecision(3) << input_vertices->at(k).tError()
	       << endl;
	  }
	}

	if(multi_vertex_dump_counter_ < 1){
	  for(unsigned int i = 0; i< ntrk; i++){
	    if( (track_vtx_assignment[i] != timing_track_vtx_input[i])
	      && (track_vtx_assignment_prob[i] > 0.01)) {
	      auto trk = tracks[timing_tracklist[i]];
	      cout << " track "  << setw(4) << i 
		   << "  z= "   << fixed << setw(8) << setprecision(4) << track_z[i] 
		   << "  +/-"  << fixed << setw(7) << setprecision(4) << track_z_error[i] 
		   << "  t = "   << fixed << setw(7) << setprecision(3) << trk.th[0] << " " << setw(7) << setprecision(3)<< trk.th[1] << " " << setw(7) << setprecision(3) << trk.th[2]
		   << "  +/-"  << fixed << setw(7) << setprecision(3) << trk.MTD_timeerror()
		   << "   changed owner " 
		   <<  "  from " << setw(4) << timing_track_vtx_input[i] << ", z,t="  
		   << fixed << setw(8) << setprecision(4) << vtx_z[timing_track_vtx_input[i]] 
		   <<  fixed << setw(7) << setprecision(3) << vtx_t[timing_track_vtx_input[i]] 
		   <<  "  to " << setw(4) << track_vtx_assignment[i]  <<", z,t="  
		 << fixed << setw(8) << setprecision(4) << vtx_z[track_vtx_assignment[i]]
		   <<  fixed << setw(7) << setprecision(3) << vtx_t[track_vtx_assignment[i]] 
		   << "   p= "  << track_vtx_assignment_prob[i];
	      if(trk.is_matched()){
		cout << "  zsim= " << fixed << setw(8) << setprecision(4) << trk.zsim() 
		     << "  tsim= " <<  fixed << setw(7) << setprecision(3) << trk.tsim();
	      }
	      if(trk.is_matched() && !(trk.is_primary())){
		cout << "  secondary";
	      }
	      cout << endl;

	    }
	  }
	}

	// refit (z)
	AdaptiveVertexFitter  theFitter(GeometricAnnealing(2.5)); 
	auto new_collection = new reco::VertexCollection;

	for(unsigned int k = 0; k < nvtx; k++){
	  // select the tracks =  timing tracks + non-timing tracks
	  vector<TransientTrack> vtxtracks;

	  // non-timing, use existing assignment
	  auto & v_input = input_vertices->at(k);
	  for(auto tk = v_input.tracks_begin(); tk != v_input.tracks_end(); tk++){
	    const reco::TrackRef trackRef = tk->castTo<reco::TrackRef>();
	    unsigned int itrk = tracks.index_from_ref(*tk);
	    if (find(timing_tracklist.begin(), timing_tracklist.end(), itrk) == timing_tracklist.end()){
	      vtxtracks.push_back(tracks[itrk].transientTrack());
	    }
	  }

	  // timing tracks, possibly re-assigned
	  for(unsigned int i = 0; i< ntrk; i++){
	    if ((track_vtx_assignment[i] == k) && (track_vtx_assignment_prob[i] > 0.4)){
	      vtxtracks.push_back(tracks[timing_tracklist[i]].transientTrack());
	    }
	  }


	  // redo the 3d-fit with the revised tracklist
	  TransientVertex v = theFitter.vertex(vtxtracks);
	  if (v.isValid()){
	    // make it a 4d vertex, use the time from the multi-vertex fit
	    auto err = v.positionError().matrix4D();
	    err(3, 3) = pow(vtx_tError[k], 2);
	    auto trkWeightMap3d = v.weightMap();  // copy the 3d-fit weights
	    v = TransientVertex(v.position(), vtx_t[k], err, 
				v.originalTracks(), v.totalChiSquared(), v.degreesOfFreedom());
	    v.weightMap(trkWeightMap3d);
	    reco::Vertex new_recv = v;
	    new_collection->push_back(new_recv);
	  }else{
	    // drop this vertex
	    cout << "mass_constrained_multi_vertex_time_from_tracks_pid : refit failed for vertex " <<  k << endl;
	  }
      	}

	delete recVtxs_[label];
	recVtxs_[label] = new_collection;
	return;
      }
      
      // converged at intermediate temperature, cool down
      if ( (tres_sq < 1e6) && ((beta < 1.))){
	beta = std::min(1., beta / cooling_factor);
      }

    }// while
 
  cout << "mass_constrained_multi_vertex_time_from_tracks_pid did not converge,   nit=" << nit << "  tres=" << scientific << tres_sq << endl;
}





/* this is the version without cluster masses, probably obsolete !!! */
void PrimaryVertexAnalyzer4PU::multi_vertex_time_from_tracks_pid(const std::string label,
							   Tracks& tracks,
							   double minquality,
							   bool verbose) {

  multi_vertex_dump_counter_++;
  // work on this collection
  auto & input_vertices = recVtxs_[label];


  // cluster properties and clusters cluster-sums
  unsigned int nvtx = input_vertices->size();
  std::vector<double> tsum(nvtx);
  std::vector<double> wsum(nvtx);
  std::vector<double> w2sum(nvtx);
  std::vector<double> vtx_t(nvtx);
  std::vector<double> vtx_tError(nvtx);
  std::vector<double> vtx_z(nvtx);     // input only



  // build a track table for the relevant tracks with timing
  std::vector<unsigned int> timing_tracklist;
  std::vector<unsigned int> timing_track_vtx_input;
  for(unsigned int k =0 ; k< nvtx; k++){
    auto & v = input_vertices->at(k);
    vtx_z[k] = v.z();
    for (auto tk = v.tracks_begin(); tk != v.tracks_end(); tk++) {
      unsigned int itrk = tracks.index_from_ref(*tk);
      if (tracks[itrk].has_timing() && (tracks[itrk].timeQuality() >= minquality)){
	if (find(timing_tracklist.begin(), timing_tracklist.end(), itrk) ==  timing_tracklist.end()){
	  timing_tracklist.push_back(itrk);
	  timing_track_vtx_input.push_back(k);
	}else{
	  // should not happen currently
	  cout << "track assigned to multiple vertices ?" << endl;
	}
      }
    }
  }


  // prepare arrays, pre-fill track z and z-error
  unsigned int ntrk = timing_tracklist.size();
  //stable_sort(timing_tracklist.begin(), timing_tracklist(end)); don't really care
  std::vector<double> Z(ntrk);
  std::vector<double> track_z_error(ntrk);
  std::vector<double> track_z(ntrk);
  std::vector<double> track_vtx_assignment_prob(ntrk);
  std::vector<unsigned int> track_vtx_assignment(ntrk);

  for(unsigned int i = 0; i < ntrk; i++){
    auto trk = tracks[timing_tracklist[i]];
    auto k = timing_track_vtx_input[i];
    track_z[i] = trk.z();
    track_z_error[i] = sqrt(pow(trk.dzError(), 2) + pow(0.0020 / tan(trk.theta()),2) + pow(0.0020,2));
    if(multi_vertex_dump_counter_ < 1){
      cout << fixed << setw(4) << setprecision(0) << i <<  ")"
	   << " k= " << fixed << setw(3) << setprecision(0) << k
	   << "  zvtx = "  << fixed << setw(8) << setprecision(4) << vtx_z[k]
	   << "  ztrk = "   << fixed << setw(8) << setprecision(4) << track_z[i] 
	   << " +/-"  << fixed << setw(7) << setprecision(4) << track_z_error[i] 
	   << "  zpull = " << fixed << setw(5) << setprecision(2) << (track_z[i]-vtx_z[k])/track_z_error[i]
	   << "     ttrk[pi/K/p] = "   << fixed << setw(7) << setprecision(3) << trk.th[0] << " " << setw(7) << setprecision(3)<< trk.th[1] << " " << setw(7) << setprecision(3) << trk.th[2]
	   << "  +/-"  << fixed << setw(7) << setprecision(3) << trk.MTD_timeerror();
      if(trk.is_matched()){
	cout << "  zsim = " <<  fixed << setw(8) << setprecision(4) << trk.zsim()
	     << "  tsim = " <<  fixed << setw(7) << setprecision(3) << trk.tsim()
	     << "  m=" << trk.get_particle_mass();
      }
      cout << endl;
    }
  }


  // constants
  double a[3] = {0.7,0.2,0.1};
  constexpr double cooling_factor = 0.5;


  // initial guess for cluster timing at high T

  // reset sums
  for(unsigned int k =0 ; k< nvtx; k++){
    tsum[k] = 0;
    wsum[k] = 0;
  }


  for(unsigned int i = 0; i < ntrk; i++){
    auto  trk = tracks[timing_tracklist[i]];
    double dt = trk.MTD_timeerror();
    double tbar = 0;
    for(unsigned int j=0; j < 3; j++){ tbar += a[j] * trk.th[j];  }

    // track assignment probability denominator (partition function)
    // at beta -> 0 this is z-only!
    double Ztrk = exp(-0.5*100.); // cut-off 
    for(unsigned int k = 0; k < nvtx; k++){
      double Ez = 0.5 * pow((track_z[i] - vtx_z[k]) / track_z_error[i], 2); 
      if( Ez < 50.) { Ztrk += exp(-Ez); } 
    }

    // use Ztrk for the the assignment probabilities now
    if(Ztrk > exp(-0.5*100.)){
      for(unsigned int k =0; k < nvtx; k++){
	double Ez = 0.5 * pow((track_z[i] - vtx_z[k]) / track_z_error[i], 2); 
	if( Ez < 50.){
	  double w = exp(-Ez) / Ztrk / (dt*dt); 
	  wsum[k] += w;
	  tsum[k] += w * tbar;
	}
      }
    }
  } // track loop

  // evaluate cluster time estimates
  for(unsigned int k=0; k< nvtx; k++){
    if(wsum[k] > 0){
      vtx_t[k] = tsum[k] / wsum[k];
      if(multi_vertex_dump_counter_ < 1){
      std::cout << "initial guess for vertex " << fixed << setw(4) << k << " at "  
		<< "  z=" << fixed << setprecision(4) << setw(8) << vtx_z[k] 
		<< "  t=" << fixed << setprecision(3) << setw(8) << vtx_t[k]
		<< endl;
      }
    }else{
      std::cout << "multi_vertex_time_from_tracks_pid : no initial guess for vertex " << k << " at "  << vtx_z[k] << endl;
      vtx_t[k] = 0.;
    }
  }


  // annealing loop

  int nit = 0;
  double beta = 1. / 256.;
  while ( (nit++) < 100)
    {

      double Z0 = exp(-0.5 * 5. * 5. - 0.5 * beta * 3. * 3. );

      // first round, determine the partition function
      for(unsigned int i = 0; i < ntrk; i++){
	auto  trk = tracks[timing_tracklist[i]];
	double dt = trk.MTD_timeerror();

	// probability denominator
	Z[i] =  Z0;
	for(unsigned int k = 0; k < nvtx; k++){
	  double Ez = 0.5 * pow((vtx_z[k] - track_z[i]) / track_z_error[i], 2); 
	  if( Ez > 50.) continue;
	  for(unsigned int j = 0; j < 3; j++){
	    double thpull =  (trk.th[j] - vtx_t[k]) / dt;
	    double Eth = 0.5 * thpull * thpull;
	    if (fabs(beta * Eth) < 100){
	      Z[i] += a[j] * exp(- beta * Eth - Ez);
	    }
	  }
	}
      }// 1st track loop

      // reset cluster sums for vertex updates in the 2nd track loop
      for(unsigned int k =0 ; k< nvtx; k++){
	tsum[k] = 0;
	wsum[k] = 0;
	w2sum[k] = 0;
      }

      // second track loop, assignment and cluster sums
      for(unsigned int i=0; i< ntrk; i++){
	track_vtx_assignment_prob[i] = 0.;
	track_vtx_assignment[i] = 0;

	auto  trk = tracks[timing_tracklist[i]];
	double dt = trk.MTD_timeerror();

	for(unsigned int k = 0; k < nvtx; k++){
	  double Ez = 0.5 * pow((vtx_z[k] - track_z[i]) / track_z_error[i], 2); 
	  if( Ez > 50.) continue;

	  double tk_vtx_prob = 0; 
	  double wsum_trk = 0;
	  for(unsigned int j = 0; j < 3; j++){
	    double thpull =  (trk.th[j] - vtx_t[k]) / dt;
	    double Eth = 0.5 * thpull * thpull;
	    if(fabs(beta * Eth) < 100){
	      double probh = a[j] * exp(- beta * Eth - Ez) / Z[i];
	      tk_vtx_prob = max(tk_vtx_prob,  probh);
	      double w = probh / (dt * dt);
	      wsum_trk += w;
	      tsum[k] += w * trk.th[j]; 
	      if( (beta>0.99) && (i==32) && ((k==0) || (k==128))){
		cout << " i=" << i << " k=" << k << "  j=" << j << "  prob="  << probh  << "   zpull=" << (vtx_z[k] - track_z[i]) / track_z_error[i] << "  tpull = " << thpull << " nit=" << nit  << endl;
	      }
	    }
	  }
	  wsum[k] += wsum_trk;
	  w2sum[k] += wsum_trk * wsum_trk * (dt * dt); 

	  if(tk_vtx_prob > track_vtx_assignment_prob[i]){
	    track_vtx_assignment_prob[i] = tk_vtx_prob;
	    track_vtx_assignment[i] = k;
	  }
	  if( (beta>0.99) && (i==32) && ((k==0) || (k==128))){
	    cout << " i=" << i << " k=" << k << "  tk_vtx_prob "  << tk_vtx_prob << endl;
	  }
	}// cluster k
	if( (beta>0.99) && (i==32)){
	  cout << " i=" << i << " track_vtx_assignment=" <<  track_vtx_assignment[i] << "  tk_vtx_assignment_prob = "  <<  track_vtx_assignment_prob[i]<< endl;
	}
      }//track i


      // update cluster times
      double tres = 0;
      for(unsigned int k = 0; k < nvtx; k++){
	if (wsum[k] < 1e-10) {
	  // essentially no timing for this cluster
	  vtx_tError[k] = 9.999;
	}else{
	  double t_new = tsum[k] / wsum[k];
	  tres += pow(t_new - vtx_t[k], 2);
	  vtx_t[k] = t_new;
	  vtx_tError[k] = sqrt(w2sum[k]) / wsum[k];
	}
      }
      
      // at final temperature and converged?
      if ((tres < 1e-8 / beta) && (beta >= 1.)){
	if(multi_vertex_dump_counter_ < 20){
	cout << "multi converged" 
	     << " iteration " <<  nit
	     << " tres = " << scientific << tres << fixed
	     << endl;
	}

	if(multi_vertex_dump_counter_ < 1){
	for(unsigned int k = 0; k < nvtx; k++){
	  cout << fixed << setw(3) << k << " ) " 
	       << "  z = " << fixed << setw(8) << setprecision(4) << vtx_z[k]
	       << "  t = " << fixed << setw(7) << setprecision(3) << vtx_t[k]
	       << "  +/- " << fixed << setw(7) << setprecision(3) << vtx_tError[k]
	       << endl;
	}
	}

	if(multi_vertex_dump_counter_ < 1){
	for(unsigned int i = 0; i< ntrk; i++){
	  if(track_vtx_assignment[i] != timing_track_vtx_input[i]){
	    auto trk = tracks[timing_tracklist[i]];
	    cout << " track "  << setw(4) << i 
		 << "  z= "   << fixed << setw(8) << setprecision(4) << track_z[i] 
		 << "  +/-"  << fixed << setw(7) << setprecision(4) << track_z_error[i] 
		 << "  t = "   << fixed << setw(7) << setprecision(3) << trk.th[0] << " " << setw(7) << setprecision(3)<< trk.th[1] << " " << setw(7) << setprecision(3) << trk.th[2]
		 << "  +/-"  << fixed << setw(7) << setprecision(3) << trk.MTD_timeerror()
		 << "   changed owner " 
		 <<  "  from " << setw(4) << timing_track_vtx_input[i] << ", z,t="  
		 << fixed << setw(8) << setprecision(4) << vtx_z[timing_track_vtx_input[i]] 
		 <<  fixed << setw(7) << setprecision(3) << vtx_t[timing_track_vtx_input[i]] 
		 <<  "  to " << setw(4) << track_vtx_assignment[i]  <<", z,t="  
		 << fixed << setw(8) << setprecision(4) << vtx_z[track_vtx_assignment[i]]
		 <<  fixed << setw(7) << setprecision(3) << vtx_t[track_vtx_assignment[i]] 
		 << "   p= "  << track_vtx_assignment_prob[i];
	    if(trk.is_matched()){
	      cout << "  zsim= " << fixed << setw(8) << setprecision(4) << trk.zsim() 
		   << "  tsim= " <<  fixed << setw(7) << setprecision(3) << trk.tsim();
	    }
	    if(trk.is_matched() && !(trk.is_primary())){
	      cout << "  secondary";
	    }
	    cout << endl;

	  }
	}
	}

	// refit
	AdaptiveVertexFitter  theFitter(GeometricAnnealing(2.5)); 
	auto new_collection = new reco::VertexCollection;

	for(unsigned int k = 0; k < nvtx; k++){
	  // select the tracks =  timing tracks + non-timing tracks
	  vector<TransientTrack> vtxtracks;

	  // non-timing, use existing assignment
	  auto & v_input = input_vertices->at(k);
	  for(auto tk = v_input.tracks_begin(); tk != v_input.tracks_end(); tk++){
	    const reco::TrackRef trackRef = tk->castTo<reco::TrackRef>();
	    unsigned int itrk = tracks.index_from_ref(*tk);
	    if (find(timing_tracklist.begin(), timing_tracklist.end(), itrk) == timing_tracklist.end()){
	      vtxtracks.push_back(tracks[itrk].transientTrack());
	    }
	  }

	  // timing tracks, possibly re-assigned
	  for(unsigned int i = 0; i< ntrk; i++){
	    if ((track_vtx_assignment[i] == k) && (track_vtx_assignment_prob[i] > 0.1)){
	      // track assigned by multi-vertex fit
	      vtxtracks.push_back(tracks[timing_tracklist[i]].transientTrack());
	    //}else if ((timing_track_vtx_input[i] == k) && (track_vtx_assignment_prob[i] == 0)){
	      // track was not re-assigned(?)
	      //  vtxtracks.push_back(tracks[timing_tracklist[i]].transientTrack());
	    }
	  }

	  // FIXME : remove the vertex time when no timing track is part of the fit
	  // fit
	  TransientVertex v = theFitter.vertex(vtxtracks);
	  if (v.isValid()){
	    // make it a 4d vertex, use the time from the multi-vertex fit
	    auto err = v.positionError().matrix4D();
	    err(3, 3) = pow(vtx_tError[k], 2);
	    auto trkWeightMap3d = v.weightMap();  // copy the 3d-fit weights
	    v = TransientVertex(v.position(), vtx_t[k], err, 
				v.originalTracks(), v.totalChiSquared(), v.degreesOfFreedom());
	    v.weightMap(trkWeightMap3d);
	    reco::Vertex new_recv = v;
	    new_collection->push_back(new_recv);
	  }else{
	    // drop this vertex
	    cout << "refit failed for vertex " <<  k << endl;
	  }
      	}

	delete recVtxs_[label];
	recVtxs_[label] = new_collection;
	return;
      }
      
      // converged at intermediate temperature, cool down
      if ( (tres < 1e6) && ((beta < 1.))){
	beta = std::min(1., beta / cooling_factor);
      }

    }// while
 
  cout << "multi_vertex_time_from_tracks_pid did not converge,   nit= " << nit << endl;
}



double PrimaryVertexAnalyzer4PU::vertex_ptmax2(const reco::Vertex& v) {
  double ptmax1 = 0;
  double ptmax2 = 0;

  for (trackit_t t = v.tracks_begin(); t != v.tracks_end(); t++) {
    if (v.trackWeight(*t) > min_trk_in_vtx_weight_) {
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

bool PrimaryVertexAnalyzer4PU::select(const MVertex& v, int level) {
  /* level
     0  !isFake  && ndof>selNdof_  (default)
     1  !isFake  && ndof>selNdof_ && prob > 0.01 
     2  !isFake  && ndof>selNdof_ && prob > 0.01 && ptmax2 > 0.4
  */
  if (v.isRecoFake())
    return false;
  if ((level == 0) && (v.ndof() >= selNdof_))
    return true;
  if ((level == 1) && (v.ndof() >= selNdof_) && (v.pxy(vertexBeamSpot_) > 0.01))
    return true;
  if ((level == 2) && (v.ndof() >= selNdof_) && (v.pxy(vertexBeamSpot_) > 0.01) && (v.ptmax2() > 0.4))
    return true;
  if ((level == 3) && (v.ndof() >= selNdof_) && (v.ptmax2() < 0.4))
    return true;
  return false;
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
							const unsigned int index,
                                                        const double deltaz,
                                                        const bool verbose)
// fill vertex quantities available from the reco vertex itself
{
  if(v->isFake()) return;
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

  if (index != NO_INDEX) {
    Fill(h, vtype + "/index", float(index));
    Fill(h, vtype + "/logndofvsindex", float(index), log(v->ndof()) / log(10.));
    Fill(h, vtype + "/ndofvsindex", float(index), v->ndof());
  }
  Fill(h, vtype + "/c2xy", c2xy);
  Fill(h, vtype + "/c2xyvsntrk", v->tracksSize(), c2xy);
  Fill(h, vtype + "/chi2", v->chi2());
  if (v->tracksSize() > 0)
    {
      Fill(h, vtype + "/chi2overntk", v->chi2() / v->tracksSize());
    }
  Fill(h, vtype + "/probxy", TMath::Prob(c2xy, 2));
  Fill(h, vtype + "/r", sqrt(dx * dx + dy * dy));
  Fill(h, vtype + "/xrecbeam", dx);
  Fill(h, vtype + "/yrecbeam", dy);
  Fill(h, vtype + "/xyrecbeam", dx, dy);
  Fill(h, vtype + "/zpullbeam", (z - vertexBeamSpot_.z0() - dzb_) / sigmaZ_);
  Fill(h, vtype + "/ndof", v->ndof());
  if((v->ndof()==0) && (vtype=="recvtx")){
    // FIXME, this should not happen
      std::cout << "fillVertexHistosNoTracks  ndof= "<< v-> ndof() << "  z=" << v->z() << "  chisq =" <<  v-> chi2() << "  index=" <<  index << std::endl;
  }
  Fill(h, vtype + "/ndofvspu", lumiPU_, v->ndof());

  if (v->ndof() > 0) {
    Fill(h, vtype + "/logndof", log(v->ndof()) / log(10.));
  }
  Fill(h, vtype + "/numtrk", v->tracksSize());
  Fill(h, vtype + "/vtxndfvsntrk", v->tracksSize(), v->ndof());
  Fill(h, vtype + "/avweight", (v->ndof() - ndof0trk_) / (2. * v->tracksSize()));
  Fill(h, vtype + "/avweightvsndof", (v->ndof() - ndof0trk_) / (2. * v->tracksSize()), v->ndof());
  Fill(h, vtype + "/avweightX", (v->ndof() - ndof0trk_) / (2. * v->tracksSize()), v->ndof()); // avweight weighted with ndof
  Fill(h, vtype + "/errx", v->xError());
  Fill(h, vtype + "/erry", v->yError());
  Fill(h, vtype + "/errz", v->zError());
  Fill(h, vtype + "/logerrz", log(v->zError()/1e-4)/log(10.));

  Fill(h, vtype + "/zvtx", z);
  if (f4D_) {
    // contains vertices that don't actually have timing information
    Fill(h, vtype + "/terrvtx", v->tError());
    if((v->tError() < 0.1) && (v->tError() > 0.0)){
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



void PrimaryVertexAnalyzer4PU::fillRecoVertexHistos(std::map<std::string, TH1*>& h,
                                                const std::string& vtype,
                                                const reco::Vertex* v,
                                                Tracks& tracks,
						const unsigned int index,
                                                const double deltaz,
                                                const bool verbose) {
  // delta z = z_this - z_other
  timer_start("fillRecoVertexHistos");
  fillVertexHistosNoTracks(h, vtype, v, index, deltaz, verbose);
  Fill(h, vtype + "/ptmax2", vertex_ptmax2(*v));
  double sumpt2 = vertex_sumpt2(*v);
  double sumpt = vertex_sumpt(*v);  // not weighted, unlike sumapt
  Fill(h, vtype + "/sumpt", sumpt);
  Fill(h, vtype + "/sumpt2", sumpt2);
  Fill(h, vtype + "/logsumpt2", log(sumpt2)/log(10.));
  Fill(h, vtype + "/logsumpt", log(sumpt)/log(10.));
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

      if (v->trackWeight(*t) > min_trk_in_vtx_weight_) {
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
  timer_stop("fillRecoVertexHistos");
}
/****************************************************************************/




/****************************************************************************/
void PrimaryVertexAnalyzer4PU::fillVertexHistos(std::map<std::string, TH1*>& h,
                                                const std::string& vtype,
                                                const MVertex & v,
                                                Tracks& tracks,
                                                const double deltaz,
                                                const bool verbose)
/****************************************************************************/
{
  if(v.isRecoFake()) return;
  // delta z = z_this - z_other
  timer_start("fillVertexHistos");
  fillVertexHistosNoTracks(h, vtype, v.recvtx, v.index(), deltaz, verbose);
  Fill(h, vtype + "/ptmax2", v.ptmax2());
  double sumpt2 = v.sumpt2();
  double sumpt = v.sumpt();  // not weighted, unlike sumapt
  Fill(h, vtype + "/sumpt", sumpt);   // code duplication, FIXME
  Fill(h, vtype + "/sumpt2", sumpt2);
  Fill(h, vtype + "/logsumpt", log(sumpt)/log(10.));
  Fill(h, vtype + "/logsumpt2", log(sumpt2)/log(10.));
  Fill(h, vtype + "/sumpt2vssumpt", sumpt, sumpt2);
  if(sumpt >0){
    Fill(h, vtype + "/sumpt2oversumpt", sumpt2/sumpt);
    Fill(h, vtype + "/sumpt2oversumptvssumpt2", sumpt2, sumpt2/sumpt);
  }
  Fill(h, vtype + "/nseltrkvtx", float(v.tracksSize()));



  int nlowt = 0;
  for( auto  tk : v.tracks ){
    fillTrackHistos(h, vtype, *tk, v.recvtx);   // select based on weight here?

    Fill(h, vtype + "/trkweight", v.trackWeight(tk));
    Fill(h, vtype + "/trkptnorm", tk->pt(), 1./  v.tracksSize());
    if ( v.trackWeight(tk) < 0.6) nlowt += 1;
  }
  Fill(h, vtype + "/numlowttrk", float(nlowt) );
  if (v.tracksSize() > 0) Fill(h, vtype + "/fraclowttrk", float(nlowt) / v.tracksSize());

  /* not ported yet , todo 
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
  double l1 = (Sxx + Syy + Ds) / (2 * Spt);x
  double l2 = (Sxx + Syy - Ds) / (2 * Spt);
  double SLT = 2 * l2 / (l1 + l2);
  Fill(h, vtype + "/sphericity", SLT);
  Fill(h, vtype + "/sphericityvsntrk", Sw, SLT);
    */
  timer_stop("fillVertexHistos");
}
/****************************************************************************/




/*************************************************************************************/
void PrimaryVertexAnalyzer4PU::fillVertexHistosMatched(std::map<std::string, TH1*>& h,
						       const std::string& vtype,
						       MVertex & v,
						       Tracks& tracks,
						       const std::vector<SimEvent>& simEvt,
						       const double deltaz,
						       const bool verbose)
/*************************************************************************************/
{
  fillVertexHistos(h, vtype, v, tracks, deltaz, verbose);

  // now fill properties wrt the matched simvertex
  assert(v.sim != NOT_MATCHED_VTX_SIM);
  const SimEvent& simevt = simEvt[v.sim];

  unsigned int ntiming = 0;
  unsigned int ntiming_qual05 = 0;
  unsigned int ntiming_qual08 = 0;
  unsigned int nseltrk = 0;
  for( auto tk : v.tracks){
    nseltrk++;
    if(f4D_){
      if (tk->has_timing()) {
	ntiming++;
	if (tk->timeQuality() > 0.5) ntiming_qual05++;
	if (tk->timeQuality() > 0.8) ntiming_qual08++;
      }
    }
  }
  // these histograms are also filled for unmatched vertices, the false flag prevents double bookings
  Fill(h, vtype + "/nseltrkvtx", float(nseltrk), simevt.is_signal(), false);
  Fill(h, vtype + "/ptmax2", v.ptmax2(), simevt.is_signal(), false);

  // note that this is redundant with some histos filled in analyzevertexcollectiontp , FIXME  consolidate

  double xsim = simevt.x;
  double ysim = simevt.y;
  double zsim = simevt.z;
  Fill(h, vtype + "/xrecsim", v.x() - xsim, simevt.is_signal());
  Fill(h, vtype + "/yrecsim", v.y() - ysim, simevt.is_signal());
  Fill(h, vtype + "/zrecsim", v.z() - zsim, simevt.is_signal());
  Fill(h, vtype + "/xrecerr", v.xError(), simevt.is_signal());  // like xerr but separated by signal/pu
  Fill(h, vtype + "/yrecerr", v.yError(), simevt.is_signal());
  Fill(h, vtype + "/zrecerr", v.zError(), simevt.is_signal());
  Fill(h, vtype + "/xrecsimHR", v.x() - xsim, simevt.is_signal());
  Fill(h, vtype + "/yrecsimHR", v.y() - ysim, simevt.is_signal());
  Fill(h, vtype + "/zrecsimHR", v.z() - zsim, simevt.is_signal());
  Fill(h, vtype + "/xrecsimpull", (v.x() - xsim) / v.xError(), simevt.is_signal());
  Fill(h, vtype + "/yrecsimpull", (v.y() - ysim) / v.yError(), simevt.is_signal());
  Fill(h, vtype + "/zrecsimpull", (v.z() - zsim) / v.zError(), simevt.is_signal());
  Fill(h, vtype + "/xrecsimpullvsntrk", float(nseltrk), (v.x() - xsim) / v.xError(), simevt.is_signal());
  Fill(h, vtype + "/yrecsimpullvsntrk", float(nseltrk), (v.y() - ysim) / v.yError(), simevt.is_signal());
  Fill(h, vtype + "/zrecsimpullvsntrk", float(nseltrk), (v.z() - zsim) / v.zError(), simevt.is_signal());
  Fill(h, vtype + "/xrecsimvsxsimbeam", xsim - vertexBeamSpot_.x(zsim), v.x() - xsim, simevt.is_signal());
  Fill(h, vtype + "/yrecsimvsysimbeam", ysim - vertexBeamSpot_.y(zsim), v.y() - ysim, simevt.is_signal());
  if (v. tracksSize() < 5){
    Fill(h, vtype + "/xrecsimvsxsimbeamntklt5", xsim - vertexBeamSpot_.x(zsim), v.x() - xsim, simevt.is_signal());
    Fill(h, vtype + "/yrecsimvsysimbeamntklt5", ysim - vertexBeamSpot_.y(zsim), v.y() - ysim, simevt.is_signal());
  }else{
    Fill(h, vtype + "/xrecsimvsxsimbeamntkge5", xsim - vertexBeamSpot_.x(zsim), v.x() - xsim, simevt.is_signal());
    Fill(h, vtype + "/yrecsimvsysimbeamntkge5", ysim - vertexBeamSpot_.y(zsim), v.y() - ysim, simevt.is_signal());
  }

  double du=0, dv=0, uError=0, vError=0, pol=0;
  if (uv_analysis(v, simevt, du, dv, uError, vError, pol)){
    if(pol < 0.1){
      Fill(h, vtype + "/urecsim",   du, simevt.is_signal());
      Fill(h, vtype + "/vrecsim",   dv, simevt.is_signal());
      Fill(h, vtype + "/urecsimHR",   du, simevt.is_signal());
      Fill(h, vtype + "/vrecsimHR",   dv, simevt.is_signal());
      Fill(h, vtype + "/urecsimpull", du / uError, simevt.is_signal());
      Fill(h, vtype + "/vrecsimpull", dv / vError, simevt.is_signal());
      Fill(h, vtype + "/urecsimpullvsntrk", float(nseltrk), du / uError, simevt.is_signal());
      Fill(h, vtype + "/vrecsimpullvsntrk", float(nseltrk), dv / vError, simevt.is_signal());
    }
  }


  Fill(h, vtype + "/xrecsimvssumptprof", v.sumpt(), std::pow(v.x() - xsim,2), simevt.is_signal());
  Fill(h, vtype + "/xrecsimvssumpt", v.sumpt(), v.x() - xsim, simevt.is_signal());
  Fill(h, vtype + "/xrecsimvsntrkprof", float(nseltrk), std::pow(v.x() - xsim,2), simevt.is_signal());
  Fill(h, vtype + "/xrecsimvsntrk", float(nseltrk), v.x() - xsim, simevt.is_signal());
  Fill(h, vtype + "/xrecsimvsnsimtrkprof", float(simevt.rtk.size()), std::pow(v.x() - xsim,2), simevt.is_signal());
  Fill(h, vtype + "/xrecsimvsnsimtrk", float(simevt.rtk.size()), v.x() - xsim, simevt.is_signal());
  Fill(h, vtype + "/xrecsimvssumptsimprof", simevt.sumpt, std::pow(v.x() - xsim,2), simevt.is_signal());
  Fill(h, vtype + "/xrecsimvssumptsim", simevt.sumpt, v.x() - xsim, simevt.is_signal());

  Fill(h, vtype + "/yrecsimvssumptprof", v.sumpt(), std::pow(v.y() - ysim,2), simevt.is_signal());
  Fill(h, vtype + "/yrecsimvssumpt", v.sumpt(), v.y() - ysim, simevt.is_signal());
  Fill(h, vtype + "/yrecsimvsntrkprof", float(nseltrk), std::pow(v.y() - ysim,2), simevt.is_signal());
  Fill(h, vtype + "/yrecsimvsntrk", float(nseltrk), v.y() - ysim, simevt.is_signal());
  Fill(h, vtype + "/yrecsimvsnsimtrkprof", float(simevt.rtk.size()), std::pow(v.y() - ysim,2), simevt.is_signal());
  Fill(h, vtype + "/yrecsimvsnsimtrk", float(simevt.rtk.size()), v.y() - ysim, simevt.is_signal());
  Fill(h, vtype + "/yrecsimvssumptsimprof", simevt.sumpt, std::pow(v.y() - ysim,2), simevt.is_signal());
  Fill(h, vtype + "/yrecsimvssumptsim", simevt.sumpt, v.y() - ysim, simevt.is_signal());

  Fill(h, vtype + "/zrecsimvssumptprof", v.sumpt(), std::pow(v.z() - zsim,2), simevt.is_signal());
  Fill(h, vtype + "/zrecsimvssumpt", v.sumpt(), v.z() - zsim, simevt.is_signal());
  Fill(h, vtype + "/zrecsimvsntrkprof", float(nseltrk), std::pow(v.z() - zsim,2), simevt.is_signal());
  Fill(h, vtype + "/zrecsimvsntrk", float(nseltrk), v.z() - zsim, simevt.is_signal());
  Fill(h, vtype + "/zrecsimvsnsimtrkprof", float(simevt.rtk.size()), std::pow(v.z() - zsim,2), simevt.is_signal());
  Fill(h, vtype + "/zrecsimvsnsimtrk", float(simevt.rtk.size()), v.z() - zsim, simevt.is_signal());
  Fill(h, vtype + "/zrecsimvssumptsimprof", simevt.sumpt, std::pow(v.z() - zsim,2), simevt.is_signal());
  Fill(h, vtype + "/zrecsimvssumptsim", simevt.sumpt, v.z() - zsim, simevt.is_signal());


  
  if (simevt.type == FROM_TRACKING_TRUTH or simevt.type == FROM_PU_SUMMARY) {

    Fill(h, vtype + "/zsim", simevt.z * simUnit_, simevt.is_signal());
    Fill(h, vtype + "/dzminsim", simevt.dzmin * simUnit_, simevt.is_signal());
    if (f4D_) Fill(h, vtype + "/tsim", simevt.t, simevt.is_signal());

    if (simevt.type == FROM_TRACKING_TRUTH) {
      Fill(h, vtype + "/nbsimtksinvtx", float(simevt.nGenTrk), simevt.is_signal());
      Fill(h, vtype + "/xsim", simevt.x * simUnit_, simevt.is_signal());
      Fill(h, vtype + "/ysim", simevt.y * simUnit_, simevt.is_signal());

      Fill(h, vtype + "/logpthatsim", std::log10(simevt.pt_hat), simevt.is_signal());
      Fill(h, vtype + "/logsumptsim", std::log10(simevt.sumpt), simevt.is_signal());
      Fill(h, vtype + "/logsumpt2sim", std::log10(simevt.sumpt2), simevt.is_signal());
    }
  }

  if (f4D_)
    {
      Fill(h, vtype + "/ntimingvtx", float(ntiming), simevt.is_signal());
      Fill(h, vtype + "/ntimingqual05vtx", float(ntiming_qual05), simevt.is_signal());
      Fill(h, vtype + "/ntimingqual08vtx", float(ntiming_qual08), simevt.is_signal());

      vertex_time_from_tracks_analysis(vertex_time_from_tracks, v, tracks, 0., simevt, "", h, vtype);
      vertex_time_from_tracks_analysis(vertex_time_from_tracks, v, tracks, 0.8, simevt, "qual", h, vtype);
      vertex_time_from_tracks_analysis(vertex_time_from_tracks_pid, v, tracks, 0. , simevt, "pid", h, vtype);
      auto qresult = vertex_time_from_tracks_analysis(vertex_time_from_tracks_pid, v, tracks, 0.8, simevt, "qual_pid", h, vtype);
      /*
      auto nresult = vertex_time_from_tracks_analysis(vertex_time_from_tracks_pid_newton, v, tracks, 0.8, simevt, "qual_pid_new", h, vtype);
      if (qresult.successful() && !nresult.successful()){
	cout << "qual successful, newton failed "<< endl;
	cout << "qual : t=" << qresult.t() << " +/- " << qresult.tError() << "   nit =  " << qresult.niteration << endl;
	nresult = vertex_time_from_tracks_analysis(vertex_time_from_tracks_pid_newton, v, tracks, 0.8, simevt, "qual_pid_new", h, vtype, true);
      }
      */
      

      double tsim = simevt.t;
      if(v.has_timing()){
	Fill(h, vtype + "/trecsim", v.t() - tsim, simevt.is_signal());
	Fill(h, vtype + "/trecerr", v.tError(), simevt.is_signal());
	Fill(h, vtype + "/trecsimpull", (v.t() - tsim) / v.tError(), simevt.is_signal());
	Fill(h, vtype + "/trecsimpullwide", (v.t() - tsim) / v.tError(), simevt.is_signal());
	Fill(h, vtype + "/trecerrvsntrkprof", float(ntiming), v.tError(), simevt.is_signal());
	Fill(h, vtype + "/trecerrvsntrk", float(ntiming), v.tError(), simevt.is_signal());
	if( v.tError() < 0.1){
	  Fill(h, vtype + "/trecsim_sigmatlt01", v.t() - tsim, simevt.is_signal());
	  Fill(h, vtype + "/trecsimpull_sigmatlt01", (v.t() - tsim) / v.tError(), simevt.is_signal());
	}
      }
    }// 4D
}
/***************filvertexhistosMatched**********************************************************************/

void PrimaryVertexAnalyzer4PU::fillTrackHistos(std::map<std::string, TH1*>& h, const std::string& ttype, MTrack& tk, const reco::Vertex* v) {
  if(not fill_track_histos_) return;
  timer_start("fillTrackHistos");

  if(tk.has_track()) fillRecoTrackHistos(h, ttype, tk.track());

  if (f4D_ && tk.has_timing())
    {
      Fill(h, ttype + "/t0trk", tk.t());
      Fill(h, ttype + "/t0errtrk", tk.dt());
      Fill(h, ttype + "/t0qualtrk", tk.timeQuality());
      if (! (v==NULL) && (v->isValid()) && (v->ndof() > 10.)  && (v->tError()>0.)  && (v->tError()<9.) ){
	double deltat = tk.t() - v->t();
	Fill(h, ttype + "/trestrk", deltat);
	Fill(h, ttype + "/tpulltrk", deltat / sqrt(pow(tk.dt(),2) + pow(v->tError(), 2)));
      }
    }

  if (v != NULL && v->isValid() && v->ndof() > 10.) {
    double zvtx = v->position().z();
    double dz2 = tk.dzError() * tk.dzError() + (pow(wx_ * cos(tk.phi()), 2) + pow(wy_ * sin(tk.phi()), 2)) / pow(tan(tk.theta()), 2); // really?
    Fill(h, ttype + "/zrestrk", tk.z() - zvtx);
    Fill(h, ttype + "/zrestrkvsphi", tk.phi(), tk.z() - zvtx);
    Fill(h, ttype + "/zrestrkvseta", tk.eta(), tk.z() - zvtx);
    Fill(h, ttype + "/zpulltrkvsphi", tk.phi(), (tk.z() - zvtx) / sqrt(dz2));
    Fill(h, ttype + "/zrestrkvsz", zvtx, tk.z() - zvtx);
    Fill(h, ttype + "/zpulltrkvsz", zvtx, (tk.z() - zvtx) / sqrt(dz2));
    Fill(h, ttype + "/zpulltrkvseta", tk.eta(), (tk.z() - zvtx) / sqrt(dz2));
    Fill(h, ttype + "/zpulltrkvsz", tk.z(), pow(tk.z() - zvtx, 2) / (pow(v->zError(), 2) + dz2));
    Fill(h, ttype + "/zpulltrk", (tk.z() - zvtx) / sqrt(pow(v->zError(), 2) + dz2));
  }

  fillTrackHistosMatched(h, ttype, tk);
  timer_stop("fillTrackHistos");
}

/* FIXME   scheduled for removal 
// basically obsolete, all of this can be done with MTracks
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
*/



void PrimaryVertexAnalyzer4PU::fillRecoTrackHistos(std::map<std::string, TH1*>& h,
                                                   const std::string& ttype,
                                                   const reco::Track& t) {
  /* reco::track version, no transient tracks here */
  Fill(h, ttype + "/eta", t.eta());
  Fill(h, ttype + "/ztrk", t.vz());
  Fill(h, ttype + "/zerrtrk", t.dzError());
  Fill(h, ttype + "/logzerrtrk", log(t.dzError()/log(10.)));
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
  Fill(h, ttype + "/ndoftrk", t.ndof());
  Fill(h, ttype + "/chi2trk", t.chi2());
  Fill(h, ttype + "/nchi2trk", t.normalizedChi2());
  Fill(h, ttype + "/nchi2trkvsz", t.vz(), t.normalizedChi2());
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

void PrimaryVertexAnalyzer4PU::fillTrackHistosMatched(std::map<std::string, TH1*>& h, const std::string& ttype, MTrack& tk) {
  // ref to MultiTrackValidator (MTV):
  // https://github.com/cms-sw/cmssw/blob/f0cac195d32cda5f4f48bfd99a751449e7624861/Validation/RecoTrack/src/MTVHistoProducerAlgoForTracker.cc#L2398
  if(not tk.is_matched()) return;

  // TrackingParticle matched to this recoTrack originates from the signal SimVertex
  auto const isSignalPV = (tk._simEvt->index == 0);

  auto const rec_pt = tk.pt();
  auto const rec_dpt = tk.ptError();
  auto const rec_dptrel = rec_dpt/rec_pt;
  auto const rec_eta = tk.eta();
  auto const rec_phi = tk.phi();
  auto const rec_z = tk.z();
  auto const rec_dz = tk.dzError();

  auto const sim_pt = tk._tpr->pt();
  auto const sim_eta = tk._tpr->eta();
  auto const sim_phi = tk._tpr->phi();
  auto const sim_z = tk.zsim();

  auto const dptDiff = (rec_pt - sim_pt);
  auto const dptDiffRel = dptDiff/rec_pt;
  auto const dptPull = dptDiff/rec_dpt;
  auto const detaDiff = (rec_eta - sim_eta);
  auto const dphiDiff = reco::deltaPhi(rec_phi, sim_phi);

  auto const dzDiff = (rec_z - sim_z);
  auto const dzPull = dzDiff/rec_dz;

  Fill(h, ttype + "/tkzrecsim", dzDiff, isSignalPV);
  Fill(h, ttype + "/tkzrecsimpull", dzPull, isSignalPV);
  Fill(h, ttype + "/tkzrecsimpullvseta", rec_eta, dzPull, isSignalPV);
  Fill(h, ttype + "/tkzrecsimpullsqvseta", rec_eta, dzPull*dzPull, isSignalPV);
  Fill(h, ttype + "/tkzrecsimvseta", rec_eta, dzDiff, isSignalPV);
  Fill(h, ttype + "/tkzrecsimvseta2d", rec_eta, dzDiff, isSignalPV);
  if(tk.pt() > trkhiptmin_) Fill(h, ttype + "/tkzrecsimvseta2dhipt", rec_eta, dzDiff, isSignalPV);
  Fill(h, ttype + "/tkzrecsimsqvseta", rec_eta, dzDiff*dzDiff, isSignalPV);
  if(tk.pt() > trkhiptmin_) Fill(h, ttype + "/tkzrecsimsqvsetahipt", rec_eta, dzDiff*dzDiff, isSignalPV);
  Fill(h, ttype + "/tkzrecsimvsz", sim_z, dzDiff, isSignalPV);
  Fill(h, ttype + "/tkzrecsimvslogpt", log(rec_pt)/log(10.), dzDiff * (rec_eta < 0. ? -1. : 1.), isSignalPV);
  Fill(h, ttype + "/tkzrecsimvsetaz", rec_eta, dzDiff * (rec_z < 0. ? -1. : 1.), isSignalPV);

  Fill(h, ttype + "/tkptrecsimrelvssimeta2d", sim_eta, dptDiffRel, isSignalPV);
  Fill(h, ttype + "/tkptrecsimpullvssimeta2d", sim_eta, dptPull, isSignalPV);
  Fill(h, ttype + "/tkdptrelrecvssimeta2d", sim_eta, rec_dptrel, isSignalPV);
  Fill(h, ttype + "/tketarecsimvssimeta2d", sim_eta, detaDiff, isSignalPV);
  Fill(h, ttype + "/tkphirecsimvssimeta2d", sim_eta, dphiDiff, isSignalPV);
  Fill(h, ttype + "/tkzrecsimvssimeta2d", sim_eta, dzDiff, isSignalPV);
  Fill(h, ttype + "/tkzrecsimpullvssimeta2d", sim_eta, dzPull, isSignalPV);
  Fill(h, ttype + "/tkdzrecvssimeta2d", sim_eta, rec_dz, isSignalPV);

  if (f4D_ and tk.has_timing()){

    auto const rec_t = tk.t();
    auto const rec_dt = tk.dt();
    auto const sim_t = tk.tsim();

    auto const dtDiff = (rec_t - sim_t);
    auto const dtPull = dtDiff/rec_dt;
    auto const dtPull2 = dtPull*dtPull;

    Fill(h, ttype + "/tktrecsim", dtDiff, isSignalPV);
    Fill(h, ttype + "/tktrecsimpull", dtPull, isSignalPV);
    Fill(h, ttype + "/tktrecsimpullwide", dtPull, isSignalPV);
    Fill(h, ttype + "/tktrecsimvseta2d", rec_eta, dtDiff, isSignalPV);
    Fill(h, ttype + "/tktrecsimpullsqvserr", rec_dt, dtPull2, isSignalPV);
    Fill(h, ttype + "/tktrecsimpullvserr", rec_dt, dtPull, isSignalPV);


    
    Fill(h, ttype + "/tkpidvsetalogpt", std::abs(rec_eta), log(rec_pt)/log(10.));
    if (tk.is_pion() || tk.is_muon()){
      Fill(h, ttype + "/tkpidpionvsetalogpt", std::abs(rec_eta), log(rec_pt)/log(10.));
    }else if(tk.is_kaon()){
	Fill(h, ttype + "/tkpidkaonvsetalogpt", std::abs(rec_eta), log(rec_pt)/log(10.));
    }else if(tk.is_proton()){
	Fill(h, ttype + "/tkpidprotonvsetalogpt", std::abs(rec_eta), log(rec_pt)/log(10.));
    }

  }
}

void PrimaryVertexAnalyzer4PU::fillTrackClusterHistos(std::map<std::string, TH1*>& h, const std::string& ttype, const reco::Track& t, const reco::Vertex* v) {
  if (!RECO_) return;

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

void PrimaryVertexAnalyzer4PU::printRecVtxs(const MVertexCollection & vcollection, std::string title) {
  int ivtx = 0;

  std::cout << std::endl << title << "  nv=" << vcollection.vtxs.size() << "" << std::endl;

  for (auto v : vcollection.vtxs){
    string vtype = " recvtx";
    if (v.isRecoFake()) {
      vtype = " \033[31mfake\033[0m  ";
    } else if (v.ndof() == -3) {
      vtype = " event   ";
    }else if (!select(v)){
      vtype = " \033[31mreject\033[0m";
    }
    std::cout << "vtx " << std::setw(3) << std::setfill(' ') << ivtx++ << vtype
	      << " #trk " << std::fixed  << std::setprecision(4) << std::setw(3) << v.tracksSize()
	      << " chi2 " << std::fixed << std::setw(5) << std::setprecision(1) << v.chi2()
	      << " sumw " << std::fixed << std::setw(6) << std::setprecision(2) << v.sumw()
	      << " ndof " << std::fixed << std::setw(6) << std::setprecision(2) << v.ndof()
	      << "  x=" << std::setw(7) << std::fixed << std::setprecision(4)
	      << v.x() << "+/-"  << std::setw(6) << v.xError()                                                // <<  std::endl
              << "  y=" << std::setw(7) << v.y() << "+/-" << std::setw(6) << v.yError()  //<< std::endl
              << "  z=" << std::setw(8) << v.z() << "+/-" << std::setw(6) << v.zError();
    if(f4D_){
      std::cout << " t=" << std::setw(7) << v.t() << "+/-" << std::setw(5) << v.tError(); 
    }
    std::cout << "  dxy= ("
              << std::setw(7) << std::fixed << std::setprecision(4) << v.x() - vertexBeamSpot_.x(v.z()) << ","
              << std::setw(7) << std::fixed << std::setprecision(4) << v.y() - vertexBeamSpot_.y(v.z()) << ")"
              << " pxy= " << std::setw(6) << std::fixed << std::setprecision(4) << v.pxy(vertexBeamSpot_)
              << " ptsum= " << std::setw(6)  << std::fixed << std::setprecision(1) << v.sumabspt()
      //<< " sumpt2= " << std::setw(6)  << std::fixed << std::setprecision(1) << v.sumpt2()
              << " ptmax2= " << std::setw(5) << std::fixed << std::setprecision(2) << v.ptmax2() << std::endl;
  }
  std::cout << std::endl;
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
    std::cout << "vtx " << std::setw(3) << std::setfill(' ') << ivtx++ << vtype
	      << " #trk " << std::fixed  << std::setprecision(4) << std::setw(3) << v->tracksSize()
	      << " chi2 " << std::fixed << std::setw(5) << std::setprecision(1) << v->chi2()
	      << " sumw " << std::fixed << std::setw(6) << std::setprecision(2) << vertex_sumw(*v)
	      << " ndof " << std::fixed << std::setw(6) << std::setprecision(2) << v->ndof()
	      << "  x=" << std::setw(7) << std::fixed << std::setprecision(4)
	      << v->x() << " +/-"  << std::setw(6) << v->xError()                                                // <<  std::endl
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
  bool lt(const std::pair<double, unsigned int>& a, const std::pair<double, unsigned int>& b) {
    return a.first < b.first;
  }
}  // namespace
/********************************************************************************************************/




/********************************************************************************************************/
std::vector<bool> PrimaryVertexAnalyzer4PU::trackClass(const reco::Track& t)
/********************************************************************************************************/
{
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
void PrimaryVertexAnalyzer4PU::testTrackParameters(Tracks& tracks)
/********************************************************************************************************/
{
  // just a bunch of sanity checks
  cout << "BEAM " 
       << " " << fixed << setprecision(5) << setw(10) << vertexBeamSpot_.x0()
       << " " << fixed << setprecision(5) << setw(10) << vertexBeamSpot_.y0()
       << " " << fixed << setprecision(5) << setw(10) << vertexBeamSpot_.z0()
       << endl;
  cout << "     " 
       << " " << fixed << setprecision(5) << setw(10) << vertexBeamSpot_.x(vertexBeamSpot_.z0())
       << " " << fixed << setprecision(5) << setw(10) << vertexBeamSpot_.y(vertexBeamSpot_.z0())
       << " " << fixed << setprecision(5) << setw(10) << vertexBeamSpot_.z0()
       << endl;

  for(auto tk : tracks){
    const reco::Track & t = tk.track();
    if(t.pt() > 0.2) continue;
    double zB = vertexBeamSpot_.z0();
    double xB = vertexBeamSpot_.x(zB);
    double yB = vertexBeamSpot_.y(zB);
    //double x0 = t.vx();
    //double y0 = t.vy();
    //double z0 = t.vz();
    // curvilinear from track
    double qoverp = t.parameter(0);
    double lambda = t.parameter(1); // lambda = pi/2 - polar angle at the given point 
    double phi = t.parameter(2);    // phi = azimuth angle at the given point (i.e. at x0,y0,z0)
    double dxy = t.parameter(3);    // dxy = -vx*sin(phi) + vy*cos(phi) 
    double dsz = t.parameter(4);    // dsz = vz*cos(lambda) - (vx*cos(phi)+vy*sin(phi))*sin(lambda)
    // eucledian from curvilinear
    double p = 1/fabs(qoverp);
    double px = p * cos(lambda) * cos(phi);
    double py = p * cos(lambda) * sin(phi);
    double pz = p * sin(lambda);
    double vx = -dxy*sin(phi) + (xB * cos(phi) + yB * sin(phi)) * cos(phi);
    double vy =  dxy*cos(phi) + (xB * cos(phi) + yB * sin(phi)) * sin(phi);
    double vz = dsz / cos(lambda) + (xB * cos(phi) + yB * sin(phi)) * tan(lambda);
    // equivalent 
    double dxyB = dxy + xB * sin(phi) - yB *cos(phi);  // exact
    double vx2 = xB - sin(phi)* dxyB;
    double vy2 = yB + cos(phi)* dxyB;
    //double dszB = dsz - zB * cos(lambda) + (xB * cos(phi) + yB * sin(phi)) * sin(lambda);  // exact?
    cout << endl;
    cout << "eucl(track) " 
	 << " p=      "  << setprecision(5) << setw(10) << t.p()
	 << " px=     "  << setprecision(5) << setw(10) << t.px()
	 << " py=     "  << setprecision(5) << setw(10) << t.py()
	 << " pz=     "  << setprecision(5) << setw(10) << t.pz()
	 << " vx=     "  << setprecision(5) << setw(10) << t.vx()
	 << " vy=     "  << setprecision(5) << setw(10) << t.vy()
	 << " vz=     "  << setprecision(5) << setw(10) << t.vz() 
	 << endl;
    cout << "eucl(test)  " 
	 << " p=      "  << setprecision(5) << setw(10) << p
	 << " px=     "  << setprecision(5) << setw(10) << px
	 << " py=     "  << setprecision(5) << setw(10) << py
	 << " pz=     "  << setprecision(5) << setw(10) << pz
	 << " vx=     "  << setprecision(5) << setw(10) << vx
	 << " vy=     "  << setprecision(5) << setw(10) << vy
	 << " vz=     "  << setprecision(5) << setw(10) << vz
	 << endl;
    cout << "eucl(test2) " 
	 << " p=      "  << setprecision(5) << setw(10) << 0
	 << " px=     "  << setprecision(5) << setw(10) << 0
	 << " py=     "  << setprecision(5) << setw(10) << 0
	 << " pz=     "  << setprecision(5) << setw(10) << 0
	 << " vx=     "  << setprecision(5) << setw(10) << vx2
	 << " vy=     "  << setprecision(5) << setw(10) << vy2
      //	 << " vz=     "  << setprecision(5) << setw(10) << vz2
         << "                                            db = " << dxyB << "   |ref-beam|=" << sqrt(pow(t.vx()-xB,2)+pow(t.vy()-yB,2))
	 << endl;
    cout << "curv(track) " 
	 << " q/p=    " << setprecision(5) << setw(10) << t.parameter(0)
	 << " lambda= " << setprecision(5) << setw(10) << t.parameter(1)
	 << " phi=    " << setprecision(5) << setw(10) << t.parameter(2)
	 << " dxy=    " << setprecision(5) << setw(10) << t.parameter(3)
	 << " dsz=    " << setprecision(5) << setw(10) << t.parameter(4)
	 << endl;

    cout << "curv(test)  " 
	 << " q/p=    " << setprecision(5) << setw(10) << 1./sqrt(pow(t.px(),2) + pow(t.py(),2) + pow(t.pz(),2))
	 << " lambda= " << setprecision(5) << setw(10) << atan2(t.pz(), t.pt())
	 << " phi=    " << setprecision(5) << setw(10) << atan2(t.py(), t.px())
	 << " dxy=    " << setprecision(5) << setw(10) << ( -t.vx() * t.py()/ t.pt() + t.vy() * t.px()/ t.pt())
	 << " dsz=    " << setprecision(5) << setw(10) << t.vz() * t.pt() / t.p() - (t.vx() * t.px() + t.vy() * t.py()) / t.pt() * t.pz() / t.p()
	 << " dsz()=  " << setprecision(5) << setw(10) << t.dsz() // vz() * theptoverp - (vx() * px() + vy() * py()) / thept * pz() * thepinv;
	 << endl;

    /*
    // linear approximation, the point of closest approach to (x1,y1) 
    double x1 = 0, y1=0;
    double spca1 = (x1-x0) * cos(phi) + (y1-y0)*sin(phi);
    double xpca_origin = x0 + spca1 * cos(phi);
    double ypca_origin = y0 + spca1 * sin(phi);
    x1 = xB;
    y1 = yB;
    spca1 = (x1-x0) * cos(phi) + (y1-y0)*sin(phi);
    double xpca_beam0 = x0 + spca1 * cos(phi);
    double ypca_beam0 = y0 + spca1 * sin(phi);
    x1 = vertexBeamSpot_.x(z0);
    y1 = vertexBeamSpot_.y(z0);
    spca1 = (x1-x0) * cos(phi) + (y1-y0)*sin(phi);
    double xpca_beamz = x0 + spca1 * cos(phi);
    double ypca_beamz = y0 + spca1 * sin(phi);
    //
    //double double db1 = xb ∗ sin(par[2]) − yb ∗ cos( par [2] );
    //
    // question 1) is vx,vy,vz the point of closest approach to the beam-spot or 'the center of CMS' or the origin or what?
    cout << "PPQ  " 
	 << " ?? "  << setprecision(5) << setw(10) << x0 + dxy *sin(phi) 
	 << " "  << setprecision(5) << setw(10) << y0 - dxy *cos(phi)  
	 << " O? "  << setprecision(5) << setw(10) <<  + dxy *sin(phi) 
	 << "   "  << setprecision(5) << setw(10) <<  - dxy *cos(phi)  
	 << " test " << fixed << setprecision(5) << setw(10) << dxy -( -x0 *sin(phi) + y0*cos(phi))
	 << " dxyE " << fixed << setprecision(5) << setw(10) << t.dxyError()
	 << " z0 " << fixed << setprecision(5) << setw(10) << z0
	 << " |O  " << fixed << setprecision(5) << setw(10) << xpca_origin << ","  << fixed << setprecision(5) << setw(10)<< ypca_origin
	 << " |B0 " << fixed << setprecision(5) << setw(10) << xpca_beam0  << ","  << fixed << setprecision(5) << setw(10)<< ypca_beam0
	 << " |BZ " << fixed << setprecision(5) << setw(10) << xpca_beamz  << ","  << fixed << setprecision(5) << setw(10)<< ypca_beamz
	 << " |   " << fixed << setprecision(5) << setw(10) << x0  << ","  << fixed << setprecision(5) << setw(10)<< y0
	 << endl;

    //double x00 = x0  - dxy *sin(phi);
    //double y00 = y0  + dxy *cos(phi);
    //double z00 = z0 - dsz / 
    */

  // <<<<<<<<<<   temporary test
  }
}


/********************************************************************************************************/
void PrimaryVertexAnalyzer4PU::analyzeTracksTP(Tracks& tracks, std::vector<SimEvent>& simEvt)
/********************************************************************************************************/
{
  //const double vertexSize = 0.1;
  const double d0CutOff = 3.;

  // now histogram all tracks
  for (unsigned int i = 0; i < tracks.size(); i++) {
    const MTrack & tk = tracks(i);

    double t_pi = 1. / (1. + exp(std::pow(tk.ip() / tk.dip(), 2) - std::pow(d0CutOff, 2)));

    if (tk.has_timing()) {
      Fill(hTrk, "ttrk_rec_all_wide", tk.t());
      Fill(hTrk, "ttrk_rec_all", tk.t());
      Fill(hTrk, "terrtrk_rec_all", tk.dt());
      Fill(hTrk, "tqualtrk_rec_all", tk.timeQuality());
    }

    if (tk.selected() && tk.has_timing()) {
      Fill(hTrk, "ttrk_rec_sel_wide", tk.t());
      Fill(hTrk, "ttrk_rec_sel", tk.t());
      Fill(hTrk, "terrtrk_rec_sel", tk.dt());
      Fill(hTrk, "tqualtrk_rec_sel", tk.timeQuality());
      if(tk.dt() < 0.1)  Fill(hTrk, "tqualtrk_sigmat01_rec_sel", tk.timeQuality());

      if (tk.is_matched()) {
	double t_pid = tk.get_t_pid();
	if (std::abs(t_pid) < 10.){
	  double tres = t_pid - tk.tsim();
	  double tpull = tres / tk.MTD_timeerror();
	  Fill(hTrk, "trestrk_selmatched_pid", tres );
	  Fill(hTrk, "tpulltrk_selmatched_pid", tpull );
	  if(tk.MTD_timeerror() < 0.1){
	    Fill(hTrk, "trestrk_sigmatlt01_selmatched_pid", tres );
	    Fill(hTrk, "tpulltrk_sigmatlt01_selmatched_pid", tpull );
	  }
	  if(tk.pt() > 1.0){
	    Fill(hTrk, "trestrk_selmatched_pid_ptgt1", tres );
	    Fill(hTrk, "tpulltrk_selmatched_pid_ptgt1", tpull );
	  }
	  if (tk.timeQuality() < 0.5){
	    Fill(hTrk, "trestrk_selmatched_pid_qlt05", tres );
	    Fill(hTrk, "tpulltrk_selmatched_pid_qlt05", tpull );
	  }else if(tk.timeQuality() < 0.8){
	    Fill(hTrk, "trestrk_selmatched_pid_q0508", tres );
	    Fill(hTrk, "tpulltrk_selmatched_pid_q0508", tpull );
	  }else{
	    Fill(hTrk, "trestrk_selmatched_pid_qgt08", tres );
	    Fill(hTrk, "tpulltrk_selmatched_pid_qgt08", tpull );
	  }
	}

        Fill(hTrk, "ttrk_rec_selmatched_wide", tk.t());
        Fill(hTrk, "ttrk_rec_selmatched", tk.t());
        Fill(hTrk, "terrtrk_rec_selmatched", tk.dt());
        Fill(hTrk, "ttrk_sim_selmatched", tk.tsim());
        Fill(hTrk, "trestrk_selmatched", tk.t() - tk.tsim());
        Fill(hTrk, "tpulltrk_selmatched", (tk.t() - tk.tsim()) / tk.dt());
	Fill(hTrk, "tqualtrk_rec_selmatched", tk.timeQuality());

	if(tk.dt() < 0.1)
	  {
	    Fill(hTrk, "tpulltrk_sigmat01_selmatched", (tk.t() - tk.tsim()) / tk.dt());
	    Fill(hTrk, "tqualtrk_sigmat01_selmatched", tk.timeQuality());
	    if(fabs(tk.eta()) <  1.4)
	      {// barrel timing layer
		Fill(hTrk, "trestrk_selmatched_barrel", tk.t() - tk.tsim());
		Fill(hTrk, "terrtrk_rec_selmatched_barrel", tk.dt() );
		Fill(hTrk, "tpulltrk_selmatched_barrel", (tk.t() - tk.tsim()) /  tk.dt() );
		Fill(hTrk, "trestrkvszrestrk_selmatched_barrel", tk.z() - tk.zsim(), tk.t() - tk.tsim());
		Fill(hTrk, "tqualtrk_rec_selmatched_barrel", tk.timeQuality());
		if(tk.pt() > trkhiptmin_){
		  Fill(hTrk, "tpulltrk_selmatched_barrel_hipt", (tk.t() - tk.tsim()) /  tk.dt() );
		  Fill(hTrk, "trestrk_selmatched_barrel_hipt", tk.t() - tk.tsim());
		  Fill(hTrk, "tqualtrk_rec_selmatched_barrel_hipt", tk.timeQuality());
		}
		if(tk.pt() < trkloptmax_){
		  Fill(hTrk, "trestrk_selmatched_barrel_lopt", tk.t() - tk.tsim());
		  Fill(hTrk, "tpulltrk_selmatched_barrel_lopt", (tk.t() - tk.tsim()) /  tk.dt() );
		  Fill(hTrk, "tqualtrk_rec_selmatched_barrel_lopt", tk.timeQuality());
		}
		
	      }
	    else if (fabs(tk.eta()) > 1.7)
	      {// endcap
		Fill(hTrk, "trestrk_selmatched_endcap", tk.t() - tk.tsim());
		Fill(hTrk, "terrtrk_rec_selmatched_endcap", tk.dt() );
		Fill(hTrk, "tpulltrk_selmatched_endcap", (tk.t() - tk.tsim()) /  tk.dt() );
		Fill(hTrk, "tqualtrk_rec_selmatched_endcap", tk.timeQuality());
		if(tk.pt() > trkhiptmin_){
		  Fill(hTrk, "trestrk_selmatched_endcap_hipt", tk.t() - tk.tsim());
		  Fill(hTrk, "tpulltrk_selmatched_endcap_hipt", (tk.t() - tk.tsim()) /  tk.dt() );
		  Fill(hTrk, "tqualtrk_rec_selmatched_endcap_hipt", tk.timeQuality());
		}
		if(tk.pt() < trkloptmax_){
		  Fill(hTrk, "trestrk_selmatched_endcap_lopt", tk.t() - tk.tsim());
		  Fill(hTrk, "tpulltrk_selmatched_endcap_lopt", (tk.t()- tk.tsim()) /  tk.dt() );
		  Fill(hTrk, "tqualtrk_rec_selmatched_endcap_lopt", tk.timeQuality());
		}
		Fill(hTrk, "trestrkvszrestrk_selmatched_endcap", tk.zres(), tk.tres());
		if(tk.eta() > 1.7)
		  {
		    Fill(hTrk, "trestrk_selmatched_fwd", tk.tres());
		  }
		else
		  {
		    Fill(hTrk, "trestrk_selmatched_bwd", tk.tres());
		  }
	      }
	  }

	// fill residuals end pulls for three regions of the quality variable
	if (tk.timeQuality() > 0.8){
	  Fill(hTrk, "trestrk_selmatched_qgt08", tk.t() - tk.tsim());
	  Fill(hTrk, "tpulltrk_selmatched_qgt08", (tk.t() - tk.tsim()) / tk.dt() );
	  if (tk.dt() < 0.1){
	    Fill(hTrk, "trestrk_sigmatlt01_selmatched_qgt08", tk.t() - tk.tsim());
	    Fill(hTrk, "tpulltrk_sigmatlt01_selmatched_qgt08", (tk.t ()- tk.tsim()) / tk.dt() );
	  }
	}else if(tk.timeQuality() < 0.5){
	  Fill(hTrk, "trestrk_selmatched_qlt05", tk.t() - tk.tsim());
	  Fill(hTrk, "tpulltrk_selmatched_qlt05", (tk.t() - tk.tsim()) / tk.dt() );
	  if (tk.dt() < 0.1){
	    Fill(hTrk, "trestrk_sigmatlt01_selmatched_qlt05", tk.t() - tk.tsim());
	    Fill(hTrk, "tpulltrk_sigmatlt01_selmatched_qlt05", (tk.t() - tk.tsim()) / tk.dt() );
	  }
	}else if((tk.timeQuality() >= 0.5) && (tk.timeQuality() <= 0.8)){
	  Fill(hTrk, "trestrk_sigmatlt01_selmatched_q0508", tk.t() - tk.tsim());
	  Fill(hTrk, "tpulltrk_sigmatlt01_selmatched_q0508", (tk.t() - tk.tsim()) / tk.dt() );
	  if (tk.dt() < 0.1){
	    Fill(hTrk, "trestrk_selmatched_q0508", tk.tres());
	    Fill(hTrk, "tpulltrk_selmatched_q0508", tk.tpull());
	  }
	}
	  
        double mtd_dt = tk.MTD_timeerror();
	if (tk.is_pion()){
	  Fill(hTrk, "trestrk_selpion", tk.tres());
	  Fill(hTrk, "tpulltrk_selpion", tk.tpull());
	  Fill(hTrk, "trestrkh_selpion", tk.th[0] - tk.tsim());
	  Fill(hTrk, "tpulltrkh_selpion", (tk.th[0] - tk.tsim()) / mtd_dt);
	  Fill(hTrk, "trestrkh_selhyp", tk.th[0] - tk.tsim());
	  Fill(hTrk, "tpulltrkh_selhyp", (tk.th[0] - tk.tsim()) / mtd_dt);
	  if(tk.dt() < 0.1) Fill(hTrk, "tpulltrk_sigmat01_selpion", tk.tpull());
	}else if (tk.is_kaon()){
	  Fill(hTrk, "trestrk_selkaon", tk.tres());
	  Fill(hTrk, "tpulltrk_selkaon", tk.tpull());
	  Fill(hTrk, "trestrkh_selkaon", tk.th[1] - tk.tsim());
	  Fill(hTrk, "tpulltrkh_selkaon", (tk.th[1] - tk.tsim()) / mtd_dt);
	  Fill(hTrk, "trestrkh_selhyp", tk.th[1] - tk.tsim());
	  Fill(hTrk, "tpulltrkh_selhyp", (tk.th[1] - tk.tsim()) / mtd_dt);
	  if(tk.dt() < 0.1) Fill(hTrk, "tpulltrk_sigmat01_selkaon", tk.tpull());
	}else if (tk.is_proton()){
	  Fill(hTrk, "trestrk_selproton", tk.tres());
	  Fill(hTrk, "tpulltrk_selproton", tk.tpull());
	  Fill(hTrk, "trestrkh_selproton", tk.th[2] - tk.tsim());
	  Fill(hTrk, "tpulltrkh_selproton", (tk.th[2] - tk.tsim()) / mtd_dt);
	  Fill(hTrk, "trestrkh_selhyp", tk.th[2] - tk.tsim());
	  Fill(hTrk, "tpulltrkh_selhyp", (tk.th[2] - tk.tsim()) / mtd_dt);
	  if(tk.dt() < 0.1) Fill(hTrk, "tpulltrk_sigmat01_selproton", tk.tpull());
	}	  
	
      } else {
        Fill(hTrk, "ttrk_rec_selunmatched_wide", tk.t());
        Fill(hTrk, "ttrk_rec_selunmatched", tk.t());
      }
    }// selected and has_timing


    Fill(hTrk, "etaprim", tk.eta());
    if(tk.selected()) Fill(hTrk, "etaprimsel", tk.eta());
    if(tk.pt() < 1.){
      Fill(hTrk, "etaprim_Pt000to001", tk.eta());
       if(tk.selected())Fill(hTrk, "etaprimsel_Pt000to001", tk.eta());
    }else if(tk.pt() < 3.){
      Fill(hTrk, "etaprim_Pt001to003", tk.eta());
       if(tk.selected())Fill(hTrk, "etaprimsel_Pt001to003", tk.eta());
    }else if(tk.pt() < 10.){
      Fill(hTrk, "etaprim_Pt003to010", tk.eta());
       if(tk.selected()) Fill(hTrk, "etaprimsel_Pt003to010", tk.eta());
    }

    // all primary tracks
    if(tk.is_primary()){
      Fill(hTrk, "d0pullprim", tk.ip() / tk.dip());
      Fill(hTrk, "zpullvsd0pullprim", std::fabs(tk.ip() / tk.dip()), std::fabs(tk.zpull()));
    }


    // plots for primary, selected tracks
    if (tk.selected() && tk.is_primary()) {
      double zpull = tk.zpull();
      Fill(hTrk, "zpulltrk_primselmatched", zpull);
      Fill(hTrk, "zpulltrkprimsel", tk.eta(), zpull);
      Fill(hTrk, "zpulltrkprimselvseta", tk.eta(), zpull);

      double logpt = log(tk.pt()) / log(10.);
      unsigned int bin = 4;
      if (tk.dzError() <= 0.01) {bin =0;}
      else if (tk.dzError() < 0.02) { bin = 1; }
      else if (tk.dzError() < 0.05) { bin = 2; }
      else if (tk.dzError() < 0.10) { bin = 3; }
      else { bin = 4;}
      const char * sbin = trkdzbin_[bin].c_str();
      

      Fill(hTrk, "zpulltrkprimselvslogpt", logpt, zpull);
      Fill(hTrk, Form("zpulltrkprimselvslogpt_%s", sbin), logpt, zpull);

      auto hitPattern = tk.track().hitPattern();
      auto nbarrel = hitPattern.pixelBarrelLayersWithMeasurement();
      if ( nbarrel < 2){
	Fill(hTrk, Form("zpulltrkprimselbpxlt2vseta_%s", sbin), tk.eta(), zpull);
      }else if (nbarrel > 2){
	Fill(hTrk, Form("zpulltrkprimselbpxgt2vseta_%s", sbin), tk.eta(), zpull);
      }

      Fill(hTrk, Form("ztailtprimselvslogpt_%s", sbin), logpt, (std::abs(zpull) > 3. ? 1. : 0.));

      if (std::abs(tk.eta()) > 1.2){
	Fill(hTrk, Form("ztailtprimselvslogpt_etahi_%s", sbin), logpt, (std::abs(zpull) > 3. ? 1. : 0.));
      }else{
	Fill(hTrk, Form("ztailtprimselvslogpt_etalo_%s", sbin), logpt, (std::abs(zpull) > 3. ? 1. : 0.));
      }

      Fill(hTrk, "tprimselvslogpteta", tk.eta(), logpt);
      Fill(hTrk, Form("tprimselvslogpteta_%s", sbin), tk.eta(), logpt);
      if(std::abs(zpull) > 3.){
	Fill(hTrk, "ztailtprimselvslogpteta", tk.eta(), logpt);
	Fill(hTrk, Form("ztailtprimselvslogpteta_%s",sbin), tk.eta(), logpt);
      }


      
      Fill(hTrk, "zrestrk_primselmatched", tk.zres());
      if (tk.has_timing()) {
        Fill(hTrk, "ttrk_rec_primselmatched", tk.t());
        Fill(hTrk, "ttrk_sim_primselmatched", tk.tsim());
        Fill(hTrk, "trestrk_primselmatched", tk.tres());
        Fill(hTrk, "tpulltrk_primselmatched", tk.tpull());
	Fill(hTrk, "zpulltrkt_primselmatched", tk.zpull());
	Fill(hTrk, "zrestrkt_primselmatched", tk.zres());
      }

      Fill(hTrk, "d0pullprimsel", tk.ip() / tk.dip());
      Fill(hTrk, "zpullvsd0pullprimsel", std::fabs(tk.ip() / tk.dip()), std::fabs(tk.zpull()));
      if (std::fabs(tk.eta()) < 1.5) {
	Fill(hTrk, "zpulltrkt_primselmatched_central", tk.zpull());
      } else if (std::fabs(tk.eta()) < 2.0) {
	Fill(hTrk, "zpulltrkt_primselmatched_inter", tk.zpull());
      } else {
        Fill(hTrk, "zpulltrkt_primselmatched_fwd", tk.zpull());
      }
      Fill(hTrk, "tpiprim", t_pi);

    }

    // secondary, not necessarily selected
    if(!tk.is_primary()){
      Fill(hTrk, "zpulltrksec", tk.zpull());
      Fill(hTrk, "d0pullsec", tk.ip() / tk.dip());
      Fill(hTrk, "zpullvsd0pullsec", std::fabs(tk.ip() / tk.dip()), std::fabs(tk.zpull()));
      Fill(hTrk, "tpisec", t_pi);
    }

  }  // track loop
}
/********************************************************************************************************/





/********************************************************************************************************/
void PrimaryVertexAnalyzer4PU::printPVTrksZT(Tracks & tracks,
                                             MVertexCollection& vcollection,
					     std::vector<SimEvent>& simEvt)
// make a printout of the tracks selected for PV reconstructions, show matching MC tracks, too
/********************************************************************************************************/
 {
  vector<MTrack> sortedTrks;

  bool show_timing = false;
  unsigned int nTimingTracks = 0;
  unsigned int nSelectedTracks = 0;

  for (auto tk : tracks){
    sortedTrks.push_back(tk);
    if (tk.has_timing()){
      nTimingTracks++;
    }
    if (tk.selected()){
      nSelectedTracks++;
    }
  }

  
  
  std::cout << "printPVTracksZT : " << nTimingTracks << " tracks out of " << tracks.size() << " have timing info"
            << std::endl;
  show_timing = (nTimingTracks > 0);
  cout << " selected " << nSelectedTracks << endl;


  stable_sort(sortedTrks.begin(), sortedTrks.end(), MTrack::lessz);

  // now dump in a format similar to the clusterizer
  cout << endl << "printPVTrksZT      run : event  =   " << run_ << " : " << event_ << endl;
  cout << "----                  z +/- dz";
  if (show_timing) {
    cout << "          t +/- dt    ";
  }
  cout << "    q1bfet-ilo:alg CSQ     ip +/-dip         pt   phi   eta";
  if ( MC_ ) {  // FIXME we need a flag for the availability of some kind of truth-matching
    cout << "  type      pdg    zvtx    zpull      match";
  }
  cout << endl;

  cout.precision(4);
  int isel = 0;
  double tz0 = -10000;

  for (unsigned int i = 0; i < sortedTrks.size(); i++) {

    auto tk = sortedTrks.at(i);
    unsigned int i0 = tk._index;  // original index in the track list before sorting (why needed?)

    unsigned int vmatch = tk.get_recv(vcollection.collection_label);
    double wmatch = 0.;
    if (vmatch != NO_RECVTX){
      wmatch = vcollection.at(vmatch).trackWeight(tk);
    }
    
    double tz = tk.z();
    double tantheta = tan(tk.theta());
    double phi = tk.phi();
    double tdz2 = pow(tk.dzError(), 2) + (pow(wx_ * cos(phi), 2) + pow(wy_ * sin(phi), 2)) / pow(tantheta, 2);
			  
    // print vertices in between tracks
    int iiv = 0;
    for (auto v : vcollection){
      if ((v.ndof() > 0) && (v.z() < tz) && (v.z() > tz0)) {
	if ((v.zError() > 0) && (v.zError() < 50.)){
	  cout << "rec [" << setw(4) << iiv << "]     " << setw(8) << fixed << setprecision(4) << v.z() << " +/-"
	       << setw(6) << v.zError();
	}else{
	  // Miniaod?
	  cout << "rec [" << setw(4) << iiv << "]     " << setw(8) << fixed << setprecision(4) << v.z() << " +/-      ";
	}	  
	
        if (show_timing) {
          if (v.has_timing()) {
	    if((v.tError() > 0) && (v.tError()<99.)){
	      cout << setw(7) << fixed << setprecision(3) << v.t() << " +/-" << setw(5) << v.tError();
	    }else{
	      cout << setw(7) << fixed << setprecision(3) << v.t() << " +/-     ";
	    }
          } else {
            cout << "        +/-     ";
          }
        }
        cout << "  ndof=" << v.ndof() << "  ptmax2=" << v.ptmax2();
        cout << "  sumpt2=" << v.sumpt2();
        if ((v.ndof() > 4) && (v.r(vertexBeamSpot_) > 0.0200)) {
          cout << "    rvtx= " << setw(8) << setprecision(4) << v.r(vertexBeamSpot_);
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
          if (simEvt[event].type == FROM_TRACKING_TRUTH) {
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

    }			  
    tz0 = tz;

    // for selected tracks print a selected track index
    if (tk.selected()) {
      cout << setw(4) << isel;
      isel++;
      //    } else if (!tt.stateAtBeamLine().isValid()) {
      //      cout << "XXXX";
    } else {
      cout << "    ";
    }

    if (vmatch != NO_RECVTX) {
      cout << "[" << setw(4) << vmatch << "]" << setw(4) << setprecision(2) << fixed << wmatch << " ";
    } else {
      cout << "           ";
    }

    if (fabs(tz) < 100) {
      cout << setw(8) << fixed << setprecision(4) << tz;
    } else {
      cout << setw(8) << fixed << setprecision(4) << 99.99 * tz / fabs(tz);
    }

    if (sqrt(tdz2) < 50) {
      cout << " +/-" << setw(6)  << fixed << setprecision(4) << sqrt(tdz2);
    } else {
      cout << " +/-      ";
    }

    if (show_timing) {
      if (tk.has_timing()){
	if(tk.dt() < 50.){
	  cout << setw(7) << fixed << setprecision(3) << tk.t() << " +/-" << setw(5) << tk.dt();
	}else{
	  cout << setw(7) << fixed << setprecision(3) << tk.t() << " +/-     ";
	}
      } else {
        cout << "                ";
      }
    }
  		  
    // track quality and hit information, see DataFormats/TrackReco/interface/HitPattern.h
    if (tk.is_highPurity()){
      cout << " *";
    } else {
      cout << "  ";
    }
    if (tk.has_validHitInFirstPixelBarrelLayer()) {
      cout << "+";
    } else {
      cout << "-";
    }
    if(tk.has_hitPattern()){
      // pixel layers with hits
      cout << setw(1) << tk.hitPattern().pixelBarrelLayersWithMeasurement();
      cout << setw(1) << tk.hitPattern().pixelEndcapLayersWithMeasurement();
      auto tepx = tk.hitPattern().pixelLayersWithMeasurement() -
	(tk.hitPattern().pixelEndcapLayersWithMeasurement() +
	 tk.hitPattern().pixelBarrelLayersWithMeasurement());  // just guessin'
      cout << setw(1) << tepx;
      // outer tracker layers with hits
      int mm = tk.hitPattern().trackerLayersWithMeasurement() - tk.hitPattern().pixelLayersWithMeasurement();
      if (mm >= 0) {
	cout << setw(1) << hex << mm << dec;
      } else {
	cout << "X";
      }
      
      // - missing hits  -[inner][outer],  2020-08-12 changed from numberOfAllHits to numberOfLostHits, seer PR #20938
      cout << "-" 
	   << setw(1) << hex << tk.hitPattern().numberOfLostHits(HitPattern::MISSING_INNER_HITS)
	   << setw(1) << hex << tk.lost()
	   << setw(1) << hex << tk.hitPattern().numberOfLostHits(HitPattern::MISSING_OUTER_HITS) << dec;
    }else{
      // no hit pattern available
      cout << "     " ;
    }

    if( tk.is_looper()){
      cout << "l";
    }else{
      cout << " ";
    }

    cout << ":" << setw(2) << dec << tk.algo();

    // track fit chi**2 / ndf
    if(tk.normalizedChi2()>0){
      cout << setw(4) << setprecision(1) << fixed << tk.normalizedChi2();
    }else{
      cout << "    ";
    }

    // xy impact paramter and error
    cout << setw(8) << setprecision(4) << tk.ip();
    if(tk.dip() < 1.){
      cout << "+/-" << setw(6) << tk.dip();
    }else{
      cout << "+/-      ";
    }

    // pt
    double relative_pterror = tk.ptError() / tk.pt();
    if (relative_pterror < 0.1) {
      cout << " " << setw(7) << setprecision(2) << tk.pt() * tk.charge();
    } else if (relative_pterror < 0.5) {
      cout << " " << setw(6) << setprecision(1) << tk.pt() * tk.charge() << "-";
    } else {
      cout << " " << setw(6) << setprecision(1) << tk.pt() * tk.charge() << "*";
    }
    // phi and eta
    cout << " " << setw(5) << setprecision(2) << tk.phi() << " " << setw(5) << setprecision(2) << tk.eta();

    // print MC info, if available
    if (MC_) {
      if ((simEvt.size() > 0) && (tracking_truth_available_)){
	//try to get by without the collection, may fail
	auto tk = tracks(i0);
	bool matched = tk.is_matched();
	TrackingParticleRef tpr = tk._tpr;
	auto trtb = tracks.ref(i0); // this will only work when we have a trackCollection (i.e. not for miniaod)
	
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
            for (unsigned int j = 0; j < simEvt.size(); j++) {
              if (simEvt[j].eventId == tpr->eventId()) {
                ez1 = simEvt[j].z;
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
            //for (View<Track>::size_type j = 0; j < tC.size(); ++j) {
            for (unsigned int j = 0; j < trkidx2tp_.size(); ++j) { // FIXME : never filled ???
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
        if (tk.is_matched()){
          if (tk.is_primary()){
            cout << " prim ";
          } else {
            cout << " sec  ";
          }
	  /* FIXME TODO, fill the equivalent of tsim[rectosim[i]] into the MTrack
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
	  */
        }
      }
    }
    cout << endl;
  }
  cout << "printPVTrksZT.end ---------------- " << endl;
}
/********************************************************************************************************/









/********************************************************************************************************/
void PrimaryVertexAnalyzer4PU::getTc(
    const std::vector<MTrack>& tracks, double& Tc, double& chsq, double& dzmax, double& dztrim) 
/********************************************************************************************************/
{
  if (tracks.size() < 2) {
    Tc = -1;
    chsq = -1;
    dzmax = -1;
    dztrim = -1;
    return;
  }

  double sumw = 0, sumwz = 0, sumww = 0, sumwwz = 0, sumwwzz = 0;
  double zmin = 1e10, zmin1 = 1e10, zmax1 = -1e10, zmax = -1e10;
  double m4 = 0, m3 = 0, m2 = 0, m1 = 0, m0 = 0;

  for(auto tk : tracks){
    double tantheta = tan(tk.theta());
    double z = tk.z();
    double dz2 = pow(tk.dzError(), 2) + pow(vertexBeamSpot_.BeamWidthX() / tantheta, 2);
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
  } else {
    chsq = 0;
    Tc = 0;
  }
  dzmax = zmax - zmin;
  dztrim = zmax1 - zmin1;  // truncated
}//getTc
/********************************************************************************************************/




/********************************************************************************************************/
bool PrimaryVertexAnalyzer4PU::select_truthMatchedTrack(const edm::RefToBase<reco::Track> track, TrackingParticleRef& tpr)const

/********************************************************************************************************/
// for a reco track select the matching tracking particle, always use this function to make sure we
// are consistent
// after calling select_truthMatchedTrack, tpr may have changed its value
// to get the TrackingParticle from the TrackingParticleRef, use ->get();
{
  if (tp_r2s_->find(track) == tp_r2s_->end()) {
    return false;
  } else {
    double f = -1e10;
    TrackingParticleRef tf;
    std::vector<std::pair<TrackingParticleRef, double>> tp = (*tp_r2s_)[track];
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
}//select_truthMatchedTrack
/********************************************************************************************************/




/********************************************************************************************************/
void PrimaryVertexAnalyzer4PU::printTruthMatchValues(edm::RefToBase<reco::Track> track)

/********************************************************************************************************/
// print the tracking truth assocation values for a reco track
{
  if (tp_r2s_->find(track) == tp_r2s_->end()) {
    return;
  }

  cout << "[" << setw(4) << setprecision(2) << fixed;
  std::vector<std::pair<TrackingParticleRef, double>> tp = (*tp_r2s_)[track];
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
}//printTruthMatchValues
/********************************************************************************************************/




/********************************************************************************************************/
std::vector<edm::RefToBase<reco::Track>> PrimaryVertexAnalyzer4PU::getTruthMatchedVertexTracks_obsolete(const reco::Vertex& v,
                                                                                               double min_weight) const
// for rec vertex v get a list of tracks for which truth matching is available
/********************************************************************************************************/
{
  std::vector<edm::RefToBase<reco::Track>> b;
  TrackingParticleRef tpr;

  for (trackit_t tv = v.tracks_begin(); tv != v.tracks_end(); tv++) {
    if (v.trackWeight(*tv) >= min_weight) {
      if (select_truthMatchedTrack(*tv, tpr)) {
        b.push_back(*tv);
      }
    }
  }
  return b;
}
/********************************************************************************************************/



/********************************************************************************************************/
std::vector<PrimaryVertexAnalyzer4PU::SimEvent> PrimaryVertexAnalyzer4PU::getSimEvents_miniaod(
											       const edm::Event & iEvent, 
											       Tracks & tracks
											       )
/********************************************************************************************************/
{
  //typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag> GlobalPoint; in DataFormats/Math/interface/Point3D.h
  Handle<GenEventVertex> genParticlesXyz0Handle;
  iEvent.getByToken(theGenParticlesXyz0Token_, genParticlesXyz0Handle);
  GenEventVertex  genXyz0 = (*genParticlesXyz0Handle.product());
  Handle<float> genParticlesT0Handle;
  iEvent.getByToken(theGenParticlesT0Token_, genParticlesT0Handle);
  float t0 = (*genParticlesT0Handle.product());
  if(verbose_){
    cout << "getSimEvents_miniaod  genXyz0.z()="<< genXyz0 <<   "  t0=" << t0 << endl;
  }

  Handle<edm::View<pat::PackedGenParticle> > genParticlesHandle;
  iEvent.getByToken(theGenParticlesToken_, genParticlesHandle);
  edm::View<pat::PackedGenParticle> gen = (*genParticlesHandle.product());
  if(verbose_){
    cout << "getSimEvents_miniaod  gen.size()="<< gen.size() << endl;
  }


  // the "important" particles
  Handle<edm::View<reco::GenParticle> > prunedGenParticlesHandle;
  iEvent.getByToken(thePrunedGenParticlesToken_, prunedGenParticlesHandle);
  edm::View<reco::GenParticle> pruned = (*prunedGenParticlesHandle.product());

  std::map<const reco::GenParticle *, unsigned int> pm;
  if(verbose_){
    cout << "getSimEvents_miniaod  pruned.size()="<< pruned.size() << endl;
    /*
    if(pruned.size()>0){
      for(unsigned int p = 0; p < pruned.size(); p++){
	const reco::GenParticle * pp = &pruned[p];
	pm[pp] = p;
	cout << "pruned " << p << " pdgid " << pruned[p].pdgId() << "  vertex = " << pruned[p].vertex() << "  " << &(pruned[p]) << endl;
      }
    }
    */
  }


  /*
  edm::Handle<GenEventInfoProduct> genEvtInfoH;
  iEvent.getByLabel( "generator", genEvtInfoH );
  cout << "qScale = " << genEvtInfoH->qScale();  
  */

  vector<SimEvent> simEvents;
  SimEvent e(0);
  e.type = FROM_WHATEVER; // partial
  e.x = genXyz0.x();
  e.y = genXyz0.y();
  e.z = genXyz0.z();
  e.t = t0;

  e.pt_hat = 0;
  for(size_t it = 0; it < gen.size(); it++){
    auto & cand = gen[it];
    e.pt_hat += cand.energy();  // not really

    // get the parent vertex (is it the parent's decay or production vertex?)
    double vx = e.x;
    double vy = e.y;
    double vz = e.z;
    bool prompt = cand.isPromptFinalState() || cand.isDirectPromptTauDecayProductFinalState();
    if(prompt){
      // use the event vertex as the production point
    }
    if((cand.numberOfMothers() > 0) && (!prompt)){
      auto parent = cand.mother(0); // politically correct, parent in "pruned"
      //std::cout << "cand  " <<  cand.pdgId() << "   parent= " <<  cand.mother(0) << "  pdgid=" << parent->pdgId() << " " << parent->vertex() << " this=" << &cand << std::endl;
      if(!((parent->vx() == 0) && (parent->vy()==0) && (parent->vz()==0))){
	vx = parent->vx();
	vy = parent->vy();
	vz = parent->vz();
	if(pow(vx-e.x,2) + pow(vy-e.y,2) + pow(vz-e.z, 2) > pow(10e-4,2)){
	  prompt = false;
	}
      }
    }

      /*
      cout << "gen " << it << "  pdg " <<  cand.pdgId() <<  "(" << parentpdgid << ") " << " charge= " <<cand.charge() << " pt=" << cand.pt()
	 <<  " phi=" << cand.phi()
	 <<  " theta=" << cand.theta()
	 << " vtx " << cand.vertex() 
	 << " z " << cand.vertex().z()             // all 0
	 << " vx="  << cand.vx()
	 << " vy="  << cand.vy()
	 << " vz="  << cand.vz()
	 << " dxy="  << cand.dxy()
	 << " dz="  << cand.dzError()                  
	 << " chi2="  << cand.vertexChi2()         // always 0
      	 << " ndau"  << cand.numberOfDaughters()  // always = 0
      	 << " nmo" << cand.numberOfMothers()      //  always =1
	 << " isPromptFinalState() " << cand.isPromptFinalState()
	 << " prompt = " << prompt
	 << endl;
      */

    if ( (fabs(cand.eta()) < etaMaxVisible_) && (cand.charge() != 0) ){
      e.sumpt += cand.pt();
      e.sumpt2 += pow(cand.pt(), 2);
      if(cand.pt() > ptMinVisible_){
	e.nGenTrk ++;
	e.pxvis += cand.px();
	e.pyvis += cand.py();
	e.ptvis += cand.pt();
	// we probably can't rely on the production vertex info
	// that makes track parameter matching kind of a gamble
	auto p = SimPart(0, prompt, vx, vy, vz, cand, fBfield_);
	e.parts.push_back(p);
      }
    }
  }
  
  simEvents.push_back(e);
  return simEvents;
}//getSimEvents_miniaod
/********************************************************************************************************/



/********************************************************************************************************/
std::vector<PrimaryVertexAnalyzer4PU::SimEvent> PrimaryVertexAnalyzer4PU::getSimEvents_simtrks
(
 const Handle<SimTrackContainer> simTrks,
 const Handle<SimVertexContainer> simVtxs,
 Tracks& tracks
 )
// SimTrack + SimVertex info for the signal vertex only, old
// assumes that all sim tracks belong to the signal vertex
/********************************************************************************************************/{
  report_counted("PrimaryVertexAnalyzer4PU::getSimEvents from SimTracks + SimVertexs",1);
  SimVertexContainer::const_iterator vsim = simVtxs->begin(); 
  vector<SimEvent> simEvt;
  SimEvent e(0); 
  e.type = FROM_HEPMC;
  e.eventId = vsim->eventId();
  e.nChTP = 0;
  e.ptvis = 0;
  e.sumpt = 0;
  e.sumpt2 = 0;
  e.pxvis = 0;
  e.pyvis = 0;
  e.sumpt2 = 0;
  e.dzmin = 999.;
  e.x = vsim->position().x() * simUnit_;
  e.y = vsim->position().y() * simUnit_;
  e.z = vsim->position().z() * simUnit_;
  e.t = vsim->position().t() * simtUnit_;

  simEvt.push_back(e);
  for (unsigned int i = 0; i < tracks.size(); i++) {
    MTrack& t = tracks[i];
    simEvt[0].trkidx.push_back(i);
    if (t.selected()) {
      simEvt[0].rtk.push_back(t);
    }

    // fill the matching info for the MTrack
    t._matched = 0;
    t._simEvt = &simEvt[0];
    t._zsim = simEvt[0].z;
    if (f4D_) {
      t._tsim = simEvt[0].t;
    }
    //if (ipdist < 10)  // originated less than 10 um from the interaction point
    //  {
        t._is_primary = true;
    //  }
  }

  return simEvt;
}//getSimEvents_simtrks
/********************************************************************************************************/






/********************************************************************************************************/
std::vector<PrimaryVertexAnalyzer4PU::SimEvent> PrimaryVertexAnalyzer4PU::getSimEvents_tp(
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
      e.type = FROM_TRACKING_TRUTH;  //full
      e.eventId = it->eventId();
      e.nChTP = 0;
      e.sumpt = 0;
      e.sumpt2 = 0;
      e.dzmin = 999.;
      const TrackingVertex* parentVertex = it->parentVertex().get();
      if (DEBUG_) {
        cout << "getSimEvents_tp: ";
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
    MTrack& tk = tracks[i];

    if (select_truthMatchedTrack(tracks.ref(i), tk._tpr)) {
      if (eventIdToEventMap.find(tk._tpr->eventId()) == eventIdToEventMap.end()) {
        if (tk._tpr->eventId().bunchCrossing() == 0) {
          cout << "Bug in getSimEvents, missing event ? " << tk._tpr->eventId().bunchCrossing() << ","
               << tk._tpr->eventId().event() << endl;
        } else if (DEBUG_) {
          cout << "track without SimEvent " << tk._tpr->eventId().bunchCrossing() << "," << tk._tpr->eventId().event()
               << endl;
        }
        //break;
        continue;
      }

      unsigned int event = eventIdToEventMap[tk._tpr->eventId()];
      simEvt[event].trkidx.push_back(i);
      const TrackingVertex* parentVertex = tk._tpr->parentVertex().get();
      double vx = parentVertex->position().x();  // problems with tpr->vz()
      double vy = parentVertex->position().y();
      double vz = parentVertex->position().z();
      double ipdist = sqrt(pow(simEvt[event].x - vx, 2) + pow(simEvt[event].y - vy, 2) + pow(simEvt[event].z - vz, 2)) *
                      1.e4;  // in um

      // fill the matching info for the MTrack
      tk._matched = event;
      tk._simEvt = &simEvt[event];
      tk._zsim = simEvt[event].z;
      if (f4D_) {
        tk._tsim = simEvt[event].t;
      }
      if (ipdist < 10)  // originated less than 10 um from the interaction point
      {
        tk._is_primary = true;
      }

      if (tk.selected()) {
        simEvt[event].trkref.push_back(tracks.ref(tk.index()));  // deprecated ?
	
        simEvt[event].rtk.push_back(tk);
        if (tk.is_primary())
          simEvt[event].rtkprim.push_back(tk);
        if ((tk.ip() / tk.dip()) < 4) // should this cut be configurable? or at least a global constant?
          simEvt[event].rtkprimsel.push_back(tk);
      }
      
    } else {
      // track not truth matched
      tk._matched = NOT_MATCHED_TK_SIM;
      tk._simEvt = NULL;
      tk._zsim = 0;
      tk._tsim = 0;
      tk._is_primary = false;
    }

  }  // end of track loop


  for (unsigned int i = 0; i < simEvt.size(); i++) {
    if (simEvt[i].rtkprim.size() > 0) {
      getTc(simEvt[i].rtkprimsel, simEvt[i].Tc, simEvt[i].chisq, simEvt[i].dzmax, simEvt[i].dztrim);
      simEvt[i].zfit = -99;
    } else {
      simEvt[i].Tc = 0;
      simEvt[i].chisq = 0;
      simEvt[i].dzmax = 0;
      simEvt[i].dztrim = 0;
      simEvt[i].zfit = -99;
      simEvt[i].tfit = -99;
    }

    auto dzmin = -1.;
    for (uint vj = 0; vj < simEvt.size(); ++vj) {
      if(i == vj) continue;
      auto const dz = std::abs(simEvt[i].z - simEvt[vj].z);
      if(dzmin < 0. or dz < dzmin) dzmin = dz;
    }
    simEvt[i].dzmin = dzmin < 0. ? 999. : dzmin;

    if (DEBUG_) {
      cout << setw(3) << i << " )   nTP=" << setw(4) << simEvt[i].tp.size() << "   z=" << setw(8) << setprecision(4)
           << fixed << simEvt[i].z << "    recTrks=" << setw(3) << simEvt[i].rtk.size() << "    recTrksPrim=" << setw(3)
           << simEvt[i].rtkprim.size() << "    allTrks=" << setw(3) << simEvt[i].trkidx.size() << endl;
    }
  }

  return simEvt;
}//getSimEvents_tp
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
}// vertexFilter
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
  auto const & tTopoH = iSetup.getHandle(trackerTopologyToken_);
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


  // select pu window if requested
  assert(puInfo.getPU_zpositions().size() == simPU_);

  if (bPuInfo && ((puInfo.getPU_zpositions().size() < nPUmin_) || (puInfo.getPU_zpositions().size() > nPUmax_))) {
    if (verbose_) {
      cout << "skipping event, out of pu window  npu=" << puInfo.getPU_zpositions().size() << endl;
    }
    return;
  }


   // get the beamspot into global "vertexBeamSpot_"
  if (!get_beamspot_data(iEvent)) {
    return;
  }


  // load vertex collections
  for (auto label : vertexCollectionLabels_){

    if (label.find("WithBS") != std::string::npos){
      selNdof_ = selNdofWithBS_;
      ndof0trk_ = 0;
    }else{
      selNdof_ = selNdofNoBS_;
      ndof0trk_ = -3.;
    }

    Handle<reco::VertexCollection> recVtxsH;
    if (iEvent.getByToken(vertexCollectionTokens_[label], recVtxsH)) {
      recVtxs_[label] = vertexFilter(recVtxsH, useVertexFilter_);
      analyzeVertexRecoCPUTime(histograms_[label], recVtxs_[label], label);
      analyzeVertexCollectionRecoNoTracks(histograms_[label], recVtxs_[label], label);
    } else {
      recVtxs_[label] = NULL;
      cout << "collection " << label << " not found " << endl;
    }
    
  }


  // load the tracks (or in case of miniaod tracks and vertices)
  Tracks tracks;
  
  //  get transient track builder and b-field
  //iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", theB_);
  theB_ = iSetup.getHandle(transientTrackBuilderToken_);
  fBfield_ = ((*theB_).field()->inTesla(GlobalPoint(0., 0., 0.))).z();


  if (MINIAOD_){

    /* miniaod tracks (candidates) have a link to a vertex
       which one is it ?
       for now assume it's the first one in our list
    */
    std::string miniaod_vertexcollection_label = *vertexCollectionLabels_.begin();

    cout << " got "  << recVtxs_[miniaod_vertexcollection_label]->size()
	 << " miniaod vertices from collection " << miniaod_vertexcollection_label << endl;

    get_miniaod_tracks(iSetup, iEvent, miniaod_vertexcollection_label, tracks);
    cout << " got mini aod tracks : " << tracks.size() << endl;

  }else if (!get_reco_and_transient_tracks(iSetup, iEvent, tracks)) {
    
    // clean up and exit if nothing is found
    for (auto label : vertexCollectionLabels_){
      delete recVtxs_[label];
    }
    return;  //FIXME some things can be done without tracks
  }


  // collect MC information
  std::vector<SimEvent> simEvt;
  if (run_ < 1000) {
    get_particle_data_table(iSetup);
      
    MC_ = get_MC_truth(iEvent, tracks, bPuInfo, puInfo, simEvt);

    if (MC_){
      if(verbose_) { cout << "MC info found" << endl;}
      fill_simvtx_histos(simEvt);
    }else{
      if (verbose_) {cout << "no MC info found" << endl;}
    }
  }

  /* FIXME re-enable later
  add_timing_to_vertex_collection(vertexCollectionLabels_[0], tracks);  // no refit
  if (frefit_) {
    //refit_recvertices_after_timing(tracks);
    mass_constrained_multi_vertex_time_from_tracks_pid(vertexCollectionLabels_[0], tracks, 0.8, true);
  }
  */


  fill_track_to_vertex_pointers(tracks);  // must do this here in case refitting rejects vertices



  // analyze the vertex collections
  
  for( auto label : vertexCollectionLabels_){

    if (recVtxs_[label] == NULL){
      report_counted("skipping empty vertex collection " + label,5);
      continue;
    }
    
    if (verbose_) {
      cout << "**** analyzing " << label << "  with size " << recVtxs_[label]->size() << endl;
    }

    // create MVertexCollections
    vertexes_[label] = MVertexCollection(label, recVtxs_[label], tracks);
    auto & vertexes = vertexes_[label];
    auto & histos = histograms_[label];

    // replace by set_ndof_globals(label);
    if (label.find("WithBS") !=std::string::npos){
      selNdof_ = selNdofWithBS_;
      ndof0trk_ = 0;
    }else{
      selNdof_ = selNdofNoBS_;
      ndof0trk_ = -3.;
    }

    timer_start("analyzeVertexCollectionReco");
    analyzeVertexCollectionReco(histos, vertexes, tracks, label);
    timer_stop("analyzeVertexCollectionReco");
    
    if (MC_) {

      timer_start("analyzeVertexCollectionSimPvNoSimTracks");
      analyzeVertexCollectionSimPvNoSimTracks(histos, vertexes, tracks, simEvt, label);
      timer_stop("analyzeVertexCollectionSimPvNoSimTracks");
      
      if(tracking_truth_available_){
	
	timer_start("tp-matching");
	tpmatch(vertexes, simEvt, tracks);    // fill the rec-to-sim matching matrix (vertex.wos[simevent])
	wos_match(vertexes, simEvt, tracks);  // make a one-to-one assignment
	timer_stop("tp-matching");

	timer_start("analyzeMergeRateTP");
	analyzeVertexMergeRateTP(histos, vertexes, tracks, simEvt, label); //recvmatch_[label]
	timer_stop("analyzeMergeRateTP");

        timer_start("analyzeVertexCollectionTP");
	analyzeVertexCollectionTP(histos, vertexes, tracks, simEvt, label);
	timer_stop("analyzeVertexCollectionTP");

	analyzeVertexCollectionZmatching(histos, vertexes, simEvt, label, zwindow_sigmas_);

	analyzeVertexCollectionPtvis(histos, vertexes, tracks, simEvt,label);
      } else {
	timer_start("analyzeVertexCollectionSimTracks");
	analyzeVertexCollectionSimTracks(histos, vertexes, tracks, simEvt, label);
	timer_stop("analyzeVertexCollectionSimTracks");
      }
    }
  }

  if (do_vertex_analysis_ && ((nCompareCollections_ < 0) || (eventcounter_ <= nCompareCollections_))) {
    cout << " comparing collections " << endl;
    compareCollections(simEvt);
  }

  // print summary info
  if ((dumpThisEvent_ && (dumpcounter_ < ndump_)) || (veryverbose_ && (eventcounter_ < ndump_)) ||
      (autoDumpCounter_-- > 0) || (forceDump_))
    {
      std::cout  << "  dumpThisEvent_ " << dumpThisEvent_
		 << "  dumpcounter "  << dumpcounter_ << " ndump " << ndump_
		 << "  veryverbose " << veryverbose_
		 << "  autoDumpCounter " << autoDumpCounter_
		 << "  forceDump " << forceDump_
		 << std::endl;
      dumpEventSummary(simEvt, tracks);
    }

  // clean up
  for ( auto label : vertexCollectionLabels_ ){
    delete recVtxs_[label];
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
void PrimaryVertexAnalyzer4PU::dumpEventSummary(std::vector<SimEvent> & simEvt, Tracks & tracks)
/***************************************************************************************/
{
// print summary info
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
    bool bsdumped = false;

    for ( auto label : vertexCollectionLabels_ ){

      if (recVtxs_[label] == NULL) continue;
      
      cout << " dumping collection " << label  << "   run:event =" << run_ <<":"<< event_ << endl;
      auto & vertexes = vertexes_[label];
      printRecVtxs( vertexes );
      
      // redo the matching for the dump (invalid if another collection has been analyzed)
      if(tracking_truth_available_){
	tpmatch(vertexes, simEvt, tracks);
	wos_match(vertexes, simEvt, tracks);
      }else{
	//crudematch(vertexes, simEvt); // TODO, just for printouts
      }

      if (matchsummaries_ > 0)
	  printMatchingSummary(vertexes, simEvt, label);
      
      //printEventSummary_notp(vertexes, tracks, simEvt, label);
      
      
      if (!trksdumped) {
	printPVTrksZT(tracks, vertexes, simEvt);
	cout << "---" << endl;
	trksdumped = true;  // only dump one track list per event
      }
      
      if (!bsdumped &&(dumpcounter_ < 2)) {
	cout << "beamspot " << vertexBeamSpot_ << endl;
	bsdumped = true;
      }
      
      if (verbose_)
      cout << endl << endl;
    }

    matchsummaries_--;

}//dumpEventSummary
/***************************************************************************************/




/***************************************************************************************/
void PrimaryVertexAnalyzer4PU::printEventSummary_notp(MVertexCollection& vtxs,
						      Tracks& tracks,
						      vector<SimEvent>& simEvt,
						      const string message)
// make a readable summary of rec and sim vertices without tracking particle information
// only using z-ordering
/***************************************************************************************/
{
  bool show_matching = true;
  cout << endl << "printEventSummary  (notp)" << endl;
  cout << "---------------------------" << endl;
  cout << "event counter = " << eventcounter_ << "   " << message << endl;
  cout << "run " << run_ << "    event " << event_ << endl;
  cout << "---------------------------" << endl;

  vector<pair<double, unsigned int>> zrecv;
  for (unsigned int idx = 0; idx < vtxs.size(); idx++) {
    zrecv.push_back(make_pair(vtxs(idx).z(), idx));
  }
  stable_sort(zrecv.begin(), zrecv.end(), lt);

  // same for simulated vertices
  vector<pair<double, unsigned int>> zsimv;
  for (unsigned int idx = 0; idx < simEvt.size(); idx++) {
    zsimv.push_back(make_pair(simEvt[idx].z, idx));
  }
  stable_sort(zsimv.begin(), zsimv.end(), lt);

  double zmatch = 0.05; // just for printing purposes, not used for anything else
  cout.precision(4);

  cout << "   sim z";
  if(f4D_) cout << ",   t   ";
  cout << "     pT^       ";
  cout << "     rec z";
  if(f4D_) cout << ",   t   ";
  if(show_matching) cout << "[  w ]";
  cout   << "                                                        #trk" << endl;

  unsigned int idxrec = 0;         // move these two pointers throught the lists
  unsigned int idxsim = 0;         // to show similar (in z) vertices side-by-side
  while ((idxrec < vtxs.size()) || (idxsim < simEvt.size())) {
    string signal = " ";
    string tag = " ";
    if ((idxsim < simEvt.size()) && (zsimv[idxsim].second == 0)) {
      signal = "*";
    }
    if ((idxrec < vtxs.size()) && (zrecv[idxrec].second == 0)) {
      tag = "*";
    }

    double ndof = 0, ptmax2 = 0, trec=0, wmatch=0, sumpt=0; // pxy=0
    if (idxrec < vtxs.size()) {
      ndof = vtxs(zrecv[idxrec].second).ndof();
      //pxy = vtxs(zrecv[idxrec].second).pxy(vertexBeamSpot_);
      ptmax2 = vtxs(zrecv[idxrec].second).ptmax2();
      sumpt = vtxs(zrecv[idxrec].second).sumpt();
      if(f4D_) trec = vtxs(zrecv[idxrec].second).t();
      wmatch = vtxs(zrecv[idxrec].second).sigwntfrac;
    }

    if ((idxrec < vtxs.size()) && (idxsim < simEvt.size()) &&
        (std::abs(zrecv[idxrec].first - zsimv[idxsim].first) <
         std::max(zmatch, 3 * vtxs(zrecv[idxrec].second).zError())) &&
        (((idxsim + 1) == simEvt.size()) ||
         (fabs(zrecv[idxrec].first - zsimv[idxsim].first) < fabs(zrecv[idxrec].first - zsimv[idxsim + 1].first))) &&
        (((idxrec + 1) == vtxs.size()) ||
         (fabs(zrecv[idxrec].first - zsimv[idxsim].first) < fabs(zrecv[idxrec + 1].first - zsimv[idxsim].first)))) {
      // rec and sim  on this line
      cout << setw(8) << fixed << setprecision(4) << zsimv[idxsim].first << signal;                      // sim z
      if(f4D_) {cout << "," << setw(7) << setprecision(4) << fixed << simEvt[zsimv[idxsim].second].t;}   // sim t
      cout << setw(6) << setprecision(1) << fixed << simEvt[zsimv[idxsim].second].pt_hat;                // sim pt_hat
      cout << "   <->    ";
      cout << setw(8) << setprecision(4) << fixed << zrecv[idxrec].first << tag;                         //rec z
      if(f4D_) {cout << "," << setw(7) << setprecision(4) << fixed << trec;}                             //rec t
      if(show_matching){
	if(wmatch > 0.01){ cout << "[" << setw(4) << setprecision(2) << wmatch << "]";}
	else{cout << "      ";}
      }
      cout << " (ndof=" << fixed << setw(6) << setprecision(1) << ndof                                   // more
	//<< ", p="  << fixed<< setw(6) << setprecision(4) << pxy 
	   << ",PT="  << fixed<< setw(6) << setprecision(1) << sumpt 
           << ", ptmax2="  << fixed<< setw(4) << setprecision(1) << ptmax2 << ")";

      if (simEvt[zsimv[idxsim].second].is_visible()) {
        cout << "            (" << fixed << setw(3) << simEvt[zsimv[idxsim].second].nGenTrk << ")";
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

    } else if (((idxrec < vtxs.size()) && (idxsim < simEvt.size()) && (zrecv[idxrec].first < zsimv[idxsim].first)) ||
               ((idxrec < vtxs.size()) && (idxsim == simEvt.size()))) {
      // rec only on this line
      cout << "                         ";                                                // no sim
      if(f4D_) cout << "        ";                                                         // no sim
      cout << setw(8) << setprecision(4) << fixed << zrecv[idxrec].first << tag;           // rec z
      if(f4D_) {cout << "," << setw(7) << setprecision(4) << fixed << trec;}               // rec t
      if(show_matching){
	if(wmatch > 0.01){ cout << "[" << setw(4) << setprecision(2) << wmatch << "]";}
	else{cout << "      ";}
      }
      cout << " (ndof=" << fixed << setw(6) << setprecision(1) << ndof;                    // more
      //cout << ", p=" << setw(6) << setprecision(4) << pxy;
      cout << ",PT="  << fixed << setw(6) << setprecision(1) << sumpt;
      cout << ", ptmax2=" << setw(4) << setprecision(1) << ptmax2 << ")";
      cout << "   fake " << endl;
      idxrec++;

    } else if (((idxrec < vtxs.size()) && (idxsim < simEvt.size()) && (zrecv[idxrec].first > zsimv[idxsim].first)) ||
               ((idxrec == vtxs.size()) && (idxsim < simEvt.size()))) {
      // sim only on this line
      cout << setw(8) << setprecision(4) <<fixed << zsimv[idxsim].first << signal;                       // sim z
      if(f4D_) {cout << "," << setw(7) << setprecision(4) << fixed << simEvt[zsimv[idxsim].second].t;}   // sim t
      cout << setw(6) << setprecision(1) << fixed << simEvt[zsimv[idxsim].second].pt_hat;                // sim pt_hat
      cout << "          ";
      cout << "          ";
      if(f4D_) {cout << "        ";}
      if (simEvt[zsimv[idxsim].second].is_visible()) {
	if(show_matching) cout << "      ";
        if (zsimv[idxsim].second == 0) {
          cout << "                                       lost signal vertex" << endl;
        } else {
          cout << "                                       lost PU  ";
          if (simEvt[zsimv[idxsim].second].is_visible()) {
            cout << "(" << setw(3) << simEvt[zsimv[idxsim].second].nGenTrk << ")";
            //   << "(" << setw(5) << setprecision(3) << simEvt[zsimv[idxsim].second].pt_hat << ")"
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
  std::cout << std::endl;
}//printEventSummary_notp
/***************************************************************************************/









/***************************************************************************************/
void PrimaryVertexAnalyzer4PU::compareCollections(vector<SimEvent>& simEvt) {
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

  // unmatched rec vertices from all collections
  for (unsigned int icol = 0; icol < ncoll; icol++) {
    auto vtxs = vertexes_[vertexCollectionLabels_[icol]];
    for (unsigned int idx = 0; idx < vtxs.size(); idx++) {
      if (vtxs(idx).isRecoFake()) continue;
      if (vtxs(idx).sim == NOT_MATCHED_VTX_SIM) {
          unsigned int* v = new unsigned int[ncoll + 1];
          for (unsigned int j = 0; j < ncoll + 1; j++) {
            v[j] = NOT_MATCHED_VTX_REC;
          }
          v[icol] = idx;
          row.push_back(make_pair(vtxs(idx).z(), v));
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
            if ((row[irow].second[i] != NOT_MATCHED_VTX_REC) && (row[irow + 1].second[i] != NOT_MATCHED_VTX_REC))
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
          if ((row[join].second[i] == NOT_MATCHED_VTX_REC) && (row[join + 1].second[i] != NOT_MATCHED_VTX_REC)) {
            row[join].second[i] = row[join + 1].second[i];
            if (i == ncoll)
              row[join].first = row[join + 1].first;
          }
        }

        //row z
        if (row[join].second[ncoll] == NOT_MATCHED_VTX_REC) {
          double zrow = 0;
          int nv = 0;
          for (unsigned int i = 0; i < ncoll; i++) {
            int iv = row[join].second[i];
            if (iv > int(vertexes_[vertexCollectionLabels_[i]].size())) {
              cout << "illegal vertex index " << iv << "    join=" << join << endl;
            }else{
	      if (iv >= 0) {
		zrow += vertexes_[vertexCollectionLabels_[i]].at(iv).z();
		nv++;
	      }
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
        v[j] = NOT_MATCHED_VTX_REC;
      }
      for (unsigned int jcol = 0; jcol < ncoll; jcol++) {
        auto vtxs = vertexes_[vertexCollectionLabels_[jcol]];
	int i = NOT_MATCHED_VTX_REC;
	for (unsigned int j = 0; j < vtxs.size(); j++) {
	  if (vertexes_[vertexCollectionLabels_[jcol]][j].sim == idx) {
	    i = j;
	  }
	}
	v[jcol] = i;
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

    if (idx0 == NOT_MATCHED_VTX_REC) {
      cout << "%                    ";
    } else {
      // sim vertex
      cout << fixed << setw(10) << setprecision(4) << z << " [" << setw(3) << simEvt[idx0].nChTP << "," << setw(3)
	   << simEvt[idx0].rtk.size() << "]";
      if (idx0 == 0) {
	cout << "*";
      } else {
	cout << " ";
      }
    }

    // count collections that  have a rec vertex for this sim vertex (for reporting)
    unsigned int nrec = 0;
    for (unsigned int jcol = 0; jcol < ncoll; jcol++) {
      if (v[jcol] != NOT_MATCHED_VTX_REC) {
        nrec++;
      }
    }
    if ((nrec > 0) && (nrec < ncoll)) {
      differences++;
    }

    for (unsigned int jcol = 0; jcol < ncoll; jcol++) {
      unsigned int idx = v[jcol];
      if (idx != NOT_MATCHED_VTX_REC) {
        auto vtxs = vertexes_[vertexCollectionLabels_[jcol]];
        if (idx < vtxs.size()){
          cout << fixed << setw(10) << setprecision(4) << vtxs(idx).z() << " (" << setw(5) << setprecision(1)
               << vtxs(idx).ndof() << ")";
        } else {
          cout << "                  ";
        }
      } else {
        if (idx0 == NOT_MATCHED_VTX_REC) {
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








/*****************************************************************************************************/
std::string PrimaryVertexAnalyzer4PU::formatMatchList(const std::map<unsigned int, double>& valuemap,
                                                      unsigned int nfield,
                                                      bool sim) 
/*****************************************************************************************************/
{
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
/*************formatMatchList**************************************************************************/






/***************************************************************************************/
void PrimaryVertexAnalyzer4PU::printMatchingSummary(MVertexCollection & vertexes,
                                                    vector<SimEvent>& simEvt,
                                                    const string message)
// dump details about the matching, tracking particles only!
/***************************************************************************************/
{
  cout << endl << "printMatchingSummary  (simEvt)" << endl;
  cout << "---------------------------" << endl;
  cout << "event counter = " << eventcounter_ << "   " << message << endl;
  cout << "run " << run_ << "    event " << event_ << endl;

  cout << "     z        t   [sim]ntrk   wnt(rec)                             z        t   (rec) ndof   wnt [sim]" << endl;
  //  dump matched and unmatched vertices here, then sort in z
  vector<pair<double, std::string>> row;
  vector<unsigned int> dumped_rec;
  
  // unmatched rec
  for (unsigned int idxrec = 0; idxrec < vertexes.size(); idxrec++) {
    if (vertexes(idxrec).isRecoFake()) continue;
    if (vertexes(idxrec).sim == NOT_MATCHED_VTX_SIM) {
      std::ostringstream s;
      s << "%" << setw(61) << " " << setw(9) << setprecision(4) << fixed << vertexes(idxrec).z() << "," << setw(7)
        << setprecision(3) << fixed << vertexes(idxrec).t() << " " << setw(1) << "(" << setw(3) << idxrec << ") "
	<< setprecision(1) << fixed << setw(5) << vertexes(idxrec).ndof() << " "
        << formatMatchList(vertexes(idxrec).wnt, 3, false);
      row.push_back(make_pair(vertexes(idxrec).z(), s.str()));
      dumped_rec.push_back(idxrec);
    }
  }

  // unmatched sim
  for (unsigned int idxsim = 0; idxsim < simEvt.size(); idxsim++) {
    if ((simEvt[idxsim].is_visible()) && (simEvt[idxsim].rec == NOT_MATCHED_VTX_REC)) {
      std::ostringstream s;
      if (idxsim == 0){
	s << "*";
      }else{
	s << " ";
      }
      s << setw(8) << setprecision(4) << fixed << simEvt[idxsim].z << "," << setw(7) << setprecision(3) << fixed
        << simEvt[idxsim].t << " " << setw(1) << "[" << setw(3) << idxsim << "]"
	<< setw(4) << simEvt[idxsim].rtk.size() << " "
        << formatMatchList(simEvt[idxsim].wnt, 3, true)
	<< "%";
      if (idxsim == 0 ) {
	s << "      signal vertex not matched   ";
      }
      row.push_back(make_pair(simEvt[idxsim].z, s.str()));
    }
  }

  // matched rec + sim
  for (unsigned int idxsim = 0; idxsim < simEvt.size(); idxsim++) {
    if (simEvt[idxsim].rec != NOT_MATCHED_VTX_REC) {
      unsigned int idxrec = simEvt[idxsim].rec;
      if( vertexes(idxrec).sim != idxsim ){ cout << "crash and burn :  " << idxrec << "," << vertexes(idxrec).sim << "," << idxsim << endl;}
      assert(vertexes(idxrec).sim == idxsim);

      std::ostringstream s;
      if (idxsim == 0){
	s << "*";
      }else{
	s << " ";
      }
      s << setw(8) << setprecision(4) << fixed << simEvt[idxsim].z << "," << setw(7) << setprecision(3) << fixed
        << simEvt[idxsim].t << " " << setw(1) << "[" << setw(3) << idxsim << "]"
	<< setw(4) << simEvt[idxsim].rtk.size() << " "
        << formatMatchList(simEvt[idxsim].wnt, 3, true) << setw(9) << setprecision(4) << fixed
        << vertexes(idxrec).z() << "," << setw(7) << setprecision(3) << fixed << vertexes(idxrec).t() << setw(1)
        << "(" << setw(4) << idxrec << ") "
	<< setprecision(1) << fixed << setw(5) << vertexes(idxrec).ndof() << " "
	<< formatMatchList(vertexes(idxrec).wnt, 3, false);

      dumped_rec.push_back(idxrec);
      row.push_back(make_pair(vertexes(idxrec).z(), s.str()));
    }
  }

  // what have we missed?
  // unmatched sim
  for (unsigned int idxrec = 0; idxrec < vertexes.size(); idxrec++) {
    if (vertexes(idxrec).isRecoFake()) continue;
    if (std::find(dumped_rec.begin(), dumped_rec.end(), idxrec) == dumped_rec.end()){
      unsigned int idxsim = vertexes(idxrec).sim;
      
      cout << "missed in printMatchingSummary  "  << idxrec <<  "  " << vertexes(idxrec).z() << "  match=" << vertexes(idxrec).sim;
      std::ostringstream s;
      s << "@" << setw(61) << " " << setw(9) << setprecision(4) << fixed << vertexes(idxrec).z() << "," << setw(7)
        << setprecision(3) << fixed << vertexes(idxrec).t() << " " << setw(1) << "(" << setw(3) << idxrec << ") "
        << formatMatchList(vertexes(idxrec).wnt, 3, false);
      row.push_back(make_pair(vertexes(idxrec).z(), s.str()));
      if(idxsim != NOT_MATCHED_VTX_SIM){
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
/*************printMatchingSummary******************************************************/







/***************************************************************************************/
void PrimaryVertexAnalyzer4PU::printEventSummary_tp(std::map<std::string, TH1*>& h,
  					    MVertexCollection & vtxs,
                                                 Tracks& tracks,
                                                 vector<SimEvent>& simEvt,
                                                 const string message)
  // make a readable summary of the vertex finding if the TrackingParticles are availabe
/***************************************************************************************/  
{
  if (simEvt.size() == 0) {
    return;
  }

  if (eventSummaryCounter_++ > nEventSummary_){
    report_counted(message + " count limit exceeded", 1);
    return;
  }
  // sort vertices in z ... for nicer printout

  vector<pair<double, unsigned int>> zrecv;
  for (unsigned int idx = 0; idx < vtxs.size(); idx++) {
    if (vtxs(idx).isRecoFake()) continue;
    zrecv.push_back(make_pair(vtxs(idx).z(), idx));
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

  cout << endl << "printEventSummary (tracking truth)" << endl;
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
    svmatch[itrec->second] = NOT_MATCHED_VTX_SIM;
    purity[itrec->second] = 0.;
    wpurity[itrec->second] = 0.;
  }

  if (zrecv.size() > 20) {

    cout << " full dump dropped because we have more than 20 vertices   nrec = " << zrecv.size() << endl;
    
    
    for (vector<pair<double, unsigned int>>::iterator itsim = zsimv.begin(); itsim != zsimv.end(); itsim++) {
      // itsim->first = z of the simvx
      // itsim->second= index of the simvtx
      if ((itsim->first < zmin) || (itsim->first > zmax))
        continue;
      SimEvent* ev = &(simEvt[itsim->second]);
      rvmatch[itsim->second] = NOT_MATCHED_VTX_REC;

      nmatch[itsim->second] = 0;   // highest number of truth-matched tracks of ev found in a recvtx
      wnmatch[itsim->second] = 0;  // highest sum of weights of truth-matched tracks of ev found in a recvtx
      //	double matchpurity=0;//,matchwpurity=0;//,matchpurity2=0;

      // compare this sim vertex to all recvertices:
      for (vector<pair<double, unsigned int>>::iterator itrec = zrecv.begin(); itrec != zrecv.end(); itrec++) {
        // itrec->first  = z coordinate of the recvtx
        // itrec->second = index of the recvtx
        unsigned int irecvtx = itrec->second;
	const auto v = & vtxs(irecvtx);

        // count tracks found in both, sim and rec
        double n = 0, wt = 0;
	for( auto tv : v->tracks ){
          if (ev->hasTrack(tv)) {
            n++;
            wt += v->trackWeight(tv);
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
      cout << setw(7) << fixed << vtxs(itrec->second).tracksSize();
      if (itrec->second == 0) {
        cout << "*";
      } else {
        cout << " ";
      }
    }
    cout << "   rec tracks" << endl;

    if(tracking_truth_available_){
    
    // truthMatchedVertexTracks[irecvtx]= number of rec tracks in that vertex for which we have truth matched simtracks
    // (not necessarily all from the same simvertex)
    map<unsigned int, int> truthMatchedVertexTracks;
    
    for (vector<pair<double, unsigned int>>::iterator itrec = zrecv.begin(); itrec != zrecv.end(); itrec++) {
      truthMatchedVertexTracks[itrec->second] =
	getTruthMatchedVertexTracks_obsolete(vtxs(itrec->second).recovertex(),min_trk_in_vtx_weight_).size();  // FIXME replace consistently
      cout << "DDDDD " << truthMatchedVertexTracks[itrec->second] << " =?=" << vtxs(itrec->second).countTruthMatchedTracks(min_trk_in_vtx_weight_);
    }
    
    cout << "                        ";
    for (vector<pair<double, unsigned int>>::iterator itrec = zrecv.begin(); itrec != zrecv.end(); itrec++) {
      if ((itrec->first < zmin) || (itrec->first > zmax))
	continue;
      cout << setw(7) << fixed << vtxs[itrec->second].sumwnt
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
      rvmatch[itsim->second] = NOT_MATCHED_VTX_REC;
      
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
        const auto v = &(vtxs(irecvtx));

        // count tracks found in both, sim and rec
        double n = 0, wt = 0;
	for( auto tv : v->tracks ){
          if (ev->hasTrack(tv)) {
            n++;
            wt += v->trackWeight(tv);
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
        } else if (ev->rec != NOT_MATCHED_VTX_REC) {
          cout << "poor,  quality " << ev->matchQuality << "  for zrec " << vtxs(ev->rec).z() << endl;
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
      if (vtxs[itrec->second].sumwnt > 0) {
        cout << setw(7) << fixed << vtxs[itrec->second].maxwnt / float(vtxs[itrec->second].sumwnt);
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

      if ((svmatch[itrec->second] != NOT_MATCHED_VTX_SIM) && (rvmatch[svmatch[itrec->second]] == itrec->second)) {
        if ((purity[itrec->second] / truthMatchedVertexTracks[itrec->second]) >= 0.5) {
          cout << "   ok   ";
        } else {
          cout << "  ugly  ";
        }
      } else {
        if (vtxs[itrec->second].split_from() >= 0) {
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
      if (vtxs[itrec->second].maxwos > 0) {
        cout << setw(8) << fixed << simEvt[vtxs[itrec->second].wosmatch].z;
      } else if (vtxs[itrec->second].maxwnt > 0) {
        cout << setw(8) << fixed << simEvt[vtxs[itrec->second].wntmatch].z;
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

    unsigned int idxrec = itrec->second;

    if ((svmatch[idxrec] == NOT_MATCHED_VTX_SIM) || (rvmatch[svmatch[idxrec]] != idxrec)) {
      cout << "zrec=" << setw(9) << fixed << setprecision(4) << itrec->first << "  ";
      if (vtxs[idxrec].maxwos > 0) {
        cout << "maxwos = " << setw(10) << fixed << setprecision(1) << vtxs[idxrec].maxwos << " / " << setw(10)
             << fixed << setprecision(1) << vtxs[idxrec].sumwos << "   zwosmatch =" << setw(9) << fixed
             << setprecision(4) << simEvt[vtxs[idxrec].wosmatch].z << "  ";
      }
      if (vtxs[idxrec].maxwnt > 0) {
        cout << "maxnt = " << setw(3) << setprecision(0) << vtxs[idxrec].maxwnt << " /" << setw(3) << setprecision(0)
             << vtxs[idxrec].sumwnt << "  z=" << setw(8) << fixed << setprecision(4)
             << simEvt[vtxs[idxrec].wntmatch].z << "  ";
      }
      cout << "  fake=" << vtxs[idxrec].is_fake();
      cout << "  qual=" << vtxs[idxrec].matchQuality;
      cout << "  split=" << vtxs[idxrec].split_from();
      cout << endl;
    }
  }

  if (!fillHistograms)
    return;

  // FIXME this code needs an overhaul
  // list problematic tracks
  /*
  for (vector<pair<double, unsigned int>>::iterator itsim = zsimv.begin(); itsim != zsimv.end(); itsim++) {
    SimEvent* ev = &(simEvt[itsim->second]);

    //for (vector<TransientTrack>::iterator te = ev->tk.begin(); te != ev->tk.end(); te++) {
    for (auto tk : ev->rtk){
      const reco::Track& RTe = *tk._trk;

      unsigned int ivassign = NOT_ASSIGNED;  // will become the index of the vertex to which a track was assigned

      for (vector<pair<double, unsigned int>>::iterator itrec = zrecv.begin(); itrec != zrecv.end(); itrec++) {
        const reco::Vertex* v = &(vertexes->at(itrec->second));

        for (trackit_t tv = v->tracks_begin(); tv != v->tracks_end(); tv++) {
          const reco::Track& RTv = *(tv->get());
          if (RTe.vz() == RTv.vz()) {
            ivassign = itrec->second;
          }
        }
      }

      double tantheta = tan(tk.theta());
      //      reco::BeamSpot beamspot = (tk.transientTrack()->stateAtBeamLine()).beamSpot();
      //double dz2 = pow(RTe.dzError(), 2) + pow(beamspot.BeamWidthX() / tantheta, 2);
      double dz2 = pow(RTe.dzError(), 2) + pow(vertexBeamSpot_.BeamWidthX() / tantheta, 2);

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
  */
}
/***************************************************************************************/








/***************************************************************************************/
void PrimaryVertexAnalyzer4PU::signalvtxmatch(MVertexCollection & vertexes, std::vector<SimEvent> & simEvt)
// if we have MC info for the signal vertex, match it to the reco vertices
// as usual, assume that the first vertex in simEvt is the signal vertex
/***************************************************************************************/
{
  if (!MC_ || simEvt.size() < 1){
    return;
  }

  auto signal = simEvt.at(0);

  vector<double> score(vertexes.size());

  for( auto v : vertexes){
    unsigned int iv = v.index();
    double dz = signal.z - v.z();
    if( fabs(dz) > 1.0){
      score[iv] = 0;
      continue;
    }
    
    score[iv] = exp(-0.5 * pow(dz / v.zError(),2));
    
    //if(f4D_ && vertexes.label.find("4D")!= std::string::npos){
    //  dt = 
    //}

    //for(auto tv : v.tracks){
    //      
    // }

  }

}
/***************************************************************************************/






/***************************************************************************************/
 void PrimaryVertexAnalyzer4PU::tpmatch(MVertexCollection & vtxs,
					std::vector<SimEvent>& simEvt,
					Tracks& tracks)
  // collects truth-matched track information for matching rec and sim vertices
  // two different matches are used to determine the dominant sim vertex in a recvertex
  // wos match  = sum of "weight-over-sigma**2"
  // ntrk match = number of tracks (implicit with weight > min_trk_in_vtx_weight_), weighted with pt (cut-off at 1Gev)
  // this information is filled into both, the MVertex objects and the SimEvents
  // a rec-to-sim match is made by identifying the the dominant sim vertex
  //   (wosmatch and wntmatch)
  // note: this is not a one-to-one match, a simvertex can dominate multiple recvtxs
  // simEvt.nwosmatch  and simEvt.nwtmatch keep count the number of recvertices dominated by a simvertex
/***************************************************************************************/
{

  // turn on debugging messages for one rec vertex and one sim vertex
  bool DEBUG = false;
  unsigned int ivdebug = 0;  // rec
  unsigned int ievdebug = 0; // sim
  if(DEBUG){ cout << "tpmatch  run:event = " << run_ << ":" << event_ << "  ivdebug=" << ivdebug << " ievdebug=" << ievdebug << endl;}


  // clear old sim->rec pointers, only valid for one vertex collection
  for (auto & ev : simEvt){
    ev.clear_matching_info();
  }

  for (unsigned int iv = 0; iv < vtxs.size(); iv++) {
    auto & v = vtxs[iv];
    v.clear_matching_info();
    
    if (v.isRecoFake()) continue;

    unsigned int iev = 0;  // simevent index
    for (auto & ev : simEvt){
      double evwos = 0;       // weight over sigma**2  of simEvt ev in the current recvtx
      double evwnt = 0;       // weighted number of tracks
      unsigned int evnt = 0;  // number of tracks from simEvt ev in the current recvtx

      assert(iev == ev.index);

      for (auto tk : v.tracks){

	double wt = v.trackWeight(tk);
	if(wt < min_trk_in_vtx_weight_) continue;  

        if (tk->is_matched() && (tk->_simEvt == &ev)) {
          double dz2_beam = pow(vertexBeamSpot_.BeamWidthX() * cos(tk->phi()) / tan(tk->theta()), 2) +
	    pow(vertexBeamSpot_.BeamWidthY() * sin(tk->phi()) / tan(tk->theta()), 2);
          double dz2 = pow(tk->dzError(), 2) + dz2_beam + pow(0.0020, 2); // added 20 um, some tracks have crazy small resolutions
          double wos = v.trackWeight(tk) / dz2;
          if (f4D_ && tk->has_timing() && (tk->dt() > 0)) {
            wos = wos / erf(tk->dt() / sigmaT_);
          }
          double wnt = v.trackWeight(tk) * min(tk->pt(), 1.0);  // (truncated-)pt-weighted track count (downweights pt < 1 GeV tracks)
          v.add_truthmatched_track(iev, wos, wnt);   // fill track(wos)  rec vtx <- sim vtx, increments sumnt and sumwos
          ev.addTrack(iv, wos, wnt);                 // fill track(wos)  sim vtx -> rec vtx
          evwos += wos;
          evwnt += wnt;
          evnt++;
	  if (DEBUG && (iev==ievdebug)){
	    cout << "tpmatch iev=" << iev << "  iv=" << iv << "    wos=" << wos << "   evwos= " << evwos << "   dz2= "  << dz2 << "  w=" <<  v.trackWeight(tk) << endl;
	  }
	  
        }
      }
      if(DEBUG && (iv == ivdebug) && (iev == ievdebug)){
	cout << "tpmatch rec " << iv << "  sim " << iev << "  evwos = " << evwos<< " v.maxwos=" << v.maxwos << "   evnt=" << evnt << endl;
      }
      // require 2 tracks for a wos-match
      if ((evwos > 0) && (evwos > v.maxwos) && (evnt > 1)) {
        v.wosmatch = iev;
        v.maxwos = evwos;
        v.maxwosnt = evnt;
      }

      // weighted track counting match, require at least one track
      if ((evnt > 0) && (evwnt > v.maxwnt)) {
        v.wntmatch = iev;
        v.maxwnt = evwnt;
      }

      iev++;
      
    } // end of simevent loop

    // now wosmatch  holds the index of the dominant simvertex for this recvtx and maxwos the wos value
    if(DEBUG && (iv==ivdebug)){
      cout << "tpmatch: recvtx " << ivdebug << ", v.maxwos=" << v.maxwos << "   v.wosmatch="  << v.wosmatch << endl;
    }

    if (v.maxwos > 0) {
      simEvt.at(v.wosmatch).nwosmatch++;  // count the recvertices dominated by a simvertex
      simEvt.at(v.wosmatch).wos_dominated_recv.push_back(iv);
      assert(iv < vtxs.size());
    }

    if (v.maxwnt > 0) {
      simEvt.at(v.wntmatch).nwntmatch++;  // count the recvertices dominated by a simvertex
    }

  }  // end of vertex loop


  // for comparison with simTrack like matching (signal vertex-only)
  double sumsigwntfrac = 0;
  for(auto & v : vtxs){
    if (v.wos.find(0) == v.wos.end()){
      v.sigwosfrac = 0;
    }else{
      assert(v.wos[0] == simEvt.at(0).wos[v.index()]);
      v.sigwosfrac = v.wos[0] / simEvt.at(0).sumwos;
    }

    if (v.wnt.find(0) == v.wnt.end()){
      v.sigwntfrac = 0;
    }else{
      assert(v.wos[0] == simEvt.at(0).wos[v.index()]);
      v.sigwntfrac = v.wnt[0] / simEvt.at(0).sumwnt;
      sumsigwntfrac += v.sigwntfrac; 
    }
  }

}//tpmatch
/*******************************************************************************************************************************/



/*******************************************************************************************************************************/
void PrimaryVertexAnalyzer4PU::wos_match(MVertexCollection& vtxs,
                               std::vector<SimEvent>& simEvt,
                               Tracks& tracks)
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
/*******************************************************************************************************************************/
{
  bool DEBUG = false;

  // reset
  if(DEBUG) {cout << "wos_match  vtxs.size() = " << vtxs.size()  << endl;}
  for (unsigned int iv = 0; iv < vtxs.size(); iv++) {
    vtxs[iv].sim = NOT_MATCHED_VTX_SIM;
    vtxs[iv].matchQuality = 0;
  }
  for (unsigned int iev = 0; iev < simEvt.size(); iev++){
    simEvt[iev].rec = NOT_MATCHED_VTX_REC;
    simEvt[iev].matchQuality = 0;
  }


  if(DEBUG){
    cout << "DEBUG: wos_match nwosmatch[0] = " << simEvt.at(0).nwosmatch << endl;
  }

  // select a rec vertex for simEvt if that rec vertex is dominated by the sim vertex
  // start with assigning sim vertices that dominate exactly one rec vertex (rank 1)
  // then repeat with sim vertices that dominate more recvertices (rank 2, 3, ...)
  // when two or more rec vertices are dominated by the same sim vertex,
  // assign the rec vertex that got more from that sim (higher wos)
  for (unsigned int rank = 1; rank < 8; rank++)
    {
      for (unsigned int iev = 0; iev < simEvt.size(); iev++)
	{
	  assert(simEvt.at(iev).wos_dominated_recv.size() == simEvt.at(iev).nwosmatch);
	  if (simEvt.at(iev).nwosmatch == 0) continue;     // doesn't dominate anything
	  if (simEvt.at(iev).nwosmatch > rank) continue;   // less ambiguous first

	  // only continue for simEvts that have not already been matched 
	  if (simEvt.at(iev).rec != NOT_MATCHED_VTX_REC) continue;
	  
	  // select a rec vertex (index iv)
	  unsigned int iv = NOT_MATCHED_VTX_REC;
	  for (unsigned int k = 0; k < simEvt.at(iev).wos_dominated_recv.size(); k++)
	    {
	      unsigned int rec = simEvt.at(iev).wos_dominated_recv.at(k); // candidate (index in vtxs)
	      if (vtxs(rec).sim != NOT_MATCHED_VTX_SIM) continue; // already matched
	      if (fabs(simEvt.at(iev).z - vtxs(rec).z()) > zWosMatchMax_) continue;// insanely far away
	      if ( (iv == NOT_MATCHED_VTX_REC) || (simEvt.at(iev).wos.at(rec) > simEvt.at(iev).wos.at(iv)) )
		{
		  iv = rec;
		}
	    }
	  // if we have found a viable candidate, make the link
	  if (iv != NOT_MATCHED_VTX_REC)
	    {
	      vtxs.at(iv).sim = iev;
	      simEvt.at(iev).rec = iv;
	      vtxs.at(iv).matchQuality = rank;
	      simEvt.at(iev).matchQuality = rank;
	      if(DEBUG) {cout << "wos_match : match made  [" << iev << "]  <-> (" << iv << ")   rank=" << rank << endl;}
	    }
	}
    }
      

  // by now we have exhausted the rec vertices that are dominated by an unmatched simvertex
  // have we?
  for(unsigned int iv = 0; iv < vtxs.size(); iv++) {
    if ((vtxs.at(iv).sim == NOT_MATCHED_VTX_SIM) && (vtxs.at(iv).maxwos > 0)){
      if (simEvt.at(vtxs.at(iv).wosmatch).rec == NOT_MATCHED_VTX_REC){
	auto ev = simEvt.at(vtxs.at(iv).wosmatch);
	cout << "wos_match :  unmatched [" << vtxs.at(iv).wosmatch << "]   dominantes unmatched  ("<< iv <<")   ????"
	     << "  rec.maxwos="  <<  vtxs.at(iv).maxwos 
	     << "  rec.maxwosnt="  <<  vtxs.at(iv).maxwosnt
	     << "  rec.sumwos=" << vtxs.at(iv).sumwos 
	     << "  sim.nwosmatch = " << ev.nwosmatch
	     << "  zrec=" << vtxs.at(iv).z()
	     << "  zsim=" << ev.z
	     <<endl;
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
	  if ((simEvt.at(sim).rec != NOT_MATCHED_VTX_REC) || (simEvt.at(sim).wos.size() == 0))
	    continue;
	
	  // ok, single simvertex sim, who is your your favorite rec vertex?
	  unsigned int rec = NOT_MATCHED_VTX_REC;
	  for (auto rv : simEvt.at(sim).wos) {
	    if ((rec == NOT_MATCHED_VTX_REC) || (rv.second > simEvt.at(sim).wos.at(rec))) {
	      rec = rv.first;
	    }
	  }

	  if (rec == NOT_MATCHED_VTX_REC){
	    cout << "wos_match : last hope failed for (" << rec << ")" << endl;
	    // second chance if wos didn't work?
	    for (auto rv : simEvt.at(sim).wnt) {
	      if ((rec == NOT_MATCHED_VTX_REC) || (rv.second > simEvt.at(sim).wnt.at(rec))) {
		rec = rv.first;
	      }
	    }
	  }
	
	if (rec == NOT_MATCHED_VTX_REC)
	  continue;  // should not happen
	if (vtxs.at(rec).sim != NOT_MATCHED_VTX_SIM)
	  continue;  // already gone
	
	// do you, recvertex rec, take this simvertex sim as your lawful wedded truthmatch?
	unsigned int rec2sim = NOT_MATCHED_VTX_SIM;
	for (auto sv : vtxs.at(rec).wos) {
	  if (simEvt.at(sv.first).rec != NOT_MATCHED_VTX_REC)
	    continue;  // already used
	  if ((rec2sim == NOT_MATCHED_VTX_SIM) || (sv.second > vtxs.at(rec).wos.at(rec2sim))) {
	    rec2sim = sv.first;
	  }
	}
	
	if (sim == rec2sim) {
	  // I do
	  vtxs.at(rec).sim = sim;
	  vtxs.at(rec).matchQuality = 8;
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
    for (unsigned int iv = 0; iv < vtxs.size(); iv++) {
      if (iv == 0) {
        dumpThisEvent_ = true;
        cout << "wos_match    recvtx = " << setw(4) << iv << "   (sim=" << vtxs[iv].sim << ")" << endl;
        cout << " splitfrom= " << vtxs[iv].split_from() << endl;
        cout << " maxwos = " << vtxs[iv].maxwos << endl;
        cout << " sumwos = " << vtxs[iv].sumwos << endl;

        for (auto w = vtxs[iv].wos.begin(); w != vtxs[iv].wos.end(); w++) {
          cout << "matching (wos)  simevent = " << setw(4) << w->first << "  zsim = " << setw(8) << fixed << setprecision(4)
               << simEvt.at(w->first).z << "  wos=" << setw(10) << w->second << endl;
        }

        for (unsigned int iev = 0; iev < simEvt.size(); iev++) {
          for (auto rv = simEvt.at(iev).wos.begin(); rv != simEvt.at(iev).wos.end(); rv++) {
            if (rv->first == iv) {
              cout << "matching (wos)  simevent " << iev << "    zsim= " << setw(8) << fixed << setprecision(4)
                   << simEvt.at(iev).z << "    zrec= " << setw(8) << fixed << setprecision(4)
                   << vtxs.at((*rv).first).z() << "  nt=" << (*rv).second << endl;
            }
          }
        }

	for (unsigned int iev = 0; iev < simEvt.size(); iev++) {
          for (auto rv = simEvt.at(iev).wnt.begin(); rv != simEvt.at(iev).wnt.end(); rv++) {
            if (rv->first == iv) {
              cout << "matching (wnt)  simevent " << iev << "    zsim= " << setw(8) << fixed << setprecision(4)
                   << simEvt.at(iev).z << "    zrec= " << setw(8) << fixed << setprecision(4)
                   << vtxs.at((*rv).first).z() << "  nt=" << (*rv).second << endl;
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
	  for(unsigned int iv=0; iv<vtxs.size(); iv++)
	    {
	      
	      if ((vtxs[iv].sim == NOT_MATCHED) && (vtxs[iv].split_from() < 0))
		{
		  double dz = simEvt.at(iev).z-vtxs.at(iv).z();
		  if( ( fabs(dz) < 0.1 )
		      ||( fabs(dz) < (2*vtxs.at(iv).zError()) )
		      )
		  {
		    vtxs[iv].matchQuality = 10;
		    vtxs[iv].sim = iev;
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
    if (iv != NOT_MATCHED_VTX_REC) {
      simEvt.at(iev).ndof = vtxs.at(iv).ndof();
      simEvt.at(iev).zrec = vtxs.at(iv).z();
    } else {
      simEvt.at(iev).ndof = 0.;
      simEvt.at(iev).zrec = 1000.;
    }
  }


  // consistency check
  for(unsigned int iev = 0; iev < simEvt.size(); iev++) {
    auto rec= simEvt.at(iev).rec;
    assert((rec == NOT_MATCHED_VTX_REC) || ((rec<vtxs.size()) && (vtxs.at(rec).sim == iev)));
  }
  for(unsigned int iv = 0; iv < vtxs.size(); iv++) {
    auto sim = vtxs.at(iv).sim;
    assert((sim == NOT_MATCHED_VTX_SIM) || ((sim < simEvt.size()) && (simEvt.at(sim).rec == iv)));
  }

} // wos_match
/******************************************************************************/








/***************************************************************************************/
void PrimaryVertexAnalyzer4PU::analyzeVertexComposition(std::map<std::string, TH1*>& h,
                                                           MVertex & v,
                                                           MVertexCollection & vtxs,
							   Tracks & tracks,
                                                           vector<SimEvent>& simEvt,
							   float npu)
// analyzes the composition of reco vertices in terms simEvts
/***************************************************************************************/
{
  // track counters per simEvt
  VertexCounter nt;     // all
  VertexCounter ntMTD;  // in MTD acceptance

  // weighted sums (aka purities)
  double sum_pt(0), sum_pt_majority(0), sum_pt_minority(0), sum_pt_unmatched(0);
  double sum_nt(0), sum_nt_majority(0), sum_nt_minority(0), sum_nt_unmatched(0);
  double sum_wt(0), sum_wt_majority(0), sum_wt_minority(0), sum_wt_unmatched(0);
  double sum_nt_timing(0), sum_nt_timing_majority(0), sum_nt_timing_minority(0), sum_nt_timing_unmatched(0);
  double sum_nt_reject(0), sum_nt_reject_minority(0), sum_nt_reject_majority(0), sum_nt_reject_unmatched(0);

  const double min_timing_quality = 0.8;  // FIXME make this configurable?

  // count tracks in the rec vertex
  for (auto tv : v.tracks){
    Float_t wt = v.trackWeight(tv);
    Float_t dz_trk_vtx = tv->z() - v.z();
    Float_t dz_trk_vtx_pull = (tv->z() - v.z())/tv->dzError();

    // alternative weight: pt (truncated at 10 GeV)
    double wpt = min(10., tv->pt());
    // 'track inside MTD acceptance' flag
    bool MTD = (tv->pt() > 0.9) && (fabs(tv->eta()) < 2.7);

    bool good_timing = tv->has_timing() && (tv->timeQuality() >= min_timing_quality);


    // (index of) the sim vertex this track has been matched to (or NOT_MATCHED_TK_SIM)
    unsigned int tk_sim = tracks.simevent_index_from_key(tv->key());
    assert(tv->_matched == tk_sim); // FIXME introduce  a function that returns the simEvent (or NOT_MATCHED_TK_SIM)
    bool correctly_assigned = (tk_sim == v.sim);
    Float_t wassign = correctly_assigned ? 1 : 0;

    Fill(h, "correctVtxFractionAllVsWeight", wt, wassign, v.is_signal()); // numerator matched and correct
    Fill(h, "vtxtrkallpurityvsputest", npu, wassign, v.is_signal());   // should be the same as vtxtrkallpurityvspu
    Fill(h, "vtxtrkpurityvsdz", dz_trk_vtx, wassign, v.is_signal());
    Fill(h, "vtxtrkpurityvsdzpull", dz_trk_vtx_pull, wassign, v.is_signal());
    Fill(h, "vtxtrkpurityvsdzerror", tv->dzError(), wassign, v.is_signal());

    if(tk_sim ==  NOT_MATCHED_TK_SIM){
	Fill(h, "unmatchedFractionVsWeight", wt, 1., v.is_signal());
	Fill(h, "trkWeightUnmatched", wt, v.is_signal());
    }else{
      Fill(h, "unmatchedFractionVsWeight", wt, 0., v.is_signal());
      if(tk_sim == v.sim){
	Fill(h, "trkWeightCorrectVtx", wt, v.is_signal());
	Fill(h, "correctVtxFractionVsWeight", wt, 1., v.is_signal()); // numerator matched and correct
      }else{ 
	Fill(h, "trkWeightIncorrectVtx", wt, v.is_signal());
	Fill(h, "correctVtxFractionVsWeight", wt, 0., v.is_signal()); // denominator only matched tracks
      }
    }

    if( wt >= min_trk_in_vtx_weight_){

      if(tk_sim != NOT_MATCHED_TK_SIM){
	nt.count(tk_sim);
	if (MTD) {
	  ntMTD.count(tk_sim);
	}
      }

      if (tk_sim == v.sim) { // majority actually stands for "track matched to the same simevent as the vertex"
	sum_pt_majority += wpt;  
	sum_nt_majority += 1.;
	sum_wt_majority += wt;
	if(good_timing) sum_nt_timing_majority +=1;
      } else if (tk_sim == NOT_MATCHED_TK_SIM) {
	sum_pt_unmatched += wpt;
	sum_nt_unmatched += 1.;
	sum_wt_unmatched += wt;
	if(good_timing) sum_nt_timing_unmatched +=1;
      } else {
	sum_pt_minority += wpt;
	sum_nt_minority += 1.;
	sum_wt_minority += wt;
	if(good_timing) sum_nt_timing_minority +=1;
      }

    }else{  // wt < min_vtxtrk_weight

      // for tracks rejected by the fitter weight cut
      // how often were they actually correctly assigned?
      if (tk_sim == v.sim) { 
	sum_nt_reject_majority += 1.;
      } else if (tk_sim == NOT_MATCHED_TK_SIM) {
	sum_nt_reject_unmatched += 1.;
      } else {
	sum_nt_reject_minority += 1.;
      }
    }
  }

  // fill histograms

  if (v.is_signal()) {
    Fill(h, "MTDTDR", sum_nt_minority);
  }

  Fill(h, "nsimevt", float(nt.nkey()), v.is_signal());
  Fill(h, "nsimevt_nt2", float(nt.nkey(2)), v.is_signal());
  Fill(h, "nsimevtMTD", float(ntMTD.nkey()), v.is_signal());

  sum_pt = sum_pt_majority + sum_pt_minority + sum_pt_unmatched;
  sum_nt = sum_nt_majority + sum_nt_minority + sum_nt_unmatched;
  sum_wt = sum_wt_majority + sum_wt_minority + sum_wt_unmatched;
  sum_nt_timing =  sum_nt_timing_majority + sum_nt_timing_minority + sum_nt_timing_unmatched;
  if (false && (sum_pt_majority == 0) && (sum_pt > 0)) {
    cout << "analyzeRecVertexComposition, vertex without majority tracks ?" << endl;
    cout << "z = " << v.z() << endl;
    cout << "nt = " << sum_nt_majority << "," << sum_nt_minority << "," << sum_nt_unmatched << endl;
    cout << "match = " << v.sim << endl;
    cout << "wosmatch = " << v.wosmatch << " maxwos =" << v.maxwos << " sumwos =" << v.sumwos << endl;
    cout << "wntmatch = " << v.wntmatch << " maxnt=" << v.maxwnt << " sumnt=" << v.sumwnt << endl;
    cout << "quality = " << v.matchQuality << endl;
    unsigned int cand = v.wosmatch;
    cout << "simEvt.at(cand).rec = " << simEvt.at(cand).rec << "    zrec=" << simEvt.at(cand).zrec << endl;
    cout << "simEvt.at(cand).nwosmatch = " << simEvt.at(cand).nwosmatch << endl;
    cout << "simEvt.at(cand).z = " << simEvt.at(cand).z << endl;
    cout << "dz = " << simEvt.at(cand).z - v.z() << " cutoff = " << zWosMatchMax_ << endl;
    cout << endl;
  }

  Fill(h, "pt_majority_frac", sum_pt_majority / sum_pt, v.is_signal());
  Fill(h, "pt_minority_frac", sum_pt_minority / sum_pt, v.is_signal());
  Fill(h, "pt_unmatched_frac", sum_pt_unmatched / sum_pt, v.is_signal());

  Fill(h, "pt_majority_frac_vsz", v.z(), sum_pt_majority / sum_pt, v.is_signal());
  Fill(h, "pt_minority_frac_vsz", v.z(), sum_pt_minority / sum_pt, v.is_signal());
  Fill(h, "pt_unmatched_frac_vsz", v.z(), sum_pt_unmatched / sum_pt, v.is_signal());

  Fill(h, "nt_majority_frac", sum_nt_majority / sum_nt, v.is_signal());
  Fill(h, "nt_minority_frac", sum_nt_minority / sum_nt, v.is_signal());
  Fill(h, "nt_unmatched_frac", sum_nt_unmatched / sum_nt, v.is_signal());

  Fill(h, "nt_majority_frac_vsz", v.z(), sum_nt_majority / sum_nt, v.is_signal());
  Fill(h, "nt_minority_frac_vsz", v.z(), sum_nt_minority / sum_nt, v.is_signal());
  Fill(h, "nt_unmatched_frac_vsz", v.z(), sum_nt_unmatched / sum_nt, v.is_signal());

  Fill(h, "wt_majority_frac", sum_wt_majority / sum_wt, v.is_signal());
  Fill(h, "wt_minority_frac", sum_wt_minority / sum_wt, v.is_signal());
  Fill(h, "wt_unmatched_frac", sum_wt_unmatched / sum_wt, v.is_signal());

  //  the denominator for purity should not include unmatched tracks because the numerator doesn't
  // note that "PU" is actually PU (with all purity 0) + fake
  double sum_nt_purity = sum_nt_majority + sum_nt_minority;
  if(sum_nt_purity > 0){
    float purity = sum_nt_majority / sum_nt_purity;
    Fill(h, "vtxtrkpurity", purity, v.is_signal());
    Fill(h, "vtxtrkpurityvsz", v.z(), purity, v.is_signal());
    Fill(h, "vtxtrkpurityvspu", npu, purity, v.is_signal());
  }

  sum_nt_reject = sum_nt_reject_majority + sum_nt_reject_minority + sum_nt_reject_unmatched;
  if((sum_nt + sum_nt_reject)> 0){// extended, i.e. including rejected (by weight) tracks
    float purity = (sum_nt_majority + sum_nt_reject_majority ) / (sum_nt + sum_nt_reject);
    Fill(h, "vtxtrkallpurity", purity, v.is_signal());
    Fill(h, "vtxtrkallpurityvsz", v.z(), purity, v.is_signal());
    Fill(h, "vtxtrkallpurityvspu", npu, purity, v.is_signal());
  }

  if(sum_nt_timing > 0){
    float timing_purity = sum_nt_timing_majority / sum_nt_timing;
    Fill(h, "vtxtrktimingpurity", timing_purity, v.is_signal());
    Fill(h, "vtxtrktimingpurityvsz", v.z(), timing_purity, v.is_signal());    // average purity vs z
    Fill(h, "vtxtrktimingpurityvspu", npu, timing_purity, v.is_signal());
  }

  // slightly different denominator for purity, not counting unmatched tracks (is that fairer?)
  float sum_nt_m = sum_nt_majority + sum_nt_minority;
  if(sum_nt_m > 0){
    float purity_m = sum_nt_majority / sum_nt_m;
    Fill(h, "vtxtrkpuritym", purity_m, v.is_signal());
    Fill(h, "vtxtrkpuritymvsz", v.z(), purity_m, v.is_signal());
    Fill(h, "vtxtrkpuritymvspu", npu, purity_m, v.is_signal());
  }

  // unbiased track assignment efficiency
  // (unbiased in contrast to the conditional assignment efficiency where the existence of the rec vertex is required)
  // here the denominator contains all tracks, even if the corresponding simEvt has no reconstructed counterpart
  // it is therefore always smaller than the biased track assignment efficiency (hmm, it's not)
  for(auto tk : tracks){
    if (tk.is_matched()){
      unsigned int ivsimrec = tk._simEvt->rec;   // index of the matched reco vertex or  NOT_MATCHED_VTX_REC
      unsigned int ivrec = tk.get_recv(vtxs.label()); // index or NO_RECVTX
      auto numtk = tk._simEvt->rtk.size(); // number of tracks in the track's sim event

      
      double w_any = 0., w_correct = 0;
      //if ((iv != NOT_MATCHED_VTX_REC ) && (select(vtxs.at(ivrec))) ){
      if ((ivrec != NO_RECVTX ) && (select(vtxs.at(ivrec))) ){
	w_any = 1.;                        // the track is attached to any recvertex
	if(ivrec == ivsimrec) w_correct = 1.;    // the track is attached to the recvertex that is matched to the simevent this track is part of
      }
      Fill(h, "utrkAssignmentEfficiency", w_any, tk._simEvt->is_signal());
      if (tk.is_primary()){
	Fill(h, "uprimtrkAssignmentEfficiency", w_any, tk._simEvt->is_signal());
	Fill(h, "uprimtrkCorrectAssignmentEfficiency", w_correct, tk._simEvt->is_signal());
      }
      Fill(h, "utrkAssignmentEfficiencyvspu", npu, w_any, tk._simEvt->is_signal());
      Fill(h, "utrkAssignmentEfficiencyvszerror", tk.dzError(), w_any, tk._simEvt->is_signal()); 
      
      Fill(h, "utrkCorrectAssignmentEfficiencyvspu", npu, w_correct, tk._simEvt->is_signal());
      Fill(h, "utrkCorrectAssignmentEfficiencyvsntrk", std::min(double(numtk), 100.5), w_correct, tk._simEvt->is_signal()); // makes no sense
      Fill(h, "utrkCorrectAssignmentEfficiencyvsdz", tk.zres(), w_correct, tk._simEvt->is_signal()); 
      Fill(h, "utrkCorrectAssignmentEfficiencyvsdzpull", tk.zpull(), w_correct, tk._simEvt->is_signal()); 
      Fill(h, "utrkCorrectAssignmentEfficiencyvszerror", tk.dzError(), w_correct, tk._simEvt->is_signal()); 
    }
  }
  
 
  
}//analyzeVertexComposition
/*********************************************************************************/


void PrimaryVertexAnalyzer4PU::analyzeVertexTrackAssociation(std::map<std::string, TH1*>& h, MVertexCollection& vtxs, Tracks& tracks, std::vector<SimEvent> const&, float const npu){
  //
  // Track-Vertex Association: Efficiency
  //
  // - denominator:
  //   - tracks matched to a TrackingParticle (TP-matched tracks)
  // - numerators:
  //   - TP-matched tracks assigned to a recoVertex
  //   - TP-matched tracks assigned to the recoVertex matched to the TrackingVertex of the matched-TrackingParticle
  //
  for(auto const& tk : tracks){
    // skip recoTracks not matched to a TrackingParticle
    if (not tk.is_matched()) continue;

    // the matched-TrackingParticle is assigned to the signal TrackingVertex (index == 0)
    auto const isSignalSimVtx = tk._simEvt->is_signal();

    // properties of the matched-TrackingParticle
    auto const tpPt = tk._tpr->pt();
    auto const tpEta = tk._tpr->eta();

    // denominator: TP-matched tracks
    Fill(h, "trkVtxAssocEffic_TPMatchedTracks_tpPt", tpPt, isSignalSimVtx);
    Fill(h, "trkVtxAssocEffic_TPMatchedTracks_tpEta", tpEta, isSignalSimVtx);
    if     (tpPt <  1.){ Fill(h, "trkVtxAssocEffic_TPMatchedTracks_tpEta_tpPt000to001", tpEta, isSignalSimVtx); }
    else if(tpPt <  3.){ Fill(h, "trkVtxAssocEffic_TPMatchedTracks_tpEta_tpPt001to003", tpEta, isSignalSimVtx); }
    else if(tpPt < 10.){ Fill(h, "trkVtxAssocEffic_TPMatchedTracks_tpEta_tpPt003to010", tpEta, isSignalSimVtx); }
    else               { Fill(h, "trkVtxAssocEffic_TPMatchedTracks_tpEta_tpPt010toINF", tpEta, isSignalSimVtx); }

    // find the recoVertex to which this recoTrack is assigned
    // (max trackWeight above 0.5)
    int bestRecoVtxIdx = -1;
    auto maxTrkWgt = 0.f;
    for (size_t vtxIdx=0; vtxIdx<vtxs.size(); ++vtxIdx){
      if(not select(vtxs[vtxIdx])) continue;

      if(not vtxs[vtxIdx].has_track_key(tk.key())) continue;

      auto const wt = vtxs[vtxIdx].trackWeight(tk);
      if(wt < min_trk_in_vtx_weight_) continue;

      if(wt > maxTrkWgt or bestRecoVtxIdx == -1){
        maxTrkWgt = wt;
        bestRecoVtxIdx = vtxIdx;
      }
    }

    // numerator #1: TP-matched recoTracks assigned to a(ny) recoVertex
    auto const belongsToOneRecoVtx = (bestRecoVtxIdx >= 0);
    if(not belongsToOneRecoVtx) continue;

    Fill(h, "trkVtxAssocEffic_TPMatchedTracksWithRecoVtx_tpPt", tpPt, isSignalSimVtx);
    Fill(h, "trkVtxAssocEffic_TPMatchedTracksWithRecoVtx_tpEta", tpEta, isSignalSimVtx);
    if     (tpPt <  1.){ Fill(h, "trkVtxAssocEffic_TPMatchedTracksWithRecoVtx_tpEta_tpPt000to001", tpEta, isSignalSimVtx); }
    else if(tpPt <  3.){ Fill(h, "trkVtxAssocEffic_TPMatchedTracksWithRecoVtx_tpEta_tpPt001to003", tpEta, isSignalSimVtx); }
    else if(tpPt < 10.){ Fill(h, "trkVtxAssocEffic_TPMatchedTracksWithRecoVtx_tpEta_tpPt003to010", tpEta, isSignalSimVtx); }
    else               { Fill(h, "trkVtxAssocEffic_TPMatchedTracksWithRecoVtx_tpEta_tpPt010toINF", tpEta, isSignalSimVtx); }

    // numerator #2: correct track-vertex association
    // - the recoVertex to which the recoTrack is assigned corresponds to
    //   the recoVertex matched to the TrackingVertex from which the matched-TrackingParticle originates

    // index of the recoVertex (if any) matched to the TrackingVertex from which the matched-TrackingParticle originates
    auto const idxOfRecoVtxOfSimVtx = tk._simEvt->rec;

    auto const belongsToCorrectRecoVtx = (uint(bestRecoVtxIdx) == idxOfRecoVtxOfSimVtx);
    if(not belongsToCorrectRecoVtx) continue;

    Fill(h, "trkVtxAssocEffic_TPMatchedTracksWithCorrectRecoVtx_tpPt", tpPt, isSignalSimVtx);
    Fill(h, "trkVtxAssocEffic_TPMatchedTracksWithCorrectRecoVtx_tpEta", tpEta, isSignalSimVtx);
    if     (tpPt <  1.){ Fill(h, "trkVtxAssocEffic_TPMatchedTracksWithCorrectRecoVtx_tpEta_tpPt000to001", tpEta, isSignalSimVtx); }
    else if(tpPt <  3.){ Fill(h, "trkVtxAssocEffic_TPMatchedTracksWithCorrectRecoVtx_tpEta_tpPt001to003", tpEta, isSignalSimVtx); }
    else if(tpPt < 10.){ Fill(h, "trkVtxAssocEffic_TPMatchedTracksWithCorrectRecoVtx_tpEta_tpPt003to010", tpEta, isSignalSimVtx); }
    else               { Fill(h, "trkVtxAssocEffic_TPMatchedTracksWithCorrectRecoVtx_tpEta_tpPt010toINF", tpEta, isSignalSimVtx); }
  }

  //
  // Track-Vertex Association: Purity
  //
  // - denominator:
  //   - recoTracks of a recoVertex
  // - numerator:
  //   - recoTracks of a recoVertex that are TP-matched,
  //     and have the correct track-vertex association
  //     (SimVertex of TP is matched to recoVertex of recoTrack)
  //
  for (auto const& v : vtxs){
    // skip recoVertexs not passing selection criteria
    if (not select(v)) continue;

    auto const isSignalVtx = v.is_signal();

    uint ntk = 0;
    uint ntk_woFakeVtxs = 0;
    uint ntk_matched = 0;
    for (auto tv : v.tracks){
      // restrict to recoTracks truly assigned to this recoVertex (trackWeight >= min)
      auto const wt = v.trackWeight(tv);
      if (wt < min_trk_in_vtx_weight_) continue;

      // denominator: recoTracks of a recoVertex
      ++ntk;

      auto const trkPt = tv->pt();
      auto const trkEta = tv->eta();

      Fill(h, "trkVtxAssocPurity_TracksOfRecoVtx_trkPt", trkPt, isSignalVtx);
      Fill(h, "trkVtxAssocPurity_TracksOfRecoVtx_trkEta", trkEta, isSignalVtx);
      if     (trkPt <  1.){ Fill(h, "trkVtxAssocPurity_TracksOfRecoVtx_trkEta_trkPt000to001", trkEta, isSignalVtx); }
      else if(trkPt <  3.){ Fill(h, "trkVtxAssocPurity_TracksOfRecoVtx_trkEta_trkPt001to003", trkEta, isSignalVtx); }
      else if(trkPt < 10.){ Fill(h, "trkVtxAssocPurity_TracksOfRecoVtx_trkEta_trkPt003to010", trkEta, isSignalVtx); }
      else                { Fill(h, "trkVtxAssocPurity_TracksOfRecoVtx_trkEta_trkPt010toINF", trkEta, isSignalVtx); }

      // restrict to recoVertexs that are matched to a SimVertex (i.e. skip fake recoVertexs)
      auto const recoVtxHasSimMatch = (v.sim != NOT_MATCHED_VTX_SIM);
      if(not recoVtxHasSimMatch) continue;

      ++ntk_woFakeVtxs;

      Fill(h, "trkVtxAssocPurity_TracksOfRecoVtxMatchedToSimVtx_trkPt", trkPt, isSignalVtx);
      Fill(h, "trkVtxAssocPurity_TracksOfRecoVtxMatchedToSimVtx_trkEta", trkEta, isSignalVtx);
      if     (trkPt <  1.){ Fill(h, "trkVtxAssocPurity_TracksOfRecoVtxMatchedToSimVtx_trkEta_trkPt000to001", trkEta, isSignalVtx); }
      else if(trkPt <  3.){ Fill(h, "trkVtxAssocPurity_TracksOfRecoVtxMatchedToSimVtx_trkEta_trkPt001to003", trkEta, isSignalVtx); }
      else if(trkPt < 10.){ Fill(h, "trkVtxAssocPurity_TracksOfRecoVtxMatchedToSimVtx_trkEta_trkPt003to010", trkEta, isSignalVtx); }
      else                { Fill(h, "trkVtxAssocPurity_TracksOfRecoVtxMatchedToSimVtx_trkEta_trkPt010toINF", trkEta, isSignalVtx); }

      // index of the SimVertex from which the TrackingParticle matched to this recoTrack originates
      // (for TP-unmatched recoTracks, this equals NOT_MATCHED_TK_SIM)
      auto const tk_sim = tracks.simevent_index_from_key(tv->key());

      // recoTrack is TP-matched and its recoVertex is matched to the SimVertex of the recoTrack's TP
      auto const hasCorrectTrkVtxAsso = (tk_sim != NOT_MATCHED_TK_SIM and tk_sim == v.sim);
      if(not hasCorrectTrkVtxAsso) continue;

      ++ntk_matched;

      Fill(h, "trkVtxAssocPurity_TracksOfRecoVtxWithCorrectSimVtxMatch_trkPt", trkPt, isSignalVtx);
      Fill(h, "trkVtxAssocPurity_TracksOfRecoVtxWithCorrectSimVtxMatch_trkEta", trkEta, isSignalVtx);
      if     (trkPt <  1.){ Fill(h, "trkVtxAssocPurity_TracksOfRecoVtxWithCorrectSimVtxMatch_trkEta_trkPt000to001", trkEta, isSignalVtx); }
      else if(trkPt <  3.){ Fill(h, "trkVtxAssocPurity_TracksOfRecoVtxWithCorrectSimVtxMatch_trkEta_trkPt001to003", trkEta, isSignalVtx); }
      else if(trkPt < 10.){ Fill(h, "trkVtxAssocPurity_TracksOfRecoVtxWithCorrectSimVtxMatch_trkEta_trkPt003to010", trkEta, isSignalVtx); }
      else                { Fill(h, "trkVtxAssocPurity_TracksOfRecoVtxWithCorrectSimVtxMatch_trkEta_trkPt010toINF", trkEta, isSignalVtx); }
    }

    if(ntk > 0){
      auto const purity = ntk_matched / float(ntk);
      Fill(h, "trkVtxAssocPurity_vs_vz", v.z(), purity, isSignalVtx);
      Fill(h, "trkVtxAssocPurity_vs_pu", npu, purity, isSignalVtx);
    }

    if(ntk_woFakeVtxs > 0){
      auto const purity = ntk_matched / float(ntk_woFakeVtxs);
      Fill(h, "trkVtxAssocPurityWithoutFakeRecoVtxs_vs_vz", v.z(), purity, isSignalVtx);
      Fill(h, "trkVtxAssocPurityWithoutFakeRecoVtxs_vs_pu", npu, purity, isSignalVtx);
    }
  }
}


/***************************************************************************************/
void PrimaryVertexAnalyzer4PU::analyzeVertexMergeRateTP(std::map<std::string, TH1*>& h,
                                                        MVertexCollection& vertexes,
                                                        Tracks & tracks,
							vector<SimEvent>& simEvt,
                                                        const string message)
/***************************************************************************************/
// definitions ?
//  * a sim vertex that was not found but would have been reconstructed if there had not been other vertices (nearby)?
//  * a rec vertex that is matched to >= 2 sim vertices (dqm) but what is the merge rate then? <> merge-fraction ?
//
{
  auto nv = vertexes.size();
  
  for (unsigned int iv = 0; iv < nv; iv++) {
    
    // only consider selected vertices
    if (!select(vertexes(iv)))
      continue;


    unsigned int nsim = 0;      // count sim vertices that make their largest contribution in this recvtx
    unsigned int nsim_any = 0;  // count sim vertices that make a contribution to this vertex
    unsigned int nsim_any_not_matched = 0;  // count sim vertices that make a contribution but are not matched to this vertex

    //cout << setw(4) << iv << "  -> [" << recvmatch[iv].wos.size() << "]  ";
    for (auto s : vertexes(iv).wos) {
      auto sim = s.first;
      nsim_any++;
      // only count if that sim vertex was unmatched or matched to this rec vertex
      if ((simEvt[sim].rec == NOT_MATCHED_VTX_REC) || (simEvt[sim].rec == iv)) {
        nsim_any_not_matched++;
        // count if this rec vertex received the largest contribution of that simvtx
        if (simEvt[sim].max_nt_vtx() == iv) {
          nsim++;
        }
      }
    }

    //    cout << " ==>  nsim=" << nsim << " nsim_any=" << nsim_any << " nsim_any_not_matched=" << nsim_any_not_matched << endl;
    
    double z = vertexes(iv).z();
    bool signal = vertexes(iv).is_signal();
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
  if ((simEvt.size()>0) && (simEvt[0].is_signal()) && (simEvt[0].is_matched()))
    {
      double zsimsignal = simEvt[0].z;
      double zrecsignal = simEvt[0].zrec;
      double recsignal =  simEvt[0].rec;

      Fill(h,"countzrecsignal",1); // denominator
      
      for (unsigned int iv = 0; iv < vertexes.size(); iv++) {
	if (!select(vertexes(iv)) || (iv == recsignal) )
	  continue;
	double dzrecsim = vertexes(iv).z() - zsimsignal;
	double dzrecrec = vertexes(iv).z() - zrecsignal;
	bool fake = (vertexes(iv).sim == NOT_MATCHED_VTX_SIM);
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
	      reportVertex(vertexes(iv).recovertex(), Form("fake vertex on top of signal at z=%f in %s", vertexes(iv).z(), message.c_str()), dump_fake_vertex_on_top_of_signal_);
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
	  Fill(h, "effvsdzsimPU_PU", dz, simEvt.at(sim2).is_matched() ? 1. : 0.);
	  Fill(h, "effvsdzsimPU_PU_hr", dz, simEvt.at(sim2).is_matched() ? 1. : 0.);
	  if(simEvt.at(sim1).is_matched())
	    {
	      Fill(h, "effvsdzsimmatchedPU_PU", dz, simEvt.at(sim2).is_matched() ? 1. : 0.);
	      Fill(h, "effvsdzsimmatchedPU_PU_hr", dz, simEvt.at(sim2).is_matched() ? 1. : 0.);
	      
	      bool sel = simEvt.at(sim2).is_matched() && select(vertexes(simEvt.at(sim2).rec)) ;
	      Fill(h, "effselvsdzsimmatchedPU_PU", dz, sel ? 1. : 0.);
	      Fill(h, "effselvsdzsimmatchedPU_PU_hr", dz, sel ? 1. : 0.);
	    }
	  else
	    {
	      Fill(h, "effvsdzsimunmatchedPU_PU", dz, simEvt.at(sim2).is_matched() ? 1. : 0.);
	      Fill(h, "effvsdzsimunmatchedPU_PU_hr", dz, simEvt.at(sim2).is_matched() ? 1. : 0.);

	      bool sel = simEvt.at(sim2).is_matched() && select(vertexes(simEvt.at(sim2).rec)) ;
	      Fill(h, "effselvsdzsimunmatchedPU_PU", dz, sel ? 1. : 0.);
	      Fill(h, "effselvsdzsimunmatchedPU_PU_hr", dz, sel ? 1. : 0.);
	    }
	}
      // efficiency vs distance to the nearest other PU sim vertex (excluding signal vertices)
      Fill(h, "effallvsdzsimminPU_PU", dzmin, simEvt.at(sim1).is_matched() ? 1. :0);
      bool sel1 = simEvt.at(sim1).is_matched() && select(vertexes(simEvt.at(sim1).rec)) ;
      Fill(h, "effselvsdzsimminPU_PU", dzmin, sel1 ? 1. :0);
    }
  
  // efficiency of pu vtx finding vs distance to the signal vertex
  if ((simEvt.size()>0) && (simEvt.at(0).is_matched()))
    {
      assert( simEvt.at(0).is_signal());
      for (unsigned int sim = 1; sim < simEvt.size(); sim++)
	{
	  double dz = simEvt.at(sim).z - simEvt.at(0).z;
	  Fill(h, "effvsdzsimsignal_PU", dz, simEvt.at(sim).is_matched() ? 1. : 0.);
	  Fill(h, "effvsdzsimsignal_PU_hr", dz, simEvt.at(sim).is_matched() ? 1. : 0.);
	  
	  bool sel = simEvt.at(sim).is_matched() && select(vertexes(simEvt.at(sim).rec)) ;
	  Fill(h, "effselvsdzsimsignal_PU", dz, sel ? 1. : 0.);
	  Fill(h, "effselvsdzsimsignal_PU_hr", dz, sel ? 1. : 0.);
	}
    }
  
}//analyzeVertexMergeRateTP
/***************************************************************************************/









/***************************************************************************************/
void PrimaryVertexAnalyzer4PU::analyzeVertexCollectionTP(std::map<std::string, TH1*>& h,
                                                         MVertexCollection& vtxs,
                                                         Tracks& tracks,
                                                         vector<SimEvent>& simEvt,
                                                         const string message)
// with track truthmatching (tp)  
/***************************************************************************************/
{
  if (verbose_) {
    cout << "analyzeVertexCollectionTP, simEvts= " << simEvt.size() << endl;
  }

  if(!tracking_truth_available_) return;

  // count pile-up
  int npu0 = simEvt.size();
  Fill(h, "npu0", npu0);
  float npu = simEvt.size();

  if (eventSummaryCounter_++ < nEventSummary_) {
    printEventSummary_tp(h, vtxs, tracks, simEvt, message);
  }

  // simtrack-like histograms, for comparison
  if(vtxs.size() > 0){
    Fill(h, "sigwosfractagged", vtxs(0).sigwosfrac);
    Fill(h, "sigwntfractagged", vtxs(0).sigwntfrac);
  }

  double sigwosfracmax = 0., sigwntfracmax = 0.;
  for(const auto & v : vtxs){
    if(v.sigwosfrac > sigwosfracmax){ sigwosfracmax = v.sigwosfrac;}
    if(v.sigwntfrac > sigwntfracmax){ sigwntfracmax = v.sigwntfrac;}
  }
  Fill(h, "sigwosfracmax", sigwosfracmax);
  Fill(h, "sigwntfracmax", sigwntfracmax);

  for(auto & v : vtxs){
    if(select(v)){
      if(v.sigwosfrac > 0) Fill(h, "sigwosfracvsdzsignal", v.z() - simEvt.at(0).z, v.sigwosfrac );
      if(v.sigwntfrac > 0) Fill(h, "sigwntfracvsdzsignal", v.z() - simEvt.at(0).z, v.sigwntfrac );
    }
  }

  // how many rec vertices have a wos value above some threshold
  const unsigned int nbin = 20;// number of bins in nrecwithsigwos_tp
  vector<unsigned int> nvc(nbin, 0.);
  for(const auto & v : vtxs){
    if (select(v) && (v.sigwosfrac > 0)){
      // count and cumulate
      unsigned int bin = v.sigwosfrac * nbin;
      if (bin > nbin){ bin = nbin;} // if sigwosfrac is == 1.0
      for(unsigned int b = 0; b < bin; b++){ nvc[b]++; }
    }
  }

  // then fill
  for(unsigned int b=0; b<nbin; b++){
    Fill(h, "nrecwithsigwos", float(b + 0.5)/nbin, float(nvc.at(b)));
  }

  unsigned int ntpfake = 0;
  unsigned int ntpreal = 0;
  unsigned int ntpsplit = 0, ntpotherfake = 0;
  unsigned int ntpfakesel = 0, ntprealsel = 0, ntpsplitsel = 0, ntpotherfakesel = 0, ntprecsel = 0, ntprec = 0;
  unsigned int ntpsplitselfromsignal = 0, ntpsplitselfrompu = 0;
  unsigned int ntpsplitfromsignal = 0, ntpsplitfrompu = 0;

  auto nv = vtxs.size();
  for (unsigned int iv = 0; iv < nv; iv++) {

    if (vtxs.at(iv).isRecoFake()) continue;
    ntprec++;

    if (select(vtxs.at(iv))) {
      ntprecsel++;

      // (signed) distance to the nearest selected other vertex
      // used for signed quantities in fillVertexHistos
      double deltaz = 1000.;
      for(unsigned int jv = 0; jv < vtxs.size(); jv++)
	{
	  if ( (jv != iv) && (fabs(vtxs.at(jv).z() - vtxs.at(iv).z()) < fabs(deltaz)) && select(vtxs.at(jv)) )
	    {
	      deltaz = vtxs.at(jv).z() - vtxs.at(iv).z();
	    }
	}
      


      if (vtxs[iv].is_real()) {
	// real selected
        ntprealsel++;
	fillVertexHistosMatched(h, "matchedvtxsel",  vtxs.at(iv), tracks, simEvt, deltaz);
    
	if(vtxs[iv].is_signal()){
	  // signal
	  for( auto tv : vtxs(iv).tracks){
	    if (tv->dzError() < 0.01) {
	      fillTrackHistos(h, "matchedsignalvtxseldriver", *tv, vtxs(iv).recvtx);
	    }
	  }
	  if (iv == 0){
	    reportEvent(Form("signal vertex selected at index [%d]", iv), false);
	  }else{
	    std::string comment = "?";
	    if (vtxs[0].is_real()){ comment = "real vertex at [0]";}
	    else if(vtxs[0].other_fake()){ comment = "other fake vertex at [0]";}
	    else if(vtxs[0].split_from() > 0){ comment = "split from " + std::to_string(vtxs[0].split_from());}
	    else if(vtxs[0].is_fake()){ comment = "fake vertex at [0]";}
	    else if(vtxs[0].is_signal()){ comment = "signal vertex at [0]?????";}
	    reportEvent(Form("signal vertex selected at index [%d]    %s", iv,comment.c_str()), false);
	  }
	}else{
	  // PU
	  if(iv == 0){
	    reportEvent(Form("PU vertex selected at index [0]"), false);
	  }
	}
	
      }else{
	
	// fake selected
        ntpfakesel++;

	fillVertexHistos(h, "fakevtxsel",  vtxs.at(iv), tracks, deltaz);
      
	if (vtxs[iv].split_from() >= 0) {
	  ntpsplitsel++;
	  fillVertexHistos(h, "splitvtxsel",  vtxs.at(iv), tracks, deltaz);
	}

	if (vtxs[iv].split_from() == 0) {
	  ntpsplitselfromsignal++;
	  fillVertexHistos(h, "splitvtxselfromsignal",  vtxs.at(iv), tracks, deltaz);
	  reportEvent(Form("split signal vertex selected at index [%d]", iv), false);

	  for (auto tv : vtxs(iv).tracks){
	    if (tv->dzError() < 0.0100) {
	      fillTrackHistos(h, "splitvtxselfromsignaldriver", *tv, vtxs(iv).recvtx);
	    }
	  }
	
	}

	if (vtxs[iv].split_from() > 0) {
	  ntpsplitselfrompu++;
	  fillVertexHistos(h, "splitvtxselfrompu",  vtxs(iv), tracks, deltaz);
	}

	if (vtxs[iv].other_fake()) {
	  ntpotherfakesel++;
	  fillVertexHistos(h, "otherfakevtxsel",  vtxs(iv), tracks, deltaz);
	}
      } // fake selected

      Fill(h, "ntpfakeselratevssimPU", npu, vtxs[iv].is_fake() ? 1. : 0.);

    } // selected vertex iv


    if (vtxs(iv).is_real()) {
      // matched, but not necessarily selected, note that matched and selected are filled above
      ntpreal++;
      fillVertexHistosMatched(h, "matchedvtx", vtxs.at(iv), tracks, simEvt);
    }

    if (vtxs(iv).is_fake()) {
      // fake (selected or not)
      ntpfake++;
      fillVertexHistos(h, "fakevtx",  vtxs.at(iv), tracks);
      if (vtxs(iv).ndof() > 50) {
        reportVertex(vtxs(iv).recovertex(), Form("big fake vertex ndof=%f , split_from=%d", 
						  vtxs(iv).ndof(), vtxs(iv).split_from()), dump_big_fakes_);
      }

      if(vtxs(iv).split_from() >=0){
	ntpsplit ++;
	if (vtxs[iv].split_from() == 0) {
	  ntpsplitfromsignal++;
	}else{
	  ntpsplitfrompu++;
	}
      }else{
	ntpotherfake ++;
      }


    }

  }  // end of recvtx loop

  if(! ((ntpfake == (ntpsplit + ntpotherfake)) && (ntpfakesel == (ntpsplitsel + ntpotherfakesel)) && ((ntpreal + ntpfake) == ntprec) && ((ntprealsel + ntpfakesel) == ntprecsel) )){
    std::cout <<" hein  ntpreal= "<< ntpreal << " ntpfake=" << ntpfake  << " ntprec="  << ntprec << "   size()=" << nv << "  nvtx=" << vtxs.nvtx() 
	      << "  ntpsplit= " << ntpsplit 
	      << "  ntpsplitsel= " << ntpsplitsel 
	      << "  ntpsplitselfromsignal= " << ntpsplitselfromsignal 
	      << "  ntpsplitselfrompu= " << ntpsplitselfrompu 
	      << "  ntpotherfake= " << ntpotherfake
	      << "  ntpotherfakesel= " << ntpotherfakesel
	      << "  npu= " << npu
	      << "  simpu_="  << simPU_
	      << std::endl;
    }

  // simpu_ = number of interaction from puInfo, some of which may not produce a vertex?
  // npu = simEvt.size()
  Fill(h, "ntpsplit", float(ntpsplit));
  Fill(h, "ntpsplitfromsignal", float(ntpsplitfromsignal));
  Fill(h, "ntpsplitfrompu", float(ntpsplitfrompu));
  Fill(h, "ntpsplitsel", float(ntpsplitsel));
  Fill(h, "ntpsplitselfromsignal", float(ntpsplitselfromsignal));
  Fill(h, "ntpsplitselfrompu", float(ntpsplitselfrompu));
  Fill(h, "ntpfake", float(ntpfake));
  Fill(h, "ntpfakesel", float(ntpfakesel));
  Fill(h, "ntpotherfake", float(ntpotherfake));
  Fill(h, "ntpotherfakesel", float(ntpotherfakesel));
  Fill(h, "ntpreal", float(ntpreal));
  Fill(h, "ntprealsel", float(ntprealsel));

  Fill(h, "ntprecselvssimPU", simPU_, float(ntprecsel));
  Fill(h, "ntpfakeselvssimPU", simPU_, float(ntpfakesel));
  Fill(h, "ntprealselvssimPU", simPU_, float(ntprealsel));
  Fill(h, "ntpsplitselvssimPU", simPU_, float(ntpsplitsel));
  Fill(h, "ntpotherfakeselvssimPU", simPU_, float(ntpotherfakesel));

  Fill(h, "ntprecvssimPU", simPU_, float(ntprec));
  Fill(h, "ntpfakevssimPU", simPU_, float(ntpfake));
  Fill(h, "ntprealvssimPU", simPU_, float(ntpreal));
  Fill(h, "ntpsplitvssimPU", simPU_, float(ntpsplit));
  Fill(h, "ntpotherfakevssimPU", simPU_, float(ntpotherfake));

  for(unsigned int iv = 0; iv < vtxs.size(); iv++){
    if( vtxs[iv].is_signal() ){
      Fill(h, "index_signal", float(iv));
    }else if( vtxs[iv].is_real() ){
      Fill(h, "index_pu", float(iv));
    }else if( vtxs[iv].split_from() == 0) {
      Fill(h, "index_splitfromsignal", float(iv));
    }else if( vtxs[iv].split_from() > 0) {
      Fill(h, "index_splitfrompu", float(iv));
    }else{
      Fill(h, "index_otherfake", float(iv));
    }
  }

  // ranking by track sumpt2:
  std::vector< std::pair<double, unsigned int> > ranking;
  for (unsigned int iv = 0; iv < vtxs.size(); iv++) { // why did it start at 1 ?
    if (select(vtxs.at(iv))){
      ranking.push_back( std::make_pair(-vtxs(iv).sumpt2(),iv) );
    }
  }

  stable_sort(ranking.begin(), ranking.end());
  for(unsigned int n = 0; n < ranking.size(); n++){
    unsigned int iv = ranking[n].second;
    if( vtxs[iv].is_signal() ){
      Fill(h, "trksumpt2rank_signal", float(n));
    }else if( vtxs[iv].is_real() ){
      Fill(h, "trksumpt2rank_pu", float(n));
    }else if( vtxs[iv].split_from() == 0) {
      Fill(h, "trksumpt2rank_splitfromsignal", float(n));
    }else if( vtxs[iv].split_from() > 0) {
      Fill(h, "trksumpt2rank_splitfrompu", float(n));
    }else{
      Fill(h, "trksumpt2rank_otherfake", float(n));
    }
  }

  // for the classification of unmatched rec- and sim-vertices
  // rec first: how many sim vertices do I need to get 50% of my wos
  for (unsigned int iv = 0; iv < vtxs.size(); iv++)
    {
      if (select(vtxs.at(iv))) {
	
	  // sort the wos values
	  std::vector<double> wos;
	  double sumwos = 0.;
	  for (auto it : vtxs[iv].wos) {
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
	  // for n==0: cumwos==vtxs[iv].sumwomaxwos, sumwos==vtxs[iv].sumwos
	  
	  n = n > 4 ? 4 : n;// fill overflow in the last bin
	  double v = log(float(sumwos))/log(10.);
	  if (vtxs[iv].is_fake()){
	    Fill(h, "wornk_unmatchedsel", v, float(n));
	    if(vtxs[iv].split_from() >= 0){
	      Fill(h, "wornk_splitsel", v, float(n));
	    }else{
	      Fill(h, "wornk_othersel", v, float(n));
	    }
	  } else {
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
      if (ev->is_matched()){
	Fill(h, "wornk_matchedsim", v, float(n));
      }else{
	Fill(h, "wornk_unmatchedsim", v, float(n));
      }
    }

  // determine the recvertex that is the signalvertex (if there is one)
  unsigned int signalv = 0;
  bool has_signal = false;
  for (unsigned int iv = 0; iv < vtxs.size(); iv++) {  // was size-1, why?
    if (select(vtxs.at(iv)) && vtxs[iv].is_real() && vtxs[iv].is_signal()) {
      signalv = iv;
      has_signal = true;
      break;
    }
  }

  // properties of vertices that are split from the signal vertex
  if (has_signal) {
    for (unsigned int jv = 0; jv < vtxs.size(); jv++) {
      if (vtxs[jv].split_from() == 0) {
        Fill(h, "ndofsignalsplit", vtxs.at(jv).ndof());

	if (select(vtxs.at(signalv)) && select(vtxs.at(jv))){
          Fill(h, "zdiffrecsignalsel", fabs(vtxs.at(signalv).z() - vtxs.at(jv).z()));
	}

        if (select(vtxs.at(jv))) {
          Fill(h, "zdiffrecselsignalsplit", fabs(vtxs.at(signalv).z() - vtxs.at(jv).z()));
          Fill(h, "zdiffrecselsignalsplit-dzbin", fabs(vtxs.at(signalv).z() - vtxs.at(jv).z()));
          // add more, e.g. sumpt, fraction of pt
        }
      }
    }
  }

  // distance histogram for pairs of selected reco vertices
  for (unsigned int iv = 0; iv < vtxs.size(); iv++) {
    if (select(vtxs.at(iv))) {
      for (unsigned int jv = 0; jv < vtxs.size(); jv++) {
        if ((!(jv == iv)) && select(vtxs.at(jv))) {
          double dz = vtxs.at(iv).z() - vtxs.at(jv).z();
          if (vtxs[jv].is_real() && vtxs[iv].is_real()) {
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
    for (unsigned int jv = 0; jv < vtxs.size(); jv++) {

      if ((!(jv == signalv)) && select(vtxs.at(jv)) && (!vtxs[jv].is_signal())) {
        double dz = vtxs.at(signalv).z() - vtxs.at(jv).z();
        if (vtxs[jv].is_real()) {
          Fill(h, "zdiffrecselsignalrealpu", std::abs(dz));
	  double dzsim = simEvt.at(vtxs[jv].sim).z - simEvt.at(0).z;
          Fill(h, "zdiffsimselsignalrealpu", std::abs(dzsim));
          if(jv == 1) Fill(h, "zdiffrecselsignalrealpuidx1", std::abs(dzsim));
        } else {
          Fill(h, "zdiffrecselsignalfake", std::abs(dz));
	  if (vtxs[jv].split_from() == 0){
	    Fill(h, "zdiffrecselsignalsplitfromsignal", std::abs(dz));
	  }else if(vtxs[jv].split_from() > 0){
	    Fill(h, "zdiffrecselsignalsplitfrompu", std::abs(dz));
	  }else{
	    Fill(h, "zdiffrecselsignalotherfake", std::abs(dz));
	  }
        }
      }

      // more or less the same, but older
      if ((!(jv == signalv)) && select(vtxs.at(jv)) && (!vtxs[jv].is_signal())) {
        double dz = vtxs.at(signalv).z() - vtxs.at(jv).z();
        if (vtxs[jv].is_real()) {
          Fill(h, "zdiffrecselsignalreal", std::abs(dz));
          Fill(h, "zdiffrecselsignalreal-dzbin", std::abs(dz));
        } else {
          Fill(h, "zdiffrecselsignalfakev2", std::abs(dz));
          Fill(h, "zdiffrecselsignalfake-dzbin", std::abs(dz));
        }
      }
    }
  }

  // histogram the neighborhood of a PU vertex
  for (unsigned int iv = 0; iv < vtxs.size(); iv++) {
    if ( select(vtxs.at(iv)) && vtxs[iv].is_real() && (!vtxs[iv].is_signal())) {
      for (unsigned int jv = 0; jv < vtxs.size(); jv++) {
        if ((!(jv == iv)) && select(vtxs.at(jv)) && (!vtxs[jv].is_signal())) {
          double dz = vtxs.at(iv).z() - vtxs.at(jv).z();
          if (vtxs[jv].is_real()) {
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
  for (unsigned int iv = 0; iv < vtxs.size(); iv++) {
    if (vtxs.at(iv).ndof() > 4.) {
      double mindistance_realreal = 1e10;

      for (unsigned int jv = iv; jv < vtxs.size(); jv++) {
        if ((!(jv == iv)) && select(vtxs.at(jv))) {
          double dz = vtxs.at(iv).z() - vtxs.at(jv).z();
          if (vtxs[iv].is_real() && vtxs[jv].is_real()) {
            Fill(h, "zdiffrecselfound", dz);
            if (fabs(dz) < fabs(mindistance_realreal)) {
              mindistance_realreal = dz;
            }
          } else if (vtxs[iv].is_fake() && vtxs[jv].is_fake()) {
            Fill(h, "zdiffrecselfakefake", dz);
          }
        }
      }

      double mindistance_fakereal = 1e10;
      double mindistance_realfake = 1e10;
      double zf = 0, zr = 0;
      for (unsigned int jv = 0; jv < vtxs.size(); jv++) {
        if ((!(jv == iv)) && select(vtxs.at(jv))) {
          double dz = vtxs.at(iv).z() - vtxs.at(jv).z();

          if (vtxs[iv].is_fake() && vtxs[jv].is_real()) {
            Fill(h, "zdiffrecselfakereal", dz);
            if (fabs(dz) < fabs(mindistance_fakereal)) {
              mindistance_fakereal = dz;
              zr = vtxs[iv].is_real() ? vtxs.at(iv).z() : vtxs.at(jv).z();
              zf = vtxs[iv].is_fake() ? vtxs.at(iv).z() : vtxs.at(jv).z();
            }
          }

          if (vtxs[iv].is_real() && vtxs[jv].is_fake() && (fabs(dz) < fabs(mindistance_realfake))) {
            mindistance_realfake = dz;
          }
        }
      }

      if (vtxs[iv].is_real()) {
        Fill(h, "zdiffmin4realreal", fabs(mindistance_realreal));
        Fill(h, "zdiffmin4realfake", fabs(mindistance_realfake));
        if (vtxs[iv].is_signal()) {
          Fill(h, "zdiffmin4signalreal", fabs(mindistance_realreal));
          Fill(h, "zdiffmin4signalfake", fabs(mindistance_realfake));
        }
      } else if (vtxs[iv].is_fake()) {
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
  if ( (nv > 1) && select(vtxs.at(0)) && select(vtxs.at(1)) ){
    double zdiff = vtxs.at(1).z() - vtxs.at(0).z();
    Fill(h, "zdiffrec10-selsel", zdiff);
    if (vtxs[0].is_signal()){
      Fill(h, "zdiffrec10-signalsel", zdiff);
      if (vtxs[1].is_real()){ Fill(h, "zdiffrec10-signalreal", zdiff);} // redundant and obsolete
      if (vtxs[1].is_real()){ Fill(h, "zdiffrec10-signalrealpu", zdiff);}
      if (vtxs[1].is_fake()){ Fill(h, "zdiffrec10-signalfake", zdiff);}
      if (vtxs[1].split_from() == 0){ 
	Fill(h, "zdiffrec10-signalsplitfromsignal", zdiff);
      }else if (vtxs[1].split_from() > 0){ 
	Fill(h, "zdiffrec10-signalsplitfrompu", zdiff);
      }else{
	Fill(h, "zdiffrec10-signalotherfake", zdiff);
      }

      if (vtxs[1].is_real() &&  (fabs(zdiff) < 0.05)) reportEvent("real pu close to signal promoted to nr 2", false);
    }
  }

  // efficiency histograms for simEvts
  for (unsigned int iev = 0; iev < simEvt.size(); iev++) {
    bool is_signal = (iev == 0);

    //simEvt matching done in tpmatch
    float wall = 0., wsel=0., wsel3sigma=0.;
    //if (simEvt.at(iev).rec != NOT_MATCHED_VTX_REC) {
    if (simEvt.at(iev).is_matched()){

      wall = 1.;

      // selected 
      if( select(vtxs.at(simEvt.at(iev).rec)) ){
	wsel = 1.;
	if (fabs(vtxs.at(simEvt.at(iev).rec).z() - simEvt.at(iev).z) < 3.* vtxs.at(simEvt.at(iev).rec).zError()){
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
    Fill(h, "effallvssimpu", simPU_, wall, is_signal);
    Fill(h, "effselvssimpu", simPU_, wsel, is_signal);
    Fill(h, "effsel3sigmavssimpu", simPU_, wsel3sigma, is_signal);
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
      
      cout << " signal vertex not found :" << message << "   run : event = " << run_ << " : " << event_ <<  " dumpflag " << dump_signal_vertex_not_tpmatched_ <<endl;
    } else if (verbose_ && (wall == 0) && (simEvt.at(iev).rtk.size() > 40)) {
       reportEvent(
		   Form("big PU vertex not found (TP) (sim=%5d   z=%8.4f    charged trks=%d)",
			iev, simEvt.at(iev).z, (int)simEvt.at(iev).rtk.size()),
		   false);
       //cout << " big PU vertex not found :" << message << "  ";
       //cout << "sim " << iev << "  z=" << simEvt.at(iev).z << " ntk=" << simEvt.at(iev).rtk.size()
       //    << "  nChTP=" << simEvt.at(iev).rtk.size() << endl;
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
    if (ev->Tc > 0) {
      Fill(h, "logTc", log(ev->Tc) / log(10.), ev == simEvt.begin());
    }

    // and the same for pairs of vertices
    for (vector<SimEvent>::iterator ev2 = ev + 1; ev2 != simEvt.end(); ev2++) {
      vector<MTrack> xt;
      if ((ev->rtkprimsel.size() > 0) && (ev2->rtkprimsel.size() > 0) &&
          (ev->rtkprimsel.size() + ev2->rtkprimsel.size()) > 1) {
        xt.insert(xt.end(), ev->rtkprimsel.begin(), ev->rtkprimsel.end());
        xt.insert(xt.end(), ev2->rtkprimsel.begin(), ev2->rtkprimsel.end());
        double xTc, xChsq, xDzmax, xDztrim;
        getTc(xt, xTc, xChsq, xDzmax, xDztrim);
        if (xTc > 0) {
          Fill(h, "xTc", xTc, ev == simEvt.begin());
          Fill(h, "logxTc", log(xTc) / log(10), ev == simEvt.begin());
          Fill(h, "xChisq", xChsq, ev == simEvt.begin());
          if (xChsq > 0) {
            Fill(h, "logxChisq", log(xChsq), ev == simEvt.begin());
          };
          Fill(h, "xdzmax", xDzmax, ev == simEvt.begin());
          Fill(h, "xdztrim", xDztrim, ev == simEvt.begin());
        }
      }
    }
  }

  // vertex pairs with >=4 tracks
  for (vector<SimEvent>::iterator ev1 = simEvt.begin(); (ev1+1) != simEvt.end(); ev1++) {

    if ((!ev1->is_matched()) || (ev1->rtk.size() < 4))   continue;

    for (vector<SimEvent>::iterator ev2 = ev1 + 1; ev2 != simEvt.end(); ev2++) {
      if ((!ev2->is_matched()) || (ev2->rtk.size() < 4))  continue;

      double deltazsim = ev2->z - ev1->z;
      Fill(h, "zdiffsimallTP", deltazsim);

      // both sim vertices of this pair were found
      Fill(h, "zdiffsimfoundTP", deltazsim);
      
      assert((ev1->rec < vtxs.size()) && (ev2->rec < vtxs.size()));
      if (select(vtxs.at(ev2->rec)) && select(vtxs.at(ev1->rec))) {
	Fill(h, "zdiffsimfoundselTP", deltazsim);
      }
      
      if (select(vtxs.at(ev2->rec)) && select(vtxs.at(ev1->rec)) && (ev1 == simEvt.begin())) {
	Fill(h, "zdiffsimfoundselSignalTP", std::abs(deltazsim));
      }
      double deltazrec = ev2->zrec - ev1->zrec;
      Fill(h, "zdiffrecvssimTP", deltazsim, deltazrec);
      Fill(h, "zdiffrecvsdsimTP", deltazsim, deltazrec - deltazsim);
      Fill(h, "zdiffrecvsdsimTPprof", deltazsim, deltazrec - deltazsim);
      
      /* nothing about mergers, yet */
    }
  }

  timer_start("TP_track_vertex_loop");
  // track properties, fake and real
  for (auto v : vtxs){
    for (auto tv : v.tracks){
      if (v.is_real()) {
        if (tv->dzError() < 0.01) {
          fillTrackHistos(h, "realvtxdriver", *tv, v.recvtx);
        }
      } else {
        if (tv->dzError() < 0.01) {
          fillTrackHistos(h, "fakevtxdriver", *tv, v.recvtx);
        }
      }
    }
  }
  timer_stop("TP_track_vertex_loop");

  for (auto v : vtxs){
    if (select(v)) {
      analyzeVertexComposition(h, v, vtxs, tracks, simEvt, npu);
    }
  }

  analyzeVertexTrackAssociation(h, vtxs, tracks, simEvt, npu);
  
  }//analyzeVertexCollectionTP
/***************************************************************************************/





/***************************************************************************************/
void PrimaryVertexAnalyzer4PU::analyzeVertexCollectionPtvis(std::map<std::string, TH1*>& h,
                                                         MVertexCollection & vtxs,
                                                         Tracks& tracks,
                                                         vector<SimEvent>& simEvt,
                                                         const string message)
/***************************************************************************************/
{
  unsigned int k = 0;
  for(; k < vtxs.size(); k++){
    if(vtxs(k).sim == 0) break;
  }
  if(k == vtxs.size()){
    report_counted("analyzeVertexCollectionPtvis : Signal vertex not found", 1);
    return;
  }
  
  auto v = vtxs(k);
  double zvtx = v.z();
  double density = simPU_ * exp(-0.5*pow((zvtx - vertexBeamSpot_.z0()) / vertexBeamSpot_.sigmaZ(),2))/sqrt(2*3.14159)/vertexBeamSpot_.sigmaZ(); 

  vector<double> ptvis = {0,0,0,0,0};
  vector<double> px =  {0,0,0,0,0};
  vector<double> py =  {0,0,0,0,0};
  

  for (unsigned int i = 0; i < tracks.size(); i++)
    {
    
      MTrack tk = tracks(i);
      double tkpx  = tk.pt() * cos(tk.phi());
      double tkpy  = tk.pt() * sin(tk.phi());
      
      if (fabs(tk.z() - zvtx) < 0.1){
	ptvis[0] += tk.pt();
	px[0] += tkpx;
	py[0] += tkpy;
      }
      if (fabs(tk.z() - zvtx) < 0.2){
	ptvis[1] += tk.pt();
	px[1] += tkpx;
	py[1] += tkpy;
      }
      double zpull = (tk.z() - zvtx) / tk.dzError();
      if (fabs(zpull) < 3.){
	ptvis[2] += tk.pt();
	px[2] += tkpx;
	py[2] += tkpy;
      }

      double w3 = 1;
      double Z4 = 0;
      for(unsigned int iv=0; iv < vtxs.size(); iv++)
	{
	  if(select(vtxs(iv))){
	    if (fabs(vtxs(iv).z() -tk.z()) <  fabs(zvtx - tk.z())){
	      w3 = 0;
	    }
	      
	    double p = (vtxs(iv).z() - tk.z()) / tk.dzError();
	    if (fabs(p) < 5.)
	      {
		Z4 +=  exp(-0.5* p * p);
	      }
	  }
	}

      if (w3 > 0){// nearest
	ptvis[3] += tk.pt();
	px[3] += tkpx;
	py[3] += tkpy;
      }
      
      if ((Z4 > 0) && (zpull < 5)){
	ptvis[4] += exp(-0.5*zpull*zpull) / Z4 * tk.pt();
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
/******************analyzeVertexCollectionPtvis*********************************************************************/


/*********************************************************************************************/
void PrimaryVertexAnalyzer4PU::analyzeVertexCollectionZmatching(std::map<std::string, TH1*>& h,
								MVertexCollection& vtxs,
								std::vector<SimEvent>& simEvts,
								const std::string message,
								const double zwindow_sigmas
								)
{
  if (verbose_) {
    cout << "analyzeVertexCollectionZmatching, simEvts= " << simEvts.size() << " vtxs= " << vtxs.size()<<  " window size = " <<zwindow_sigmas  <<  endl;
  }

  unsigned int nsim = simEvts.size();
  unsigned int nvtx = vtxs.size();
  
  std::vector<unsigned int> nsimmatch(nvtx);  // count sim vertices within <zwindow_sigmas> sigma of the reconstructed position
  std::vector<unsigned int> nrecmatch(nsim);  // count rec vertices with reconstructed position within  <zwindow_sigmas> sigma 
  // repeat with the additional requirement of 2 common truth matched tracks
  std::vector<unsigned int> nsimmatch_tp(nvtx);
  std::vector<unsigned int> nrecmatch_tp(nsim);
  // and once more with a modified z-matching
  std::vector<unsigned int> nsimmatch_c(nvtx);
  std::vector<unsigned int> nrecmatch_c(nsim);


  for(unsigned int k=0; k < nvtx; k++){

    const auto v = vtxs.at(k);
    if (v.isRecoFake()) continue;

    unsigned int nearest_sim = NOT_MATCHED_VTX_REC;

    for(unsigned int i=0; i < nsim; i++){

      if (std::abs(simEvts.at(i).z - v.z()) < (zwindow_sigmas * v.zError())){
	//std::cout << "zmatching " << message << "  [" << i << "]" << simEvts.at(i).z << "   " << v.z() << "+/-" << v.zError() << " (" << k << ")" << endl;
	nrecmatch[i]++;
	nsimmatch[k]++;
	//nmatchall++;
	if(nearest_sim == NOT_MATCHED_VTX_REC){
	  nearest_sim = i;
	}else{
	  if(std::abs(simEvts.at(i).z - v.z()) < std::abs(v.z()-simEvts.at(nearest_sim).z)){
	    nearest_sim = i;
	  }
	}
      }

      // same but requiring two matched tracks
      if ((std::abs(simEvts.at(i).z - v.z()) < (zwindow_sigmas * v.zError())) && (simEvts.at(i).countVertexTracks(v, 0.2) > 1)){
	nrecmatch_tp[i]++;
	nsimmatch_tp[k]++;
      }

      // modified z-matching
      double dzmax = std::min(0.1, std::max(0.0100, zwindow_sigmas * v.zError()));
      if (std::abs(simEvts.at(i).z - v.z()) < dzmax){
	nrecmatch_c[i]++;
	nsimmatch_c[k]++;
      }
      
    }// sim events
  }// rec vtxs



  // fill histograms
  // 
  for(unsigned int i=0; i < nsim; i++){
    Fill(h, "zmatcheffvspu", nsim, nrecmatch[i] > 0 ? 1. : 0.);
    Fill(h, "zmatchambigvspu", nsim, nrecmatch[i] > 1 ? 1. : 0.);
    Fill(h, "zmatchnrecmatch", nrecmatch[i]);

    Fill(h, "zcmatcheffvspu", nsim, nrecmatch_c[i] > 0 ? 1. : 0.);
    Fill(h, "zcmatchambigvspu", nsim, nrecmatch_c[i] > 1 ? 1. : 0.);
    Fill(h, "zcmatchnrecmatch", nrecmatch_c[i]);

    Fill(h, "ztpmatcheffvspu", nsim, nrecmatch_tp[i] > 0 ? 1. : 0.);
    Fill(h, "ztpmatchambigvspu", nsim, nrecmatch_tp[i] > 1 ? 1. : 0.);
    Fill(h, "ztpmatchnrecmatch", nrecmatch_tp[i]);
  }

  // fake here means : no sim vertex within <zwindow_sigmas> x sigma_z
  for(unsigned int k=0; k < nvtx; k++){
    Fill(h, "zmatchfakevspu", nsim, nsimmatch[k] == 0 ? 1. : 0);
    Fill(h, "zcmatchfakevspu", nsim, nsimmatch_c[k] == 0 ? 1. : 0);
    Fill(h, "ztpmatchfakevspu", nsim, nsimmatch_tp[k] == 0 ? 1. : 0);
  }

  for(unsigned int k=0; k < nvtx; k++){
    Fill(h, "zmatchnsimmatch", nsimmatch[k]);
    Fill(h, "zcmatchnsimmatch", nsimmatch_c[k]);
    Fill(h, "ztpmatchnsimmatch", nsimmatch_tp[k]);
  }

  
  auto rand = TRandom();
  // the "random efficiency"
  for(unsigned int i = 0; i < 1000; i++){
    double ztoy = rand.Gaus(0, sigmaZ_); // beam-spot
    unsigned  int nmatch = 0;
    for(unsigned int k=0; k < nvtx; k++){
      const auto v = vtxs.at(k);
      if (v.isRecoFake()) continue;
      if (std::abs(ztoy - v.z()) < (zwindow_sigmas * v.zError())) nmatch++;
    }
    Fill(h, "zrandomeffvspu", nsim, nmatch > 0 ? 1.: 0.);
  }
  // the "random fake rate"
  for(unsigned int k = 0; k < nvtx; k++){
    const auto v = vtxs.at(k);
    if (v.isRecoFake()) continue;
    double ztoy = rand.Gaus(0, sigmaZ_);
    unsigned  int nmatch = 0;
    for(unsigned int i=0; i < nsim; i++){
      if ((std::abs(simEvts.at(i).z - ztoy) < (zwindow_sigmas * v.zError()))) nmatch++;
    }
    Fill(h, "zrandomfakevspu", nsim, nmatch == 0 ? 1.: 0.);
  }



    
  // go one step beyond looking at z-windows :  maximally greedy matching
  std::vector<std::pair<unsigned int, unsigned int>> recsim;    // based on z-distance only
  std::vector<std::pair<unsigned int, unsigned int>> recsim_c;    // based on z-distance, truncated
  std::vector<std::pair<unsigned int, unsigned int>> recsim_tp;  // additionally require at least one truth matched track
  // and the same for selected vertices
  std::vector<std::pair<unsigned int, unsigned int>> recselsim;    // based on z-distance only
  std::vector<std::pair<unsigned int, unsigned int>> recselsim_c;    // based on z-distance, truncated
  std::vector<std::pair<unsigned int, unsigned int>> recselsim_tp;  // additionally require at least one truth matched track

  unsigned int nvtxrec = 0, nvtxsel = 0;
  for(unsigned int k = 0; k < nvtx; k++){
    const auto v = vtxs.at(k);
    if( v.isRecoFake() ) continue;
    
    nvtxrec ++;
    if(select(v)) nvtxsel++;


    for(unsigned int i = 0; i < nsim; i++){

      // z-distance alone
      if (std::abs(simEvts.at(i).z - v.z()) < (zwindow_sigmas * v.zError())){
	recsim.emplace_back(k,i);
	if(select(v)) recselsim.emplace_back(k,i);
      }

      // truncated z-distance, allow at least 100 um, do not allow more than 1 mm
      double dzmax = std::min(0.1, std::max(0.0100, zwindow_sigmas * v.zError()));
      if (std::abs(simEvts.at(i).z - v.z()) < dzmax){
	recsim_c.emplace_back(k,i);
	if(select(v)) recselsim_c.emplace_back(k,i);
      }

      // z-distance + at least two tracks with weight > 0.5
      if ((std::abs(simEvts.at(i).z - v.z()) < (zwindow_sigmas * v.zError())) && (simEvts.at(i).countVertexTracks(v, 0.5) > 1)){
	recsim_tp.emplace_back(k,i);
	if(select(v)) recselsim_tp.emplace_back(k,i);
      }
    }
  }

  if( nsim > 0){
    int max_match = 0, max_match_c = 0,max_match_tp = 0;
    FFA(nvtxrec, nsim, recsim, max_match);
    Fill(h, "FFAzmatcheffvspu", nsim, float(max_match) / nsim);

    FFA(nvtxrec, nsim, recsim_c, max_match_c);
    Fill(h, "FFAzcmatcheffvspu", nsim, float(max_match_c) / nsim);

    FFA(nvtxrec, nsim, recsim_tp, max_match_tp);
    Fill(h, "FFAztpmatcheffvspu", nsim, float(max_match_tp) / nsim);

    if(nvtxrec > 0){
      Fill(h, "FFAzmatchfakevspu", nsim, float(nvtxrec - max_match) / nvtxrec);
      Fill(h, "FFAzcmatchfakevspu", nsim, float(nvtxrec - max_match_c) / nvtxrec);
      Fill(h, "FFAztpmatchfakevspu", nsim, float(nvtxrec - max_match_tp) / nvtxrec);
      //should the denominator contain only vertices with tp matchted tracks?
    }

    // include selection
    max_match = 0; 
    max_match_c = 0;
    max_match_tp = 0;
    FFA(nvtxsel, nsim, recselsim, max_match);
    Fill(h, "FFAzmatchseleffvspu", nsim, float(max_match) / nsim);

    FFA(nvtxsel, nsim, recselsim_c, max_match_c);
    Fill(h, "FFAzcmatchseleffvspu", nsim, float(max_match_c) / nsim);

    FFA(nvtxsel, nsim, recselsim_tp, max_match_tp);
    Fill(h, "FFAztpmatchseleffvspu", nsim, float(max_match_tp) / nsim);

    if(nvtxsel > 0){
      Fill(h, "FFAzmatchselfakevspu", nsim, float(nvtxsel - max_match) / nvtxsel);
      Fill(h, "FFAzcmatchselfakevspu", nsim, float(nvtxsel - max_match_c) / nvtxsel);
      Fill(h, "FFAztpmatchselfakevspu", nsim, float(nvtxsel - max_match_tp) / nvtxsel);
    }

  }
  
}
/*****************************analyzeVertexCollectionZmatching**********************************/


/***************************************************************************************
*  for samples with rec tracks, but without sim tracks (for PU)                        *
***************************************************************************************/
void PrimaryVertexAnalyzer4PU::analyzeVertexCollectionSimPvNoSimTracks(std::map<std::string, TH1*>& h,
                                                                       MVertexCollection& vtxs,
                                                                       Tracks& tracks,
                                                                       std::vector<SimEvent>& simEvts,
                                                                       const std::string message) {
  if(verbose_){
    cout << "PrimaryVertexAnalyzer4PU::analyzeVertexCollectionSimPvNoSimTracks   " << message;
    cout << "       size = " << simEvts.size() << endl;
  }
  if (simEvts.size()==0)return;

  auto signal = simEvts[0];

  /* all reco vs signal (sim)*/
  for(const auto & v : vtxs){
    double dzrecsim =  v.z() - signal.z;
    Fill(h, "dzrecsimSignal", dzrecsim);
    Fill(h, "dzrecsimSignalHR", dzrecsim);
    Fill(h, "dzrecsimSignalVHR", dzrecsim);
  }

  /* tagged vs signal */
  double dztag = 1000.;
  bool signal_tagged = false;
  if(vtxs.size() > 0){
    dztag = vtxs[0].z() - signal.z;
    signal_tagged = fabs(dztag) < 0.2;
    if(fabs(dztag) > 0.2){
      reportVertex(vtxs[0].recovertex(), Form("tagged vertex at z=%f8.4 far from signal at z=%f8.4 ",vtxs[0].z(), signal.z), true);
    }

    Fill(h, "dzrecsimtaggedsignalwide", dztag);
    Fill(h, "dzrecsimtaggedsignalnarrow", dztag);
    Fill(h, "dzrecsimtaggedsignalnarrowovfl", fabs(dztag)<0.1 ? fabs(dztag) : 0.0999);
    Fill(h, "tag01mmvspu", float(simPU_), fabs(dztag) < 0.01 ? 1. : 0.);
    Fill(h, "tag02mmvspu", float(simPU_), fabs(dztag) < 0.02 ? 1. : 0.);
    Fill(h, "tag05mmvspu", float(simPU_), fabs(dztag) < 0.05 ? 1. : 0.);
    Fill(h, "tag1mmvspu", float(simPU_), fabs(dztag) < 0.1 ? 1. : 0.);
    Fill(h, "tag2mmvspu", float(simPU_), fabs(dztag) < 0.2 ? 1. : 0.);
    Fill(h, "tag5mmvspu", float(simPU_), fabs(dztag) < 0.5 ? 1. : 0.);
  }


  /* neighbourhood of the signal vertex  */
  unsigned int nrec5mm=0, nrec2mm=0, nrec1mm=0, nrec05mm=0, nrec02mm=0, nrec01mm=0;
  for( auto & v : vtxs){
    if (select(v)){
	double dz = v.z() - signal.z;
	if (fabs(dz) < 0.5) nrec5mm++;
	if (fabs(dz) < 0.2) nrec2mm++;
	if (fabs(dz) < 0.1) nrec1mm++;
	if (fabs(dz) < 0.05) nrec05mm++;
	if (fabs(dz) < 0.02) nrec02mm++;
	if (fabs(dz) < 0.01) nrec01mm++;
      }
  }   
 
  unsigned int nsim5mm=0, nsim2mm=0, nsim1mm=0, nsim05mm=0, nsim02mm=0, nsim01mm=0;
  for(unsigned int i = 1; i < simEvts.size(); i++){
    if (simEvts[i].is_visible()){
      double dz = simEvts[i].z - signal.z;
	if (fabs(dz) < 0.5) nsim5mm++;
	if (fabs(dz) < 0.2) nsim2mm++;
	if (fabs(dz) < 0.1) nsim1mm++;
	if (fabs(dz) < 0.05) nsim05mm++;
	if (fabs(dz) < 0.02) nsim02mm++;
	if (fabs(dz) < 0.01) nsim01mm++;
    }      
  }
  
  string tag = "_tagged";
  if(!signal_tagged){
    tag="_nottagged";
  }
  Fill(h, "nrecwithin5mm"+tag, double(nrec5mm));
  Fill(h, "nsimwithin5mm"+tag, double(nsim5mm));
  Fill(h, "nrecwithin1mm"+tag, double(nrec1mm));
  Fill(h, "nsimwithin1mm"+tag, double(nsim1mm));
  Fill(h, "nrecvsnsimwithin5mm"+tag, double(nsim5mm), double(nrec5mm));
  Fill(h, "nrecvsnsimwithin2mm"+tag, double(nsim2mm), double(nrec2mm));
  Fill(h, "nrecvsnsimwithin1mm"+tag, double(nsim1mm), double(nrec1mm));
  Fill(h, "nrecvsnsimwithin05mm"+tag, double(nsim05mm), double(nrec05mm));
  Fill(h, "nrecvsnsimwithin02mm"+tag, double(nsim02mm), double(nrec02mm));
  Fill(h, "nrecvsnsimwithin01mm"+tag, double(nsim01mm), double(nrec01mm));


  // index (=ranking) vs z-distance-to-signal ranking
  vector<pair<double, unsigned int> >dzsig;
  for( auto & v : vtxs){
    if(select(v)){
      dzsig.push_back(make_pair(fabs(v.z()-signal.z), v.index()));
    }
  }

  stable_sort(dzsig.begin(), dzsig.end(), lt);
  
  for(unsigned int i = 0; i < dzsig.size(); i++){
    float ptrank = dzsig[i].second < 10 ? float(dzsig[i].second) : 9.9999;
    float dzrank = dzsig[i].first < 10 ? float(dzsig[i].first) : 9.9999;
    //cout << "dz-sorted " << i << " ) " << "  dz=" << dzsig[i].first << "  index=" << dzsig[i].second  << "  ptrank " << int(ptrank) << endl;
    
    Fill(h,"indexvsdistancerank", dzrank, ptrank);
    if(signal_tagged){
      Fill(h,"indexvsdistancerank_tagged", dzrank, ptrank);
    }else{
      Fill(h,"indexvsdistancerank_nottagged", dzrank, ptrank);
    }
  }


  for( auto & v : vtxs){

    if (!select(v, 0))  continue;
    double ptmax2 = v.ptmax2();

    // for this reconstructed primary find the distance to the nearest sim vertex
    double dzmin = 1e100;
    for (unsigned int idx = 1; idx < simEvts.size(); idx++) {
      double dz = simEvts[idx].z - v.z();
      if (std::abs(dz) < std::abs(dzmin)) {
        dzmin = dz;
      }

      // also fill zrec-zsim histos for all pairs
      Fill(h, "dzrecsim", dz);
      Fill(h, "dzrecsimHR", dz);
      if (ptmax2 < 0.4) {
        Fill(h, "dzrecsimptmax2lt04", dz);
      } else {
        Fill(h, "dzrecsimptmax2gt04", dz);
      }
    }

    // closest sim vertex to a given rec-vertex
    Fill(h, "dzrecsimmin", dzmin);
    if (ptmax2 < 0.4) {
      Fill(h, "dzrecsimminptmax2lt04", dzmin);
    } else {
      Fill(h, "dzrecsimminptmax2gt04", dzmin);
    }

    if (std::abs(dzmin) < 0.1) {
      Fill(h, "zmatched01", v.z());
      Fill(h, "ptmax2matched01", ptmax2);
    }
    if (std::abs(dzmin) < 0.2) {
      Fill(h, "zmatched02", v.z());
      Fill(h, "ptmax2matched02", ptmax2);
    }

    // fake?, no real vertex within 1 (or 2 ) mm
    if (std::abs(dzmin) > 0.1) {
      Fill(h, "zunmatched01", v.z());
      Fill(h, "ptmax2unmatched01", ptmax2);
      if (ptmax2 < 0.4) {
        Fill(h, "zunmatched01ptmax2lt04", v.z());
      }
    }
    if (std::abs(dzmin) > 0.2) {
      Fill(h, "zunmatched02", v.z());
      Fill(h, "ptmax2unmatched02", ptmax2);
      if (ptmax2 < 0.4) {
        Fill(h, "zunmatched02ptmax2lt04", v.z());
      }
    }
  }

  // fill track property histograms for selected rec tracks that are far away from any simulated vertex
  for(auto tk : tracks){
    if (tk.selected()) {
      double z0 = tk.z();
      double dzmin = 1e10;
      for (unsigned int idx = 1; idx < simEvts.size(); idx++) {
        double dz = simEvts[idx].z - z0;
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
void PrimaryVertexAnalyzer4PU::analyzeVertexRecoCPUTime(std::map<std::string, TH1*>& h,
                                                        const reco::VertexCollection* recVtxs,
                                                        const std::string message)
/***************************************************************************************/
{
  double tclu = 0;
  double tfit = 0;
  double ttime = 0;
  int nsel = 0;
  bool found_timing_info = false;


  for (reco::VertexCollection::const_iterator v = recVtxs->begin(); v != recVtxs->end(); ++v) {
    if (v->isFake()) {
      tclu = v->covariance(iX, iX);
      tfit = v->covariance(iY, iY);
      ttime = v->covariance(iZ, iZ);
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
    Fill(h, "cputime/ttimevsSimPU", simPU_, ttime);
    Fill(h, "cputime/tcluvsLPU", lumiPU_, tclu);
    Fill(h, "cputime/tfitvsLPU", lumiPU_, tfit);
    Fill(h, "cputime/ttimevsLPU", lumiPU_, ttime);
    Fill(h, "cputime/tcluvsnsel", nsel, tclu);
    Fill(h, "cputime/tfitvsnsel", nsel, tfit);
    Fill(h, "cputime/ttimevsnsel", nsel, ttime);
  }
}
/***************************************************************************************/





//******* the following code does not require MC and will/should work for data, requires tracks **********
/***************************************************************************************/
void PrimaryVertexAnalyzer4PU::analyzeVertexCollectionReco(std::map<std::string, TH1*>& h,
                                                           MVertexCollection& vtxs,
                                                           Tracks& tracks,
                                                           const std::string message)
/***************************************************************************************/
{
  cout << "PrimaryVertexAnalyzer4PU::analyzeVertexCollectionReco (MVertexCollection)" << message << endl;
  int nrectrks = tracks.size();

  // repeat some counting here from analyzeVertexCollectionNotrack, needed for eventclassification and correlation
  int nrecvtx = 0;
  int nselvtx = 0;
  int nvtxselgt1sigmaz = 0;

  for (unsigned int iv = 0; iv < vtxs.size(); iv++) {
    MVertex & v = vtxs[iv];
    if (!(v.isRecoFake()) && (v.ndof() > 0)) {
      nrecvtx++;
    }
    
    fillVertexHistos(h, "recvtx", v, tracks);

    if(iv == 0){fillVertexHistos(h, "taggedvtx", v, tracks);}

    if (select(v)) {
      nselvtx++;
      fillVertexHistos(h, "selectedvtx", v, tracks);
      if (std::abs(v.z() - vertexBeamSpot_.z0() - dzb_) > sigmaZ_) {
        nvtxselgt1sigmaz++;
      }
    }
  }

  Fill(h, "lPUbybx", (float)bunchCrossing_, lumiPU_);
  Fill(h, "nselvtxbybx", (float)bunchCrossing_, float(nselvtx));
  Fill(h, "nselvtxbyLS", (float)luminosityBlock_, float(nselvtx));
  Fill(h, "nselvtx", float(nselvtx));
  Fill(h, "LPU", lumiPU_);
  Fill(h, "lBX", 1.e3 * instBXLumi_);
  Fill(h, "nlPU", lumiPU_);
  Fill(h, "lPULS", lsglb_, lumiPU_);
  Fill(h, "sigmazLS", lsglb_, sigmaZ_);
  Fill(h, "zbeamLS", lsglb_, vertexBeamSpot_.z0());
  Fill(h, "LPUvsLPU", lumiPU_, lumiPU_);
  Fill(h, "LPUvsavgLPU", avglumiPU_, lumiPU_);
  Fill(h, "avgLPUvsLPU", lumiPU_, avglumiPU_);
  Fill(h, "nselvtxvsLPU", lumiPU_, float(nselvtx));
  Fill(h, "nselvtxvssimPU", simPU_, float(nselvtx));
  if(MC_){
    Fill(h, "simPU", float(simPU_));
  }  

  // distance to first vertex and ndof for different positions in the vertex list
  if( (vtxs.size() > 1) && (select(vtxs.at(0))) ){

    double z0 = vtxs.at(0).z();
    Fill(h, "ndof-vtx0", vtxs.at(0).ndof());
    Fill(h, "logsumpt2-vtx0", log(vtxs.at(0).sumpt2()));
    Fill(h, "logndof-vtx0", log(vtxs.at(0).ndof())/log(10.));

    for (unsigned int iv = 1; iv < vtxs.size(); iv++) {
      const MVertex & v = vtxs.at(iv);
      if(select(v)){
	if(iv < 10){
	    Fill(h, Form("zdiffrec-vtx%d", iv), v.z() - z0);
	    Fill(h, Form("ndof-vtx%d", iv), v.ndof());
	    Fill(h, Form("logndof-vtx%d", iv), log(v.ndof())/log(10.));
	    Fill(h, Form("logsumpt2-vtx%d", iv), log(v.sumpt2()));
	}else{
	  Fill(h, "zdiffrec-vtxgt9", v.z() - z0);
	  Fill(h, "ndof-vtxgt9", v.ndof());
	  Fill(h, "logndof-vtxgt9", log(v.ndof())/log(10.));
	  Fill(h, "logsumpt2-vtxgt9", log(v.sumpt2()));
	}
      }
    }
  }

  // FIXME  move the following  to a separate analysis method ? e.g. analyzeTracksWithVertices ?

  bool is_hiPU = lumiPU_ > 50;
  //bool is_tail =  (lumiPU_ > 40) && (nselvtx > (0.8*lumiPU_ + 3*sqrt(0.8*lumiPU_)));
  bool is_tail = (lumiPU_ > 40) && (nvtxselgt1sigmaz > (0.24 * lumiPU_ + 3 * sqrt(0.24 * lumiPU_)));

  // -----------------  reconstructed tracks  ------------------------
  // the list of selected tracks can only be correct if the selection parameterset  and track collection
  // is the same that was used for the reconstruction

  auto const* pv0_ptr = vtxs.size() > 0 and vtxs[0].recvtx->isValid() ? vtxs[0].recvtx : NULL;

  for (unsigned int i = 0; i < tracks.size(); i++) {
    MTrack tk = tracks(i);

    fillTrackHistos(h, "trkall", tk, pv0_ptr);

    if (MC_)
    // histogram the fraction of truth-matched tracks vs z, requires TP and should not be here!!!
    {
      if (tracks(i).is_matched()) {
        Fill(h, "matchedallfractionvsz", tk.z(), 1.);
        Fill(h, "unmatchedallfractionvsz", tk.z(), 0.);
      } else {
        Fill(h, "unmatchedallfractionvsz", tk.z(), 1.);
        Fill(h, "matchedallfractionvsz", tk.z(), 0.);
      }
    }

    if (MC_ && f4D_ && tk.is_matched() && tk.has_timing())
      {
	if ((tk.dt() >0) && (tk.dt() < 0.1) && (fabs(tk.tres()) > 5. * tk.dt()))
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

  for( auto tk : tracks){
    
    if ( tk.selected() ){

      fillTrackHistos(h, "trksel", tk, pv0_ptr);
      if((tk.has_timing()) && (tk.timeQuality()>0.8)) {fillTrackHistos(h, "trkseltiming", tk, pv0_ptr);}
      // high pt track category
      if (abs(tk.pt()) > trkhiptmin_){
	  fillTrackHistos(h, "trkselhipt", tk);
	  if(abs(tk.eta())< trkcentraletamax_){
	    fillTrackHistos(h, "trkselhiptcentral", tk);
	  }else if(abs(tk.eta())> trkhiptmin_){
	    fillTrackHistos(h, "trkselhiptfwd", tk);
	  }
      }
      // low pt track categoy
      if (abs(tk.pt()) < trkloptmax_){
	  fillTrackHistos(h, "trksellopt", tk);
	  if(abs(tk.eta()) < trkcentraletamax_){
	    fillTrackHistos(h, "trkselloptcentral", tk);
	  }else if(abs(tk.eta())> trkloptmax_){
	    fillTrackHistos(h, "trkselloptfwd", tk);
	  }
      }
	
      nseltrks++;
      if (tk.pt() < 0.4)
        nseltrksptlt04++;
      if (tk.pt() < 0.3)
        nseltrksptlt03++;
      if (tk.pt() < 0.2)
        nseltrksptlt02++;

      auto nmiss = tk.lost_inner_hits();  
      if(nmiss == 0){
	fillTrackHistos(h, "missing_inner0", tk, pv0_ptr);
      }else if(nmiss == 1){
	fillTrackHistos(h, "missing_inner1", tk, pv0_ptr);
      }else if (nmiss == 2){
	fillTrackHistos(h, "missing_inner2", tk, pv0_ptr);
      }

      // don't have barrel layers in miniaod
      if((tk.has_hitPattern()) && (tk.has_track())){
	auto nbarrel = tk.track().hitPattern().pixelBarrelLayersWithMeasurement();
	if(nbarrel <2){
	  fillTrackHistos(h, "nbarrel_lt2", tk, pv0_ptr);
	}else if(nbarrel ==2){
	  fillTrackHistos(h, "nbarrel_eq2", tk, pv0_ptr);
	}      
      }
  
      if((abs(tk.z()) < 4) && (tk.dzError() < 0.01) && (abs(tk.eta()) > 1.4) && (abs(tk.eta()) < 1.8)){
	fillTrackHistos(h, "highetadriver", tk, pv0_ptr);
      }
      if((abs(tk.z()) < 4) && (tk.dzError() < 0.01)){
	fillTrackHistos(h, "alletadriver", tk, pv0_ptr);
      }

      if (MC_) {
        if (tk.is_matched()) {
          fillTrackHistos(h, "seltpmatched", tk);
          if(tk.pt() > trkhiptmin_) fillTrackHistos(h, "seltpmatchedhipt", tk);
	  if (tk._simEvt == 0){
            fillTrackHistos(h, "seltpmatchedSignal", tk);
	  }else{
            fillTrackHistos(h, "seltpmatchedPU", tk);
	  }
          Fill(h, "matchedselfractionvsz", tk.z(), 1.);
          Fill(h, "unmatchedselfractionvsz", tk.z(), 0.);
        } else {
          Fill(h, "matchedselfractionvsz", tk.z(), 0.);
          Fill(h, "unmatchedselfractionvsz", tk.z(), 1.);
          fillTrackHistos(h, "seltpunmatched", tk);
        }
      } // (M_C)

      // is the track tk (not) part of any vertex?
      if (tk.get_recv(vtxs.label()) == NO_RECVTX){
        fillTrackHistos(h, "sellost", tk);
      }
      if (tk.get_weight(vtxs.label()) > min_trk_in_vtx_weight_){
        fillTrackHistos(h, "selused", tk);
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
  if ((nrecvtx == 0) || (vtxs[0].isRecoFake())) {
    Fill(h, "nrectrk0vtx", nrectrks);
    Fill(h, "nseltrk0vtx", nseltrks);
  }

  for (unsigned int iv = 0; iv < vtxs.size(); iv++) {
    MVertex & v = vtxs[iv];

    // count vertex tracks
    double npt1 = 0, ntrkwgt05 = 0;
    for (auto tk : v.tracks){
      double wt = v.trackWeight(tk);
      if (wt > min_trk_in_vtx_weight_) {
        ntrkwgt05++;
        Fill(h, "trackwtgt05vsz", v.z(), wt);
      }
      if (tk->pt() > 1.0)
        npt1++;
      Fill(h, "trackwtvsz", v.z(), wt);
    }

    float ntrk = v.tracksSize();
    if (v.ndof() > 4) {
      Fill(h, "ntrkpt1vsz", v.z(), npt1);
      Fill(h, "ntrkwgt05vsz", v.z(), ntrkwgt05);
      Fill(h, "ftrkwgt05vsz", v.z(), ntrkwgt05 / ntrk);
    }

    Fill(h, "nbtksinvtx", ntrk);
    if (instBXLumi_ > 0) {
      Fill(h, "nbtksinvtxvsL", 1e3 * instBXLumi_, ntrk);
      Fill(h, "nbtksinvtxvsLPU", lumiPU_, ntrk);
      if (v.ndof() > 4)
        Fill(h, "nbtksinvtx4vsL", 1e3 * instBXLumi_, ntrk);
      if (v.ndof() > 4)
        Fill(h, "nbtksinvtx4vsLPU", lumiPU_, ntrk);
    }
    Fill(h, "nbtksinvtx2", ntrk);

    double sumw = v.sumw();
    if (v.ndof() > 4) {
      Fill(h, "sumwvsz", v.z(), sumw);
      Fill(h, "sumntkvsz", v.z(), (double) ntrk);
      Fill(h, "vtxndfoverntk", v.ndof() / ntrk);
      Fill(h, "vtxndf2overntk", (v.ndof() + 2) / ntrk);
      Fill(h, "sumwoverntkvsz", v.z(), sumw / ntrk);
      Fill(h, "sumwoverntkvsz4", v.z(), sumw / ntrk);
      if (ntrkwgt05 > 0) {
        Fill(h, "sumwoverntkwgt05vsz", v.z(), sumw / ntrkwgt05);
      }
    }
    if ((v.ndof() > 4) && (v.ndof() < 10)) {
      Fill(h, "sumwoverntkvszlo", v.z(), sumw / ntrk);
    }
    if ((v.ndof() > 4) && (v.ndof() < 10)) {
      Fill(h, "sumwoverntkvszlo", v.z(), sumw / ntrk);
    }
    if (v.ndof() > 20) {
      Fill(h, "sumwoverntkvszhi", v.z(), sumw / ntrk);
    }
    Fill(h, "ntrkvsz", v.z(), ntrk);
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
  // some quick z-distributins, partly redundant with fillVertexHistos
  for(unsigned int idx = 0; idx < recVtxs->size(); idx++) {
    const auto & v = recVtxs->at(idx);
    if(v.isFake()) continue;
    double zvtx = v.z();
    Fill(h, "zvtx_rec", zvtx);
    if(select(v, 0)) Fill(h, "zvtx_sel0", zvtx);
    if(select(v, 1)) Fill(h, "zvtx_sel1", zvtx);
    if(select(v, 2)) Fill(h, "zvtx_sel2", zvtx);
    if(select(v, 3)) Fill(h, "zvtx_sel3", zvtx);
  }

  // make a z-sorted list for gap histogramming and neighbour searches
  std::vector<pair<double, unsigned int>> zrecv;

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

    // this is the test for the "empty" vertices, remove it
    for (unsigned int idx = 0; idx < zrecv.size(); idx++) {
      if ( (recVtxs->at(idx).ndof() > 4.09) && (recVtxs->at(idx).ndof() < 4.11) 
	   && (recVtxs->at(idx).tracksSize() == 0)
	   && (recVtxs->at(idx).chi2() > 0.9)
	   && (recVtxs->at(idx).chi2() < 1.1)){
	Fill(h, "indexempty", idx);
      }
    }
  }


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
  if (zrecv.size() > 0 and !(recVtxs->begin()->isFake())) {
    ztag = recVtxs->begin()->position().z();
    Fill(h, "vtxz_tag", ztag);
  }

  int nrec = 0, nrecsel = 0;
  for (reco::VertexCollection::const_iterator v = recVtxs->begin(); v != recVtxs->end(); ++v) {
    if (!(v->isFake()) && v->ndof() > 0) {
      nrec++;
      if (select(*v)) nrecsel++;
    }
  }
  Fill(h, "nrecvtx", nrec);

  if (instBXLumi_ > 0) {
    Fill(h, "nrecvtxvsLPU", lumiPU_, float(nrec));
    Fill(h, "nrecvtxvsLPUprof", lumiPU_, float(nrec));
    Fill(h, "nrecvtx4vsLPU", lumiPU_, float(nrecsel));
    Fill(h, "nrecvtx4vsLPUprof", lumiPU_, float(nrecsel));
    Fill(h, "nrecvtx4vsavgLPU", avglumiPU_, float(nrecsel));
    Fill(h, "nrecvtx4vsavgLPUprof", avglumiPU_, float(nrecsel));
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

    if (v->ndof() > 4) {
      Fill(h, "ndofvsz4", v->position().z(), v->ndof());

      Fill(h, "vtxndofvsLS", lsglb_, v->ndof());
      if (fabs(v->position().z() - vertexBeamSpot_.z0() - dzb_) > sigmaZ_) {
        Fill(h, "vtxndofzgt1vsLS", lsglb_, v->ndof());
      }
      if (v == recVtxs->begin()) {
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
          if (nrecsel == 2)
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
        }

        if ((v->ndof() > 7) && (v1->ndof() > 7)) {
          Fill(h, "zdiffvsz7", z1 - z0, zbar);
          Fill(h, "zdiffrec7", z1 - z0);
          Fill(h, "zdiffvszp7", z1 - z0, zbarp);
        }

        if ((v->ndof() > 8) && (v1->ndof() > 8)) {
          Fill(h, "zdiffvsz8", z1 - z0, zbar);
          Fill(h, "zdiffrec8", z1 - z0);
          Fill(h, "zdiffvszp8", z1 - z0, zbarp);
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
    /* moved
    Fill(h, "LPUvsLPU", lumiPU_, lumiPU_);
    Fill(h, "LPUvsavgLPU", avglumiPU_, lumiPU_);
    Fill(h, "avgLPUvsLPU", lumiPU_, avglumiPU_);
    Fill(h, "nselvtxvsLPU", lumiPU_, nvtxsel);
    Fill(h, "nselvtxvssimPU", simPU_, nvtxsel);
    */
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


    Fill(h, "nselvtxvsavgLPU", avglumiPU_, nvtxsel);
    Fill(h, "nselvtx1sigmaz", nvtxsel1sigmaz);
    Fill(h, "nselvtx2sigmaz", nvtxsel2sigmaz);

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

    Fill(h, "nvtxgt1vslt1sigmaz", float(nvtxsel1sigmaz), float(nvtxselgt1sigmaz));

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
