
// -*- C++ -*-
//
// Package:    MyPrimaryVertexAnalyzer4PU
// Class:      MyPrimaryVertexAnalyzer4PU
//
/**\class PrimaryVertexAnalyzer4PU PrimaryVertexAnalyzer4PU.cc Validation/RecoVertex/src/PrimaryVertexAnalyzer4PU.cc

 Description: primary vertex analyzer for events with pile-up

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Wolfram Erdmann

#define NOT_MATCHED 666666
#define NOT_ASSIGNED 666667
#define NO_RECVTX 999999

#define DO_BEAMSPOT_ANALYSIS false
#define DO_ZVSZ_ANALYSIS false
#define DO_DENSITY_ANALYSIS false
#define DO_ZDIFFVSZ_ANALYSIS false

// system include files
#include <memory>
#include <string>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

//generator level
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"
#include "HepMC/GenParticle.h"

// vertex stuff
/////#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

// simulated vertices,..., add <use name=SimDataFormats/Vertex> and <../Track>
#include <SimDataFormats/Vertex/interface/SimVertex.h>
#include <SimDataFormats/Vertex/interface/SimVertexContainer.h>
#include <SimDataFormats/Track/interface/SimTrack.h>
#include <SimDataFormats/TrackingHit/interface/PSimHitContainer.h>
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

#include <SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h>
#include <SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h>
#include <SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h>

#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"

#include "SimDataFormats/EncodedEventId/interface/EncodedEventId.h"


#include "DataFormats/Luminosity/interface/LumiInfo.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"// obsolete??
#include "DataFormats/Luminosity/interface/LumiSummary.h"// obsolete??
#include "DataFormats/Scalers/interface/LumiScalers.h"// obsolete??


// AOD
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//Track et al
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

// Pileup
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// vertex associator
#include "SimDataFormats/Associations/interface/VertexToTrackingVertexAssociator.h"

// for pixel cluster stuff
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/Common/interface/DetSetVector.h"
// TCDS
#include <DataFormats/TCDS/interface/TCDSRecord.h>

// Root
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TObjString.h>
#include <TString.h>

// timing
#include <chrono>

#include "RecoVertex/PrimaryVertexProducer/interface/TrackFilterForPVFinding.h"

typedef reco::Vertex::trackRef_iterator trackit_t;

// class declaration
class PrimaryVertexAnalyzer4PU : public edm::EDAnalyzer {
  typedef math::XYZTLorentzVector LorentzVector;
  typedef reco::TrackBase::ParameterVector ParameterVector;

  struct SimPart {
    ParameterVector par;
    int type;      // 0 = primary
    double zdcap;  // z@dca' (closest approach to the beam
    double ddcap;
    double zvtx;  // z of the production vertex
    double xvtx;  // x of the production vertex
    double yvtx;  // y of the production vertex
    int pdg;      // particle pdg id
    int rec;
    int simpvidx;  // index of the primary vertex
    double ldec;   // distance of the production vertex from the primary vertex
    double eta;    // convenience
    double phi;
  };

  // auxiliary class holding simulated primary vertices
  class simPrimaryVertex {
  public:
    simPrimaryVertex(double x1, double y1, double z1, double t1) : x(x1), y(y1), z(z1), t(t1), ptsq(0), nGenTrk(0) {
      ptot.setPx(0);
      ptot.setPy(0);
      ptot.setPz(0);
      ptot.setE(0);
      p4 = LorentzVector(0, 0, 0, 0);
      type = 0;
      sumpT = 0;
      pt_hat = 0;
      is_visible = false;
      nTrk = 0;
      nTrkPrim = 0;
      nTrkSec = 0;
    };
    int type;  // 0=not defined, 1=full,  2 = from PileupSummaryInfo
    double x, y, z, t;
    HepMC::FourVector ptot;
    LorentzVector p4;
    double ptsq;
    double sumpT;
    int nGenTrk;
    int nTrk, nTrkPrim, nTrkSec;  // charged particles inside acceptance (--> getSimPVs)
    float pt_hat;
    bool is_visible;

    //int nMatchedTracks;  --> EXTRAS
    //int event;
    EncodedEventId eventId;
    std::vector<int> finalstateParticles;
    std::vector<int> simTrackIndex;
    std::vector<int> matchedRecTrackIndex;
    std::vector<int> genVertex;
    std::vector<reco::Track> reconstructedTracks;
    const reco::Vertex* recVtx;

    // temporary entry, to be removed when "mergedrate" issue is resolved
    double closest_vertex_distance_z;
  };

  class RecoTrack;  // forward declaration

  
  // auxiliary class holding simulated events
  class SimEvent {
  public:
    SimEvent(unsigned int idx) {
      index = idx;
      type = 0;
      nChTP = 0;
      ptvis = 0;
      pxvis = 0;
      pyvis = 0;
      z = -99;
      sumpt2rec = 0.;
      sumpt2 = 0;
      sumpt = 0;
      Tc = -1;
      dzmax = 0;
      dztrim = 0;
      chisq = 0;
      trkidx.clear();

      nwosmatch = 0;
      nwntmatch = 0;
      wnt.clear();  // formerly known as recvnt
      wos.clear();

      wos_dominated_recv.clear();    // list of wos dominated rec vertices
      matchQuality = 0;

      rec = NOT_MATCHED;  // index of a matched rec vertex, if any
      ndof = 0;           // ndof of the matched rec vertex
      zrec = 1000.;
    };

    unsigned int index;  // =index in the SimEvent list
    EncodedEventId eventId;
    int type;
    // 0=not filled,
    //1=full (e.g. from TrackingParticles),
    //2=partially filled (from PileUpSummary)

    double x, y, z, t;
    double xfit, yfit, zfit, tfit;
    int nChTP;
    double ptvis, pxvis, pyvis;

    std::vector<const TrackingParticle*> tp;

    std::vector<reco::TransientTrack> tk;              // deprecated, use rtk instead
    std::vector<reco::TransientTrack> tkprim;          // deprecated, use rktprim insted
    std::vector<reco::TransientTrack> tkprimsel;       // deprecated, use rtkprimsel instead
    std::vector<edm::RefToBase<reco::Track> > trkref;  // deprecated, don't use at all

    std::vector<RecoTrack> rtk;              // all selected RecoTracks matched to this simevent
    std::vector<RecoTrack> rtkprim;          // all selected RecoTracks matched to this simevent that come from within 5 microns of the simvertex
    std::vector<RecoTrack> rtkprimsel;       // all selected RecoTracks matched to this simevent that pass the ip/sigma(ip) < 4 cut (rec wrt beam) ?????
    std::vector<unsigned int> trkidx;        // list of RecoTrack indices, used anywher?
    
    double Tc, chisq, dzmax, dztrim, m4m2;  // filled by getSimEvents via "getTc"
    double sumpt2, sumpt;                   // filled by getSimEvents
    double sumpt2rec;                       // who fills this?
    
    // rec vertex matching                 // this block should be obsolete
    int nmatch, nmatch2;
    double zmatchn, zmatchn2;
    double pmatchn, pmatchn2;
    double wmatch;
    double zmatchw;

    unsigned int nwosmatch;  // number of recvertices dominated by this simevt (by wos)
    unsigned int nwntmatch;  // number of recvertices dominated by this simevt  (by nt)
    std::vector<unsigned int> wos_dominated_recv;  // list of dominated recv (by wos, size==nwosmatch)

    // number of tracks in recvtx with index i (VertexCollection->at(i))
    std::map<unsigned int, int> ntInRecVi;

    std::map<unsigned int, double> wnt;     // weighted number of tracks in recvtx (by index)
    std::map<unsigned int, double> wos;     // sum of wos in recvtx (by index) // oops -> this was int before 04-22
    double sumwos;                          // sum of wos in any recvtx
    double sumwt;                           // sum of weighted tracks


    // index of matched rec vertex or NOT_MATCHED
    unsigned int rec;
    unsigned int matchQuality;
    double ndof;
    double zrec;

    bool matched() const { return (rec != NOT_MATCHED); }
    bool is_signal() const { return (index == 0); }

    void addTrack(unsigned int irecv, double twos, double wt) {
      sumwt += wt;
      if (wnt.find(irecv) == wnt.end()) {
        wnt[irecv] = wt;
      } else {
        wnt[irecv] += wt;
      }

      sumwos += twos;
      if (wos.find(irecv) == wos.end()) {
        wos[irecv] = twos;
      } else {
        wos[irecv] += twos;
      }
    };

    bool hasRecoTrack(const edm::RefToBase<reco::Track>& t) {
      for (unsigned int i = 0; i < trkref.size(); i++) {
        if (t.key() == trkref[i].key()) {
          return true;
        }
      }
      return false;
    }

    int countVertexTracks(const reco::Vertex& v, double min_weight = 0.5) {
      // count the number of tracks from this simevent in recvertex v
      int n = 0;
      for (trackit_t t = v.tracks_begin(); t != v.tracks_end(); t++) {
        if (v.trackWeight(*t) >= min_weight) {
          for (unsigned int i = 0; i < trkref.size(); i++) {
            if (t->key() == trkref[i].key()) {
              n++;
              break;
            }
          }
        }
      }
      return n;
    }

    int countVertexTracksMTD(const reco::Vertex& v, double min_weight = 0.5) {
      int n = 0;
      for (trackit_t t = v.tracks_begin(); t != v.tracks_end(); t++) {
        if ((v.trackWeight(*t) >= min_weight) && ((**t).pt() > 0.9) && (fabs((**t).eta()) < 2.7)) {
          for (unsigned int i = 0; i < trkref.size(); i++) {
            if (t->key() == trkref[i].key()) {
              n++;
              break;
            }
          }
        }
      }
      return n;
    }

    unsigned int max_wos_vtx() {
      double maxwos = 0.;
      unsigned int vtx = NOT_MATCHED;
      for (auto it : wos) {
        if (it.second > maxwos) {
          maxwos = it.second;
          vtx = it.first;
        }
      }
      return vtx;  // this is often, but not necessarily the matched vertex
    }

    unsigned int max_nt_vtx() {
      double maxnt = 0.;
      unsigned int vtx = NOT_MATCHED;
      for (auto it : wnt) {
        if (it.second > maxnt) {
          maxnt = it.second;
          vtx = it.first;
        }
      }
      return vtx;  // this is often, but not necessarily the matched vertex
    }

    bool was_found(double min_ndof = 0.) { return (rec != NOT_MATCHED) && (ndof > min_ndof); }

    void clear_matching_info(){
      rec = NOT_MATCHED;
      matchQuality = 0;

      wos.clear();
      sumwos = 0.;
      nwosmatch = 0;
      wos_dominated_recv.clear();

      wnt.clear();
      nwntmatch = 0;
    }
    
    bool operator==(const SimEvent& rh) const { return (index == rh.index); }
  };

  /* helper class holding recvertex -> simvertex matching information */
  // FIXME  create a RecoVertex class that holds the recvtx, the RSmatch and
  // auxilary information so we don't need to do the same things over and over
  class RSmatch {
  public:
    RSmatch() {
      wos.clear();
      wnt.clear();
      wosmatch = NOT_MATCHED;
      wntmatch = NOT_MATCHED;
      sumwos = 0;
      sumwnt = 0;
      maxwos = 0.;
      maxwnt = 0;
      maxwosnt = 0;

      matchQuality = 0;
      sim = NOT_MATCHED;
    }

    void addTrack(unsigned int iev, double twos, double twt) {
      sumwnt += twt;
      if (wnt.find(iev) == wnt.end()) {
        wnt[iev] = twt;
      } else {
        wnt[iev] += twt;
      }

      sumwos += twos;
      if (wos.find(iev) == wos.end()) {
        wos[iev] = twos;
      } else {
        wos[iev] += twos;
      }
    }

    bool is_real() { return (matchQuality > 0) && (matchQuality < 99); }

    bool is_fake() { return (matchQuality <= 0) || (matchQuality >= 99); }

    bool is_signal() { return (sim == 0); }

    int split_from() {
      if (is_real())
        return -1;
      if ((maxwos > 0) && (maxwos > 0.3 * sumwos))
        return wosmatch;
      return -1;
    }
    bool other_fake() { return (is_fake() & (split_from() < 0)); }

    std::map<unsigned int, double> wos;  // simevent -> wos
    std::map<unsigned int, double> wnt;  // simevent -> weighted number of truth matched tracks
    unsigned int wosmatch;               // index of the simevent providing the largest contribution to wos
    unsigned int wntmatch;               // index of the simevent providing the highest number of tracks
    double sumwos;                       // total sum of wos of all truth matched tracks
    double sumwnt;                       // total weighted number of truth matchted tracks
    double maxwos;                       // largest wos sum from one sim event (wosmatch)
    double maxwnt;                       // largest weighted  number of tracks from one sim event (ntmatch)
    int maxwosnt;                        // number of tracks from the simevt with highest wos
    unsigned int sim;                    // best match  (NO_MATCH if not matched)
    unsigned int matchQuality;           // quality flag
  };

  
  /* RecoVertexExtras, container for additional recvertex info
     
   */
  class RecoVertexExtras {
  public:
    RecoVertexExtras(unsigned int idx, const reco::Vertex* recvtx) {
      index = idx;
      vtx = recvtx;
      deltaz = 1000.;
    }
    unsigned int index;
    const reco::Vertex* vtx;
    double deltaz; // signed distance to the neares selected recvertex (1000 if there is no other)
  };


  
  /* RecoTrack is helper class for collecting rectrk related information */
  class RecoTrack {
  public:
    RecoTrack(unsigned int idx,
              const reco::Track* rectrk,
              reco::TransientTrack* transienttrk,
              unsigned int trkkey,
              bool f4D) {
      index = idx;
      trk = rectrk;
      tt = transienttrk;
      key = trkkey;

      z = (tt->stateAtBeamLine().trackStateAtPCA()).position().z();
      dz = trk->dzError();
      pt = trk->pt();
      eta = trk->eta();
      phi = (tt->stateAtBeamLine().trackStateAtPCA()).momentum().phi();
      theta = (tt->stateAtBeamLine().trackStateAtPCA()).momentum().theta();
      Measurement1D atIP = tt->stateAtBeamLine().transverseImpactParameter();  // error contains beamspot
      ip = atIP.value();
      dip = atIP.error();

      selected = false;
      has_timing = false;
      t = 0;
      dt = 1e10;
      timeQuality = -1.;  // will stay -1. for no timing info
      MTD_pathlength = 0;
      MTD_time = -1.;
      MTD_timeerror = 1e10;
      MTD_momentum = 0;
      if (f4D) {
        t = tt->timeExt();
        dt = tt->dtErrorExt();
        if ((dt > 0) && (dt < 0.5) && (!((t==0) && (std::abs(dt-0.35)<1e-5)))) {
          has_timing = true;
	  timeQuality = 1.; // temporary, overridden externally for now
        }
      }

      recv_index = NO_RECVTX;  // deprecated
      recv.clear();

      // filled later by getSimEvents
      matched = false;
      simEvt = NULL;
      zsim = 0;
      tsim = 0;
      is_primary = false;
    }

    unsigned int get_recv(std::string& vtxcollection) {
      if (recv.find(vtxcollection) == recv.end()) {
        return NO_RECVTX;
      } else {
        return recv[vtxcollection];
      }
    }

    bool is_pion(){
      return (matched && (!tpr.isNull()) && (abs(tpr->pdgId()) == 211));
    }
    
    bool is_kaon(){
      return (matched && (!tpr.isNull()) && (abs(tpr->pdgId()) == 321));
    }
    bool is_proton(){
      return (matched && (!tpr.isNull()) && (abs(tpr->pdgId()) == 2212));
    }
    double get_particle_mass(){
      if (matched &&  (!tpr.isNull())){
	return tpr->mass();
      }else{
	return 0;
      }
    }
    double get_t_pid(double mass =0){ //mass corrected reconstructed time at the beamline
      if (mass == 0 ){ // use the true mass (if known)
	mass = get_particle_mass();
      }
      if ((mass > 0) && (MTD_pathlength > 0) ){
	double gammasq= 1. + MTD_momentum * MTD_momentum / (mass * mass);
	double v = 2.99792458e1 * std::sqrt(1. - 1. / gammasq);  // cm / ns
	return MTD_time - MTD_pathlength / v;
      }else{
	return -100.;
      }
    }
    
    unsigned int index;
    unsigned int key;
    const reco::Track* trk;  // FIXME get rid of this
    reco::TransientTrack* tt;

    bool selected;  // result of the track filter

    // convenience
    double z;
    double dz;
    double pt, eta, phi, ip, dip, theta;
    bool has_timing;
    double t;
    double dt;
    // the following variables are filled in PrimaryVertexAnalyzer4PU::get_reco_and_transient_tracks
    double timeQuality;
    double MTD_pathlength, MTD_time, MTD_timeerror, MTD_momentum;
    double th[3];  // track time for particle hypotheses : 0=pion, 1=kaon, 2=proton

    // filled later if available
    unsigned int recv_index;  // deprecated
    std::map<std::string, unsigned int> recv;

    // MC truth related
    bool matched;
    SimEvent* simEvt;
    //TrackingParticle * tp;
    TrackingParticleRef tpr;
    // convenience
    double zsim;
    double tsim;
    bool is_primary;
  };

  /* collect information on reco tracks in one place
     basically a vector of RecoTracks with some extra info and pointers
   */

  class Tracks {
  public:
    edm::Handle<edm::View<reco::Track> > trackCollectionH;
    std::vector<RecoTrack> recotracks;
    // some pointers
    std::map<unsigned int, unsigned int> key2idx;  // get the recotrack from a reftobase key
    // e.g.  i = tracks.key2idx[vt->key()]   for track_it vt

    Tracks() {
      recotracks.clear();
      key2idx.clear();
    }

    RecoTrack& operator()(int i) { return recotracks.at(i); }

    edm::RefToBase<reco::Track> ref(unsigned int i) { return edm::RefToBase(trackCollectionH, i); }

    RecoTrack& from_key(const unsigned int key){
      // make the map if needed
      if (key2idx.size() != recotracks.size()) {
        for (unsigned int i = 0; i < recotracks.size(); i++) {
          key2idx[recotracks[i].key] = i;
        }
      }
      return recotracks[key2idx[key]];
    }

    RecoTrack& from_ref(const edm::RefToBase<reco::Track> ref) { return from_key(ref.key()); }

    unsigned int simevent_index_from_key(const unsigned int key){
      RecoTrack& tk = from_key(key);
      if (tk.matched){
	return tk.simEvt->index;
      }else{
	return NOT_MATCHED;
      }
    }
    // for convenience behave like a vector of RecoTrack
    unsigned int size() { return recotracks.size(); }
    void clear() { recotracks.clear(); }
    void push_back(RecoTrack t) { recotracks.push_back(t); }
  };


  
  /* helper class for histogramming/counting keys*/
  class VertexCounter {
  public:
    VertexCounter() { n.clear(); }

    void count(int key) {
      if (key == NOT_MATCHED) {
        unmatched += 1;
      } else {
        matched += 1;
        if (n.count(key) == 0) {
          n[key] = 1;
        } else {
          n[key]++;
        }
      }
    }

    int nkey() { return n.size(); }

    int nkey(int threshold) {
      int keycount = 0;
      for (auto it = n.begin(); it != n.end(); ++it) {
        if (it->second >= threshold) {
          keycount++;
        }
      }
      return keycount;
    }

    // member variables
    int matched = 0;
    int unmatched = 0;
    std::map<unsigned int, int> n;
  };

public:
  explicit PrimaryVertexAnalyzer4PU(const edm::ParameterSet&);
  ~PrimaryVertexAnalyzer4PU();

  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  /* EXTRAS
  void analyzePixelClusters(const edm::Event&, const edm::EventSetup&);
  */

  virtual void beginJob();
  virtual void endJob();

  void report_counted(const std::string msg, const int max_count);
  void report_counted(const std::string msg, const std::string msg2, const int max_count);

  void analyzeTracksTP(Tracks& tracks, std::vector<SimEvent>& simEvt);

private:
  void analyzeTracksTP(edm::Handle<edm::View<reco::Track> > trackCollectionH,
                       const reco::VertexCollection* recVtxs,
                       std::vector<SimPart>& tsim,
                       std::vector<simPrimaryVertex>& simpv,
                       std::vector<SimEvent>& simEvt);

  void printPVTrksZT(edm::Handle<edm::View<reco::Track> > trackCollectionH,
                     const reco::VertexCollection* recVtxs,
                     std::vector<SimPart>& tsim,
                     std::vector<simPrimaryVertex>& simpv,
                     std::vector<SimEvent>& simEvt);

  void reportEvent(const char*, bool dumpvertex = false);
  void reportVertex(const reco::Vertex&, const char*, bool);

  int* supf(std::vector<SimPart>& simtrks, const edm::View<reco::Track>&);
  int* supfv(std::vector<SimPart>& simtrks, const std::vector<reco::Track>&);

  static bool match(const ParameterVector& a, const ParameterVector& b);
  std::vector<SimPart> getSimTrkParameters(edm::Handle<edm::SimTrackContainer>& simTrks,
                                           edm::Handle<edm::SimVertexContainer>& simVtcs,
                                           double simUnit = 1.0);
  std::vector<SimPart> getSimTrkParameters(const edm::Handle<reco::GenParticleCollection>);
  void getTc(const std::vector<reco::TransientTrack>&, double&, double&, double&, double&, double&);
  void getTc(const std::vector<RecoTrack>&, double&, double&, double&, double&, double&);

  double vertex_pxy(const reco::Vertex&);
  double vertex_aptsum(const reco::Vertex&);
  double vertex_r(const reco::Vertex&);
  double vertex_ptmax2(const reco::Vertex&);
  double vertex_yum(const reco::Vertex&);
  double vertex_maxfrac(const reco::Vertex&);
  double vertex_sumpt2(const reco::Vertex&);
  double vertex_sumpt(const reco::Vertex&);
  bool vertex_time_from_tracks(const reco::Vertex&, Tracks& tracks, double minquality, double& t, double& tError);
  bool vertex_time_from_tracks_pid(const reco::Vertex&, Tracks& tracks, double minquality, double& t, double& tError);

  bool select(const reco::Vertex&, const int level = 0);


  void addn(std::map<std::string, TH1*>& h, TH1* hist) {
    // add a histogram in a subdirectory and use the subdirectory name in the map key
    h[gDirectory->GetName() + std::string("/") + hist->GetName()] = hist;
    hist->StatOverflows(kTRUE);
    hist->SetDirectory(gDirectory);
  }
  
  void add(std::map<std::string, TH1*>& h, TH1* hist) {
    // add a histogram
    h[hist->GetName()] = hist;
    hist->StatOverflows(kTRUE);
    hist->SetDirectory(gDirectory);
  }
  
  void add(std::map<std::string, TH1*>& h, TH1* hist, const std::string& vtype) {
    // add a histogram
    std::string hname(hist->GetName());
    std::string htitle(hist->GetTitle());
    TH1* histn = (TH1*)hist->Clone((hname + "_" + vtype).c_str());
    histn->SetTitle((htitle + " (" + vtype + ")").c_str());
    h[histn->GetName()] = histn;
    histn->StatOverflows(kTRUE);
    histn->SetDirectory(gDirectory);
    delete hist;
  }

  void addSP(std::map<std::string, TH1*>& h, TH1* hist) {

    // add a histogram + two versions (Signal and PU)
    std::string name = hist->GetName();
    std::string title = hist->GetTitle();
    //std::cout << "addSP  " << name << " gDirectory = " <<  gDirectory->GetName() << std::endl;
    h[hist->GetName()] = hist;
    hist->StatOverflows(kTRUE);
    hist->SetDirectory(gDirectory);
    TH1* hS = (TH1*)hist->Clone((name + "Signal").c_str());
    TH1* hP = (TH1*)hist->Clone((name + "PU").c_str());
    hS->SetTitle((title + "(Signal)").c_str());
    hP->SetTitle((title + "(PU)").c_str());
    h[name + "Signal"] = hS;
    h[name + "PU"] = hP;
    hP->SetDirectory(gDirectory);
    hS->SetDirectory(gDirectory);
    hP->StatOverflows(kTRUE);
    hS->StatOverflows(kTRUE);
  }
  
  void addnSP(std::map<std::string, TH1*>& h, TH1* hist) {
    // add a histogram + two versions (Signal and PU) and use the subdirectory name in the map key
    std::string name(hist->GetName());
    std::string title(hist->GetTitle());
    addn(h, hist);
    TH1* hS = (TH1*)hist->Clone((name + "Signal").c_str());
    TH1* hP = (TH1*)hist->Clone((name + "PU").c_str());
    hS->SetTitle((title + "(Signal)").c_str());
    hP->SetTitle((title + "(PU)").c_str());
    addn(h, hS);
    addn(h, hP);
  }

  void addSP_obsolete(std::map<std::string, TH1*>& h, TH1* hist, const std::string& type) {
    std::string hname(hist->GetName());
    std::string htitle(hist->GetTitle());
    TH1* histn = (TH1*)hist->Clone((hname + "_" + type).c_str());
    histn->SetTitle((htitle + " (" + type + ")").c_str());
    addSP(h, histn);
    delete hist;
  }

  std::string fillmsg(const std::string name, const double value1, const double value2=-12345){
    std::stringstream s;
    s << " " << name << " with " << value1;
    if (value2 != -12345){
      s<< "," << value2;
    }
    return s.str();
  }
  
  void Fill(std::map<std::string, TH1*>& h, std::string s, double x) {
    if (h.count(s) == 0) {
      report_counted("PrimaryVertexAnalyzer4PU::Trying to fill non-existing 1d Histogram",
		     fillmsg(s, x), 1000);
      return;
    }
    h[s]->Fill(x);
  }

  void Fill(std::map<std::string, TH1*>& h, std::string s, double x, double y) {
    //    cout << "Fill2 " << s << endl;
    if (h.count(s) == 0) {
      std::cout << "Trying to fill non-existing 2d Histogram named " << s << " with " << x << "," << y << std::endl;
      return;
    }
    h[s]->Fill(x, y);
  }

  void Fillw(std::map<std::string, TH1*>& h, std::string s, double x, double y, double w) {
    if (h.count(s) == 0) {
      std::cout << "Trying to fill non-existing 2d weighted Histogram named " << s << std::endl;
      return;
    }
    (static_cast<TH2F*>(h[s]))->Fill(x, y, w);
  }

  void Fill(std::map<std::string, TH1*>& h, std::string s, double x, bool signal, bool fill_all = true) {
    if (h.count(s) == 0) {
      std::cout << "Trying to fill non-existing Histogram named " << s << std::endl;
      return;
    }

    if(fill_all){
      h[s]->Fill(x);
    }
    
    if (signal) {
      if (h.count(s + "Signal") == 0) {
        std::cout << "Trying to fill non-existing Histogram named " << s + "Signal" << std::endl;
        return;
      }
      h[s + "Signal"]->Fill(x);
    } else {
      if (h.count(s + "PU") == 0) {
        std::cout << "Trying to fill non-existing Histogram named " << s + "PU" << std::endl;
        return;
      }
      h[s + "PU"]->Fill(x);
    }
  }

  void Fill(std::map<std::string, TH1*>& h, std::string s, double x, double y, bool signal) {
    if (h.count(s) == 0) {
      std::cout << "Trying to fill non-existing Histogram named " << s << std::endl;
      return;
    }

    h[s]->Fill(x, y);
    if (signal) {
      if (h.count(s + "Signal") == 0) {
        std::cout << "Trying to fill non-existing Histogram named " << s + "Signal" << std::endl;
        return;
      }
      h[s + "Signal"]->Fill(x, y);
    } else {
      if (h.count(s + "PU") == 0) {
        std::cout << "Trying to fill non-existing Histogram named " << s + "PU" << std::endl;
        return;
      }
      h[s + "PU"]->Fill(x, y);
    }
  }

  void Fill(std::map<std::string, TH1*>& h, std::string s, bool yesno, bool signal) {
    if (yesno) {
      Fill(h, s, 1., signal);
    } else {
      Fill(h, s, 0., signal);
    }
  }

  void Cumulate(std::map<std::string, TH1*>& h, std::string s) {
    if (h.count(s) == 0) {
      std::cout << "Trying to cumulate non-existing Histogram named " << s << std::endl;
      return;
    }

    if ((h[s]->GetEntries() == 0) || (h[s]->Integral() <= 0)) {
      //std::cout << "DEBUG : Cumulate called with empty histogram " << h->GetTitle() << std::endl;
      return;
    }
    try {
      h[s]->ComputeIntegral();
      Double_t* integral = h[s]->GetIntegral();
      h[s]->SetContent(integral);
    } catch (...) {
      std::cout << "DEBUG : an error occurred cumulating  " << h[s]->GetTitle() << std::endl;
    }
  }

  void Cumulate(TH1* h) {
    if ((h->GetEntries() == 0) || (h->Integral() <= 0)) {
      //std::cout << "DEBUG : Cumulate called with empty histogram " << h->GetTitle() << std::endl;
      return;
    }
    try {
      h->ComputeIntegral();
      Double_t* integral = h->GetIntegral();
      h->SetContent(integral);
    } catch (...) {
      std::cout << "DEBUG : an error occurred cumulating  " << h->GetTitle() << std::endl;
    }
  }

  std::map<std::string, TH1*> bookVertexHistograms(TDirectory * dir);
  void bookTrackHistograms(const char * directory_name);
  void bookSimPVHistograms(const char * directory_name);
  void bookEventHistograms(const char * directory_name);
  //void bookClusterHistograms(); --> EXTRAS

  void get_luminosity_infos(const edm::Event& iEvent);
  void get_particle_data_table(const edm::EventSetup&);
  bool get_beamspot_data(const edm::Event&);
  bool get_reco_and_transient_tracks(const edm::EventSetup&, const edm::Event&, Tracks&);

  double muvtx(double z);

  bool matchVertex(const simPrimaryVertex& vsim, const reco::Vertex& vrec);
  bool isResonance(const HepMC::GenParticle* p);
  bool isFinalstateParticle(const HepMC::GenParticle* p);
  bool isCharged(const HepMC::GenParticle* p);
  void fillVertexHistosNoTracks(std::map<std::string, TH1*>& h,
                                const std::string& vtype,
                                const reco::Vertex* v = NULL,
				const int index = -1, 
                                const double deltaz = 0,
                                const bool verbose = false);
  void fillVertexHistos(std::map<std::string, TH1*>& h,
                        const std::string& vtype,
                        const reco::Vertex* v,
                        Tracks& tracks,
			const int index = -1,
                        const double deltaz = 0,
                        const bool verbose = false);
  void fillVertexHistosMatched(std::map<std::string, TH1*>& h,
			       const std::string& vtype,
			       const reco::Vertex* v,
			       Tracks& tracks,
			       const RSmatch& rs,
			       const std::vector<SimEvent>& simEvt,
                               const int index = -1,
			       const double deltaz = 0,
			       const bool verbose = false);
  
  void fillTrackHistos(std::map<std::string, TH1*>& h,
                       const std::string& ttype,
                       RecoTrack& tk,
                       const reco::Vertex* v = NULL);
  void fillTrackHistosMatched(std::map<std::string, TH1*>& h,
                       const std::string& ttype,
                       RecoTrack& tk);
 void fillTransientTrackHistos(std::map<std::string, TH1*>& h,
                                const std::string& ttype,
                                const reco::TransientTrack* tt,
                                const reco::Vertex* v = NULL);
  void fillRecoTrackHistos(std::map<std::string, TH1*>& h, const std::string& ttype, const reco::Track& t);
  void fillTrackClusterHistos(std::map<std::string, TH1*>& h,
                              const std::string& ttype,
                              const reco::Track& t,
                              const reco::Vertex* v = NULL);
  void dumpHitInfo(const reco::Track& t);
  void printRecTrks(const edm::View<reco::Track>& recTrks);
  void printRecVtxs(const reco::VertexCollection* recVtxs, std::string title = "Reconstructed Vertices");
  void printSimVtxs(const edm::Handle<edm::SimVertexContainer> simVtxs);
  void printSimTrks(const edm::Handle<edm::SimTrackContainer> simVtrks);
  bool getPuInfo(const edm::Event& iEvent, PileupSummaryInfo& puInfo);
  std::vector<simPrimaryVertex> getSimPVs(const edm::Handle<edm::HepMCProduct> evtMC);
  std::vector<simPrimaryVertex> getSimPVs(const edm::Handle<reco::GenParticleCollection>);
  std::vector<simPrimaryVertex> getSimPVs(const edm::Handle<edm::SimVertexContainer> simVtxs,
                                          const edm::Handle<edm::SimTrackContainer> simTrks);
  std::vector<simPrimaryVertex> getSimPVs(const edm::Handle<TrackingVertexCollection>);
  std::vector<simPrimaryVertex> getSimPVs(const edm::Handle<TrackingVertexCollection>,
                                          std::vector<PrimaryVertexAnalyzer4PU::simPrimaryVertex>&);

  Int_t getAssociatedRecoTrackIndex(const edm::Handle<reco::TrackCollection>& recTrks, TrackingParticleRef tpr);
  bool truthMatchedTrack(edm::RefToBase<reco::Track>, TrackingParticleRef&);
  std::vector<edm::RefToBase<reco::Track> > getTruthMatchedVertexTracks(const reco::Vertex&, double min_weight = 0.5);
  void printTruthMatchValues(edm::RefToBase<reco::Track> track);

  bool get_MC_truth(const edm::Event& iEvent,
                    Tracks& tracks,
                    bool bPuInfo,
                    PileupSummaryInfo& puInfo,
                    std::vector<SimEvent>& simEvt,
                    std::vector<simPrimaryVertex>& simpv,
                    std::vector<SimPart>& tsim);

  void fill_simvtx_histos(std::vector<simPrimaryVertex>& simpv);

  std::vector<PrimaryVertexAnalyzer4PU::SimEvent> getSimEvents(edm::Handle<TrackingParticleCollection>,
                                                               Tracks& tracks);

  std::vector<PrimaryVertexAnalyzer4PU::SimEvent> getSimEvents_notp(
								    const edm::Handle<edm::SimTrackContainer> simTrks,
								    const edm::Handle<edm::SimVertexContainer> simVtxs,
								    Tracks& tracks);

  void analyzeVertexRecoCPUTime(std::map<std::string, TH1*>& h,
                                const reco::VertexCollection* recVtxs,
                                const std::string message = "");
  void analyzeVertexCollectionRecoNoTracks(std::map<std::string, TH1*>& h,
                                           const reco::VertexCollection* recVtxs,
                                           const std::string message = "");

  void analyzeVertexCollectionReco(std::map<std::string, TH1*>& h,
                                   const reco::VertexCollection* recVtxs,
                                   Tracks& tracks,
                                   const std::string message = "");

  void analyzeVertexCollectionSimPvNoSimTracks(std::map<std::string, TH1*>& h,
                                               const reco::VertexCollection* recVtxs,
                                               Tracks& tracks,
                                               std::vector<simPrimaryVertex>& simpv,
                                               const std::string message = "");

  void analyzeVertexCollectionSimPv(std::map<std::string, TH1*>& h,
                                    const reco::VertexCollection* recVtxs,
                                    Tracks& tracks,
                                    std::vector<simPrimaryVertex>& simpv,
                                    std::vector<PrimaryVertexAnalyzer4PU::SimPart>& tsim,
                                    const std::string message = "");

  void analyzeVertexCollectionDQMMC(std::map<std::string, TH1*>& h,
                                    const reco::VertexCollection* recVtxs,
                                    std::vector<simPrimaryVertex>& simpv,
                                    std::vector<SimPart>& tsim,
                                    const std::string message);

  void analyzeVertexMergeRateTP(std::map<std::string, TH1*>& h,
                                const reco::VertexCollection* recVtxs,
                                std::vector<SimEvent>& simEvt,
                                std::vector<RSmatch>& recvmatch,
                                const std::string message = "");

  void analyzeVertexCollectionTP(std::map<std::string, TH1*>& h,
                                 const reco::VertexCollection* recVtxs,
                                 Tracks& tracks,
                                 std::vector<SimEvent>& simEvt,
                                 std::vector<RSmatch>& recvmatch,
                                 const std::string message = "");

  void analyzeVertexCollectionPtvis(std::map<std::string, TH1*>& h,
				    const reco::VertexCollection* recVtxs,
				    Tracks& tracks,
				    std::vector<SimEvent>& simEvt,
				    std::vector<RSmatch>& recvmatch,
				    const std::string message = "");

  void analyzeRecVertexComposition(std::map<std::string, TH1*>& h,
                                   const reco::Vertex& v,
				   Tracks& tracks,
                                   RSmatch& rs,
                                   std::vector<SimEvent>& simEvt,
				   float npu);

  std::vector<RSmatch> tpmatch(const reco::VertexCollection* recVtxs,
                               std::vector<SimEvent>& simEvt,
                               Tracks& tracks);
  
  void wos_match(std::vector<RSmatch> & recvmatch,
                               const reco::VertexCollection* recVtxs,
                               std::vector<SimEvent>& simEvt,
                               Tracks& tracks);

  void nwt_match(std::vector<RSmatch> & recvmatch,
                               const reco::VertexCollection* recVtxs,
                               std::vector<SimEvent>& simEvt,
                               Tracks& tracks);

  void ihmatch(const reco::VertexCollection* recVtxs, std::vector<SimEvent>& simEvt, std::vector<RSmatch>& recvmatch);

  std::string formatMatchList(std::map<unsigned int, double>&, unsigned int nfield, bool sim);

  void printMatchingSummary(const reco::VertexCollection* recVtxs,
                            std::vector<SimEvent>& simEvt,
                            std::vector<RSmatch>& recvmatch,
                            const std::string message);

  void printEventSummary(std::map<std::string, TH1*>& h,
                         const reco::VertexCollection* recVtxs,
                         Tracks& tracks,
                         std::vector<SimEvent>& simEvt,
                         std::vector<RSmatch>& recvmatch,
                         const std::string message);

  void printEventSummary(std::map<std::string, TH1*>& h,
                         const reco::VertexCollection* recVtxs,
                         Tracks& tracks,
                         std::vector<simPrimaryVertex>& simpv,
                         const std::string message);

  reco::VertexCollection* vertexFilter(edm::Handle<reco::VertexCollection>, bool filter);

  void compareCollections(std::vector<SimEvent>& simEvt, std::vector<simPrimaryVertex>& simpv);

  void history(const edm::Handle<edm::View<reco::Track> >& tracks, const size_t trackindex = 10000);
  std::string particleString(int) const;
  std::string vertexString(TrackingParticleRefVector, TrackingParticleRefVector) const;

  std::string vertexString(HepMC::GenVertex::particles_in_const_iterator,
                           HepMC::GenVertex::particles_in_const_iterator,
                           HepMC::GenVertex::particles_out_const_iterator,
                           HepMC::GenVertex::particles_out_const_iterator) const;

  std::vector<bool> trackClass(const reco::Track&);

  
  void set_ndof_globals(std::string & vertexcollection){
    // ndof with and without beam constraint:
    //                                 0    1   2   3   4   5 
    // no BS    ndof = 2 * nt - 3     -3   -1   1   3   5   7
    // withBS   ndof = 2 * nt - 1     -1    1   3   5   7
    // with an adaptive fitter, replace nt by <w>
    // hence <w> = (ndof - ndof0trk_) / 2.
    if (vertexcollection.find("WithBS") != std::string::npos){
      selNdof_ = selNdofWithBS_;
      ndof0trk_ = -1.;
    }else{
      selNdof_ = selNdofNoBS_;
      ndof0trk_ = -3.;
    }
  }

  //timers
  std::map<std::string, std::chrono::time_point<std::chrono::steady_clock> > timer_start_;
  std::map<std::string, double> timers_;
  void inline timer_start(const std:: string label){timer_start_[label] = std::chrono::steady_clock::now();}
  void timer_stop(const std::string & label){
    auto stop = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - timer_start_[label]);
    timers_[label] += duration.count();
  }

  /*  --> EXTRAS
  bool getBXLumi_from_files();

  int time_of_last_reset(const edm::Event& iEvent);

  bool load_lumisection_lumi_from_file();
 uint32_t orbitOfLastReset_;
  int tres_;
  */

  // ----------member data ---------------------------

  // in the order of initialization
  bool verbose_;
  bool veryverbose_;
  bool do_pixel_cluster_analysis_;
  bool do_pixel_cluster_analysis_all_layers_;
  bool do_vertex_analysis_;
  bool doMatching_;
  bool dumpThisEvent_;
  bool forceDump_;
  int analyzeLS_;  //
  bool fill_track_histos_;
  bool DEBUG_;
  int eventcounter_;
  int dumpcounter_;
  int emptyeventcounter_;
  int autoDumpCounter_;
  int ndump_;
  int ndump_tracks_;
  long long int run_;
  int luminosityBlock_;
  long long int event_;
  int bunchCrossing_;
  uint32_t orbitNumber_;
  double fBfield_;
  double simUnit_;
  double simtUnit_;
  double zmatch_;
  int eventSummaryCounter_;
  int nEventSummary_;
  double zWosMatchMax_;  // cut-off for wos matching, relax for exotica

  int nEventNsel_;           // counts events entering certain histograms
  unsigned int nbindzcorr_;  // number of bins for correlation counting

  double sigmaZoverride_;
  double sigmaZ_;
  double sigmaT_;
  unsigned int nPUmin_;
  unsigned int nPUmax_;
  bool useVertexFilter_;
  int bxFilter_;
  double wxy2_, wx_, wy_;

  std::string outputFile_;  // output file
  TrackFilterForPVFinding theTrackFilter;
  std::string recoTrackProducer_;
  std::string trackAssociatorLabel_;
  double trackAssociatorMin_;
  TObjString* info_;
  TObjString* build_;
  std::vector<std::string> vtxSample_;  // make this a a vector to keep cfg compatibility with PrimaryVertexAnalyzer
  TFile* rootFile_;
  //  TDirectory * root_dir_; // temporary for backward compatibility of "add"

  edm::InputTag simG4_;

  edm::ESHandle<ParticleDataTable> pdt_;
  math::XYZPoint myBeamSpot;
  // local counters

  int nCompareCollections_;

  int nfake_;
  int npair_;
  int currentLS_;
  bool MC_;

  unsigned int minNumberOfRecTrks_;
  int minNumberOfSelTrks_;

  std::vector<std::string> vertexCollectionLabels_;
  std::map<std::string, edm::EDGetTokenT<reco::VertexCollection> > vertexCollectionTokens_;
  std::map<std::string, std::map<std::string, TH1*> > histograms_;
  std::map<std::string, reco::VertexCollection*> recVtxs_;
  std::map<std::string, std::vector<RSmatch> > recvmatch_;
  int nihmatch_;
  int matchsummaries_;

  std::map<std::string, TH1*> hsimPV;
  std::map<std::string, TH1*> hTrk;
  std::map<std::string, TH1*> hEvt;
  std::map<std::string, TH1*> hClu;

  unsigned int reset_nbin_;  //
  double reset_period_;
  unsigned int max_LS_;  //

  std::map<unsigned int, TrackingParticleRef> trkidx2tp_;  // reco::track index    --> tracking particle
  std::map<unsigned int, unsigned int> trkidx2simevt_;
  std::map<unsigned int, unsigned int> trkkey2simevt_;
  std::map<unsigned int, unsigned int> trkkey2recvtx_;

  std::map<unsigned int, reco::TransientTrack*> trkkey2ttrk_obsolete;

  reco::BeamSpot vertexBeamSpot_;
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle_;
  edm::ESHandle<TransientTrackBuilder> theB_;

  std::vector<reco::TransientTrack> t_tks_;  // FIXME  --> deprecated

  std::map<std::string, std::pair<unsigned int,unsigned int> > counted_messages_;


  // control verbosity for some report types
  static constexpr bool dump_signal_vertex_not_tpmatched_ = true;
  static constexpr bool dump_fake_vertex_on_top_of_signal_ = false;
  static constexpr bool dump_big_fakes_ = false;
  
  bool f4D_;
  bool fTrackTime_;

  bool RECO_;
  double instBXLumi_;
  double avginstBXLumi_;

  //const
  const TrackerTopology* tTopo_;

  double lumiHistoRange_;
  double lumiPUHistoRange_;
  double sigma_pp_;  // [ub]  pp-cross section for PU

  edm::EDGetTokenT<TCDSRecord> tcdsrecord_;

  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > vecPileupSummaryInfoToken_;

  edm::EDGetTokenT<reco::GenParticleCollection> genParticleCollection_Token_;
  edm::EDGetTokenT<reco::TrackCollection> recoTrackCollectionToken_;
  // see RecoVertex/Configuration/python/RecoVertex_phase2_timing_cff.py
  edm::EDGetTokenT<edm::ValueMap<float> > trkTimesToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > trkTimeResosToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > trkTimeQualityToken_;
  double trkTimeQualityThreshold_;
  edm::EDGetTokenT<edm::ValueMap<float> > MTD_pathlength_Token_;
  edm::EDGetTokenT<edm::ValueMap<float> > MTD_time_Token_;
  edm::EDGetTokenT<edm::ValueMap<float> > MTD_timeerror_Token_ ;
  edm::EDGetTokenT<edm::ValueMap<float> > MTD_momentum_Token_;
  
  edm::EDGetTokenT<reco::BeamSpot> recoBeamSpotToken_;
  edm::EDGetTokenT<edm::View<reco::Track> > edmView_recoTrack_Token_;
  edm::EDGetTokenT<edm::SimVertexContainer> edmSimVertexContainerToken_;
  edm::EDGetTokenT<edm::SimTrackContainer> edmSimTrackContainerToken_;
  edm::EDGetTokenT<edm::HepMCProduct> edmHepMCProductToken_;
  edm::EDGetTokenT<reco::TrackToTrackingParticleAssociator> recoTrackToTrackingParticleAssociatorToken_;

  const reco::RecoToSimCollection* r2s_;
  bool tracking_truth_available_;

  std::vector<edm::EDGetTokenT<edm::View<reco::Vertex> > > reco_vertex_view_tokens_;  // FIXME used ?
  std::vector<edm::InputTag> reco_vertex_collections_;

  std::vector<edm::EDGetTokenT<edm::PSimHitContainer> > hitCollectionTokens_;
  std::vector<std::string> hitCollections_;

  edm::EDGetTokenT<TrackingParticleCollection> trackingParticleCollectionToken_;
  edm::EDGetTokenT<TrackingVertexCollection> trackingVertexCollectionToken_;
  edm::EDGetTokenT<reco::SimToRecoCollection> simToRecoAssociationToken_;
  edm::EDGetTokenT<reco::RecoToSimCollection> recoToSimAssociationToken_;
  edm::EDGetTokenT<reco::VertexToTrackingVertexAssociator> vertexAssociatorToken_;
  edm::EDGetTokenT<LumiDetails> lumiDetailsToken_;
  edm::EDGetTokenT<LumiSummary> lumiSummaryToken_;
  edm::EDGetTokenT<LumiScalersCollection> lumiScalersToken_;
  edm::EDGetTokenT<LumiInfo> lumiInfoToken_;

  edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster> > pixelClusters_;
  bool l1garbled_;

  static const int iX = 0;
  static const int iY = 1;
  static const int iZ = 2;

  static const bool dcEffByCore_ = false;
  std::vector<double> lumibinmax_;

  double etaMaxVisible_;
  double ptMinVisible_;
  int numTrkHitsVisible_;
  double selNdof_;
  double selNdofWithBS_;
  double selNdofNoBS_;
  double ndof0trk_;

  static const int nzbins_ = 40;
  double zbinmax_;
  double dxb_;
  double dyb_;
  double dzb_;
  double wxb_;
  double wyb_;
  double wzb_;

  std::map<unsigned int, unsigned int> lsindex_;
  unsigned int nls_glb_;
  unsigned int lsglb_;
  std::map<unsigned int, std::vector<float> > lumidata_;
  std::map<unsigned int, std::vector<float> > bxlumidata_;  // obsolete
  std::map<unsigned int, std::vector<float> > bxpu_;
  /* EXTRA
  TH2S* hxpu_;
  TH1F* hapu_;
  TFile* tfx_;
  std::string era_;
  int bxdatarun_;
  */
  double lumiPU_;
  double avglumiPU_;
  bool firstBXinTrain_;
  double simPU_;

  std::vector<std::string> reports_;

  std::vector<std::string> trkdzbin_{"dz000-100","dz100-200","dz200-500","dz500-1000","dzgt1000"};

  std::vector<double> dzbins_{0.0050,
                              0.0150,
                              0.0250,
                              0.0350,
                              0.0450,
                              0.0550,
                              0.0650,
                              0.0750,
                              0.0850,
                              0.1000,
                              0.2000,
                              0.3000,
                              0.4000,
                              0.5000,
                              1.000,
                              2.000,
                              3.000};


};
