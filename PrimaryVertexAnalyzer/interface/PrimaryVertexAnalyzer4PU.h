// -*- C++ -*-
//
// Package:    MyPrimaryVertexAnalyzer4PU
// Class:      MyPrimaryVertexAnalyzer4PU
//
/**\class PrimaryVertexAnalyzer4PU PrimaryVertexAnalyzer4PU.cc Validation/RecoVertex/src/PrimaryVertexAnalyzer4PU.cc

 Description: primary vertex analyzer for events with pile-up

*/
//
// Original Author:  Wolfram Erdmann

//#define NOT_MATCHED 666666
#define NOT_MATCHED_VTX_REC 666665
#define NOT_MATCHED_VTX_SIM 666667
#define NOT_MATCHED_TK_SIM 666668 
#define NOT_MATCHED_TK_REC 666664 
#define NOT_ASSIGNED 666669
#define NO_RECVTX 999999
#define NO_KEY 999998
#define NO_INDEX 999997
#define FROM_TRACKING_TRUTH 999996
#define FROM_PU_SUMMARY 999995
#define FROM_WHATEVER 999994
#define FROM_HEPMC 999993
#define NOT_FILLED 999992

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
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
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

// MINIAOD
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

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
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TObjString.h>
#include <TString.h>

// timing
#include <chrono>

#include "RecoVertex/PrimaryVertexProducer/interface/TrackFilterForPVFinding.h"

typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag> GenEventVertex;

typedef reco::Vertex::trackRef_iterator trackit_t;

// class declaration
class PrimaryVertexAnalyzer4PU : public edm::one::EDAnalyzer<> {
  typedef math::XYZTLorentzVector LorentzVector;
  typedef reco::TrackBase::ParameterVector ParameterVector;


  
public:
  class SimPart {
  public:
    SimPart(){
      // only for backwards compatibility, eventually get rid of this constructor
      simpvidx=0;
      rec=NOT_MATCHED_TK_REC;
      type = 0;
      pdg = 0;
      xvtx = 0;
      yvtx = 0;
      zvtx = 0;
      tvtx = 0;
      charge = 0;
      pt = 0;
      phi = 0;
      eta = 0;
      ldec = 0;
      ddcap = 0;
      zdcap = 0;
      par = ParameterVector(0,0,0,0,0);
    };
    SimPart(int evt0, int type0, int pdgCode, double x1, double y1, double z1, double t1, double px1, double py1, double pz1, double BfieldT){
      // x0,y0,z0 should be the production point relative to the beam spot coordinates
      simpvidx = evt0;
      type = type0;
      pdg = pdgCode;
      xvtx = x1;
      yvtx = y1;
      zvtx = z1;
      tvtx = t1;
      set_charge(pdg);
      set_trkpar_c(charge, x1, y1, z1, px1, py1, pz1, BfieldT);
      rec = NOT_MATCHED_TK_REC;
      ldec = 0;
      ddcap = 0;
      zdcap = 0;
    };
    
    SimPart(unsigned int evt0, bool is_prompt, double x0, double y0, double z0, const pat::PackedGenParticle & cand, double BfieldT){
      // x0,y0,z0 is assumed to be the production point  (not the track reference point)
      simpvidx = evt0;
      if (is_prompt){
	type = 0;
      }else{
	type = 1;
      }
      pdg = cand.pdgId();
      xvtx = x0;
      yvtx = y0;
      zvtx = z0;
      tvtx = 0;
      ldec = 0;
      ddcap = 0;
      zdcap = 0;
      set_charge(pdg);
      if(charge != cand.charge()) {
	std::cout << "charge mismatch for pdg " << cand.pdgId() <<  " " << charge <<  "<> " << cand.charge() << std::endl;
      }
      set_trkpar_p(charge, x0, y0, z0,  cand.pt(), cand.phi(), cand.theta(), BfieldT);
      pt = cand.pt();
      phi = cand.phi();
      eta = cand.eta();
      rec = NOT_MATCHED_TK_REC;
    };

    
    void set_charge(int pdgCode){
      charge = 0;
      if ((pdgCode == 11) || (pdgCode == 13) || (pdgCode == 15) || (pdgCode == -211) || (pdgCode == -2212) ||
	  (pdgCode == -321) || (pdgCode == -3222) || (pdgCode == 3312) || (pdgCode == 3112)) {
	charge = -1;
      } else if ((pdgCode == -11) || (pdgCode == -13) || (pdgCode == -15) || (pdgCode == 211) || (pdgCode == 2212) ||
		 (pdgCode == 321) || (pdgCode == 3222) || (pdgCode == -3312) || (pdgCode == -3112)) {
	charge = 1;
      }else{
	std::cout << "SimPart: unknown pdg  " << pdgCode << std::endl;
	charge = 0;
      }
    };
    
    
    void set_trkpar_c(double Q, double x1, double y1,double z1, double px,double py,double pz, double BfieldT){
      pt = sqrt(px*px + py*py);
      double cosphi = px / pt;
      double sinphi = py / pt;
      double ptot = sqrt(px*px + py*py + pz*pz);

      double kappa = -Q * 0.002998 * BfieldT / pt;
      double D0 = x1 * sinphi - y1 * cosphi - 0.5 * kappa * (x1 * x1 + y1 * y1);
      double q = sqrt(1. - 2. * kappa * D0);
      double s0 = (x1 * cosphi + y1 * sinphi) / q;
      double s1;
      if (fabs(kappa * s0) > 0.001) {
	s1 = asin(kappa * s0) / kappa;
      } else {
	double ks02 = (kappa * s0) * (kappa * s0);
	s1 = s0 * (1. + ks02 / 6. + 3. / 40. * ks02 * ks02 + 5. / 112. * pow(ks02, 3));
      }
      par[reco::TrackBase::i_qoverp] = Q / ptot;
      par[reco::TrackBase::i_lambda] = M_PI / 2. - atan2(pt, pz);
      par[reco::TrackBase::i_phi] = phi - asin(kappa * s0);
      par[reco::TrackBase::i_dxy] = -2. * D0 / (1. + q);
      par[reco::TrackBase::i_dsz] = z1 * pt / ptot - s1 * pz / ptot;

      std::cout << std::endl;
      std::cout  << " full phi " << phi - asin(kappa * s0) << " linear " << phi << std::endl;
      std::cout  << " full dxy " << -2. * D0 / (1. + q) << " linear " <<   -(x1 * sinphi - y1 * cosphi) << std::endl;
      std::cout  << " full dszi " <<  z1 * pt / ptot - s1 * pz / ptot << " linear " <<  z1 * pt / ptot - (x1 * cosphi + y1 * sinphi) * pz / ptot << std::endl;
      /*
      // linear approximation
      par[reco::TrackBase::i_qoverp] = Q / ptot;
      par[reco::TrackBase::i_lambda] = M_PI / 2. - atan2(pt, pz);
      par[reco::TrackBase::i_phi] = phi;
      par[reco::TrackBase::i_dxy] = -(x1 * sinphi - y1 * cosphi);
      par[reco::TrackBase::i_dsz] = z1 * pt / ptot - (x1 * cosphi + y1 * sinphi) * pz / ptot;
      */

    }
    
    void set_trkpar_p(double Q, double x1, double y1,double z1, double pt1, double phi1, double theta1, double BfieldT){
      // assuming that phi1 is given at (x1,y1,z1)
      // note that the reference point for reco::Track is the point of closest approach to the beam,
      double kappa = -Q * 0.002998 * BfieldT / pt1;
      double D0 = x1 * sin(phi1) - y1 * cos(phi1) - 0.5 * kappa * (x1 * x1 + y1 * y1);
      double q = sqrt(1. - 2. * kappa * D0);
      double s0 = (x1 * cos(phi1) + y1 * sin(phi1)) / q;
      double s1;
      if (fabs(kappa * s0) > 0.001) {
	s1 = asin(kappa * s0) / kappa;
      } else {
	double ks02 = (kappa * s0) * (kappa * s0);
	s1 = s0 * (1. + ks02 / 6. + 3. / 40. * ks02 * ks02 + 5. / 112. * pow(ks02, 3));
      }
      double ptot = pt1 / sin(theta1);
      par[reco::TrackBase::i_qoverp] = Q / ptot;
      par[reco::TrackBase::i_lambda] = M_PI / 2. - theta1;
      par[reco::TrackBase::i_phi] = phi1 - asin(kappa * s0);
      par[reco::TrackBase::i_dxy] = -2. * D0 / (1. + q);
      par[reco::TrackBase::i_dsz] = z1 * sin(theta1) - s1 * cos(theta1);
    }

    /* put this into a separate function, doesn't seem to be used, anyway
      // now get zpca  (get perigee wrt beam)
      double x10 = x0 - vertexBeamSpot_.x(z0);
      double y10 = y0 - vertexBeamSpot_.y(z0);
      
      D0 = x10 * sin(p.phi()) - y10 * cos(p.phi()) - 0.5 * kappa * (x10 * x10 + y10 * y10);
      q = sqrt(1. - 2. * kappa * D0);
      s0 = (x1 * cos(p.phi()) + y1 * sin(p.phi())) / q;
      if (fabs(kappa * s0) > 0.001) {
	s1 = asin(kappa * s0) / kappa;
      } else {
	double ks02 = (kappa * s0) * (kappa * s0);
	s1 = s0 * (1. + ks02 / 6. + 3. / 40. * ks02 * ks02 + 5. / 112. * pow(ks02, 3));
      }
      ddcap = -2. * D0 / (1. + q);
      zdcap = z0 - s1 / tan(p.theta());
          sp.zvtx = z0;
          sp.xvtx = x0;
          sp.yvtx = y0;

          sp.phi = p.phi();
          sp.eta = p.eta();
      
    }
    */
    
    ParameterVector par;
    int type;      // 0 = primary
    double zvtx;  // z of the production vertex (not necessarily the primary)
    double xvtx;  // x of the production vertex
    double yvtx;  // y of the production vertex
    double tvtx;  // t of the production vertex (not available for miniaod)
    double charge;
    int pdg;      // particle pdg id
    unsigned int rec;       // index of the matched MTrack
    unsigned int simpvidx;  // index of the primary vertex
    /* defined but never used */
    double zdcap;  // z@dca' (closest approach to the beam)
    double ddcap;
    double ldec;   // distance of the production vertex from the primary vertex
    // convenience, may or may not be filled
    double eta;
    double phi;
    double pt; 
  };

  

  // forward declarations
  class MTrack;
  class MVertex;
  class Tracks;

  // a class holding simulated events
  class SimEvent {
  public:
    SimEvent(unsigned int idx) {
      index = idx;
      type = 0;
      x = 0;          // simulated vertex position
      y = 0;
      z = -99;
      t = 0;
      dzmin = 999.;
      nChTP = 0;     // number of charged particles
      ptvis = 0;
      pxvis = 0;
      pyvis = 0;
      pt_hat = 0;
      nGenTrk = 0;
      parts.clear();
      sumpt2 = 0;
      sumpt = 0;
      Tc = -1;
      dzmax = 0;
      dztrim = 0;
      chisq = 0;

      // track truth-matching
      trkidx.clear();         // indices of truth-matched tracks

      // MVertex matching
      clear_matching_info();

    };

    unsigned int index;      // =index in the SimEvent list
    EncodedEventId eventId;  // 
    int type;
    // 0 = not filled,
    // 1 = full (e.g. from TrackingParticles),
    // 2 = partially filled (from PileUpSummary)

    double x, y, z, t;
    double dzmin;
    double xfit, yfit, zfit, tfit;
    int nChTP;
    double ptvis, pxvis, pyvis;
    int nGenTrk;
    double pt_hat;
    std::vector<SimPart> parts;           // list of simulated particles that might have tracks

    std::vector<const TrackingParticle*> tp;
    std::vector<edm::RefToBase<reco::Track> > trkref;  // hm, can we get rid of this?

    std::vector<MTrack> rtk;              // all selected MTracks matched to this simevent
    std::vector<MTrack> rtkprim;          // all selected MTracks matched to this simevent that come from within 5 microns of the simvertex
    std::vector<MTrack> rtkprimsel;       // all selected MTracks matched to this simevent that pass the ip/sigma(ip) < 4 cut (rec wrt beam) ?????
    std::vector<unsigned int> trkidx;        // list of MTrack indices, used anywhere?

    double Tc, chisq, dzmax, dztrim;        // filled by getSimEvents via "getTc"
    double sumpt2, sumpt;                   // filled by getSimEvents
    
    
    unsigned int nwosmatch;  // number of recvertices dominated by this simevt (by wos)
    unsigned int nwntmatch;  // number of recvertices dominated by this simevt  (by nt)
    std::vector<unsigned int> wos_dominated_recv;  // list of dominated recv (by wos, size==nwosmatch)

    // number of tracks in recvtx with index i (VertexCollection->at(i))
    std::map<unsigned int, int> ntInRecVi;

    std::map<unsigned int, double> wnt;     // weighted number of tracks in recvtx (by index)
    std::map<unsigned int, double> wos;     // sum of wos in recvtx (by index) // oops -> this was int before 04-22
    double sumwos;                          // sum of wos in any recvtx
    double sumwnt;                          // sum of weighted tracks



    unsigned int rec;      // index of matched rec vertex or NOT_MATCHED_VTX_REC
    std::map<std::string, unsigned int> _recv;  // per collection, better switch to this

    unsigned int matchQuality;
    double zrec;
    double ndof;

    bool is_matched() const { return (rec != NOT_MATCHED_VTX_REC); }
    bool is_signal() const { return (index == 0); }


    void addTrack(unsigned int irecv, double twos, double twnt) {
      sumwnt += twnt;
      if (wnt.find(irecv) == wnt.end()) {
        wnt[irecv] = twnt;
      } else {
        wnt[irecv] += twnt;
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

    bool hasTrack(const MTrack * t) {
      for (unsigned int i = 0; i < rtk.size(); i++) {
        if (t->key() == rtk[i].key()) {
          return true;
        }
      }
      return false;
    }

    int countVertexTracks(const reco::Vertex& v, double min_weight = 0.5) {// see also min_trk_in_vtx_weight_
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

    int countVertexTracks(const MVertex & v, double min_weight = 0.5) {
      return countVertexTracks(v.recovertex(), min_weight);
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
      unsigned int vtx = NOT_MATCHED_VTX_REC;
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
      unsigned int vtx = NOT_MATCHED_VTX_REC;
      for (auto it : wnt) {
        if (it.second > maxnt) {
          maxnt = it.second;
          vtx = it.first;
        }
      }
      return vtx;  // this is often, but not necessarily the matched vertex
    }

    bool was_found(double min_ndof = 0.) { return (rec != NOT_MATCHED_VTX_REC) && (ndof > min_ndof); } //hmmm, never used?

    void clear_matching_info(){
      rec = NOT_MATCHED_VTX_REC;
      ndof = 0;           // ndof of the matched rec vertex (deprecated, not always filled)
      zrec = 1000.;       // ditto
      matchQuality = 0;

      wos.clear();
      sumwos = 0.;
      nwosmatch = 0;              // FIXME, this should be the size of  wos_dominated_recv, make it a function
      wos_dominated_recv.clear();

      wnt.clear();
      sumwnt = 0.;
      nwntmatch = 0;
    }
    
    bool operator==(const SimEvent& rh) const { return (index == rh.index); }


    bool is_visible(){
      if(type == NOT_FILLED){
	return false;
      }else if(type == FROM_TRACKING_TRUTH){
	return nChTP > 1;
      }else if(type == FROM_PU_SUMMARY){
	return pt_hat > 0.7; // tune to get similary visible vertex counts as with tracking particles, 0.5 is too low, 0.8 is slightly too high
      }else if(type == FROM_HEPMC){
	return nGenTrk > 2;
      }else if(type == FROM_WHATEVER){
	return true;
      }
      std::cout << "SimEvent.is_visible() : undefined type " << type << " encountered" << std::endl;
      return false;
    }

  };
  // SimEvent




  class MVertex{
    // container for a reco::Vertex  + extras
    // unified access for various formats including miniaod
  public:
    void init(){
      // reco part
      recvtx = NULL;
      tracks.clear();
      weights.clear();
      _index = NO_INDEX;

      // truth-matching part (formerly know as RSmatch)
      clear_matching_info();
    }

    MVertex(const reco::Vertex* v){
      init();
      recvtx = v;
      _index = NO_INDEX;
    };

    MVertex(const reco::Vertex* v, unsigned int idx){
      init();
      recvtx = v;      
      _index = idx;
    };

    double sumw() const {
      double s = 0;
      for(auto tk : tracks){
	s += trackWeight(tk);
      }
      return s;
    };

    double sumpt() const {
      double s = 0;
      for(auto tk : tracks){
	s += tk->pt();
      }
      return s;
    };

    double sumpt2() const {
      double s = 0;
      for(auto tk : tracks){
	s += pow(tk->pt(), 2);
      }
      return s;
    };

    double ptmax2() const {
      double ptmax_1 = 0;
      double ptmax_2 = 0;
      for(auto tk : tracks){
	if (tk->pt() > ptmax_1){
	  ptmax_2 = ptmax_1;
	  ptmax_1 = tk->pt();
	}else if (tk->pt() > ptmax_2){
	  ptmax_2 = tk->pt();
	}
      }
      return ptmax_2;
    };

    const double sumabspt() const {
      double aptsum = 0.;
      //double waptsum = 0.;
      //double wsum = 0.;
      for(auto tk : tracks){
	aptsum += tk->pt();
	//double w = trackWeight(tk);
	//wsum += w;
	//waptsum += tk->pt * w;
      }
      return aptsum;
    }

    int index() const{ return _index;};

    const double c2xy(const reco::BeamSpot & beam) const {
      double vxx = recvtx->covariance(iX, iX) + pow(beam.BeamWidthX(), 2);
      double vyy = recvtx->covariance(iY, iY) + pow(beam.BeamWidthY(), 2);

      double vxy = recvtx->covariance(iX, iY);
      double dx = recvtx->x() - beam.x(recvtx->z());
      double dy = recvtx->y() - beam.y(recvtx->z());
      double D = vxx * vyy - vxy * vxy;
      return pow(dx, 2) * vyy / D + pow(dy, 2) * vxx / D - 2 * dx * dy * vxy / D;
    };
    
    const double pxy(const reco::BeamSpot & beam) const {
      return TMath::Prob(c2xy(beam), 2);
    };
    
    double r(const reco::BeamSpot & beam) const {
      double z = recvtx->z();
      double dx = recvtx->x() - beam.x(z);
      double dy = recvtx->y() - beam.y(z);
      return sqrt(dx * dx + dy * dy);
    };
    
    // some pass-through functions to the underlying reco::Vertex for convenience
    double ndof() const { return recvtx->ndof();};
    bool isRecoFake() const { return recvtx->isFake();};  //bool isFake() const { return (chi2_ == 0 && ndof_ == 0 && tracks_.empty()); }
    double x() const { return recvtx->position().x();};
    double y() const { return recvtx->position().y();};
    double z() const { return recvtx->position().z();};
    double t() const { return recvtx->t();};
    double xError() const { return recvtx->xError();};
    double yError() const { return recvtx->yError();};
    double zError() const { return recvtx->zError();};
    double tError() const { return recvtx->tError();};
    double xyCorrelation() const { return recvtx->covariance(0,1) / sqrt(recvtx->covariance(0,0) * recvtx->covariance(1,1));};
    double covariance(const int i, const int j) const { return recvtx->covariance(i,j);};
    double chi2() const {return recvtx->chi2();};

    const reco::Vertex & recovertex() const{  // get a reference of the reco::Vertex
      if(recvtx == NULL){
	std::cout << "MVertex:recovertex()  invalid reco vertex requested" << std::endl;
	throw std::runtime_error("MVertex:recovertex()  invalid reco vertex requested");
      }else{
	return *recvtx;
      }
    }

    void add_track(MTrack *tk, float weight){
      tracks.push_back(tk);
      weights[tk->key()] = weight;
    }

    bool has_track_key(uint const tkey) const { return weights.count(tkey); }

    unsigned int tracksSize() const {return tracks.size();};
    float trackWeight (MTrack const * tk)const { return weights.at(tk->key());};
    float trackWeight (MTrack const & tk)const { return weights.at(tk.key());};

    unsigned int countTruthMatchedTracks(const double min_weight = 0) const{ // replace getTruthMatchedVertexTracks
      if(isRecoFake()) return 0;
      unsigned int ntk = 0;
      for(auto tk : tracks){
	if( (trackWeight(tk) > min_weight) && (tk->is_matched()) ){ntk++;}
      }
      return ntk;
    }

    bool has_timing() const {return (t() != 0) && (tError()>0.) && (tError()<9.) ;};// FIXME improve


    const reco::Vertex * recvtx;   // pointer to the underlying reco::Vertex
    unsigned int _index;           // position in the reco::Vertex collection

    std::vector<MTrack *> tracks;    // these are pointers!, be sure the source track list doesn't move 
    std::map<unsigned int, float> weights;  // assignment weight from fit



    // truth-matching
    void add_truthmatched_track(unsigned int iev, double twos, double twt) { // formerly know as addTrack
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

    bool is_real() const{ return (matchQuality > 0) && (matchQuality < 99); }

    bool is_fake() const{ return (matchQuality <= 0) || (matchQuality >= 99); }  // not to be confused with isRecoFake

    bool is_signal() const{ return (sim == 0); }

    int split_from() const{
      if (is_real())
        return -1;
      if ((maxwos > 0) && (maxwos > 0.3 * sumwos))
        return wosmatch;
      return -1;
    }
    bool other_fake() const{ return (is_fake() & (split_from() < 0)); }

    void clear_matching_info(){
      sim = NOT_MATCHED_VTX_SIM;
      wos.clear();
      wnt.clear();
      wosmatch = NOT_MATCHED_VTX_SIM;
      wntmatch = NOT_MATCHED_VTX_SIM;
      sumwos = 0;
      sumwnt = 0;
      maxwos = 0.;
      maxwnt = 0;
      maxwosnt = 0;
      matchQuality = 0;
      sigwosfrac = 0.;
      sigwntfrac = 0.;
    }
    
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
    double sigwosfrac;                   // for simtrack matching : fraction of the signal wos
    double sigwntfrac;                   // for simtrack matching : fraction of the weighted signal tracks
  };



  class MVertexCollection{
    // container for all Vertexs of a collection, 
  public:

    MVertexCollection(){
      collection_label = "";
      vtxs.clear();
    };

    MVertexCollection(const std::string & label, const reco::VertexCollection* recVtxs, Tracks & tracks){
      collection_label = label;
      vtxs.clear();
      if(recVtxs == NULL){
	std::cout << "invalid collection " << label << std::endl;
      }else{
	
	for(unsigned iv = 0; iv < recVtxs->size(); iv++){
	  vtxs.push_back(MVertex( &recVtxs->at(iv), iv));
	}

	// vertex track lists, should also work for miniaod
	for(unsigned int i = 0; i < tracks.size(); i++){
	  auto iv = tracks(i).get_recv(collection_label);
	  if (iv != NO_RECVTX){
	    assert(iv <= vtxs.size());
	    //std::cout << " adding track " << i << " to vertex " << iv << "  with weight " <<  tracks[i].get_weight(collection_label) << std::endl;
	    vtxs[iv].add_track(&tracks[i], tracks(i).get_weight(collection_label)); 
	  }
	}
      }
    }

    std::string label() const{ return collection_label;}
    MVertex& at(const unsigned int i){ return vtxs.at(i); }
    MVertex& operator[] (const unsigned int i){ return vtxs.at(i); }
    const MVertex& operator() (const unsigned int i){ return vtxs.at(i); }const 
    unsigned int size() const{return vtxs.size();}
    unsigned int nvtx() const{  // correct for a possible unphysical dummy vertex at the end of the list
      if(vtxs.size()==0){
	return 0;
      }else{
	if( vtxs.at(vtxs.size()-1).isRecoFake() ){
	  return vtxs.size()-1;
        }else{
          return vtxs.size();
       }
      }
    }
    std::vector<MVertex>::iterator begin(){return vtxs.begin();}
    std::vector<MVertex>::iterator end(){return vtxs.end();}
   
    std::string collection_label;
    std::vector<MVertex> vtxs;
  };







  /* MTrack is a helper class for collecting track related information and provide a common interface for reco::Track and minoad candidate based tracks */
  class MTrack {
  public:

    // fill some default values
    void init(){
      _from_collection = false;  
      _selected = false;
      _has_transientTrack = true;
      _lost_inner_hits = 0;
      _has_validHitInFirstPixelBarrelLayer = false;
      _has_hitPattern = 0;
      _has_timing = false;
      _is_highPurity = false;
      _algo = reco::TrackBase::TrackAlgorithm::undefAlgorithm;
      _ptError = 1.;         

      _t = 0;
      _dt = 1e10;
      _timeQuality = -1.;  // will stay -1. for no timing info
      _MTD_pathlength = 0;
      _MTD_time = -1.;
      _MTD_timeerror = 1e10;
      _MTD_momentum = 0;

      _recv.clear();    // for reco this is filled in get_reco_and_transient_tracks, for miniaod in the MVertex constructor
      _weight.clear();

      // filled later by getSimEvents
      _matched = NOT_MATCHED_TK_SIM;
      _simEvt = NULL;
      _zsim = 0;
      _tsim = 0;
      _is_primary = false;
    }


    // constructor for tracks from a collection
    MTrack(unsigned int idx,
              const reco::Track* trk,
       	      const reco::TransientTrack ttr,  // passed by value, makes a local copy, doesn't it?
              unsigned int trkkey,
              bool f4D) {
      init();
      _from_collection = true;
      _index = idx;            // index into the original collection
      _trk = trk;              // pointer into the collection  (only guaranteed to be valid when _from_collection == true)
      _transientTrack = ttr;   // transient track
      _has_transientTrack = true; // don't use the pointer to the track, it may become invalid when the MTrack is moved around (vectors do that)
      _key = trkkey;           // this really seems to always be the same as the index, is it guaranteed?
      _has_hitPattern = true; // well, if this is RECO
      _z = ( _transientTrack.stateAtBeamLine().trackStateAtPCA()).position().z();
      _dz = trk->dzError();
      _pt = trk->pt();
      _eta = trk->eta();
      _phi = (_transientTrack.stateAtBeamLine().trackStateAtPCA()).momentum().phi();
      _theta = (_transientTrack.stateAtBeamLine().trackStateAtPCA()).momentum().theta();
      Measurement1D atIP = _transientTrack.stateAtBeamLine().transverseImpactParameter();  // error contains beamspot
      _ip = atIP.value();
      _dip = atIP.error();
      _ptError = trk->ptError();
      _charge = trk->charge();
      _normalizedChi2 = _transientTrack.normalizedChi2();

      _lost_inner_hits = trk->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS); 
      _is_highPurity = trk->quality(reco::TrackBase::highPurity);
      _has_validHitInFirstPixelBarrelLayer = trk->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::PixelBarrel, 1);
  
      if (f4D) {
        _t = _transientTrack.timeExt();
        _dt = _transientTrack.dtErrorExt();
        if ((_dt > 0) && (_dt < 0.5)){// ?? && (!((_t==0) && (std::abs(_dt) <-0.35)<1e-5)))) {
          _has_timing = true;
	  _timeQuality = 1.; // temporary, overridden externally for now in get_reco_and_transient_tracks
        }
      }
    }



    // constructor for tracks from miniaod packed candidates (miniaod)
    MTrack(unsigned int idx,
	   pat::PackedCandidate cand,
	   unsigned int unique_id,
	   edm::ESHandle<TransientTrackBuilder> & builder,
	   const reco::BeamSpot & beamspot,
	   const std::string & vertexcollection_label,
	   bool f4D){

      init();
      _from_collection = false;
      _index = idx;
      _key = unique_id;
      _has_hitPattern = false;
      _normalizedChi2 = -1;
      _charge = cand.charge(); // hmm, is the mc truth based?
      
      // map to the (reco) primary vertex,
      // the inverse mapping is created later with the MVertexCollection
      unsigned int pv = NO_RECVTX;
      float weight = 0;
      if (cand.pvAssociationQuality() == pat::PackedCandidate::UsedInFitTight){
	pv = cand.vertexRef().key();
	weight = 0.8;  // > 0.5
      }else if (cand.pvAssociationQuality() == pat::PackedCandidate::UsedInFitLoose){
	pv = cand.vertexRef().key();
	weight = 0.1; // < 0.5
      }
      if ((_recv.find(vertexcollection_label) == _recv.end() || _weight[vertexcollection_label] < weight) ){
	_recv[vertexcollection_label] = pv;
	_weight[vertexcollection_label] = weight;
      }
      
      // timing
      if (f4D) {
        _t = cand.time();
        _dt = cand.timeError();
        if ((_dt > 0) && (_dt < 0.5) && (!((_t==0) && (std::abs(_dt-0.35)<1e-5)))) {
          _has_timing = true;
	  _timeQuality = 1.; // not available in miniod
        }
      }

      // fill mor info, depending or whether or not we have trackDetails
      // what is trk = cand.bestTrack() ??
      if(cand.hasTrackDetails()){
	
	// we have a pseudo track and a transient track, fill as ususal
	_pseudoTrack = cand.pseudoTrack();
	_trk = &_pseudoTrack;  // invalid when copies are made, don't use unless you kow what you are doing
	_transientTrack = (*builder).build( _pseudoTrack );
	_transientTrack.setBeamSpot(beamspot);
	_has_transientTrack = true;

	_normalizedChi2 = _transientTrack.normalizedChi2();
	_z = (_transientTrack.stateAtBeamLine().trackStateAtPCA()).position().z();
	_dz = _trk->dzError();
	_pt = _trk->pt();
	_eta = _trk->eta();
	_phi = (_transientTrack.stateAtBeamLine().trackStateAtPCA()).momentum().phi();
	_theta = (_transientTrack.stateAtBeamLine().trackStateAtPCA()).momentum().theta();
	Measurement1D atIP = _transientTrack.stateAtBeamLine().transverseImpactParameter();
	_ip = atIP.value();
	_dip = atIP.error();
	// limited hit info
	_has_validHitInFirstPixelBarrelLayer = false;
	if (cand.lostInnerHits() == pat::PackedCandidate::validHitInFirstPixelBarrelLayer){
	  _has_validHitInFirstPixelBarrelLayer = true;
	  _lost_inner_hits = 0; //?
	}else if (cand.lostInnerHits() == pat::PackedCandidate::noLostInnerHits){
	  _lost_inner_hits = 0;
	}else if (cand.lostInnerHits() == pat::PackedCandidate::oneLostInnerHit){
	  _lost_inner_hits =1 ;
	}else if(cand.lostInnerHits() == pat::PackedCandidate::moreLostInnerHits){
	  _lost_inner_hits =2 ;
	}
	//cand.pixelLayersWithMeasurement();
	//cand.trackerLayersWithMeasurement();
	_is_highPurity = cand.trackHighPurity();
	
      }else{
	
	// candidate without details
	//
	//auto pv = cand.vertex(); ??
	_trk = NULL;
	_z = cand.vz() + cand.dzAssociatedPV();
	_pt = cand.pt();
	_eta = cand.eta();
	_phi = cand.phi();
	_theta = cand.theta();
	//const auto beam = Point(vertexBeamSpot_.x(z), vertexBeamSpot_.y(z), z);
	//ip = cand.dxy(beam);
	_ip = 0.;
	//?? only if we have a covariance matrix?
	_dz = 100.; //cand.dzError();
	_dip= 100.; //cand.dxyError(); 
	_is_highPurity = cand.trackHighPurity();
	_has_transientTrack = false;

      }

    }

    
    // trivial accessor methods
    double t () const {return _t;};
    double dt () const {return _dt;};
    double z () const {return _z;};
    double dz () const {return _dz;}; // obsolete
    double dzError () const {return _dz;};
    double pt () const {return _pt;}
    double px () const {return _pt * cos(_phi);}
    double py () const {return _pt * sin(_phi);}
    double eta () const {return _eta;}
    double phi () const {return _phi;}
    double ip () const {return _ip;}
    double dip () const {return _dip;}
    double theta () const {return _theta;}
    double normalizedChi2() const{return _normalizedChi2;}
    double ptError () const {return _ptError;}
    double charge () const {return _charge;}
    double timeQuality () const {return _timeQuality;};
    double MTD_pathlength () const {return _MTD_pathlength;};
    double MTD_time () const {return _MTD_time;};
    double MTD_timeerror () const {return _MTD_timeerror;};
    double MTD_momentum () const {return _MTD_momentum;};
    reco::TrackBase::TrackAlgorithm algo () const {return _algo;}
    bool has_track () const { return !(_trk == NULL); }
    bool has_transienttrack () const { return _has_transientTrack; }
    const reco::TransientTrack & transientTrack() { return _transientTrack;}
    const reco::Track & track() const { if (_from_collection){return *_trk;}else{return _pseudoTrack;}};
    bool has_hitPattern () const {return _has_hitPattern;}
    bool selected () const {return _selected;}
    bool has_timing () const {return _has_timing;}
    bool is_highPurity () const {return _is_highPurity;}
    bool has_validHitInFirstPixelBarrelLayer () const {return _has_validHitInFirstPixelBarrelLayer;}
    unsigned int lost_inner_hits () const {return _lost_inner_hits;}
    unsigned int key () const {return _key;}
    unsigned int index () const {return _index;}
    

    // returns the index of the associated MVertex in it's collection
    unsigned int get_recv(const std::string& vtxcollection) const  {
      if (_recv.find(vtxcollection) == _recv.end()) {
        return NO_RECVTX;
      } else {
        return _recv.at(vtxcollection);
      }
    }

     // returns the weight of the associated MVertex in it's collection
    double get_weight(const std::string& vtxcollection) const  {
      if (_recv.find(vtxcollection) == _recv.end()) {
        return 0.;
      } else {
        return _weight.at(vtxcollection);
      }
    }

   // the following fields require MC truth
    bool is_matched() const {return _matched != NOT_MATCHED_TK_SIM;}
    bool is_primary () const {return _is_primary;}
    double zsim () const {return _zsim;}
    double tsim () const {return _tsim;}
    double tres () const {return _t - _tsim;}
    double zres () const {return _z - _zsim;}
    double tpull () const {return (_t - _tsim) / _dt;}
    double zpull () const {return (_z - _zsim) / _dz;}

    bool is_muon () const {
      return (is_matched() && (!_tpr.isNull()) && (abs(_tpr->pdgId()) == 13));
    }
    
    bool is_pion () const {
      return (is_matched() && (!_tpr.isNull()) && (abs(_tpr->pdgId()) == 211));
    }
    
    bool is_kaon () const {
      return (is_matched() && (!_tpr.isNull()) && (abs(_tpr->pdgId()) == 321));
    }
    bool is_proton () const {
      return (is_matched() && (!_tpr.isNull()) && (abs(_tpr->pdgId()) == 2212));
    }
    double get_particle_mass () const {
      if (is_matched() &&  (!_tpr.isNull())){
	return _tpr->mass();
      }else{
	return 0;
      }
    }
    double get_t_pid(double mass = 0) const{ //mass corrected reconstructed time at the beamline
      if (mass == 0 ){ // use the true mass (if known)
	mass = get_particle_mass();
      }
      if ((mass > 0) && (_MTD_pathlength > 0) ){
	double gammasq= 1. + _MTD_momentum * _MTD_momentum / (mass * mass);
	double v = 2.99792458e1 * std::sqrt(1. - 1. / gammasq);  // cm / ns
	return _MTD_time - _MTD_pathlength / v;
      }else{
	return -100.;
      }
    }

    bool is_looper () const {
      return false; // FIXME not implemented yet
    }

    // for z-sorting
    static bool lessz(const MTrack & tk1, const MTrack & tk2){
      return tk1._z < tk2._z;
    }

    reco::HitPattern hitPattern() const{
      if (_has_hitPattern){
	return track().hitPattern();
      }else{
	std::cout << "non-existing hitPattern requested" << std::endl;
	//reco::HitPattern dummy = reco::HitPattern();
	//return dummy;
	return reco::HitPattern();
      }
    }
      
    unsigned int lost() const{//hits
      if (_has_hitPattern){
	return track().lost();
      }else{
	return 0;
      }
    }

    bool _from_collection;               // true if the track comes from a reco trackcollection (then _trk points there)
    unsigned int _index;                 // index in the source collection (if applicable)
    unsigned int _key;                   // unique identifier  (usually the index)
  private:
    const reco::Track* _trk;             // points to the underlying reco::Track (or pseudoTrack), may be NULL, better not access directly
  public:
    reco::Track _pseudoTrack;            // used for miniaod only, _trk points to it in that case (and _from_collection is false)
    reco::TransientTrack _transientTrack;// local copy of the transient track (if it exists)
    bool _has_transientTrack;            // tells us if the transient track exists
    
    bool _selected;  // result of the track filter

    // convenience
    double _z;
    double _dz;
    double _pt, _eta, _phi, _ip, _dip, _theta;
    double _ptError;
    double _charge;
    bool _has_timing;
    double _t;
    double _dt;
    double _normalizedChi2;
    reco::TrackBase::TrackAlgorithm _algo;
    
    // the following variables are filled in PrimaryVertexAnalyzer4PU::get_reco_and_transient_tracks
    double _timeQuality;
    double _MTD_pathlength, _MTD_time, _MTD_timeerror, _MTD_momentum;
    double th[3];  // track time for particle hypotheses : 0=pion, 1=kaon, 2=proton
    
    bool _is_highPurity;
    // subset of hit info available in miniaod
    bool _has_hitPattern;
    unsigned int _lost_inner_hits;
    bool _has_validHitInFirstPixelBarrelLayer;
    
    // filled later if available : track appears int the track list of this reco vertex
    std::map<std::string, unsigned int> _recv;
    std::map<std::string, float> _weight;       // vertex.trackWeight(track)
    
    // MC truth related
    unsigned int _matched;  // idx of the matched simEvent or NOT_MATCHED_TK_SIM
    SimEvent* _simEvt;
    //TrackingParticle * tp;
    TrackingParticleRef _tpr;
    // convenience
    double _zsim;
    double _tsim;
    bool _is_primary;
  };


  /* collect information on reco tracks in one place
     basically a vector of MTracks with some extra info and pointers
  */
  class Tracks {
  public:
    edm::Handle<edm::View<reco::Track> > trackCollectionH;
    std::vector<MTrack> recotracks;
    std::map<unsigned int, unsigned int> key2idx;  // get the recotrack from a reftobase key?

    Tracks() {
      clear();
    }

    //const MTrack& operator() (const unsigned int i) { return recotracks.at(i); }const 
    const MTrack& operator() (const unsigned int i) const { return recotracks.at(i); }
    MTrack& operator[](const unsigned int i) { return recotracks[i]; }

    edm::RefToBase<reco::Track> ref(unsigned int i) { return edm::RefToBase(trackCollectionH, i); }

    // for tracks from a collection : access by the key of the original track
    MTrack& from_key(const unsigned int key){
      auto index = index_from_key(key);
      assert(recotracks[index].key() == key);
      return recotracks[index];
    }

    unsigned int index_from_key(const unsigned int key){
      // make the map if needed
      if (key2idx.size() != recotracks.size()) {
	key2idx.clear();
        for (unsigned int i = 0; i < recotracks.size(); i++) {
          key2idx[recotracks[i].key()] = i;
        }
      }
      assert(recotracks[key2idx[key]].key() == key);
      return key2idx[key];
    }

    MTrack& from_ref(const edm::RefToBase<reco::Track> ref) { return from_key(ref.key()); }
    unsigned int  index_from_ref(const edm::RefToBase<reco::Track> ref) { return index_from_key(ref.key()); }

    unsigned int simevent_index_from_key(const unsigned int key) {
      MTrack& tk = from_key(key);
      if (tk.is_matched()){
	return tk._simEvt->index;
      }else{
	return NOT_MATCHED_TK_SIM;
      }
    }

    // for convenience behave like a vector of MTrack
    unsigned int size() const { return recotracks.size(); }
    void clear() { 
      recotracks.clear();
      key2idx.clear();
    }
    void reserve(unsigned int n) { 
      recotracks.reserve(n);
    }
    void push_back(MTrack t) { recotracks.push_back(t); }// careful, does not automatically keep other arrays synchronized
    std::vector<MTrack>::const_iterator begin() const{return recotracks.begin();}
    std::vector<MTrack>::const_iterator end()const{return recotracks.end();}
  };


  
  /* helper class for histogramming/counting keys*/
  class VertexCounter {
  public:
    VertexCounter() { n.clear(); }

    void count(int key) {
      if (key == NOT_MATCHED_TK_SIM) {
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
  void dumpEventSummary(std::vector<SimEvent>&, Tracks & );
  
  void analyzeTracksTP(Tracks& tracks, std::vector<SimEvent>& simEvt);
  void testTrackParameters(Tracks& tracks);

private:

  void printPVTrksZT(Tracks & tracks,
                     MVertexCollection & recVtxs,
                     std::vector<SimEvent>& simEvt);

  void reportEvent(const char*, bool dumpvertex = false);
  void reportVertex(const reco::Vertex&, const char*, bool);

  std::vector<unsigned int> supfT(std::vector<SimPart>& simtrks, const Tracks &);

  static bool match(const ParameterVector& a, const ParameterVector& b);
  std::vector<SimPart> getSimTrkParameters(edm::Handle<edm::SimTrackContainer>& simTrks,
                                           edm::Handle<edm::SimVertexContainer>& simVtcs,
                                           double simUnit = 1.0);
  std::vector<SimPart> getSimTrkParameters(const edm::Handle<reco::GenParticleCollection>);
  void getTc(const std::vector<MTrack>&, double&, double&, double&, double&);

  double vertex_pxy(const reco::Vertex&);
  double vertex_sumw(const reco::Vertex&);
  double vertex_aptsum(const reco::Vertex&);
  double vertex_r(const reco::Vertex&);
  double vertex_ptmax2(const reco::Vertex&);
  //double vertex_yum(const reco::Vertex&);
  double vertex_maxfrac(const reco::Vertex&);
  double vertex_sumpt2(const reco::Vertex&);
  double vertex_sumpt(const reco::Vertex&);

  class Vertex_time_result {
  public:
    unsigned int status; 
    double tvtx, tvtxError;
    unsigned int niteration;
    std::vector<float> tk_tweight;
    double Tc;
    Vertex_time_result(){
      status = 1;
      tvtx = 0;
      tvtxError = 0;
      niteration = 0;
      tk_tweight.clear();
      Tc = 0;
    }
    void success(const double t, const double tError, const unsigned int iterations, const double Tclast=0.){ 
      tvtx = t;
      tvtxError = tError;
      niteration = iterations;
      status = 0;
      Tc = Tclast;
    }
    const double t(){ return tvtx;};
    const double tError(){ return tvtxError;};
    const bool successful(){ return status==0;}
    const bool no_status(){ return status==1;}
    const bool no_tracks_with_timing(){ return status==2;}
    const bool not_converged(){ return status==3;}
  };

    
  static Vertex_time_result vertex_time_from_tracks(const reco::Vertex&, Tracks& tracks, double minquality, bool verbose=false);
  static Vertex_time_result vertex_time_from_tracks_pid(const reco::Vertex&, Tracks& tracks, double minquality, bool verbose=false);
  static Vertex_time_result vertex_time_from_tracks_pid_newton(const reco::Vertex&, Tracks& tracks, double minquality, bool verbose=false);

  PrimaryVertexAnalyzer4PU::Vertex_time_result vertex_time_from_tracks_analysis(PrimaryVertexAnalyzer4PU::Vertex_time_result time_method(const reco::Vertex&,
                                                       Tracks& ,
						       double ,
						        bool),
					const MVertex & v,
					Tracks& tracks,
					double minquality,
					const SimEvent &,
					const std::string label, 
					std::map<std::string, TH1*>& h,
					const std::string vtype,
					bool verbose=false);


  void multi_vertex_time_from_tracks_pid(const std::string label, Tracks& tracks, double minquality, bool verbose=false);
  void mass_constrained_multi_vertex_time_from_tracks_pid(const std::string label, Tracks& tracks, double minquality, bool verbose=false);


  bool select(const reco::Vertex&, const int level = 0);
  bool select(const MVertex&, const int level = 0);

  bool uv_analysis(const MVertex & vtx, const SimEvent &, double & du, double & dv, double & uError, double & vError, double &pol);

  void addn(std::map<std::string, TH1*>& h, TH1* hist) {
    // add a histogram in a subdirectory and use the subdirectory name in the map key
    auto key = gDirectory->GetName() + std::string("/") + hist->GetName();
    if (h.find(key) != h.end()){
      std::cout << "addn: Warning ! adding already existing histogram " << key << std::endl;
    }
    h[gDirectory->GetName() + std::string("/") + hist->GetName()] = hist;
    hist->StatOverflows(kTRUE);
    hist->SetDirectory(gDirectory);
  }

  void add(std::map<std::string, TH1*>& h, TH1* hist) {
    // add a histogram
    if (h.find(hist->GetName()) != h.end()){
      std::cout << "warning ! adding already existing histogram " << hist->GetName() << std::endl;
    }
    h[hist->GetName()] = hist;
    hist->StatOverflows(kTRUE);
    hist->SetDirectory(gDirectory);
  }
  
  void add(std::map<std::string, TH1*>& h, TH1* hist, const std::string& vtype) {
    // add a histogram
    if (h.find(hist->GetName()) != h.end()){
      std::cout << "warning ! adding already existing histogram " << hist->GetName() << std::endl;
    }
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
    if (h.find(name + "Signal") != h.end()){
      std::cout << "addSP: warning ! adding already existing histogram " << name << std::endl;
    }
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
  bool get_miniaod_tracks(const edm::EventSetup&, const edm::Event&, const std::string &, Tracks&);
  void add_timing_to_vertex_collection(const std::string & label, Tracks& tracks);
  void refit_recvertices_after_timing(Tracks&, double min_tk_vtx_weight=0.5);
  void fill_track_to_vertex_pointers(Tracks&);

  double muvtx(double z);

  bool isResonance(const HepMC::GenParticle* p);
  bool isFinalstateParticle(const HepMC::GenParticle* p);
  bool isCharged(const HepMC::GenParticle* p);
  
  void fillVertexHistos(std::map<std::string, TH1*>& h,
                        const std::string& vtype,
                        const MVertex & v,
                        Tracks& tracks,
                        const double deltaz = 0,
                        const bool verbose = false);

  void fillVertexHistosNoTracks(std::map<std::string, TH1*>& h,
                                const std::string& vtype,
                                const reco::Vertex* v = NULL,
				const unsigned int index = NO_INDEX, 
                                const double deltaz = 0,
                                const bool verbose = false);
  
  void fillRecoVertexHistos(std::map<std::string, TH1*>& h,
                        const std::string& vtype,
                        const reco::Vertex* v,
                        Tracks& tracks,
			const unsigned int index = NO_INDEX,
                        const double deltaz = 0,
                        const bool verbose = false);
  
  
  void fillVertexHistosMatched(std::map<std::string, TH1*>& h,
			       const std::string& vtype,
			       MVertex & v,
			       Tracks& tracks,
			       const std::vector<SimEvent>& simEvt,
			       const double deltaz = 0,
			       const bool verbose = false);
  
  void fillTrackHistos(std::map<std::string, TH1*>& h,
                       const std::string& ttype,
                       MTrack& tk,
                       const reco::Vertex* v = NULL);
  void fillTrackHistosMatched(std::map<std::string, TH1*>& h,
                       const std::string& ttype,
                       MTrack& tk);
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
  void printRecVtxs(const MVertexCollection& , std::string title = "Reconstructed Vertices");
  void printSimVtxs(const edm::Handle<edm::SimVertexContainer> simVtxs);
  void printSimTrks(const edm::Handle<edm::SimTrackContainer> simVtrks);
  bool getPuInfo(const edm::Event& iEvent, PileupSummaryInfo& puInfo);


  Int_t getAssociatedRecoTrackIndex(const edm::Handle<reco::TrackCollection>& recTrks, TrackingParticleRef tpr);
  bool select_truthMatchedTrack(const edm::RefToBase<reco::Track>, TrackingParticleRef&)const;
  std::vector<edm::RefToBase<reco::Track> > getTruthMatchedVertexTracks_obsolete(const reco::Vertex&, double min_weight = 0.5) const;
  void printTruthMatchValues(edm::RefToBase<reco::Track> track);


  bool get_MC_truth(const edm::Event& iEvent,
                    Tracks& tracks,
                    bool bPuInfo,
                    PileupSummaryInfo& puInfo,
                    std::vector<SimEvent>& simEvt);

  void fill_simvtx_histos(std::vector<SimEvent>& simEvts);

  void getSimEvents_pu(PileupSummaryInfo& puInfo, std::vector<SimEvent>& simEvents);

  std::vector<PrimaryVertexAnalyzer4PU::SimEvent> getSimEvents_tp(edm::Handle<TrackingParticleCollection>,
                                                               Tracks& tracks);

  std::vector<PrimaryVertexAnalyzer4PU::SimEvent> getSimEvents_simtrks(
								    const edm::Handle<edm::SimTrackContainer> simTrks,
								    const edm::Handle<edm::SimVertexContainer> simVtxs,
								    Tracks& tracks);

  std::vector<PrimaryVertexAnalyzer4PU::SimEvent> getSimEvents_miniaod(
								       const edm::Event &, 
								       Tracks&);

void analyzeVertexCollectionZmatching(std::map<std::string, TH1*>& h,
				      MVertexCollection& vtxs,
				      std::vector<SimEvent>& simEvts,
				      const std::string message,
				      const double zwindow_sigmas = 3.0);
  
  void analyzeVertexRecoCPUTime(std::map<std::string, TH1*>& h,
                                const reco::VertexCollection* recVtxs,
                                const std::string message = "");
  void analyzeVertexCollectionRecoNoTracks(std::map<std::string, TH1*>& h,
                                           const reco::VertexCollection* recVtxs,
                                           const std::string message = "");

  void analyzeVertexCollectionReco(std::map<std::string, TH1*>& h,
                                   MVertexCollection& recVtxs,
                                   Tracks& tracks,  // do I even need this?
                                   const std::string message = "");

  void analyzeVertexCollectionSimTracks(std::map<std::string, TH1*>& h,
					MVertexCollection& vtxs,
					Tracks& tracks,
					std::vector<SimEvent>& ,
					const std::string message = "");

  void analyzeVertexCollectionSimPvNoSimTracks(std::map<std::string, TH1*>& h,
                                               MVertexCollection& vtxs,
                                               Tracks& tracks,
                                               std::vector<SimEvent>& ,
                                               const std::string message = "");


  void analyzeVertexMergeRateTP(std::map<std::string, TH1*>& h,
                                MVertexCollection& recVtxs,
				Tracks & tracks, 
				std::vector<SimEvent>& simEvt,
                                const std::string message = "");


  void analyzeVertexCollectionTP(std::map<std::string, TH1*>& h,
                                 MVertexCollection& recVtxs,
                                 Tracks& tracks,
                                 std::vector<SimEvent>& simEvt,
                                 const std::string message = "");

  void analyzeVertexCollectionPtvis(std::map<std::string, TH1*>& h,
				    MVertexCollection & vtxs,
				    Tracks& tracks,
				    std::vector<SimEvent>& simEvt,
				    const std::string message = "");

  void analyzeVertexTrackAssociation(std::map<std::string, TH1*>& h, MVertexCollection& vtxs, Tracks& tracks, std::vector<SimEvent> const& simEvt, float const npu);

  void analyzeVertexComposition(std::map<std::string, TH1*>& h,
                                   MVertex & v,
			           MVertexCollection & vtxs,
				   Tracks& tracks,
                                   std::vector<SimEvent>& simEvt,
				   float npu);
  void signalvtxmatch(MVertexCollection &, std::vector<SimEvent> &);

  void tpmatch(MVertexCollection& vtxs,
	       std::vector<SimEvent>& simEvt,
	       Tracks& tracks);
  
  void wos_match(MVertexCollection& recVtxs,
		 std::vector<SimEvent>& simEvt,
		 Tracks& tracks);

  std::string formatMatchList(const std::map<unsigned int, double>&, unsigned int nfield, bool sim);


  void printMatchingSummary(MVertexCollection& recVtxs,
                            std::vector<SimEvent>& simEvt,
                            const std::string message);

  void printEventSummary_tp(std::map<std::string, TH1*>& h,
			    MVertexCollection&,
			    Tracks& tracks,
			    std::vector<SimEvent>& simEvt,
			    const std::string message);

  void printEventSummary_notp(
                         MVertexCollection & vtxs,
                         Tracks& tracks,
                         std::vector<SimEvent>& simEvt,
                         const std::string message);

  reco::VertexCollection* vertexFilter(edm::Handle<reco::VertexCollection>, bool filter);

  void compareCollections(std::vector<SimEvent>& simEvt);

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
    // withBS   ndof = 2 * nt          0    2   4   6   8  10
    // with an adaptive fitter, replace nt by <w>
    // hence <w> = (ndof - ndof0trk_) / 2.
    if (vertexcollection.find("WithBS") != std::string::npos){
      selNdof_ = selNdofWithBS_;
      ndof0trk_ = 0;
    }else{
      selNdof_ = selNdofNoBS_;
      ndof0trk_ = -3.;
    }
  }

  //timers
  class PVTimer{
  public:
    PVTimer(){duration=0; counter=0;};
    double duration;
    unsigned int counter;
  };
  
  std::map<std::string, std::chrono::time_point<std::chrono::steady_clock> > timer_start_;
  std::map<std::string, double> timers_;
  std::map<std::string, PVTimer> pvtimers_;

  void inline timer_start(const std:: string label){
    timer_start_[label] = std::chrono::steady_clock::now();
  }
  void timer_stop(const std::string & label){
    auto stop = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - timer_start_[label]);
    timers_[label] += duration.count();
    auto & timer = pvtimers_[label];
    timer.duration += duration.count();
    timer.counter ++;
  }

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
  double zwindow_sigmas_;
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
  int matchsummaries_;



  /* global containers */
  std::vector<std::string> vertexCollectionLabels_;
  std::map<std::string, edm::EDGetTokenT<reco::VertexCollection> > vertexCollectionTokens_;
  std::map<std::string, std::map<std::string, TH1*> > histograms_;

  std::map<std::string, reco::VertexCollection*> recVtxs_;     // pointers to the reco::Vertex collections
  std::map<std::string, MVertexCollection > vertexes_;         // MVertices, parallel to recVtxs,  hold track- and mc truth-links if available

  std::map<std::string, TH1*> hsimPV;
  std::map<std::string, TH1*> hTrk;
  std::map<std::string, TH1*> hEvt;
  std::map<std::string, TH1*> hClu;

  unsigned int reset_nbin_;  //
  double reset_period_;
  unsigned int max_LS_;  //

  std::map<unsigned int, TrackingParticleRef> trkidx2tp_;  // reco::track index    --> tracking particle  (FIXME, seems to be unused or at least not filled)

  reco::BeamSpot vertexBeamSpot_;
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle_;
  edm::ESHandle<TransientTrackBuilder> theB_;

  std::map<std::string, std::pair<unsigned int,unsigned int> > counted_messages_;


  // control verbosity for some report types
  static constexpr bool dump_signal_vertex_not_tpmatched_ = true;
  static constexpr bool dump_fake_vertex_on_top_of_signal_ = false;
  static constexpr bool dump_big_fakes_ = false;
  
  bool f4D_;
  bool frefit_;
  bool fTrackTime_;

  bool RECO_;
  bool MINIAOD_;
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

  const reco::RecoToSimCollection* tp_r2s_;
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

  //Kyril's names, for miniaod
  edm::EDGetTokenT<edm::View<pat::PackedCandidate> > theTracksToken_;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate> > theLostTracksToken_;
  edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > theGenParticlesToken_;
  //edm::EDGetTokenT<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag> > theGenParticlesXyz0Token_;
  edm::EDGetTokenT<GenEventVertex> theGenParticlesXyz0Token_;
  edm::EDGetTokenT<float> theGenParticlesT0Token_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > thePrunedGenParticlesToken_;
  //edm::EDGetTokenT<edm::GenEventInfoProduct> GenEventInfoProductToken_;

  edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster> > pixelClusters_;
  bool l1garbled_;

  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideHowToGetDataFromES
  //edm::EDGetTokenT<edm::ESHandle<ParticleDataTable> > pdtToken_;
  const edm::ESGetToken<HepPDT::ParticleDataTable, edm::DefaultRecord> pdtToken_;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> transientTrackBuilderToken_;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> trackerTopologyToken_;

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
  double min_trk_in_vtx_weight_;
  double ndof0trk_;
  double trkcentraletamax_;
  double trkloptmax_;
  double trkhiptmin_;

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
