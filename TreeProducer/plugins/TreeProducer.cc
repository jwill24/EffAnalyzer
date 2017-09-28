// -*- C++ -*-
//
// Package:    EffAnalyzer/TreeProducer
// Class:      TreeProducer
// 
/**\class TreeProducer TreeProducer.cc EffAnalyzer/TreeProducer/plugins/TreeProducer.cc

*/
//
// Original Author:  Justin Williams
//         Created:  Thu, 24 Aug 2017 08:35:26 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/Common/interface/DetSetVector.h"

#include "DataFormats/CTPPSDetId/interface/TotemRPDetId.h"
#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"
#include "DataFormats/CTPPSReco/interface/TotemRPCluster.h"
#include "DataFormats/CTPPSReco/interface/TotemRPRecHit.h"
#include "DataFormats/CTPPSReco/interface/TotemRPUVPattern.h"
#include "DataFormats/CTPPSDigi/interface/CTPPSDiamondDigi.h"

#include "EffAnalyzer/TreeProducer/interface/HLTMatcher.h"
#include "EffAnalyzer/TreeProducer/interface/SelectionUtils.h"
#include "EffAnalyzer/TreeProducer/interface/XiInterpolator.h"
#include "EffAnalyzer/TreeProducer/interface/FillNumberLUTHandler.h"
#include "EffAnalyzer/TreeProducer/interface/AlignmentLUTHandler.h"

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

using namespace std;

//
// class declaration
//
#define MAX_PROTON_TRK 20
#define MAX_PHOTON 10

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class TreeProducer : public edm::one::EDAnalyzer<edm::one::SharedResources>  
{
public:
  explicit TreeProducer(const edm::ParameterSet&);
  ~TreeProducer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  
  // ----------member data ---------------------------
  
  void clearTree();
  void analyzeTriggers( const edm::Event&, const edm::EventSetup& );
  
  edm::EDGetTokenT<edm::DetSetVector<TotemRPUVPattern> > totemRPUVPatternsToken_;
  edm::EDGetTokenT<edm::DetSetVector<TotemRPCluster> > totemRPClustersToken_;
  edm::EDGetTokenT<edm::DetSetVector<TotemRPLocalTrack> > totemRPTracksToken_;
  edm::EDGetTokenT<edm::DetSetVector<TotemRPRecHit> > totemRPRecHitsToken_;  
  edm::EDGetTokenT<std::vector<reco::Photon> > photonsToken_;
  
  bool isData_;    
  double sqrtS_;
  
  std::string filename_;
  double maxGenLevelDR_;
  
  bool useXiInterp_;
  std::unique_ptr<ProtonUtils::XiInterpolator> xiInterp_;
  
  std::unique_ptr<CTPPSAlCa::FillNumberLUTHandler> fillLUTHandler_;
  std::unique_ptr<CTPPSAlCa::AlignmentLUTHandler> alignmentLUTHandler_;
  std::unique_ptr<edm::LumiReWeighting> lumiReweighter_;
    
  TFile* file_;
  TTree* tree_;
  
  
  //---tree components---  
  unsigned int fBX, fFill, fRun, fLumiSection;
  unsigned long long fEventNum;
  
  
  //Track variables
  unsigned int fProtonTrackNum;
  float fProtonTrackX[MAX_PROTON_TRK], fProtonTrackY[MAX_PROTON_TRK];
  float fProtonTrackXCorr[MAX_PROTON_TRK], fProtonTrackYCorr[MAX_PROTON_TRK];
  float fProtonTrackXi[MAX_PROTON_TRK], fProtonTrackXiError[MAX_PROTON_TRK];
  float fProtonTrackChi2[MAX_PROTON_TRK], fProtonTrackNormChi2[MAX_PROTON_TRK];
  unsigned int fProtonTrackSide[MAX_PROTON_TRK], fProtonTrackPot[MAX_PROTON_TRK], fProtonTrackLinkNF[MAX_PROTON_TRK];
  float fProtonTrackMinLinkDist[MAX_PROTON_TRK];
  
  //RecHit variables
  unsigned int fRHNum;
  double_t fRecHitPos[MAX_PROTON_TRK];
  double_t fRecHitSig[MAX_PROTON_TRK];

  //Cluster variables
  unsigned int fClusterNum;
  uint16_t fClusterStripBegin[MAX_PROTON_TRK];
  unsigned int fClusterStripEnd[MAX_PROTON_TRK];
  unsigned int fClusterNumStrips[MAX_PROTON_TRK];
  float fClusterCenterStripPosition[MAX_PROTON_TRK];
  
  //UV pattern variables
  unsigned int fUVPNum;  
  double  fUVPatternA[MAX_PROTON_TRK];  
  double fUVPatternB[MAX_PROTON_TRK];
  double fUVPatternW[MAX_PROTON_TRK];

  //Photon Variables
  float fPhotonPT[MAX_PROTON_TRK];
  int fPhotonNum;
  
  float fPileupWeight;
  
};

//
// constructors and destructor
//
TreeProducer::TreeProducer( const edm::ParameterSet& iConfig ) :
  totemRPUVPatternsToken_ ( consumes<edm::DetSetVector<TotemRPUVPattern> >     ( iConfig.getParameter<edm::InputTag>( "totemRPUVPatternsLabel") ) ),
  totemRPClustersToken_   ( consumes< edm::DetSetVector<TotemRPCluster> >      ( iConfig.getParameter<edm::InputTag> ( "totemRPClustersLabel") ) ),
  totemRPTracksToken_     ( consumes< edm::DetSetVector<TotemRPLocalTrack> >   ( iConfig.getParameter<edm::InputTag>( "totemRPTracksLabel") ) ),
  totemRPRecHitsToken_    ( consumes< edm::DetSetVector<TotemRPRecHit> >       ( iConfig.getParameter<edm::InputTag>( "totemRPRecHitsLabel") ) ),
  photonsToken_           (consumes<vector<reco::Photon> >                     ( iConfig.getParameter<edm::InputTag>("photonsLabel"))),
  sqrtS_              ( iConfig.getParameter<double>     ( "sqrtS" ) ),
  filename_           ( iConfig.getParameter<std::string>( "outputFilename" ) ),
  useXiInterp_        ( iConfig.getParameter<bool>       ( "useXiInterpolation" ) ),
  file_( 0 ), tree_( 0 )
  
  
{  
  xiInterp_ = std::make_unique<ProtonUtils::XiInterpolator>();
  if ( useXiInterp_ ) {
    xiInterp_->loadInterpolationGraphs( iConfig.getParameter<edm::FileInPath>( "xiInterpolationFile" ).fullPath().c_str() );
  }
  alignmentLUTHandler_ = std::make_unique<CTPPSAlCa::AlignmentLUTHandler>( iConfig.getParameter<edm::FileInPath>( "alignmentLUTFile" ).fullPath().c_str() );
  
  
  file_ = new TFile( filename_.c_str(), "recreate" );
  file_->cd();
  
  tree_ = new TTree( "ntp", "mb ntuple" );
  
}


TreeProducer::~TreeProducer()
  
{
  
  if ( file_ ) {
    file_->Write();
    file_->Close();
    delete file_;
  }
  if ( tree_ ) delete tree_;
  
  
}


//
// member functions
//


void
TreeProducer::clearTree()
{
  fBX = fRun = fLumiSection = fEventNum = 0;
  
  fProtonTrackNum = 0;
  fRHNum = 0;
  fClusterNum =0;
  fUVPNum =0;
  for ( unsigned int i=0; i<MAX_PROTON_TRK; i++ ) {
    fProtonTrackX[i] = fProtonTrackY[i] = -1.;
    fProtonTrackXCorr[i] = fProtonTrackYCorr[i] = -1.;
    fProtonTrackXi[i] = fProtonTrackXiError[i] = -1.;
    fProtonTrackChi2[i] = fProtonTrackNormChi2[i] = -1.;
    fProtonTrackSide[i] = 2; //invalid
    fProtonTrackPot[i] = 0;
    fProtonTrackLinkNF[i] = 999;
    fProtonTrackMinLinkDist[i] = -1.;
    
    fRecHitPos[i] = -1;
    fRecHitSig[i] = -1;
    fClusterStripBegin[i] = -1;
    fClusterStripEnd[i] = -1;
    fClusterNumStrips[i] = -1;
    fClusterCenterStripPosition[i] = -1;
    fUVPatternA[i] = -1;
    fUVPatternB[i] = -1;
    fUVPatternW[i] = -1;
  }
  

  fPhotonNum = 0;
  for ( unsigned int i=0; i<MAX_PHOTON; i++ ) {
    fPhotonPT[i] = -1.;
  }
  
  fPileupWeight = 1.;
}



// ------------ method called for each event  ------------


void
TreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

#include <iostream> // for debugging purposes

  clearTree();
  
  // Run and BX information
  fBX = iEvent.bunchCrossing();
  fRun = iEvent.id().run();
  fLumiSection = iEvent.luminosityBlock();
  fEventNum = iEvent.id().event();
  
  // get the fill number from the run id <-> fill number LUT
  if ( fillLUTHandler_ ) {
    fFill = fillLUTHandler_->getFillNumber( iEvent.id().run() );
  }
  else fFill = 1;
  
  
  //----- photon information as a test -----
  
  edm::Handle< std::vector<reco::Photon> > photons;
  iEvent.getByToken( photonsToken_, photons );
  
  fPhotonNum=0;
  for ( unsigned int i=0; i<photons->size() && fPhotonNum<MAX_PHOTON; i++ ) {
    const reco::Photon &photon = (*photons)[i];
    fPhotonPT[fPhotonNum] = photon.pt();
    fPhotonNum++;
  }

  //----- forward RP information -----


  // Clusters Information
  edm::Handle< edm::DetSetVector<TotemRPCluster> > digCluster;
  iEvent.getByToken( totemRPClustersToken_, digCluster );
  fClusterNum=0;
  for (  edm::DetSetVector<TotemRPCluster>::const_iterator it = digCluster->begin(); it != digCluster->end(); it++ ) {
    const TotemRPDetId detId(  it->detId() );
    for (edm::DetSet<TotemRPCluster>::const_iterator dit = it->begin(); dit != it->end(); ++dit) {
      fClusterStripBegin[fClusterNum] = dit->getStripBegin();
      fClusterStripEnd[fClusterNum] = dit->getStripEnd();
      fClusterNumStrips[fClusterNum] = dit->getStripEnd();
      fClusterCenterStripPosition[fClusterNum] = dit->getCenterStripPosition();
      fProtonTrackNum++;
    }
  }
  
    
  // UVPatterns Information
  edm::Handle< edm::DetSetVector<TotemRPUVPattern> > patterns;
  iEvent.getByToken( totemRPUVPatternsToken_, patterns );    
  fUVPNum=0;
  for ( edm::DetSetVector<TotemRPUVPattern>::const_iterator it = patterns->begin(); it != patterns->end(); it++ ) {
    const TotemRPDetId detId( it->detId() );
    for(edm::DetSet<TotemRPUVPattern>::const_iterator dit = it->begin(); dit != it->end(); ++dit) {
      fUVPatternA[fUVPNum] = dit->getA();
      fUVPatternB[fUVPNum] = dit->getB();
      fUVPatternW[fUVPNum] = dit->getW();
      fProtonTrackNum++;
    }
  }

  //RecHits Information
  edm::Handle< edm::DetSetVector<TotemRPRecHit> > rechits;
  iEvent.getByToken( totemRPRecHitsToken_, rechits );
  fRHNum=0;
  for ( edm::DetSetVector<TotemRPRecHit>::const_iterator it = rechits->begin(); it != rechits->end(); ++it ) {
    const TotemRPDetId detId( it->detId() );
    for(edm::DetSet<TotemRPRecHit>::const_iterator dit = it->begin(); dit != it->end(); ++dit ) {
      fRecHitPos[fRHNum] = dit->getPosition();
      fRecHitSig[fRHNum] = dit->getSigma();
      fProtonTrackNum++;
    }
  }
  
  
  // Tracks Information 
  
  edm::Handle< edm::DetSetVector<TotemRPLocalTrack> > rpLocalTracks;
  iEvent.getByToken( totemRPTracksToken_, rpLocalTracks );
  
  const CTPPSAlCa::RPAlignmentConstants align = alignmentLUTHandler_->getAlignmentConstants( fFill ); // fill-based alignment corrections
  xiInterp_->setAlignmentConstants( align );
  xiInterp_->setCalibrationConstants( fRun ); // run-based calibration parameters
  
  typedef std::pair<unsigned int, const TotemRPLocalTrack&> localtrack_t; // RP id -> local track object
  std::map<unsigned int,localtrack_t> map_near, map_far; // local track id in the tree -> localtrack_t object
  fProtonTrackNum = 0;
  for ( edm::DetSetVector<TotemRPLocalTrack>::const_iterator it=rpLocalTracks->begin(); it!=rpLocalTracks->end(); it++ ) {
    const TotemRPDetId detid( it->detId() );
    const unsigned short side = detid.arm(),
      pot = detid.rp();
    const CTPPSAlCa::RPAlignmentConstants::Quantities align_quant = align.quantities( it->detId() );
    
    for ( edm::DetSet<TotemRPLocalTrack>::const_iterator trk=it->begin(); trk!=it->end(); trk++ ) {
      if ( !trk->isValid() ) { continue; }
      
      fProtonTrackX[fProtonTrackNum] = trk->getX0() * 1.e-3; // store in m
      fProtonTrackY[fProtonTrackNum] = trk->getY0() * 1.e-3; // store in m
      fProtonTrackXCorr[fProtonTrackNum] = ( trk->getX0() + align_quant.x ) * 1.e-3; // store in m
      fProtonTrackYCorr[fProtonTrackNum] = ( trk->getY0() - align_quant.y ) * 1.e-3; // store in m
      fProtonTrackSide[fProtonTrackNum] = side; // 0 = left (45) ; 1 = right (56)
      fProtonTrackPot[fProtonTrackNum] = pot; // 2 = 210m ; 3 = 220m
      
      // x-to-xi interpolation
      float xi = -1., err_xi = -1.;
      if ( useXiInterp_ ) { xiInterp_->computeXiSpline( detid, *trk, xi, err_xi ); }
      else                { xiInterp_->computeXiLinear( detid, *trk, xi, err_xi ); }
      fProtonTrackXi[fProtonTrackNum] = xi;
      fProtonTrackXiError[fProtonTrackNum] = err_xi;
      fProtonTrackChi2[fProtonTrackNum] = trk->getChiSquared();
      fProtonTrackNormChi2[fProtonTrackNum] = trk->getChiSquaredOverNDF();
      
      switch ( pot ) {
      case 2: { map_near.insert( std::make_pair( fProtonTrackNum, localtrack_t( it->detId(), *trk ) ) ); } break;
      case 3: {  map_far.insert( std::make_pair( fProtonTrackNum, localtrack_t( it->detId(), *trk ) ) ); } break;
      }
      
      fProtonTrackNum++;
    }
  }



  /*   
  // second loop to associate near-far tracks
  for ( std::map<unsigned int,localtrack_t>::const_iterator it_n=map_near.begin(); it_n!=map_near.end(); it_n++ ) {
  float min_dist = 999.999;
      unsigned int cand = 999;
      for ( std::map<unsigned int,localtrack_t>::const_iterator it_f=map_far.begin(); it_f!=map_far.end(); it_f++ ) {
      const float dist = ProtonUtils::tracksDistance( align, it_n->second, it_f->second ); // in cm
      if ( dist<0. ) continue; // skip the comparison if opposite side tracks
      if ( dist<min_dist ) {
      min_dist = dist;
      cand = it_f->first;
      }
      }
      if ( cand!=999 ) { // near-far match found
      fProtonTrackLinkNF[it_n->first] = cand;
      fProtonTrackLinkNF[cand] = it_n->first;
      fProtonTrackMinLinkDist[it_n->first] = fProtonTrackMinLinkDist[cand] = min_dist;
      }
      }
  */

  
  //    std::cout << "Photon pT:  " << fPhotonPT << std::endl;
  //    std::cout << ">> Filling the tree" << std::endl;
  
  
  tree_->Fill();

  
}




// ------------ method called once each job just before starting event loop  ------------
void 
TreeProducer::beginJob()
{

  tree_->Branch( "run_id", &fRun, "run_id/i");
  tree_->Branch( "fill_number", &fFill, "fill_number/i");
  tree_->Branch( "lumisection", &fLumiSection, "lumisection/i");
  tree_->Branch( "bunch_crossing", &fBX, "bunch_crossing/i");
  tree_->Branch( "event_number", &fEventNum, "event_number/l");

  tree_->Branch( "photon_num", &fPhotonNum, "photon_num/i");
  tree_->Branch( "photon_pT", &fPhotonPT, "photon_pT/i");

  tree_->Branch( "cluster_strip_begin", &fClusterStripBegin, "cluster_strip_begin/s" );
  tree_->Branch( "cluster_strip_end", &fClusterStripEnd, "cluster_strip_end/s" );
  tree_->Branch( "cluster_num_strips", &fClusterNumStrips, "cluster_num_strips/s" );
  tree_->Branch( "cluster_center_strip_pos", &fClusterCenterStripPosition, "cluster_center_strip_pos/i" );

  tree_->Branch( "uv_a", &fUVPatternA, "uv_a/s" );
  tree_->Branch( "uv_b", &fUVPatternB, "uv_b/s" );
  tree_->Branch( "uv_w", &fUVPatternW, "uv_w/s" );

  tree_->Branch( "rh_pos", &fRecHitPos, "rh_pos/F" );
  tree_->Branch( "rh_sigma", &fRecHitSig, "rh_sigma/F" );

  tree_->Branch( "num_proton_track", &fProtonTrackNum, "num_proton_track/i" );
  tree_->Branch( "proton_track_x", fProtonTrackX, "proton_track_x[num_proton_track]/F" );
  tree_->Branch( "proton_track_y", fProtonTrackY, "proton_track_y[num_proton_track]/F" );
  tree_->Branch( "proton_track_x_corr", fProtonTrackXCorr, "proton_track_x_corr[num_proton_track]/F" );
  tree_->Branch( "proton_track_y_corr", fProtonTrackYCorr, "proton_track_y_corr[num_proton_track]/F" );
  tree_->Branch( "proton_track_xi", fProtonTrackXi, "proton_track_xi[num_proton_track]/F" );
  tree_->Branch( "proton_track_xi_error", fProtonTrackXiError, "proton_track_xi_error[num_proton_track]/F" );
  tree_->Branch( "proton_track_side", fProtonTrackSide, "proton_track_side[num_proton_track]/i" );
  tree_->Branch( "proton_track_chi2", fProtonTrackChi2, "proton_track_chi2[num_proton_track]/F" );
  tree_->Branch( "proton_track_normchi2", fProtonTrackNormChi2, "proton_track_normchi2[num_proton_track]/F" );
  tree_->Branch( "proton_track_pot", fProtonTrackPot, "proton_track_pot[num_proton_track]/i" );
  tree_->Branch( "proton_track_link_nearfar", fProtonTrackLinkNF, "proton_track_link_nearfar[num_proton_track]/i" );
  tree_->Branch( "proton_track_link_mindist", fProtonTrackMinLinkDist, "proton_track_link_mindist[num_proton_track]/F" );

}

// ------------ method called once each job just after ending the event loop  ------------


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TreeProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TreeProducer);
