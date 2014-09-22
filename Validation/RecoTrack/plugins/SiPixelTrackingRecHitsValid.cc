
// File: SiPixelTrackingRecHitsValid.cc
// Authors:  Xingtao Huang (Puerto Rico Univ.)
//           Gavril Giurgiu (JHU)
// Creation Date: Oct. 2006

#include <memory>
#include <string>
#include <iostream>
#include <TMath.h>
#include "Validation/RecoTrack/interface/SiPixelTrackingRecHitsValid.h"

#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/GluedGeomDet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/DetId/interface/DetId.h" 
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"

#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyBuilder.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "DataFormats/SiPixelDetId/interface/PixelBarrelName.h"
#include "DataFormats/SiPixelDetId/interface/PixelEndcapName.h"

#include "DQM/SiPixelCommon/interface/SiPixelFolderOrganizer.h"
#include "DQM/SiPixelCommon/src/SiPixelFolderOrganizer.cc"


#include <TTree.h>
#include <TFile.h>

using namespace std;

// End job: write and close the ntuple file
void SiPixelTrackingRecHitsValid::endJob() 
{
  if(debugNtuple_.size()!=0){
  tfile_->Write();
  tfile_->Close();
  }
}

void SiPixelTrackingRecHitsValid::beginJob()
{
  std::cout<<"beginning...."<<std::endl;

  if(debugNtuple_.size()!=0){
    tfile_ = new TFile (debugNtuple_.c_str() , "RECREATE");
    
    t_ = new TTree("Ntuple", "Ntuple");
    int bufsize = 64000;
  
    t_->Branch("subdetId", &subdetId, "subdetId/I", bufsize);
  
    t_->Branch("layer" , &layer , "layer/I" , bufsize);
    t_->Branch("ladder", &ladder, "ladder/I", bufsize);
    t_->Branch("mod"   , &mod   , "mod/I"   , bufsize);
    t_->Branch("side"  , &side  , "side/I"  , bufsize);
    t_->Branch("disk"  , &disk  , "disk/I"  , bufsize);
    t_->Branch("blade" , &blade , "blade/I" , bufsize);
    t_->Branch("panel" , &panel , "panel/I" , bufsize);
    t_->Branch("plaq"  , &plaq  , "plaq/I"  , bufsize);
    
    t_->Branch("rechitx"    , &rechitx    , "rechitx/F"    , bufsize);
    t_->Branch("rechity"    , &rechity    , "rechity/F"    , bufsize);
    t_->Branch("rechitz"    , &rechitz    , "rechitz/F"    , bufsize);
    t_->Branch("rechiterrx" , &rechiterrx , "rechiterrx/F" , bufsize);
    t_->Branch("rechiterry" , &rechiterry , "rechiterry/F" , bufsize);
    t_->Branch("rechitresx" , &rechitresx , "rechitresx/F" , bufsize);
    t_->Branch("rechitresy" , &rechitresy , "rechitresy/F" , bufsize);
    t_->Branch("rechitpullx", &rechitpullx, "rechitpullx/F", bufsize);
    t_->Branch("rechitpully", &rechitpully, "rechitpully/F", bufsize);
    
    t_->Branch("npix"  , &npix  , "npix/I"  , bufsize);
    t_->Branch("nxpix" , &nxpix , "nxpix/I" , bufsize);
    t_->Branch("nypix" , &nypix , "nypix/I" , bufsize);
    t_->Branch("charge", &charge, "charge/F", bufsize);
    
    t_->Branch("alpha", &alpha, "alpha/F", bufsize);
    t_->Branch("beta" , &beta , "beta/F" , bufsize);
    
    t_->Branch("phi", &phi, "phi/F", bufsize);
    t_->Branch("eta", &eta, "eta/F", bufsize);
    
    t_->Branch("half"   , &half   , "half/I"   , bufsize);
    t_->Branch("flipped", &flipped, "flipped/I", bufsize);
    
    t_->Branch("simhitx", &simhitx, "simhitx/F", bufsize);
    t_->Branch("simhity", &simhity, "simhity/F", bufsize);

    t_->Branch("nsimhit", &nsimhit, "nsimhit/I", bufsize);
    t_->Branch("pidhit" , &pidhit , "pidhit/I" , bufsize);
    
    t_->Branch("evt", &evt, "evt/I", bufsize);
    t_->Branch("run", &run, "run/I", bufsize);
  }
  
}

SiPixelTrackingRecHitsValid::SiPixelTrackingRecHitsValid(const edm::ParameterSet& ps):conf_(ps), dbe_(0), tfile_(0), t_(0)
{
  std::cout<<"beginning....in constructor"<<std::endl;
  //Read config file
  MTCCtrack_ = ps.getParameter<bool>("MTCCtrack");
  runStandalone = ps.getParameter<bool>("runStandalone");
  outputFile_ = ps.getUntrackedParameter<std::string>("outputFile", "pixeltrackingrechitshisto.root");
  siPixelRecHitCollectionToken_ = consumes<SiPixelRecHitCollection>( edm::InputTag( "siPixelRecHits" ) );
  recoTrackCollectionToken_ = consumes<reco::TrackCollection>( edm::InputTag( ps.getUntrackedParameter<std::string>( "src" ) ) );
  builderName_ = ps.getParameter<std::string>("TTRHBuilder");   
  checkType_ = ps.getParameter<bool>("checkType");
  genType_ = ps.getParameter<int>("genType");
  debugNtuple_=ps.getUntrackedParameter<std::string>("debugNtuple", "SiPixelTrackingRecHitsValid_Ntuple.root");


  //GET PARAMETERS--------------------TH2


  edm::ParameterSet ParametersResXvsAlphaBarrel = conf_.getParameter<edm::ParameterSet>("TH2ResXvsAlphaBarrel");
  edm::ParameterSet ParametersResYvsAlphaBarrel = conf_.getParameter<edm::ParameterSet>("TH2ResYvsAlphaBarrel");
  edm::ParameterSet ParametersResXvsBetaBarrel = conf_.getParameter<edm::ParameterSet>("TH2ResXvsBetaBarrel");
  edm::ParameterSet ParametersResYvsBetaBarrel = conf_.getParameter<edm::ParameterSet>("TH2ResYvsBetaBarrel");
  edm::ParameterSet ParametersPullXvsAlphaBarrel = conf_.getParameter<edm::ParameterSet>("TH2PullXvsAlphaBarrel");
  edm::ParameterSet ParametersPullYvsAlphaBarrel = conf_.getParameter<edm::ParameterSet>("TH2PullYvsAlphaBarrel");
  edm::ParameterSet ParametersPullXvsBetaBarrel = conf_.getParameter<edm::ParameterSet>("TH2PullXvsBetaBarrel");
  edm::ParameterSet ParametersPullYvsBetaBarrel = conf_.getParameter<edm::ParameterSet>("TH2PullYvsBetaBarrel");
  edm::ParameterSet ParametersPullXvsPhiBarrel = conf_.getParameter<edm::ParameterSet>("TH2PullXvsPhiBarrel");
  edm::ParameterSet ParametersPullYvsPhiBarrel = conf_.getParameter<edm::ParameterSet>("TH2PullYvsPhiBarrel");
  edm::ParameterSet ParametersPullXvsEtaBarrel = conf_.getParameter<edm::ParameterSet>("TH2PullXvsEtaBarrel");
  edm::ParameterSet ParametersPullYvsEtaBarrel = conf_.getParameter<edm::ParameterSet>("TH2PullYvsEtaBarrel");

  edm::ParameterSet ParametersWPullXvsAlphaBarrel = conf_.getParameter<edm::ParameterSet>("TH2WPullXvsAlphaBarrel");
  edm::ParameterSet ParametersWPullYvsAlphaBarrel = conf_.getParameter<edm::ParameterSet>("TH2WPullYvsAlphaBarrel");
  edm::ParameterSet ParametersWPullXvsBetaBarrel = conf_.getParameter<edm::ParameterSet>("TH2WPullXvsBetaBarrel");
  edm::ParameterSet ParametersWPullYvsBetaBarrel = conf_.getParameter<edm::ParameterSet>("TH2WPullYvsBetaBarrel");

  edm::ParameterSet ParametersResXvsAlphaZmPanel1 = conf_.getParameter<edm::ParameterSet>("TH2ResXvsAlphaZmPanel1");
  edm::ParameterSet ParametersResYvsAlphaZmPanel1 = conf_.getParameter<edm::ParameterSet>("TH2ResYvsAlphaZmPanel1");
  edm::ParameterSet ParametersResXvsBetaZmPanel1 = conf_.getParameter<edm::ParameterSet>("TH2ResXvsBetaZmPanel1");
  edm::ParameterSet ParametersResYvsBetaZmPanel1 = conf_.getParameter<edm::ParameterSet>("TH2ResYvsBetaZmPanel1");

  edm::ParameterSet ParametersPullXvsAlphaZmPanel1 = conf_.getParameter<edm::ParameterSet>("TH2PullXvsAlphaZmPanel1");
  edm::ParameterSet ParametersPullYvsAlphaZmPanel1 = conf_.getParameter<edm::ParameterSet>("TH2PullYvsAlphaZmPanel1");
  edm::ParameterSet ParametersPullXvsBetaZmPanel1 = conf_.getParameter<edm::ParameterSet>("TH2PullXvsBetaZmPanel1");
  edm::ParameterSet ParametersPullYvsBetaZmPanel1 = conf_.getParameter<edm::ParameterSet>("TH2PullYvsBetaZmPanel1");

  edm::ParameterSet ParametersPullXvsPhiZmPanel1 = conf_.getParameter<edm::ParameterSet>("TH2PullXvsPhiZmPanel1");
  edm::ParameterSet ParametersPullYvsPhiZmPanel1 = conf_.getParameter<edm::ParameterSet>("TH2PullYvsPhiZmPanel1");
  edm::ParameterSet ParametersPullXvsEtaZmPanel1 = conf_.getParameter<edm::ParameterSet>("TH2PullXvsEtaZmPanel1");
  edm::ParameterSet ParametersPullYvsEtaZmPanel1 = conf_.getParameter<edm::ParameterSet>("TH2PullYvsEtaZmPanel1");

  edm::ParameterSet ParametersWPullXvsAlphaZmPanel1 = conf_.getParameter<edm::ParameterSet>("TH2WPullXvsAlphaZmPanel1");
  edm::ParameterSet ParametersWPullYvsAlphaZmPanel1 = conf_.getParameter<edm::ParameterSet>("TH2WPullYvsAlphaZmPanel1");
  edm::ParameterSet ParametersWPullXvsBetaZmPanel1 = conf_.getParameter<edm::ParameterSet>("TH2WPullXvsBetaZmPanel1");
  edm::ParameterSet ParametersWPullYvsBetaZmPanel1 = conf_.getParameter<edm::ParameterSet>("TH2WPullYvsBetaZmPanel1");

  edm::ParameterSet ParametersResXvsAlphaZpPanel1 = conf_.getParameter<edm::ParameterSet>("TH2ResXvsAlphaZpPanel1");
  edm::ParameterSet ParametersResYvsAlphaZpPanel1 = conf_.getParameter<edm::ParameterSet>("TH2ResYvsAlphaZpPanel1");
  edm::ParameterSet ParametersResXvsBetaZpPanel1 = conf_.getParameter<edm::ParameterSet>("TH2ResXvsBetaZpPanel1");
  edm::ParameterSet ParametersResYvsBetaZpPanel1 = conf_.getParameter<edm::ParameterSet>("TH2ResYvsBetaZpPanel1");

  edm::ParameterSet ParametersPullXvsAlphaZpPanel1 = conf_.getParameter<edm::ParameterSet>("TH2PullXvsAlphaZpPanel1");
  edm::ParameterSet ParametersPullYvsAlphaZpPanel1 = conf_.getParameter<edm::ParameterSet>("TH2PullYvsAlphaZpPanel1");
  edm::ParameterSet ParametersPullXvsBetaZpPanel1 = conf_.getParameter<edm::ParameterSet>("TH2PullXvsBetaZpPanel1");
  edm::ParameterSet ParametersPullYvsBetaZpPanel1 = conf_.getParameter<edm::ParameterSet>("TH2PullYvsBetaZpPanel1");
  edm::ParameterSet ParametersPullXvsPhiZpPanel1 = conf_.getParameter<edm::ParameterSet>("TH2PullXvsPhiZpPanel1");
  edm::ParameterSet ParametersPullYvsPhiZpPanel1 = conf_.getParameter<edm::ParameterSet>("TH2PullYvsPhiZpPanel1");
  edm::ParameterSet ParametersPullXvsEtaZpPanel1 = conf_.getParameter<edm::ParameterSet>("TH2PullXvsEtaZpPanel1");
  edm::ParameterSet ParametersPullYvsEtaZpPanel1 = conf_.getParameter<edm::ParameterSet>("TH2PullYvsEtaZpPanel1");

  edm::ParameterSet ParametersWPullXvsAlphaZpPanel1 = conf_.getParameter<edm::ParameterSet>("TH2WPullXvsAlphaZpPanel1");
  edm::ParameterSet ParametersWPullYvsAlphaZpPanel1 = conf_.getParameter<edm::ParameterSet>("TH2WPullYvsAlphaZpPanel1");
  edm::ParameterSet ParametersWPullXvsBetaZpPanel1 = conf_.getParameter<edm::ParameterSet>("TH2WPullXvsBetaZpPanel1");
  edm::ParameterSet ParametersWPullYvsBetaZpPanel1 = conf_.getParameter<edm::ParameterSet>("TH2WPullYvsBetaZpPanel1");

  edm::ParameterSet ParametersResXvsAlphaZmPanel2 = conf_.getParameter<edm::ParameterSet>("TH2ResXvsAlphaZmPanel2");
  edm::ParameterSet ParametersResYvsAlphaZmPanel2 = conf_.getParameter<edm::ParameterSet>("TH2ResYvsAlphaZmPanel2");
  edm::ParameterSet ParametersResXvsBetaZmPanel2 = conf_.getParameter<edm::ParameterSet>("TH2ResXvsBetaZmPanel2");
  edm::ParameterSet ParametersResYvsBetaZmPanel2 = conf_.getParameter<edm::ParameterSet>("TH2ResYvsBetaZmPanel2");
  edm::ParameterSet ParametersPullXvsAlphaZmPanel2 = conf_.getParameter<edm::ParameterSet>("TH2PullXvsAlphaZmPanel2");
  edm::ParameterSet ParametersPullYvsAlphaZmPanel2 = conf_.getParameter<edm::ParameterSet>("TH2PullYvsAlphaZmPanel2");
  edm::ParameterSet ParametersPullXvsBetaZmPanel2 = conf_.getParameter<edm::ParameterSet>("TH2PullXvsBetaZmPanel2");
  edm::ParameterSet ParametersPullYvsBetaZmPanel2 = conf_.getParameter<edm::ParameterSet>("TH2PullYvsBetaZmPanel2");
  edm::ParameterSet ParametersPullXvsPhiZmPanel2 = conf_.getParameter<edm::ParameterSet>("TH2PullXvsPhiZmPanel2");
  edm::ParameterSet ParametersPullYvsPhiZmPanel2 = conf_.getParameter<edm::ParameterSet>("TH2PullYvsPhiZmPanel2");
  edm::ParameterSet ParametersPullXvsEtaZmPanel2 = conf_.getParameter<edm::ParameterSet>("TH2PullXvsEtaZmPanel2");
  edm::ParameterSet ParametersPullYvsEtaZmPanel2 = conf_.getParameter<edm::ParameterSet>("TH2PullYvsEtaZmPanel2");


  edm::ParameterSet ParametersResXvsAlphaZpPanel2 = conf_.getParameter<edm::ParameterSet>("TH2ResXvsAlphaZpPanel2");
  edm::ParameterSet ParametersResYvsAlphaZpPanel2 = conf_.getParameter<edm::ParameterSet>("TH2ResYvsAlphaZpPanel2");
  edm::ParameterSet ParametersResXvsBetaZpPanel2 = conf_.getParameter<edm::ParameterSet>("TH2ResXvsBetaZpPanel2");
  edm::ParameterSet ParametersResYvsBetaZpPanel2 = conf_.getParameter<edm::ParameterSet>("TH2ResYvsBetaZpPanel2");
  edm::ParameterSet ParametersPullXvsAlphaZpPanel2 = conf_.getParameter<edm::ParameterSet>("TH2PullXvsAlphaZpPanel2");
  edm::ParameterSet ParametersPullYvsAlphaZpPanel2 = conf_.getParameter<edm::ParameterSet>("TH2PullYvsAlphaZpPanel2");
  edm::ParameterSet ParametersPullXvsBetaZpPanel2 = conf_.getParameter<edm::ParameterSet>("TH2PullXvsBetaZpPanel2");
  edm::ParameterSet ParametersPullYvsBetaZpPanel2 = conf_.getParameter<edm::ParameterSet>("TH2PullYvsBetaZpPanel2");
  edm::ParameterSet ParametersPullXvsPhiZpPanel2 = conf_.getParameter<edm::ParameterSet>("TH2PullXvsPhiZpPanel2");
  edm::ParameterSet ParametersPullYvsPhiZpPanel2 = conf_.getParameter<edm::ParameterSet>("TH2PullYvsPhiZpPanel2");
  edm::ParameterSet ParametersPullXvsEtaZpPanel2 = conf_.getParameter<edm::ParameterSet>("TH2PullXvsEtaZpPanel2");
  edm::ParameterSet ParametersPullYvsEtaZpPanel2 = conf_.getParameter<edm::ParameterSet>("TH2PullYvsEtaZpPanel2");


  edm::ParameterSet ParametersWPullXvsAlphaZmPanel2 = conf_.getParameter<edm::ParameterSet>("TH2WPullXvsAlphaZmPanel2");
  edm::ParameterSet ParametersWPullYvsAlphaZmPanel2 = conf_.getParameter<edm::ParameterSet>("TH2WPullYvsAlphaZmPanel2");
  edm::ParameterSet ParametersWPullXvsBetaZmPanel2 = conf_.getParameter<edm::ParameterSet>("TH2WPullXvsBetaZmPanel2");
  edm::ParameterSet ParametersWPullYvsBetaZmPanel2 = conf_.getParameter<edm::ParameterSet>("TH2WPullYvsBetaZmPanel2");

  edm::ParameterSet ParametersWPullXvsAlphaZpPanel2 = conf_.getParameter<edm::ParameterSet>("TH2WPullXvsAlphaZpPanel2");
  edm::ParameterSet ParametersWPullYvsAlphaZpPanel2 = conf_.getParameter<edm::ParameterSet>("TH2WPullYvsAlphaZpPanel2");
  edm::ParameterSet ParametersWPullXvsBetaZpPanel2 = conf_.getParameter<edm::ParameterSet>("TH2WPullXvsBetaZpPanel2");
  edm::ParameterSet ParametersWPullYvsBetaZpPanel2 = conf_.getParameter<edm::ParameterSet>("TH2WPullYvsBetaZpPanel2");


  ///////////////////////////////////// DEFINE BOOL


  switchResXvsAlphaBarrelFlippedLaddersLayer = ParametersResXvsAlphaBarrel.getParameter<bool>("switchon_FlippedLaddersLayer");
  switchResXvsAlphaBarrelNonFlippedLaddersLayer = ParametersResXvsAlphaBarrel.getParameter<bool>("switchon_NonFlippedLaddersLayer");
  switchResXvsAlphaBarrelTotal = ParametersResXvsAlphaBarrel.getParameter<bool>("switchon_total");
  switchResXvsAlphaBarrelFlippedLaddersTotal = ParametersResXvsAlphaBarrel.getParameter<bool>("switch_FlippedLaddersLayerTotal");
  switchResXvsAlphaBarrelNonFlippedLaddersTotal = ParametersResXvsAlphaBarrel.getParameter<bool>("switch_NonFlippedLaddersTotal");
  switchResXvsAlphaBarrelLayerModule = ParametersResXvsAlphaBarrel.getParameter<bool>("switch_LayerModule");

  switchResYvsAlphaBarrelFlippedLaddersLayer = ParametersResYvsAlphaBarrel.getParameter<bool>("switchon_FlippedLaddersLayer");
  switchResYvsAlphaBarrelNonFlippedLaddersLayer = ParametersResYvsAlphaBarrel.getParameter<bool>("switchon_NonFlippedLaddersLayer");
  switchResYvsAlphaBarrelTotal = ParametersResYvsAlphaBarrel.getParameter<bool>("switchon_total");
  switchResYvsAlphaBarrelFlippedLaddersTotal = ParametersResYvsAlphaBarrel.getParameter<bool>("switch_FlippedLaddersLayerTotal");
  switchResYvsAlphaBarrelNonFlippedLaddersTotal = ParametersResYvsAlphaBarrel.getParameter<bool>("switch_NonFlippedLaddersTotal");
  switchResYvsAlphaBarrelLayerModule = ParametersResYvsAlphaBarrel.getParameter<bool>("switch_LayerModule");

  switchResXvsBetaBarrelFlippedLaddersLayer = ParametersResXvsBetaBarrel.getParameter<bool>("switchon_FlippedLaddersLayer");
  switchResXvsBetaBarrelNonFlippedLaddersLayer = ParametersResXvsBetaBarrel.getParameter<bool>("switchon_NonFlippedLaddersLayer");
  switchResXvsBetaBarrelTotal = ParametersResXvsBetaBarrel.getParameter<bool>("switchon_total");
  switchResXvsBetaBarrelFlippedLaddersTotal = ParametersResXvsBetaBarrel.getParameter<bool>("switch_FlippedLaddersLayerTotal");
  switchResXvsBetaBarrelNonFlippedLaddersTotal = ParametersResXvsBetaBarrel.getParameter<bool>("switch_NonFlippedLaddersTotal");
  switchResXvsBetaBarrelLayerModule = ParametersResXvsBetaBarrel.getParameter<bool>("switch_LayerModule");

  switchResYvsBetaBarrelFlippedLaddersLayer = ParametersResYvsBetaBarrel.getParameter<bool>("switchon_FlippedLaddersLayer");
  switchResYvsBetaBarrelNonFlippedLaddersLayer = ParametersResYvsBetaBarrel.getParameter<bool>("switchon_NonFlippedLaddersLayer");
  switchResYvsBetaBarrelTotal = ParametersResYvsBetaBarrel.getParameter<bool>("switchon_total");
  switchResYvsBetaBarrelFlippedLaddersTotal = ParametersResYvsBetaBarrel.getParameter<bool>("switch_FlippedLaddersLayerTotal");
  switchResYvsBetaBarrelNonFlippedLaddersTotal = ParametersResYvsBetaBarrel.getParameter<bool>("switch_NonFlippedLaddersTotal");
  switchResYvsBetaBarrelLayerModule = ParametersResYvsBetaBarrel.getParameter<bool>("switch_LayerModule");

  switchPullXvsAlphaBarrelTotal = ParametersPullXvsAlphaBarrel.getParameter<bool>("switchon_total");
  switchPullXvsAlphaBarrelFlippedLaddersTotal = ParametersPullXvsAlphaBarrel.getParameter<bool>("switchon_FlippedLaddersTotal");
  switchPullXvsAlphaBarrelNonFlippedLaddersTotal = ParametersPullXvsAlphaBarrel.getParameter<bool>("switchon_NonFlippedLaddersTotal");
  switchPullXvsAlphaBarrelLayerModule = ParametersPullXvsAlphaBarrel.getParameter<bool>("switchon_LayerModule");
  switchPullYvsAlphaBarrelTotal = ParametersPullYvsAlphaBarrel.getParameter<bool>("switchon_total");
  switchPullYvsAlphaBarrelFlippedLaddersTotal = ParametersPullYvsAlphaBarrel.getParameter<bool>("switchon_FlippedLaddersTotal");
  switchPullYvsAlphaBarrelNonFlippedLaddersTotal = ParametersPullYvsAlphaBarrel.getParameter<bool>("switchon_NonFlippedLaddersTotal");
  switchPullYvsAlphaBarrelLayerModule = ParametersPullYvsAlphaBarrel.getParameter<bool>("switchon_LayerModule");

  switchPullXvsBetaBarrelTotal = ParametersPullXvsBetaBarrel.getParameter<bool>("switchon_total");
  switchPullXvsBetaBarrelFlippedLaddersTotal = ParametersPullXvsBetaBarrel.getParameter<bool>("switchon_FlippedLaddersTotal");
  switchPullXvsBetaBarrelNonFlippedLaddersTotal = ParametersPullXvsBetaBarrel.getParameter<bool>("switchon_NonFlippedLaddersTotal");
  switchPullXvsBetaBarrelLayerModule = ParametersPullXvsBetaBarrel.getParameter<bool>("switchon_LayerModule");
  switchPullYvsBetaBarrelTotal = ParametersPullYvsBetaBarrel.getParameter<bool>("switchon_total");
  switchPullYvsBetaBarrelFlippedLaddersTotal = ParametersPullYvsBetaBarrel.getParameter<bool>("switchon_FlippedLaddersTotal");
  switchPullYvsBetaBarrelNonFlippedLaddersTotal = ParametersPullYvsBetaBarrel.getParameter<bool>("switchon_NonFlippedLaddersTotal");
  switchPullYvsBetaBarrelLayerModule = ParametersPullYvsBetaBarrel.getParameter<bool>("switchon_LayerModule");

  switchPullXvsPhiBarrelTotal = ParametersPullXvsPhiBarrel.getParameter<bool>("switchon_total");
  switchPullXvsPhiBarrelFlippedLaddersTotal = ParametersPullXvsPhiBarrel.getParameter<bool>("switchon_FlippedLaddersTotal");
  switchPullXvsPhiBarrelNonFlippedLaddersTotal = ParametersPullXvsPhiBarrel.getParameter<bool>("switchon_NonFlippedLaddersTotal");
  switchPullXvsPhiBarrelLayerModule = ParametersPullXvsPhiBarrel.getParameter<bool>("switchon_LayerModule");
  switchPullYvsPhiBarrelTotal = ParametersPullYvsPhiBarrel.getParameter<bool>("switchon_total");
  switchPullYvsPhiBarrelFlippedLaddersTotal = ParametersPullYvsPhiBarrel.getParameter<bool>("switchon_FlippedLaddersTotal");
  switchPullYvsPhiBarrelNonFlippedLaddersTotal = ParametersPullYvsPhiBarrel.getParameter<bool>("switchon_NonFlippedLaddersTotal");
  switchPullYvsPhiBarrelLayerModule = ParametersPullYvsPhiBarrel.getParameter<bool>("switchon_LayerModule");

  switchPullXvsEtaBarrelTotal = ParametersPullXvsEtaBarrel.getParameter<bool>("switchon_total");
  switchPullXvsEtaBarrelFlippedLaddersTotal = ParametersPullXvsEtaBarrel.getParameter<bool>("switchon_FlippedLaddersTotal");
  switchPullXvsEtaBarrelNonFlippedLaddersTotal = ParametersPullXvsEtaBarrel.getParameter<bool>("switchon_NonFlippedLaddersTotal");
  switchPullXvsEtaBarrelLayerModule = ParametersPullXvsEtaBarrel.getParameter<bool>("switchon_LayerModule");
  switchPullYvsEtaBarrelTotal = ParametersPullYvsEtaBarrel.getParameter<bool>("switchon_total");
  switchPullYvsEtaBarrelFlippedLaddersTotal = ParametersPullYvsEtaBarrel.getParameter<bool>("switchon_FlippedLaddersTotal");
  switchPullYvsEtaBarrelNonFlippedLaddersTotal = ParametersPullYvsEtaBarrel.getParameter<bool>("switchon_NonFlippedLaddersTotal");
  switchPullYvsEtaBarrelLayerModule = ParametersPullYvsEtaBarrel.getParameter<bool>("switchon_LayerModule");

  switchWPullXvsAlphaBarrelFlippedLaddersTotal = ParametersWPullXvsAlphaBarrel.getParameter<bool>("switchon_FlippedLaddersTotal");
  switchWPullXvsAlphaBarrelNonFlippedLaddersTotal = ParametersWPullXvsAlphaBarrel.getParameter<bool>("switchon_NonFlippedLaddersTotal");
  switchWPullYvsAlphaBarrelFlippedLaddersTotal = ParametersWPullYvsAlphaBarrel.getParameter<bool>("switchon_FlippedLaddersTotal");
  switchWPullYvsAlphaBarrelNonFlippedLaddersTotal = ParametersWPullYvsAlphaBarrel.getParameter<bool>("switchon_NonFlippedLaddersTotal");
  switchWPullXvsBetaBarrelFlippedLaddersTotal = ParametersWPullXvsBetaBarrel.getParameter<bool>("switchon_FlippedLaddersTotal");
  switchWPullXvsBetaBarrelNonFlippedLaddersTotal = ParametersWPullXvsBetaBarrel.getParameter<bool>("switchon_NonFlippedLaddersTotal");
  switchWPullYvsBetaBarrelFlippedLaddersTotal = ParametersWPullYvsBetaBarrel.getParameter<bool>("switchon_FlippedLaddersTotal");
  switchWPullYvsBetaBarrelNonFlippedLaddersTotal = ParametersWPullYvsBetaBarrel.getParameter<bool>("switchon_NonFlippedLaddersTotal");

  switchResXvsAlphaZmPanel1Total = ParametersResXvsAlphaZmPanel1.getParameter<bool>("switchon_total");
  switchResXvsAlphaZmPanel1DiskPlaquette = ParametersResXvsAlphaZmPanel1.getParameter<bool>("switchon_DiskPlaquette");
  switchResYvsAlphaZmPanel1Total = ParametersResYvsAlphaZmPanel1.getParameter<bool>("switchon_total");
  switchResYvsAlphaZmPanel1DiskPlaquette = ParametersResYvsAlphaZmPanel1.getParameter<bool>("switchon_DiskPlaquette");

  switchResXvsBetaZmPanel1Total = ParametersResXvsBetaZmPanel1.getParameter<bool>("switchon_total");
  switchResXvsBetaZmPanel1DiskPlaquette = ParametersResXvsBetaZmPanel1.getParameter<bool>("switchon_DiskPlaquette");
  switchResYvsBetaZmPanel1Total = ParametersResYvsBetaZmPanel1.getParameter<bool>("switchon_total");
  switchResYvsBetaZmPanel1DiskPlaquette = ParametersResYvsBetaZmPanel1.getParameter<bool>("switchon_DiskPlaquette");

  switchPullXvsAlphaZmPanel1Total = ParametersPullXvsAlphaZmPanel1.getParameter<bool>("switchon_total");
  switchPullXvsAlphaZmPanel1DiskPlaquette = ParametersPullXvsAlphaZmPanel1.getParameter<bool>("switchon_DiskPlaquette");
  switchPullYvsAlphaZmPanel1Total = ParametersPullYvsAlphaZmPanel1.getParameter<bool>("switchon_total");
  switchPullYvsAlphaZmPanel1DiskPlaquette = ParametersPullYvsAlphaZmPanel1.getParameter<bool>("switchon_DiskPlaquette");

  switchPullXvsBetaZmPanel1Total = ParametersPullXvsBetaZmPanel1.getParameter<bool>("switchon_total");
  switchPullXvsBetaZmPanel1DiskPlaquette = ParametersPullXvsBetaZmPanel1.getParameter<bool>("switchon_DiskPlaquette");
  switchPullYvsBetaZmPanel1Total = ParametersPullYvsBetaZmPanel1.getParameter<bool>("switchon_total");
  switchPullYvsBetaZmPanel1DiskPlaquette = ParametersPullYvsBetaZmPanel1.getParameter<bool>("switchon_DiskPlaquette");

  switchPullXvsPhiZmPanel1Total = ParametersPullXvsPhiZmPanel1.getParameter<bool>("switchon_total");
  switchPullXvsPhiZmPanel1DiskPlaquette = ParametersPullXvsPhiZmPanel1.getParameter<bool>("switchon_DiskPlaquette");
  switchPullYvsPhiZmPanel1Total = ParametersPullYvsPhiZmPanel1.getParameter<bool>("switchon_total");
  switchPullYvsPhiZmPanel1DiskPlaquette = ParametersPullYvsPhiZmPanel1.getParameter<bool>("switchon_DiskPlaquette");

  switchPullXvsEtaZmPanel1Total = ParametersPullXvsEtaZmPanel1.getParameter<bool>("switchon_total");
  switchPullXvsEtaZmPanel1DiskPlaquette = ParametersPullXvsEtaZmPanel1.getParameter<bool>("switchon_DiskPlaquette");
  switchPullYvsEtaZmPanel1Total = ParametersPullYvsEtaZmPanel1.getParameter<bool>("switchon_total");
  switchPullYvsEtaZmPanel1DiskPlaquette = ParametersPullYvsEtaZmPanel1.getParameter<bool>("switchon_DiskPlaquette");

  switchWPullXvsAlphaZmPanel1Total = ParametersWPullXvsAlphaZmPanel1.getParameter<bool>("switchon_total");
  switchWPullYvsAlphaZmPanel1Total = ParametersWPullYvsAlphaZmPanel1.getParameter<bool>("switchon_total");
  switchWPullXvsBetaZmPanel1Total = ParametersWPullXvsBetaZmPanel1.getParameter<bool>("switchon_total");
  switchWPullYvsBetaZmPanel1Total = ParametersWPullYvsBetaZmPanel1.getParameter<bool>("switchon_total");

  switchResXvsAlphaZpPanel1Total = ParametersResXvsAlphaZpPanel1.getParameter<bool>("switchon_total");
  switchResXvsAlphaZpPanel1DiskPlaquette = ParametersResXvsAlphaZpPanel1.getParameter<bool>("switchon_DiskPlaquette");
  switchResYvsAlphaZpPanel1Total = ParametersResYvsAlphaZpPanel1.getParameter<bool>("switchon_total");
  switchResYvsAlphaZpPanel1DiskPlaquette = ParametersResYvsAlphaZpPanel1.getParameter<bool>("switchon_DiskPlaquette");

  switchResXvsBetaZpPanel1Total = ParametersResXvsBetaZpPanel1.getParameter<bool>("switchon_total");
  switchResXvsBetaZpPanel1DiskPlaquette = ParametersResXvsBetaZpPanel1.getParameter<bool>("switchon_DiskPlaquette");
  switchResYvsBetaZpPanel1Total = ParametersResYvsBetaZpPanel1.getParameter<bool>("switchon_total");
  switchResYvsBetaZpPanel1DiskPlaquette = ParametersResYvsBetaZpPanel1.getParameter<bool>("switchon_DiskPlaquette");

  switchPullXvsAlphaZpPanel1Total = ParametersPullXvsAlphaZpPanel1.getParameter<bool>("switchon_total");
  switchPullXvsAlphaZpPanel1DiskPlaquette = ParametersPullXvsAlphaZpPanel1.getParameter<bool>("switchon_DiskPlaquette");
  switchPullYvsAlphaZpPanel1Total = ParametersPullYvsAlphaZpPanel1.getParameter<bool>("switchon_total");
  switchPullYvsAlphaZpPanel1DiskPlaquette = ParametersPullYvsAlphaZpPanel1.getParameter<bool>("switchon_DiskPlaquette");

  switchPullXvsBetaZpPanel1Total = ParametersPullXvsBetaZpPanel1.getParameter<bool>("switchon_total");
  switchPullXvsBetaZpPanel1DiskPlaquette = ParametersPullXvsBetaZpPanel1.getParameter<bool>("switchon_DiskPlaquette");
  switchPullYvsBetaZpPanel1Total = ParametersPullYvsBetaZpPanel1.getParameter<bool>("switchon_total");
  switchPullYvsBetaZpPanel1DiskPlaquette = ParametersPullYvsBetaZpPanel1.getParameter<bool>("switchon_DiskPlaquette");

  switchPullXvsPhiZpPanel1Total = ParametersPullXvsPhiZpPanel1.getParameter<bool>("switchon_total");
  switchPullXvsPhiZpPanel1DiskPlaquette = ParametersPullXvsPhiZpPanel1.getParameter<bool>("switchon_DiskPlaquette");
  switchPullYvsPhiZpPanel1Total = ParametersPullYvsPhiZpPanel1.getParameter<bool>("switchon_total");
  switchPullYvsPhiZpPanel1DiskPlaquette = ParametersPullYvsPhiZpPanel1.getParameter<bool>("switchon_DiskPlaquette");

  switchPullXvsEtaZpPanel1Total = ParametersPullXvsEtaZpPanel1.getParameter<bool>("switchon_total");
  switchPullXvsEtaZpPanel1DiskPlaquette = ParametersPullXvsEtaZpPanel1.getParameter<bool>("switchon_DiskPlaquette");
  switchPullYvsEtaZpPanel1Total = ParametersPullYvsEtaZpPanel1.getParameter<bool>("switchon_total");
  switchPullYvsEtaZpPanel1DiskPlaquette = ParametersPullYvsEtaZpPanel1.getParameter<bool>("switchon_DiskPlaquette");

  switchWPullXvsAlphaZpPanel1Total = ParametersWPullXvsAlphaZpPanel1.getParameter<bool>("switchon_total");
  switchWPullYvsAlphaZpPanel1Total = ParametersWPullYvsAlphaZpPanel1.getParameter<bool>("switchon_total");
  switchWPullXvsBetaZpPanel1Total = ParametersWPullXvsBetaZpPanel1.getParameter<bool>("switchon_total");
  switchWPullYvsBetaZpPanel1Total = ParametersWPullYvsBetaZpPanel1.getParameter<bool>("switchon_total");

  //Zm----Panel2
  switchResXvsAlphaZmPanel2Total = ParametersResXvsAlphaZmPanel2.getParameter<bool>("switchon_total");
  switchResXvsAlphaZmPanel2DiskPlaquette = ParametersResXvsAlphaZmPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchResYvsAlphaZmPanel2Total = ParametersResYvsAlphaZmPanel2.getParameter<bool>("switchon_total");
  switchResYvsAlphaZmPanel2DiskPlaquette = ParametersResYvsAlphaZmPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchResXvsBetaZmPanel2Total = ParametersResXvsBetaZmPanel2.getParameter<bool>("switchon_total");
  switchResXvsBetaZmPanel2DiskPlaquette = ParametersResXvsBetaZmPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchResYvsBetaZmPanel2Total = ParametersResYvsBetaZmPanel2.getParameter<bool>("switchon_total");
  switchResYvsBetaZmPanel2DiskPlaquette = ParametersResYvsBetaZmPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchPullXvsAlphaZmPanel2Total = ParametersPullXvsAlphaZmPanel2.getParameter<bool>("switchon_total");
  switchPullXvsAlphaZmPanel2DiskPlaquette = ParametersPullXvsAlphaZmPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchPullYvsAlphaZmPanel2Total = ParametersPullYvsAlphaZmPanel2.getParameter<bool>("switchon_total");
  switchPullYvsAlphaZmPanel2DiskPlaquette = ParametersPullYvsAlphaZmPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchPullXvsBetaZmPanel2Total = ParametersPullXvsBetaZmPanel2.getParameter<bool>("switchon_total");
  switchPullXvsBetaZmPanel2DiskPlaquette = ParametersPullXvsBetaZmPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchPullYvsBetaZmPanel2Total = ParametersPullYvsBetaZmPanel2.getParameter<bool>("switchon_total");
  switchPullYvsBetaZmPanel2DiskPlaquette = ParametersPullYvsBetaZmPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchPullXvsPhiZmPanel2Total = ParametersPullXvsPhiZmPanel2.getParameter<bool>("switchon_total");
  switchPullXvsPhiZmPanel2DiskPlaquette = ParametersPullXvsPhiZmPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchPullYvsPhiZmPanel2Total = ParametersPullYvsPhiZmPanel2.getParameter<bool>("switchon_total");
  switchPullYvsPhiZmPanel2DiskPlaquette = ParametersPullYvsPhiZmPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchPullXvsEtaZmPanel2Total = ParametersPullXvsEtaZmPanel2.getParameter<bool>("switchon_total");
  switchPullXvsEtaZmPanel2DiskPlaquette = ParametersPullXvsEtaZmPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchPullYvsEtaZmPanel2Total = ParametersPullYvsEtaZmPanel2.getParameter<bool>("switchon_total");
  switchPullYvsEtaZmPanel2DiskPlaquette = ParametersPullYvsEtaZmPanel2.getParameter<bool>("switchon_DiskPlaquette");

  //Zp---Panel2
  switchResXvsAlphaZpPanel2Total = ParametersResXvsAlphaZpPanel2.getParameter<bool>("switchon_total");
  switchResXvsAlphaZpPanel2DiskPlaquette = ParametersResXvsAlphaZpPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchResYvsAlphaZpPanel2Total = ParametersResYvsAlphaZpPanel2.getParameter<bool>("switchon_total");
  switchResYvsAlphaZpPanel2DiskPlaquette = ParametersResYvsAlphaZpPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchResXvsBetaZpPanel2Total = ParametersResXvsBetaZpPanel2.getParameter<bool>("switchon_total");
  switchResXvsBetaZpPanel2DiskPlaquette = ParametersResXvsBetaZpPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchResYvsBetaZpPanel2Total = ParametersResYvsBetaZpPanel2.getParameter<bool>("switchon_total");
  switchResYvsBetaZpPanel2DiskPlaquette = ParametersResYvsBetaZpPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchPullXvsAlphaZpPanel2Total = ParametersPullXvsAlphaZpPanel2.getParameter<bool>("switchon_total");
  switchPullXvsAlphaZpPanel2DiskPlaquette = ParametersPullXvsAlphaZpPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchPullYvsAlphaZpPanel2Total = ParametersPullYvsAlphaZpPanel2.getParameter<bool>("switchon_total");
  switchPullYvsAlphaZpPanel2DiskPlaquette = ParametersPullYvsAlphaZpPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchPullXvsBetaZpPanel2Total = ParametersPullXvsBetaZpPanel2.getParameter<bool>("switchon_total");
  switchPullXvsBetaZpPanel2DiskPlaquette = ParametersPullXvsBetaZpPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchPullYvsBetaZpPanel2Total = ParametersPullYvsBetaZpPanel2.getParameter<bool>("switchon_total");
  switchPullYvsBetaZpPanel2DiskPlaquette = ParametersPullYvsBetaZpPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchPullXvsPhiZpPanel2Total = ParametersPullXvsPhiZpPanel2.getParameter<bool>("switchon_total");
  switchPullXvsPhiZpPanel2DiskPlaquette = ParametersPullXvsPhiZpPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchPullYvsPhiZpPanel2Total = ParametersPullYvsPhiZpPanel2.getParameter<bool>("switchon_total");
  switchPullYvsPhiZpPanel2DiskPlaquette = ParametersPullYvsPhiZpPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchPullXvsEtaZpPanel2Total = ParametersPullXvsEtaZpPanel2.getParameter<bool>("switchon_total");
  switchPullXvsEtaZpPanel2DiskPlaquette = ParametersPullXvsEtaZpPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchPullYvsEtaZpPanel2Total = ParametersPullYvsEtaZpPanel2.getParameter<bool>("switchon_total");
  switchPullYvsEtaZpPanel2DiskPlaquette = ParametersPullYvsEtaZpPanel2.getParameter<bool>("switchon_DiskPlaquette");


  switchWPullXvsAlphaZmPanel2Total = ParametersWPullXvsAlphaZmPanel2.getParameter<bool>("switchon_total");
  switchWPullYvsAlphaZmPanel2Total = ParametersWPullYvsAlphaZmPanel2.getParameter<bool>("switchon_total");
  switchWPullXvsBetaZmPanel2Total = ParametersWPullXvsAlphaZmPanel2.getParameter<bool>("switchon_total");
  switchWPullYvsBetaZmPanel2Total = ParametersWPullYvsAlphaZmPanel2.getParameter<bool>("switchon_total");


  switchWPullXvsAlphaZpPanel2Total = ParametersWPullXvsAlphaZpPanel2.getParameter<bool>("switchon_total");
  switchWPullYvsAlphaZpPanel2Total = ParametersWPullYvsAlphaZpPanel2.getParameter<bool>("switchon_total");
  switchWPullXvsBetaZpPanel2Total = ParametersWPullXvsAlphaZpPanel2.getParameter<bool>("switchon_total");
  switchWPullYvsBetaZpPanel2Total = ParametersWPullYvsAlphaZpPanel2.getParameter<bool>("switchon_total");

  ////////////////

  ////////////////---------------------GET PARAMETERS FOR HISTOGRAMS------>  TH1
 
   
  edm::ParameterSet ParametersPosxBarrel = conf_.getParameter<edm::ParameterSet>("TH1PosxBarrel");
  edm::ParameterSet ParametersPosyBarrel = conf_.getParameter<edm::ParameterSet>("TH1PosyBarrel");
  edm::ParameterSet ParametersErrxBarrel = conf_.getParameter<edm::ParameterSet>("TH1ErrxBarrel");
  edm::ParameterSet ParametersErryBarrel = conf_.getParameter<edm::ParameterSet>("TH1ErryBarrel");
  edm::ParameterSet ParametersResxBarrel = conf_.getParameter<edm::ParameterSet>("TH1ResxBarrel");
  edm::ParameterSet ParametersResyBarrel = conf_.getParameter<edm::ParameterSet>("TH1ResyBarrel");
  edm::ParameterSet ParametersPullxBarrel = conf_.getParameter<edm::ParameterSet>("TH1PullxBarrel");
  edm::ParameterSet ParametersPullyBarrel = conf_.getParameter<edm::ParameterSet>("TH1PullyBarrel");
  edm::ParameterSet ParametersNpixBarrel = conf_.getParameter<edm::ParameterSet>("TH1NpixBarrel");
  edm::ParameterSet ParametersNxpixBarrel = conf_.getParameter<edm::ParameterSet>("TH1NxpixBarrel");
  edm::ParameterSet ParametersNypixBarrel = conf_.getParameter<edm::ParameterSet>("TH1NypixBarrel");
  edm::ParameterSet ParametersChargeBarrel = conf_.getParameter<edm::ParameterSet>("TH1ChargeBarrel");

  edm::ParameterSet ParametersPosxZmPanel1 = conf_.getParameter<edm::ParameterSet>("TH1PosxZmPanel1");
  edm::ParameterSet ParametersPosyZmPanel1 = conf_.getParameter<edm::ParameterSet>("TH1PosyZmPanel1");
  edm::ParameterSet ParametersErrxZmPanel1 = conf_.getParameter<edm::ParameterSet>("TH1ErrxZmPanel1");
  edm::ParameterSet ParametersErryZmPanel1 = conf_.getParameter<edm::ParameterSet>("TH1ErryZmPanel1");
  edm::ParameterSet ParametersResxZmPanel1 = conf_.getParameter<edm::ParameterSet>("TH1ResxZmPanel1");
  edm::ParameterSet ParametersResyZmPanel1 = conf_.getParameter<edm::ParameterSet>("TH1ResyZmPanel1");
  edm::ParameterSet ParametersPullxZmPanel1 = conf_.getParameter<edm::ParameterSet>("TH1PullxZmPanel1");
  edm::ParameterSet ParametersPullyZmPanel1 = conf_.getParameter<edm::ParameterSet>("TH1PullyZmPanel1");
  edm::ParameterSet ParametersNpixZmPanel1 = conf_.getParameter<edm::ParameterSet>("TH1NpixZmPanel1");
  edm::ParameterSet ParametersNxpixZmPanel1 = conf_.getParameter<edm::ParameterSet>("TH1NxpixZmPanel1");
  edm::ParameterSet ParametersNypixZmPanel1 = conf_.getParameter<edm::ParameterSet>("TH1NypixZmPanel1");
  edm::ParameterSet ParametersChargeZmPanel1 = conf_.getParameter<edm::ParameterSet>("TH1ChargeZmPanel1");

  edm::ParameterSet ParametersPosxZpPanel1 = conf_.getParameter<edm::ParameterSet>("TH1PosxZpPanel1");
  edm::ParameterSet ParametersPosyZpPanel1 = conf_.getParameter<edm::ParameterSet>("TH1PosyZpPanel1");
  edm::ParameterSet ParametersErrxZpPanel1 = conf_.getParameter<edm::ParameterSet>("TH1ErrxZpPanel1");
  edm::ParameterSet ParametersErryZpPanel1 = conf_.getParameter<edm::ParameterSet>("TH1ErryZpPanel1");
  edm::ParameterSet ParametersResxZpPanel1 = conf_.getParameter<edm::ParameterSet>("TH1ResxZpPanel1");
  edm::ParameterSet ParametersResyZpPanel1 = conf_.getParameter<edm::ParameterSet>("TH1ResyZpPanel1");
  edm::ParameterSet ParametersPullxZpPanel1 = conf_.getParameter<edm::ParameterSet>("TH1PullxZpPanel1");
  edm::ParameterSet ParametersPullyZpPanel1 = conf_.getParameter<edm::ParameterSet>("TH1PullyZpPanel1");
  edm::ParameterSet ParametersNpixZpPanel1 = conf_.getParameter<edm::ParameterSet>("TH1NpixZpPanel1");
  edm::ParameterSet ParametersNxpixZpPanel1 = conf_.getParameter<edm::ParameterSet>("TH1NxpixZpPanel1");
  edm::ParameterSet ParametersNypixZpPanel1 = conf_.getParameter<edm::ParameterSet>("TH1NypixZpPanel1");
  edm::ParameterSet ParametersChargeZpPanel1 = conf_.getParameter<edm::ParameterSet>("TH1ChargeZpPanel1");

  edm::ParameterSet ParametersPosxZmPanel2 = conf_.getParameter<edm::ParameterSet>("TH1PosxZmPanel2");
  edm::ParameterSet ParametersPosyZmPanel2 = conf_.getParameter<edm::ParameterSet>("TH1PosyZmPanel2");
  edm::ParameterSet ParametersErrxZmPanel2 = conf_.getParameter<edm::ParameterSet>("TH1ErrxZmPanel2");
  edm::ParameterSet ParametersErryZmPanel2 = conf_.getParameter<edm::ParameterSet>("TH1ErryZmPanel2");
  edm::ParameterSet ParametersResxZmPanel2 = conf_.getParameter<edm::ParameterSet>("TH1ResxZmPanel2");
  edm::ParameterSet ParametersResyZmPanel2 = conf_.getParameter<edm::ParameterSet>("TH1ResyZmPanel2");
  edm::ParameterSet ParametersPullxZmPanel2 = conf_.getParameter<edm::ParameterSet>("TH1PullxZmPanel2");
  edm::ParameterSet ParametersPullyZmPanel2 = conf_.getParameter<edm::ParameterSet>("TH1PullyZmPanel2");
  edm::ParameterSet ParametersNpixZmPanel2 = conf_.getParameter<edm::ParameterSet>("TH1NpixZmPanel2");
  edm::ParameterSet ParametersNxpixZmPanel2 = conf_.getParameter<edm::ParameterSet>("TH1NxpixZmPanel2");
  edm::ParameterSet ParametersNypixZmPanel2 = conf_.getParameter<edm::ParameterSet>("TH1NypixZmPanel2");
  edm::ParameterSet ParametersChargeZmPanel2 = conf_.getParameter<edm::ParameterSet>("TH1ChargeZmPanel2");


  edm::ParameterSet ParametersPosxZpPanel2 = conf_.getParameter<edm::ParameterSet>("TH1PosxZpPanel2");
  edm::ParameterSet ParametersPosyZpPanel2 = conf_.getParameter<edm::ParameterSet>("TH1PosyZpPanel2");
  edm::ParameterSet ParametersErrxZpPanel2 = conf_.getParameter<edm::ParameterSet>("TH1ErrxZpPanel2");
  edm::ParameterSet ParametersErryZpPanel2 = conf_.getParameter<edm::ParameterSet>("TH1ErryZpPanel2");
  edm::ParameterSet ParametersResxZpPanel2 = conf_.getParameter<edm::ParameterSet>("TH1ResxZpPanel2");
  edm::ParameterSet ParametersResyZpPanel2 = conf_.getParameter<edm::ParameterSet>("TH1ResyZpPanel2");
  edm::ParameterSet ParametersPullxZpPanel2 = conf_.getParameter<edm::ParameterSet>("TH1PullxZpPanel2");
  edm::ParameterSet ParametersPullyZpPanel2 = conf_.getParameter<edm::ParameterSet>("TH1PullyZpPanel2");
  edm::ParameterSet ParametersNpixZpPanel2 = conf_.getParameter<edm::ParameterSet>("TH1NpixZpPanel2");
  edm::ParameterSet ParametersNxpixZpPanel2 = conf_.getParameter<edm::ParameterSet>("TH1NxpixZpPanel2");
  edm::ParameterSet ParametersNypixZpPanel2 = conf_.getParameter<edm::ParameterSet>("TH1NypixZpPanel2");
  edm::ParameterSet ParametersChargeZpPanel2 = conf_.getParameter<edm::ParameterSet>("TH1ChargeZpPanel2");



  edm::ParameterSet ParametersPosxBarrel_all_hits = conf_.getParameter<edm::ParameterSet>("TH1PosxBarrel_all_hits");
  edm::ParameterSet ParametersPosyBarrel_all_hits = conf_.getParameter<edm::ParameterSet>("TH1PosyBarrel_all_hits");


  edm::ParameterSet ParametersPosxZmPanel1_all_hits = conf_.getParameter<edm::ParameterSet>("TH1PosxZmPanel1_all_hits");
  edm::ParameterSet ParametersPosyZmPanel1_all_hits = conf_.getParameter<edm::ParameterSet>("TH1PosyZmPanel1_all_hits");
  edm::ParameterSet ParametersPosxZmPanel2_all_hits = conf_.getParameter<edm::ParameterSet>("TH1PosxZmPanel2_all_hits");
  edm::ParameterSet ParametersPosyZmPanel2_all_hits = conf_.getParameter<edm::ParameterSet>("TH1PosyZmPanel2_all_hits");

  edm::ParameterSet ParametersPosxZpPanel1_all_hits = conf_.getParameter<edm::ParameterSet>("TH1PosxZpPanel1_all_hits");
  edm::ParameterSet ParametersPosyZpPanel1_all_hits = conf_.getParameter<edm::ParameterSet>("TH1PosyZpPanel1_all_hits");
  edm::ParameterSet ParametersPosxZpPanel2_all_hits = conf_.getParameter<edm::ParameterSet>("TH1PosxZpPanel2_all_hits");
  edm::ParameterSet ParametersPosyZpPanel2_all_hits = conf_.getParameter<edm::ParameterSet>("TH1PosyZpPanel2_all_hits");


  edm::ParameterSet ParametersTracksPerEvent = conf_.getParameter<edm::ParameterSet>("TH1TracksPerEvent");
  edm::ParameterSet ParametersPixRecHitsPerTrack = conf_.getParameter<edm::ParameterSet>("TH1PixRecHitsPerTrack");



  ///////////////////////////////////////////////   DEFINE BOOL
  //////////////////////////////////////////////

  switchPosxBarrelTotal = ParametersPosxBarrel.getParameter<bool>("switchon_total");
  switchPosxBarrelHalfModuleTotal = ParametersPosxBarrel.getParameter<bool>("switchon_HalfModuleTotal");
  switchPosxBarrelFullModuleTotal = ParametersPosxBarrel.getParameter<bool>("switchon_FullModuleTotal");
  switchPosxBarrelFlippedLaddersTotal = ParametersPosxBarrel.getParameter<bool>("switchon_FlippedLaddersTotal");
  switchPosxBarrelNonFlippedLaddersTotal = ParametersPosxBarrel.getParameter<bool>("switchon_NonFlippedLaddersTotal");
  switchPosxBarrelLayerModule = ParametersPosxBarrel.getParameter<bool>("switchon_LayerModule");

  switchPosyBarrelTotal = ParametersPosyBarrel.getParameter<bool>("switchon_total");
  switchPosyBarrelHalfModuleTotal = ParametersPosyBarrel.getParameter<bool>("switchon_HalfModuleTotal");
  switchPosyBarrelFullModuleTotal = ParametersPosyBarrel.getParameter<bool>("switchon_FullModuleTotal");
  switchPosyBarrelFlippedLaddersTotal = ParametersPosyBarrel.getParameter<bool>("switchon_FlippedLaddersTotal");
  switchPosyBarrelNonFlippedLaddersTotal = ParametersPosyBarrel.getParameter<bool>("switchon_NonFlippedLaddersTotal");
  switchPosyBarrelLayerModule = ParametersPosyBarrel.getParameter<bool>("switchon_LayerModule");

  switchErrxBarrelTotal = ParametersErrxBarrel.getParameter<bool>("switchon_total");
  switchErrxBarrelLayerModule = ParametersErrxBarrel.getParameter<bool>("switchon_LayerModule");
  switchErryBarrelTotal = ParametersErryBarrel.getParameter<bool>("switchon_total");
  switchErryBarrelLayerModule = ParametersErryBarrel.getParameter<bool>("switchon_LayerModule");

  switchResxBarrelTotal = ParametersResxBarrel.getParameter<bool>("switchon_total");
  switchResxBarrelLayerModule = ParametersResxBarrel.getParameter<bool>("switchon_LayerModule");
  switchResxBarrelLayer = ParametersResxBarrel.getParameter<bool>("switchon_Layer");
  switchResyBarrelTotal = ParametersResyBarrel.getParameter<bool>("switchon_total");
  switchResyBarrelLayerModule = ParametersResyBarrel.getParameter<bool>("switchon_LayerModule");
  switchResyBarrelLayer = ParametersResyBarrel.getParameter<bool>("switchon_Layer");

  switchPullxBarrelTotal = ParametersPullxBarrel.getParameter<bool>("switchon_total");
  switchPullxBarrelLayerModule = ParametersPullxBarrel.getParameter<bool>("switchon_LayerModule");
  switchPullxBarrelLayer = ParametersPullxBarrel.getParameter<bool>("switchon_Layer");
  switchPullyBarrelTotal = ParametersPullyBarrel.getParameter<bool>("switchon_total");
  switchPullyBarrelLayerModule = ParametersPullyBarrel.getParameter<bool>("switchon_LayerModule");
  switchPullyBarrelLayer = ParametersPullyBarrel.getParameter<bool>("switchon_Layer");

  switchNpixBarrelTotal = ParametersNpixBarrel.getParameter<bool>("switchon_total");
  switchNpixBarrelLayerModule = ParametersNpixBarrel.getParameter<bool>("switchon_LayerModule");
  switchNxpixBarrelTotal = ParametersNxpixBarrel.getParameter<bool>("switchon_total");
  switchNxpixBarrelLayerModule = ParametersNxpixBarrel.getParameter<bool>("switchon_LayerModule");
  switchNypixBarrelTotal = ParametersNypixBarrel.getParameter<bool>("switchon_total");
  switchNypixBarrelLayerModule = ParametersNypixBarrel.getParameter<bool>("switchon_LayerModule");

  switchChargeBarrelTotal = ParametersChargeBarrel.getParameter<bool>("switchon_total");
  switchChargeBarrelLayerModule = ParametersChargeBarrel.getParameter<bool>("switchon_LayerModule");


  switchPosxZmPanel1Total = ParametersPosxZmPanel1.getParameter<bool>("switchon_total");
  switchPosxZmPanel1DiskPlaquette = ParametersPosxZmPanel1.getParameter<bool>("switchon_DiskPlaquette");
  switchPosyZmPanel1Total = ParametersPosyZmPanel1.getParameter<bool>("switchon_total");
  switchPosyZmPanel1DiskPlaquette = ParametersPosyZmPanel1.getParameter<bool>("switchon_DiskPlaquette");
  switchErrxZmPanel1Total = ParametersErrxZmPanel1.getParameter<bool>("switchon_total");
  switchErrxZmPanel1DiskPlaquette = ParametersErrxZmPanel1.getParameter<bool>("switchon_DiskPlaquette");
  switchErryZmPanel1Total = ParametersErryZmPanel1.getParameter<bool>("switchon_total");
  switchErryZmPanel1DiskPlaquette = ParametersErryZmPanel1.getParameter<bool>("switchon_DiskPlaquette");
  switchResxZmPanel1Total = ParametersResxZmPanel1.getParameter<bool>("switchon_total");
  switchResxZmPanel1DiskPlaquette = ParametersResxZmPanel1.getParameter<bool>("switchon_DiskPlaquette");
  switchResyZmPanel1Total = ParametersResyZmPanel1.getParameter<bool>("switchon_total");
  switchResyZmPanel1DiskPlaquette = ParametersResyZmPanel1.getParameter<bool>("switchon_DiskPlaquette");
  switchPullxZmPanel1Total = ParametersPullxZmPanel1.getParameter<bool>("switchon_total");
  switchPullxZmPanel1DiskPlaquette = ParametersPullxZmPanel1.getParameter<bool>("switchon_DiskPlaquette");
  switchPullyZmPanel1Total = ParametersPullyZmPanel1.getParameter<bool>("switchon_total");
  switchPullyZmPanel1DiskPlaquette = ParametersPullyZmPanel1.getParameter<bool>("switchon_DiskPlaquette");
  switchNpixZmPanel1Total = ParametersNpixZmPanel1.getParameter<bool>("switchon_total");
  switchNpixZmPanel1DiskPlaquette = ParametersNpixZmPanel1.getParameter<bool>("switchon_DiskPlaquette");
  switchNxpixZmPanel1Total = ParametersNxpixZmPanel1.getParameter<bool>("switchon_total");
  switchNxpixZmPanel1DiskPlaquette = ParametersNxpixZmPanel1.getParameter<bool>("switchon_DiskPlaquette");
  switchNypixZmPanel1Total = ParametersNypixZmPanel1.getParameter<bool>("switchon_total");
  switchNypixZmPanel1DiskPlaquette = ParametersNypixZmPanel1.getParameter<bool>("switchon_DiskPlaquette");
  switchChargeZmPanel1Total = ParametersChargeZmPanel1.getParameter<bool>("switchon_total");
  switchChargeZmPanel1DiskPlaquette = ParametersChargeZmPanel1.getParameter<bool>("switchon_DiskPlaquette");

  switchPosxZpPanel1Total = ParametersPosxZpPanel1.getParameter<bool>("switchon_total");
  switchPosxZpPanel1DiskPlaquette = ParametersPosxZpPanel1.getParameter<bool>("switchon_DiskPlaquette");
  switchPosyZpPanel1Total = ParametersPosyZpPanel1.getParameter<bool>("switchon_total");
  switchPosyZpPanel1DiskPlaquette = ParametersPosyZpPanel1.getParameter<bool>("switchon_DiskPlaquette");

  switchErrxZpPanel1Total = ParametersErrxZpPanel1.getParameter<bool>("switchon_total");
  switchErrxZpPanel1DiskPlaquette = ParametersErrxZpPanel1.getParameter<bool>("switchon_DiskPlaquette");
  switchErryZpPanel1Total = ParametersErryZpPanel1.getParameter<bool>("switchon_total");
  switchErryZpPanel1DiskPlaquette = ParametersErryZpPanel1.getParameter<bool>("switchon_DiskPlaquette");

  switchResxZpPanel1Total = ParametersResxZpPanel1.getParameter<bool>("switchon_total");
  switchResxZpPanel1DiskPlaquette = ParametersResxZpPanel1.getParameter<bool>("switchon_DiskPlaquette");
  switchResyZpPanel1Total = ParametersResyZpPanel1.getParameter<bool>("switchon_total");
  switchResyZpPanel1DiskPlaquette = ParametersResyZpPanel1.getParameter<bool>("switchon_DiskPlaquette");

  switchPullxZpPanel1Total = ParametersPullxZpPanel1.getParameter<bool>("switchon_total");
  switchPullxZpPanel1DiskPlaquette = ParametersPullxZpPanel1.getParameter<bool>("switchon_DiskPlaquette");
  switchPullyZpPanel1Total = ParametersPullyZpPanel1.getParameter<bool>("switchon_total");
  switchPullyZpPanel1DiskPlaquette = ParametersPullyZpPanel1.getParameter<bool>("switchon_DiskPlaquette");

  switchNpixZpPanel1Total = ParametersNpixZpPanel1.getParameter<bool>("switchon_total");
  switchNpixZpPanel1DiskPlaquette = ParametersNpixZpPanel1.getParameter<bool>("switchon_DiskPlaquette");

  switchNxpixZpPanel1Total = ParametersNxpixZpPanel1.getParameter<bool>("switchon_total");
  switchNxpixZpPanel1DiskPlaquette = ParametersNxpixZpPanel1.getParameter<bool>("switchon_DiskPlaquette");
  switchNypixZpPanel1Total = ParametersNypixZpPanel1.getParameter<bool>("switchon_total");
  switchNypixZpPanel1DiskPlaquette = ParametersNypixZpPanel1.getParameter<bool>("switchon_DiskPlaquette");

  switchChargeZpPanel1Total = ParametersChargeZpPanel1.getParameter<bool>("switchon_total");
  switchChargeZpPanel1DiskPlaquette = ParametersChargeZpPanel1.getParameter<bool>("switchon_DiskPlaquette");

  //Zm-Panel2
  switchPosxZmPanel2Total = ParametersPosxZmPanel2.getParameter<bool>("switchon_total");
  switchPosxZmPanel2DiskPlaquette = ParametersPosxZmPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchPosyZmPanel2Total = ParametersPosyZmPanel2.getParameter<bool>("switchon_total");
  switchPosyZmPanel2DiskPlaquette = ParametersPosyZmPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchErrxZmPanel2Total = ParametersErrxZmPanel2.getParameter<bool>("switchon_total");
  switchErrxZmPanel2DiskPlaquette = ParametersErrxZmPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchErryZmPanel2Total = ParametersErryZmPanel2.getParameter<bool>("switchon_total");
  switchErryZmPanel2DiskPlaquette = ParametersErryZmPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchResxZmPanel2Total = ParametersResxZmPanel2.getParameter<bool>("switchon_total");
  switchResxZmPanel2DiskPlaquette = ParametersResxZmPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchResyZmPanel2Total = ParametersResyZmPanel2.getParameter<bool>("switchon_total");
  switchResyZmPanel2DiskPlaquette = ParametersResyZmPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchPullxZmPanel2Total = ParametersPullxZmPanel2.getParameter<bool>("switchon_total");
  switchPullxZmPanel2DiskPlaquette = ParametersPullxZmPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchPullyZmPanel2Total = ParametersPullyZmPanel2.getParameter<bool>("switchon_total");
  switchPullyZmPanel2DiskPlaquette = ParametersPullyZmPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchNpixZmPanel2Total = ParametersNpixZmPanel2.getParameter<bool>("switchon_total");
  switchNpixZmPanel2DiskPlaquette = ParametersNpixZmPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchNxpixZmPanel2Total = ParametersNxpixZmPanel2.getParameter<bool>("switchon_total");
  switchNxpixZmPanel2DiskPlaquette = ParametersNxpixZmPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchNypixZmPanel2Total = ParametersNypixZmPanel2.getParameter<bool>("switchon_total");
  switchNypixZmPanel2DiskPlaquette = ParametersNypixZmPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchChargeZmPanel2Total = ParametersChargeZmPanel2.getParameter<bool>("switchon_total");
  switchChargeZmPanel2DiskPlaquette = ParametersChargeZmPanel2.getParameter<bool>("switchon_DiskPlaquette");

  //Zp-Panel2
  switchPosxZpPanel2Total = ParametersPosxZpPanel2.getParameter<bool>("switchon_total");
  switchPosxZpPanel2DiskPlaquette = ParametersPosxZpPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchPosyZpPanel2Total = ParametersPosyZpPanel2.getParameter<bool>("switchon_total");
  switchPosyZpPanel2DiskPlaquette = ParametersPosyZpPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchErrxZpPanel2Total = ParametersErrxZpPanel2.getParameter<bool>("switchon_total");
  switchErrxZpPanel2DiskPlaquette = ParametersErrxZpPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchErryZpPanel2Total = ParametersErryZpPanel2.getParameter<bool>("switchon_total");
  switchErryZpPanel2DiskPlaquette = ParametersErryZpPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchResxZpPanel2Total = ParametersResxZpPanel2.getParameter<bool>("switchon_total");
  switchResxZpPanel2DiskPlaquette = ParametersResxZpPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchResyZpPanel2Total = ParametersResyZpPanel2.getParameter<bool>("switchon_total");
  switchResyZpPanel2DiskPlaquette = ParametersResyZpPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchPullxZpPanel2Total = ParametersPullxZpPanel2.getParameter<bool>("switchon_total");
  switchPullxZpPanel2DiskPlaquette = ParametersPullxZpPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchPullyZpPanel2Total = ParametersPullyZpPanel2.getParameter<bool>("switchon_total");
  switchPullyZpPanel2DiskPlaquette = ParametersPullyZpPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchNpixZpPanel2Total = ParametersNpixZpPanel2.getParameter<bool>("switchon_total");
  switchNpixZpPanel2DiskPlaquette = ParametersNpixZpPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchNxpixZpPanel2Total = ParametersNxpixZpPanel2.getParameter<bool>("switchon_total");
  switchNxpixZpPanel2DiskPlaquette = ParametersNxpixZpPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchNypixZpPanel2Total = ParametersNypixZpPanel2.getParameter<bool>("switchon_total");
  switchNypixZpPanel2DiskPlaquette = ParametersNypixZpPanel2.getParameter<bool>("switchon_DiskPlaquette");

  switchChargeZpPanel2Total = ParametersChargeZpPanel2.getParameter<bool>("switchon_total");
  switchChargeZpPanel2DiskPlaquette = ParametersChargeZpPanel2.getParameter<bool>("switchon_DiskPlaquette");


  switchPosxBarrel_all_hits = ParametersPosxBarrel_all_hits.getParameter<bool>("switchon_total");
  switchPosyBarrel_all_hits = ParametersPosyBarrel_all_hits.getParameter<bool>("switchon_total");

  switchPosxZmPanel1_all_hits = ParametersPosxZmPanel1_all_hits.getParameter<bool>("switchon_total");
  switchPosyZmPanel1_all_hits = ParametersPosyZmPanel1_all_hits.getParameter<bool>("switchon_total");
  switchPosxZmPanel2_all_hits = ParametersPosxZmPanel2_all_hits.getParameter<bool>("switchon_total");
  switchPosyZmPanel2_all_hits = ParametersPosyZmPanel2_all_hits.getParameter<bool>("switchon_total");

  switchPosxZpPanel1_all_hits = ParametersPosxZpPanel1_all_hits.getParameter<bool>("switchon_total");
  switchPosyZpPanel1_all_hits = ParametersPosyZpPanel1_all_hits.getParameter<bool>("switchon_total");
  switchPosxZpPanel2_all_hits = ParametersPosxZpPanel2_all_hits.getParameter<bool>("switchon_total");
  switchPosyZpPanel2_all_hits = ParametersPosyZpPanel2_all_hits.getParameter<bool>("switchon_total");


  switchTracksPerEvent = ParametersTracksPerEvent.getParameter<bool>("switchon_total");
  switchPixRecHitsPerTrack = ParametersPixRecHitsPerTrack.getParameter<bool>("switchon_total");


  ////////////////
  ////////////////



}


MonitorElement* SiPixelTrackingRecHitsValid::bookME2D(DQMStore::IBooker &ibooker, const char* ParameterSetLabel, const char* HistoName, const char* HistoTitle, const char* ST) {

  Parameters = conf_.getParameter<edm::ParameterSet>(ParameterSetLabel);

  return ibooker.bookProfile(HistoName,HistoTitle,
                        Parameters.getParameter<int32_t>("Nbinx"),
                        Parameters.getParameter<double>("xmin"),
                        Parameters.getParameter<double>("xmax"),
                        Parameters.getParameter<int32_t>("Nbiny"),
                        Parameters.getParameter<double>("ymin"),
                        Parameters.getParameter<double>("ymax"),
                        ST
                       );
  return 0;

}


MonitorElement* SiPixelTrackingRecHitsValid::bookME1D(DQMStore::IBooker &ibooker, const char* ParameterSetLabel, const char* HistoName, const char* HistoTitle) {

  Parameters = conf_.getParameter<edm::ParameterSet>(ParameterSetLabel);

  return ibooker.book1D( HistoName, HistoTitle,
                         Parameters.getParameter<int32_t>("Nbinx"),
                         Parameters.getParameter<double>("xmin"),
                         Parameters.getParameter<double>("xmax")
                       );
  return 0;
}


void SiPixelTrackingRecHitsValid::createTotalMEs(DQMStore::IBooker &ibooker, std::string label){

    
    std::string HistoName;
    std::string TitleName;
    std::string ST;
    
    if (switchPosxBarrelTotal){

           HistoName="mePosxBarrel";
           TitleName="mePosxBarrel";

           totalMEs->mePosxBarrel=bookME1D(ibooker, "TH1PosxBarrel",HistoName.c_str(), TitleName.c_str());

    }  

    if (switchPosyBarrelTotal){

           HistoName="mePosyBarrel";
           TitleName="mePosyBarrel";

           totalMEs->mePosxBarrel=bookME1D(ibooker, "TH1PosyBarrel",HistoName.c_str(), TitleName.c_str());

    }

    if (switchErrxBarrelTotal){

           HistoName="meErrxBarrel";
           TitleName="meErrxBarrel";

           totalMEs->meErrxBarrel=bookME1D(ibooker, "TH1ErrxBarrel",HistoName.c_str(), TitleName.c_str());

    }

    if (switchErryBarrelTotal){

           HistoName="meErryBarrel";
           TitleName="meErryBarrel";

           totalMEs->meErrxBarrel=bookME1D(ibooker, "TH1ErryBarrel",HistoName.c_str(), TitleName.c_str());

    }

     if (switchResxBarrelTotal){

           HistoName="meResxBarrel";
           TitleName="meResxBarrel";

           totalMEs->meResxBarrel=bookME1D(ibooker, "TH1ResxBarrel",HistoName.c_str(), TitleName.c_str());

     }
   
     if (switchResyBarrelTotal){

           HistoName="meResyBarrel";
           TitleName="meResyBarrel";

           totalMEs->meResyBarrel=bookME1D(ibooker, "TH1ResyBarrel",HistoName.c_str(), TitleName.c_str());

     }

     if (switchPullxBarrelTotal){

           HistoName="mePullxBarrel";
           TitleName="mePullxBarrel";

           totalMEs->mePullxBarrel=bookME1D(ibooker, "TH1PullxBarrel",HistoName.c_str(), TitleName.c_str());

     }

     if (switchPullyBarrelTotal){

           HistoName="mePullyBarrel";
           TitleName="mePullyBarrel";

           totalMEs->mePullyBarrel=bookME1D(ibooker, "TH1PullyBarrel",HistoName.c_str(), TitleName.c_str());

     }

      if (switchNpixBarrelTotal){

           HistoName="meNpixBarrel";
           TitleName="meNpixBarrel";

           totalMEs->meNpixBarrel=bookME1D(ibooker, "TH1NpixBarrel",HistoName.c_str(), TitleName.c_str());

     }
   
     if (switchNxpixBarrelTotal){

           HistoName="meNxpixBarrel";
           TitleName="meNxpixBarrel";

           totalMEs->meNxpixBarrel=bookME1D(ibooker, "TH1NxpixBarrel",HistoName.c_str(), TitleName.c_str());

     }

     if (switchNypixBarrelTotal){

           HistoName="meNypixBarrel";
           TitleName="meNypixBarrel";

           totalMEs->meNypixBarrel=bookME1D(ibooker, "TH1NypixBarrel",HistoName.c_str(), TitleName.c_str());

     }

     if (switchChargeBarrelTotal){

           HistoName="meChargeBarrel";
           TitleName="meChargeBarrel";

           totalMEs->meNypixBarrel=bookME1D(ibooker, "TH1ChargeBarrel",HistoName.c_str(), TitleName.c_str());

     }

     if (switchResXvsAlphaBarrelTotal){

           HistoName="meResXvsAlphaBarrel";
           TitleName="meResXvsAlphaBarrel";
           ST="";

           totalMEs->meResXvsAlphaBarrel=bookME2D(ibooker, "TH2ResXvsAlphaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchResYvsAlphaBarrelTotal){

           HistoName="meResYvsAlphaBarrel";
           TitleName="meResYvsAlphaBarrel";
           ST="";

           totalMEs->meResYvsAlphaBarrel=bookME2D(ibooker, "TH2ResYvsAlphaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchResXvsBetaBarrelTotal){

           HistoName="meResXvsBetaBarrel";
           TitleName="meResXvsBetaBarrel";
           ST="";

           totalMEs->meResXvsBetaBarrel=bookME2D(ibooker, "TH2ResXvsBetaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchResYvsBetaBarrelTotal){

           HistoName="meResYvsBetaBarrel";
           TitleName="meResYvsBetaBarrel";
           ST="";

           totalMEs->meResYvsBetaBarrel=bookME2D(ibooker, "TH2ResYvsBetaBarrelTotal", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }


     if (switchPullXvsAlphaBarrelTotal){

           HistoName="mePullXvsAlphaBarrel";
           TitleName="mePullXvsAlphaBarrel";
           ST="";

           totalMEs->mePullXvsAlphaBarrel=bookME2D(ibooker, "TH2PullXvsAlphaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchPullYvsAlphaBarrelTotal){

           HistoName="mePullYvsAlphaBarrel";
           TitleName="mePullYvsAlphaBarrel";
           ST="";

           totalMEs->mePullYvsAlphaBarrel=bookME2D(ibooker, "TH2PullYvsAlphaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchPullXvsBetaBarrelTotal){

           HistoName="mePullXvsBetaBarrel";
           TitleName="mePullXvsBetaBarrel";
           ST="";

           totalMEs->mePullXvsBetaBarrel=bookME2D(ibooker, "TH2PullXvsBetaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchPullYvsBetaBarrelTotal){

           HistoName="mePullYvsBetaBarrel";
           TitleName="mePullYvsBetaBarrel";
           ST="";

           totalMEs->mePullYvsBetaBarrel=bookME2D(ibooker, "TH2PullYvsBetaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

    
     if (switchPullXvsPhiBarrelTotal){

           HistoName="mePullXvsPhiBarrel";
           TitleName="mePullXvsPhiBarrel";
           ST="";

           totalMEs->mePullXvsPhiBarrel=bookME2D(ibooker, "TH2PullXvsPhiBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchPullYvsPhiBarrelTotal){

           HistoName="mePullYvsPhiBarrel";
           TitleName="mePullYvsPhiBarrel";
           ST="";

           totalMEs->mePullYvsPhiBarrel=bookME2D(ibooker, "TH2PullYvsPhiBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }


     if (switchPullXvsEtaBarrelTotal){

           HistoName="mePullXvsEtaBarrel";
           TitleName="mePullXvsEtaBarrel";
           ST="";

           totalMEs->mePullXvsEtaBarrel=bookME2D(ibooker, "TH2PullXvsEtaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchPullYvsEtaBarrelTotal){

           HistoName="mePullYvsEtaBarrel";
           TitleName="mePullYvsEtaBarrel";
           ST="";

           totalMEs->mePullYvsEtaBarrel=bookME2D(ibooker, "TH2PullYvsEtaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }


     if (switchPosxBarrelHalfModuleTotal){

           HistoName="mePosxBarrelHalfModule";
           TitleName="mePosxBarrelHalfModule";

           totalMEs->mePosxBarrelHalfModule=bookME1D(ibooker, "TH1PosxBarrel", HistoName.c_str(), TitleName.c_str());

     }


     if (switchPosxBarrelFullModuleTotal){

           HistoName="mePosxBarrelFullModule";
           TitleName="mePosxBarrelFullModule";

           totalMEs->mePosxBarrelFullModule=bookME1D(ibooker, "TH1PosxBarrel", HistoName.c_str(), TitleName.c_str());

     }

     if (switchPosxBarrelFlippedLaddersTotal){

           HistoName="mePosxBarrelFlippedLadders";
           TitleName="mePosxBarrelFlippedLadders";

           totalMEs->mePosxBarrelFlippedLadders=bookME1D(ibooker, "TH1PosxBarrel", HistoName.c_str(), TitleName.c_str());

     }

     if (switchPosxBarrelNonFlippedLaddersTotal){

           HistoName="mePosxBarrelNonFlippedLadders";
           TitleName="mePosxBarrelNonFlippedLadders";

           totalMEs->mePosxBarrelNonFlippedLadders=bookME1D(ibooker, "TH1PosxBarrel", HistoName.c_str(), TitleName.c_str());

     }

     if (switchPosxBarrelHalfModuleTotal){

           HistoName="mePosyBarrelHalfModule";
           TitleName="mePosyBarrelHalfModule";

           totalMEs->mePosyBarrelHalfModule=bookME1D(ibooker, "TH1PosyBarrel", HistoName.c_str(), TitleName.c_str());

     }

     if (switchPosyBarrelFullModuleTotal){

           HistoName="mePosyBarrelFullModule";
           TitleName="mePosyBarrelFullModule";

           totalMEs->mePosyBarrelFullModule=bookME1D(ibooker, "TH1PosyBarrel", HistoName.c_str(), TitleName.c_str());

     }
    
     if (switchPosyBarrelFlippedLaddersTotal){

           HistoName="mePosyBarrelFlippedLadders";
           TitleName="mePosyBarrelFlippedLadders";

           totalMEs->mePosyBarrelFlippedLadders=bookME1D(ibooker, "TH1PosyBarrel", HistoName.c_str(), TitleName.c_str());

     }

     if (switchPosyBarrelNonFlippedLaddersTotal){

           HistoName="mePosyBarrelNonFlippedLadders";
           TitleName="mePosyBarrelNonFlippedLadders";

           totalMEs->mePosyBarrelNonFlippedLadders=bookME1D(ibooker, "TH1PosyBarrel", HistoName.c_str(), TitleName.c_str());

     }

     if (switchResXvsAlphaBarrelFlippedLaddersTotal){

           HistoName="meResXvsAlphaBarrelFlippedLadders";
           TitleName="meResXvsAlphaBarrelFlippedLadders";

           totalMEs->meResXvsAlphaBarrelFlippedLadders=bookME2D(ibooker, "TH2ResXvsAlphaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

      if (switchResYvsAlphaBarrelFlippedLaddersTotal){

           HistoName="meResYvsAlphaBarrelFlippedLadders";
           TitleName="meResYvsAlphaBarrelFlippedLadders";

           totalMEs->meResYvsAlphaBarrelFlippedLadders=bookME2D(ibooker, "TH2ResYvsAlphaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

      if (switchResXvsBetaBarrelFlippedLaddersTotal){

           HistoName="meResXvsBetaBarrelFlippedLadders";
           TitleName="meResXvsBetaBarrelFlippedLadders";

           totalMEs->meResXvsBetaBarrelFlippedLadders=bookME2D(ibooker, "TH2ResXvsBetaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchResYvsBetaBarrelFlippedLaddersTotal){

           HistoName="meResYvsBetaBarrelFlippedLadders";
           TitleName="meResYvsBetaBarrelFlippedLadders";

           totalMEs->meResYvsBetaBarrelFlippedLadders=bookME2D(ibooker, "TH2ResYvsBetaBarrelFlippedLaddersLayerTotal", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchPullXvsAlphaBarrelFlippedLaddersTotal){

           HistoName="mePullXvsAlphaBarrelFlippedLadders";
           TitleName="mePullXvsAlphaBarrelFlippedLadders";

           totalMEs->mePullXvsAlphaBarrelFlippedLadders=bookME2D(ibooker, "TH2PullXvsAlphaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchPullYvsAlphaBarrelFlippedLaddersTotal){

           HistoName="mePullYvsAlphaBarrelFlippedLadders";
           TitleName="mePullYvsAlphaBarrelFlippedLadders";

           totalMEs->mePullYvsAlphaBarrelFlippedLadders=bookME2D(ibooker, "TH2PullYvsAlphaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchPullXvsBetaBarrelFlippedLaddersTotal){

           HistoName="mePullXvsBetaBarrelFlippedLadders";
           TitleName="mePullXvsBetaBarrelFlippedLadders";

           totalMEs->mePullXvsBetaBarrelFlippedLadders=bookME2D(ibooker, "TH2PullXvsBetaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchPullYvsBetaBarrelFlippedLaddersTotal){

           HistoName="mePullYvsBetaBarrelFlippedLadders";
           TitleName="mePullYvsBetaBarrelFlippedLadders";

           totalMEs->mePullYvsBetaBarrelFlippedLadders=bookME2D(ibooker, "TH2PullYvsBetaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchPullXvsPhiBarrelFlippedLaddersTotal){

           HistoName="mePullXvsPhiBarrelFlippedLadders";
           TitleName="mePullXvsPhiBarrelFlippedLadders";

           totalMEs->mePullXvsPhiBarrelFlippedLadders=bookME2D(ibooker, "TH2PullXvsPhiBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchPullYvsPhiBarrelFlippedLaddersTotal){

           HistoName="mePullYvsPhiBarrelFlippedLadders";
           TitleName="mePullYvsPhiBarrelFlippedLadders";

           totalMEs->mePullYvsPhiBarrelFlippedLadders=bookME2D(ibooker, "TH2PullYvsPhiBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchPullXvsEtaBarrelFlippedLaddersTotal){

           HistoName="mePullXvsEtaBarrelFlippedLadders";
           TitleName="mePullXvsEtaBarrelFlippedLadders";

           totalMEs->mePullXvsEtaBarrelFlippedLadders=bookME2D(ibooker, "TH2PullXvsEtaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchPullYvsEtaBarrelFlippedLaddersTotal){

           HistoName="mePullYvsEtaBarrelFlippedLadders";
           TitleName="mePullYvsEtaBarrelFlippedLadders";

           totalMEs->mePullYvsEtaBarrelFlippedLadders=bookME2D(ibooker, "TH2PullYvsEtaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

  
     if (switchWPullXvsAlphaBarrelFlippedLaddersTotal){

           HistoName="meWPullXvsAlphaBarrelFlippedLadders";
           TitleName="meWPullXvsAlphaBarrelFlippedLadders";

           totalMEs->meWPullXvsAlphaBarrelFlippedLadders=bookME2D(ibooker, "TH2WPullXvsAlphaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchWPullYvsAlphaBarrelFlippedLaddersTotal){

           HistoName="meWPullYvsAlphaBarrelFlippedLadders";
           TitleName="meWPullYvsAlphaBarrelFlippedLadders";

           totalMEs->meWPullYvsAlphaBarrelFlippedLadders=bookME2D(ibooker, "TH2WPullYvsAlphaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchWPullXvsBetaBarrelFlippedLaddersTotal){

           HistoName="meWPullXvsBetaBarrelFlippedLadders";
           TitleName="meWPullXvsBetaBarrelFlippedLadders";

           totalMEs->meWPullXvsBetaBarrelFlippedLadders=bookME2D(ibooker, "TH2WPullXvsBetaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchWPullYvsBetaBarrelFlippedLaddersTotal){

           HistoName="meWPullYvsBetaBarrelFlippedLadders";
           TitleName="meWPullYvsBetaBarrelFlippedLadders";

           totalMEs->meWPullYvsBetaBarrelFlippedLadders=bookME2D(ibooker, "TH2WPullYvsBetaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchResXvsAlphaBarrelNonFlippedLaddersTotal){

           HistoName="meResXvsAlphaBarrelNonFlippedLadders";
           TitleName="meResXvsAlphaBarrelNonFlippedLadders";

           totalMEs->meResXvsAlphaBarrelNonFlippedLadders=bookME2D(ibooker, "TH2ResXvsAlphaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchResYvsAlphaBarrelNonFlippedLaddersTotal){

           HistoName="meResYvsAlphaBarrelNonFlippedLadders";
           TitleName="meResYvsAlphaBarrelNonFlippedLadders";

           totalMEs->meResYvsAlphaBarrelNonFlippedLadders=bookME2D(ibooker, "TH2ResYvsAlphaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

      if (switchResXvsBetaBarrelNonFlippedLaddersTotal){

           HistoName="meResXvsBetaBarrelNonFlippedLadders";
           TitleName="meResXvsBetaBarrelNonFlippedLadders";

           totalMEs->meResXvsBetaBarrelNonFlippedLadders=bookME2D(ibooker, "TH2ResXvsBetaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchResYvsBetaBarrelNonFlippedLaddersTotal){

           HistoName="meResYvsBetaBarrelNonFlippedLadders";
           TitleName="meResYvsBetaBarrelNonFlippedLadders";

           totalMEs->meResYvsBetaBarrelNonFlippedLadders=bookME2D(ibooker, "TH2ResYvsBetaBarrelNonFlippedLaddersTotal", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchPullXvsAlphaBarrelNonFlippedLaddersTotal){

           HistoName="mePullXvsAlphaBarrelNonFlippedLadders";
           TitleName="mePullXvsAlphaBarrelNonFlippedLadders";

           totalMEs->mePullXvsAlphaBarrelNonFlippedLadders=bookME2D(ibooker, "TH2PullXvsAlphaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchPullYvsAlphaBarrelNonFlippedLaddersTotal){

           HistoName="mePullYvsAlphaBarrelNonFlippedLadders";
           TitleName="mePullYvsAlphaBarrelNonFlippedLadders";

           totalMEs->mePullYvsAlphaBarrelNonFlippedLadders=bookME2D(ibooker, "TH2PullYvsAlphaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchPullXvsBetaBarrelNonFlippedLaddersTotal){

           HistoName="mePullXvsBetaBarrelNonFlippedLadders";
           TitleName="mePullXvsBetaBarrelNonFlippedLadders";

           totalMEs->mePullXvsBetaBarrelNonFlippedLadders=bookME2D(ibooker, "TH2PullXvsBetaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchPullYvsBetaBarrelNonFlippedLaddersTotal){

           HistoName="mePullYvsBetaBarrelNonFlippedLadders";
           TitleName="mePullYvsBetaBarrelNonFlippedLadders";

           totalMEs->mePullYvsBetaBarrelNonFlippedLadders=bookME2D(ibooker, "TH2PullYvsBetaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchPullXvsPhiBarrelNonFlippedLaddersTotal){

           HistoName="mePullXvsPhiBarrelNonFlippedLadders";
           TitleName="mePullXvsPhiBarrelNonFlippedLadders";

           totalMEs->mePullXvsPhiBarrelNonFlippedLadders=bookME2D(ibooker, "TH2PullXvsPhiBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchPullYvsPhiBarrelNonFlippedLaddersTotal){

           HistoName="mePullYvsPhiBarrelNonFlippedLadders";
           TitleName="mePullYvsPhiBarrelNonFlippedLadders";

           totalMEs->mePullYvsPhiBarrelNonFlippedLadders=bookME2D(ibooker, "TH2PullYvsPhiBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchPullXvsEtaBarrelNonFlippedLaddersTotal){

           HistoName="mePullXvsEtaBarrelNonFlippedLadders";
           TitleName="mePullXvsEtaBarrelNonFlippedLadders";

           totalMEs->mePullXvsEtaBarrelNonFlippedLadders=bookME2D(ibooker, "TH2PullXvsEtaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchPullYvsEtaBarrelNonFlippedLaddersTotal){

           HistoName="mePullYvsEtaBarrelNonFlippedLadders";
           TitleName="mePullYvsEtaBarrelNonFlippedLadders";

           totalMEs->mePullYvsEtaBarrelNonFlippedLadders=bookME2D(ibooker, "TH2PullYvsEtaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }
 
     if (switchWPullXvsAlphaBarrelNonFlippedLaddersTotal){

           HistoName="meWPullXvsAlphaBarrelNonFlippedLadders";
           TitleName="meWPullXvsAlphaBarrelNonFlippedLadders";

           totalMEs->meWPullXvsAlphaBarrelNonFlippedLadders=bookME2D(ibooker, "TH2WPullXvsAlphaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchWPullYvsAlphaBarrelNonFlippedLaddersTotal){

           HistoName="meWPullYvsAlphaBarrelNonFlippedLadders";
           TitleName="meWPullYvsAlphaBarrelNonFlippedLadders";

           totalMEs->meWPullYvsAlphaBarrelNonFlippedLadders=bookME2D(ibooker, "TH2WPullYvsAlphaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchWPullXvsBetaBarrelNonFlippedLaddersTotal){

           HistoName="meWPullXvsBetaBarrelNonFlippedLadders";
           TitleName="meWPullXvsBetaBarrelNonFlippedLadders";

           totalMEs->meWPullXvsBetaBarrelNonFlippedLadders=bookME2D(ibooker, "TH2WPullXvsBetaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchWPullYvsBetaBarrelNonFlippedLaddersTotal){

           HistoName="meWPullYvsBetaBarrelNonFlippedLadders";
           TitleName="meWPullYvsBetaBarrelNonFlippedLadders";

           totalMEs->meWPullYvsBetaBarrelNonFlippedLadders=bookME2D(ibooker, "TH2WPullYvsBetaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }


     if (switchPosxZmPanel1Total){

           HistoName="mePosxZmPanel1";
           TitleName="mePosxZmPanel1";

           totalMEs->mePosxZmPanel1 = bookME1D(ibooker, "TH1PosxZmPanel1", HistoName.c_str(), TitleName.c_str());

     }

     if (switchPosyZmPanel1Total){

           HistoName="mePosyZmPanel1";
           TitleName="mePosyZmPanel1";

           totalMEs->mePosyZmPanel1 = bookME1D(ibooker, "TH1PosyZmPanel1", HistoName.c_str(), TitleName.c_str());

     }

     if (switchErrxZmPanel1Total){

           HistoName="meErrxZmPanel1";
           TitleName="meErrxZmPanel1";

           totalMEs->meErrxZmPanel1 = bookME1D(ibooker, "TH1ErrxZmPanel1", HistoName.c_str(), TitleName.c_str());

     }

     if (switchErryZmPanel1Total){

           HistoName="meErryZmPanel1";
           TitleName="meErryZmPanel1";

           totalMEs->meErryZmPanel1 = bookME1D(ibooker, "TH1ErryZmPanel1", HistoName.c_str(), TitleName.c_str());

     }


     if (switchResxZmPanel1Total){

           HistoName="meResxZmPanel1";
           TitleName="meResxZmPanel1";

           totalMEs->meResxZmPanel1 = bookME1D(ibooker, "TH1ResxZmPanel1", HistoName.c_str(), TitleName.c_str());

     }

     if (switchResyZmPanel1Total){

           HistoName="meResyZmPanel1";
           TitleName="meResyZmPanel1";

           totalMEs->meResyZmPanel1 = bookME1D(ibooker, "TH1ResyZmPanel1", HistoName.c_str(), TitleName.c_str());

     }
    
     if (switchPullxZmPanel1Total){

           HistoName="mePullxZmPanel1";
           TitleName="mePullxZmPanel1";

           totalMEs->mePullxZmPanel1 = bookME1D(ibooker, "TH1PullxZmPanel1", HistoName.c_str(), TitleName.c_str());

     }

     if (switchPullyZmPanel1Total){

           HistoName="mePullyZmPanel1";
           TitleName="mePullyZmPanel1";

           totalMEs->mePullyZmPanel1 = bookME1D(ibooker, "TH1PullyZmPanel1", HistoName.c_str(), TitleName.c_str());

     }


      if (switchNpixZmPanel1Total){

           HistoName="meNpixZmPanel1";
           TitleName="meNpixZmPanel1";

           totalMEs->meNpixZmPanel1 = bookME1D(ibooker, "TH1NpixZmPanel1", HistoName.c_str(), TitleName.c_str());

     }

     if (switchNxpixZmPanel1Total){

           HistoName="meNxpixZmPanel1";
           TitleName="meNxpixZmPanel1";

           totalMEs->meNxpixZmPanel1 = bookME1D(ibooker, "TH1NxpixZmPanel1", HistoName.c_str(), TitleName.c_str());

     }

     if (switchNypixZmPanel1Total){

           HistoName="meNypixZmPanel1";
           TitleName="meNypixZmPanel1";

           totalMEs->meNypixZmPanel1 = bookME1D(ibooker, "TH1NypixZmPanel1", HistoName.c_str(), TitleName.c_str());

     }

     if (switchChargeZmPanel1Total){

           HistoName="meChargeZmPanel1";
           TitleName="meChargeZmPanel1";

           totalMEs->meChargeZmPanel1 = bookME1D(ibooker, "TH1ChargeZmPanel1", HistoName.c_str(), TitleName.c_str());

     }

      if (switchResXvsAlphaZmPanel1Total){

           HistoName="meResXvsAlphaZmPanel1";
           TitleName="meResXvsAlphaZmPanel1";

           totalMEs->meResXvsAlphaZmPanel1 = bookME2D(ibooker, "TH2ResXvsAlphaZmPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

      if (switchResYvsAlphaZmPanel1Total){

           HistoName="meResYvsAlphaZmPanel1";
           TitleName="meResYvsAlphaZmPanel1";

           totalMEs->meResYvsAlphaZmPanel1 = bookME2D(ibooker, "TH2ResYvsAlphaZmPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchResXvsBetaZmPanel1Total){

           HistoName="meResXvsBetaZmPanel1";
           TitleName="meResXvsBetaZmPanel1";

           totalMEs->meResXvsBetaZmPanel1 = bookME2D(ibooker, "TH2ResXvsBetaZmPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchResYvsBetaZmPanel1Total){

           HistoName="meResYvsBetaZmPanel1";
           TitleName="meResYvsBetaZmPanel1";

           totalMEs->meResYvsBetaZmPanel1 = bookME2D(ibooker, "TH2ResYvsBetaZmPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }
   
      if (switchPullXvsAlphaZmPanel1Total){

           HistoName="mePullXvsAlphaZmPanel1";
           TitleName="mePullXvsAlphaZmPanel1";

           totalMEs->mePullXvsAlphaZmPanel1 = bookME2D(ibooker, "TH2PullXvsAlphaZmPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

      if (switchPullYvsAlphaZmPanel1Total){

           HistoName="mePullYvsAlphaZmPanel1";
           TitleName="mePullYvsAlphaZmPanel1";

           totalMEs->mePullYvsAlphaZmPanel1 = bookME2D(ibooker, "TH2PullYvsAlphaZmPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

      if (switchPullXvsBetaZmPanel1Total){

           HistoName="mePullXvsBetaZmPanel1";
           TitleName="mePullXvsBetaZmPanel1";

           totalMEs->mePullXvsBetaZmPanel1 = bookME2D(ibooker, "TH2PullXvsBetaZmPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

      if (switchPullYvsBetaZmPanel1Total){

           HistoName="mePullYvsBetaZmPanel1";
           TitleName="mePullYvsBetaZmPanel1";

           totalMEs->mePullYvsBetaZmPanel1 = bookME2D(ibooker, "TH2PullYvsBetaZmPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

      if (switchPullXvsPhiZmPanel1Total){

           HistoName="mePullXvsPhiZmPanel1";
           TitleName="mePullXvsPhiZmPanel1";

           totalMEs->mePullXvsPhiZmPanel1 = bookME2D(ibooker, "TH2PullXvsPhiZmPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

      if (switchPullYvsPhiZmPanel1Total){

           HistoName="mePullYvsPhiZmPanel1";
           TitleName="mePullYvsPhiZmPanel1";

           totalMEs->mePullYvsPhiZmPanel1 = bookME2D(ibooker, "TH2PullYvsPhiZmPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

      if (switchPullXvsEtaZmPanel1Total){

           HistoName="mePullXvsEtaZmPanel1";
           TitleName="mePullXvsEtaZmPanel1";

           totalMEs->mePullXvsEtaZmPanel1 = bookME2D(ibooker, "TH2PullXvsEtaZmPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }


      if (switchPullYvsEtaZmPanel1Total){

           HistoName="mePullYvsEtaZmPanel1";
           TitleName="mePullYvsEtaZmPanel1";

           totalMEs->mePullYvsEtaZmPanel1 = bookME2D(ibooker, "TH2PullYvsEtaZmPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }


     if (switchWPullXvsAlphaZmPanel1Total){

           HistoName="meWPullXvsAlphaZmPanel1";
           TitleName="meWPullXvsAlphaZmPanel1";

           totalMEs->meWPullXvsAlphaZmPanel1 = bookME2D(ibooker, "TH2WPullXvsAlphaZmPanel1Total", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

      if (switchWPullYvsAlphaZmPanel1Total){

           HistoName="meWPullYvsAlphaZmPanel1";
           TitleName="meWPullYvsAlphaZmPanel1";

           totalMEs->meWPullYvsAlphaZmPanel1 = bookME2D(ibooker, "TH2WPullYvsAlphaZmPanel1Total", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

      if (switchWPullXvsBetaZmPanel1Total){

           HistoName="meWPullXvsBetaZmPanel1";
           TitleName="meWPullXvsBetaZmPanel1";

           totalMEs->meWPullXvsBetaZmPanel1 = bookME2D(ibooker, "TH2WPullXvsBetaZmPanel1Total", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

      if (switchWPullYvsBetaZmPanel1Total){

           HistoName="meWPullYvsBetaZmPanel1";
           TitleName="meWPullYvsBetaZmPanel1";

           totalMEs->meWPullYvsBetaZmPanel1 = bookME2D(ibooker, "TH2WPullYvsBetaZmPanel1Total", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     //Total-Zp-Panel
     if (switchPosxZpPanel1Total){

           HistoName="mePosxZpPanel1";
           TitleName="mePosxZpPanel1";

           totalMEs->mePosxZpPanel1 = bookME1D(ibooker, "TH1PosxZpPanel1", HistoName.c_str(), TitleName.c_str());

     }

     if (switchPosyZpPanel1Total){

           HistoName="mePosyZpPanel1";
           TitleName="mePosyZpPanel1";

           totalMEs->mePosyZpPanel1 = bookME1D(ibooker, "TH1PosyZpPanel1", HistoName.c_str(), TitleName.c_str());

     }

     if (switchErrxZpPanel1Total){

           HistoName="meErrxZpPanel1";
           TitleName="meErrxZpPanel1";

           totalMEs->meErrxZpPanel1 = bookME1D(ibooker, "TH1ErrxZpPanel1", HistoName.c_str(), TitleName.c_str());

     }

     if (switchErryZpPanel1Total){

           HistoName="meErryZpPanel1";
           TitleName="meErryZpPanel1";

           totalMEs->meErryZpPanel1 = bookME1D(ibooker, "TH1ErryZpPanel1", HistoName.c_str(), TitleName.c_str());

     }


     if (switchResxZpPanel1Total){

           HistoName="meResxZpPanel1";
           TitleName="meResxZpPanel1";

           totalMEs->meResxZpPanel1 = bookME1D(ibooker, "TH1ResxZpPanel1", HistoName.c_str(), TitleName.c_str());

     }

     if (switchResyZpPanel1Total){

           HistoName="meResyZpPanel1";
           TitleName="meResyZpPanel1";

           totalMEs->meResyZpPanel1 = bookME1D(ibooker, "TH1ResyZpPanel1", HistoName.c_str(), TitleName.c_str());

     }
    
     if (switchPullxZpPanel1Total){

           HistoName="mePullxZpPanel1";
           TitleName="mePullxZpPanel1";

           totalMEs->mePullxZpPanel1 = bookME1D(ibooker, "TH1PullxZpPanel1", HistoName.c_str(), TitleName.c_str());

     }

     if (switchPullyZpPanel1Total){

           HistoName="mePullyZpPanel1";
           TitleName="mePullyZpPanel1";

           totalMEs->mePullyZpPanel1 = bookME1D(ibooker, "TH1PullyZpPanel1", HistoName.c_str(), TitleName.c_str());

     }


      if (switchNpixZpPanel1Total){

           HistoName="meNpixZpPanel1";
           TitleName="meNpixZpPanel1";

           totalMEs->meNpixZpPanel1 = bookME1D(ibooker, "TH1NpixZpPanel1", HistoName.c_str(), TitleName.c_str());

     }

     if (switchNxpixZpPanel1Total){

           HistoName="meNxpixZpPanel1";
           TitleName="meNxpixZpPanel1";

           totalMEs->meNxpixZpPanel1 = bookME1D(ibooker, "TH1NxpixZpPanel1", HistoName.c_str(), TitleName.c_str());

     }

     if (switchNypixZpPanel1Total){

           HistoName="meNypixZpPanel1";
           TitleName="meNypixZpPanel1";

           totalMEs->meNypixZpPanel1 = bookME1D(ibooker, "TH1NypixZpPanel1", HistoName.c_str(), TitleName.c_str());

     }

     if (switchChargeZpPanel1Total){

           HistoName="meChargeZpPanel1";
           TitleName="meChargeZpPanel1";

           totalMEs->meChargeZpPanel1 = bookME1D(ibooker, "TH1ChargeZpPanel1", HistoName.c_str(), TitleName.c_str());

     }

      if (switchResXvsAlphaZpPanel1Total){

           HistoName="meResXvsAlphaZpPanel1";
           TitleName="meResXvsAlphaZpPanel1";
           ST="";

           totalMEs->meResXvsAlphaZpPanel1 = bookME2D(ibooker, "TH2ResXvsAlphaZpPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

      if (switchResYvsAlphaZpPanel1Total){

           HistoName="meResYvsAlphaZpPanel1";
           TitleName="meResYvsAlphaZpPanel1";
           ST="";

           totalMEs->meResYvsAlphaZpPanel1 = bookME2D(ibooker, "TH2ResYvsAlphaZpPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchResXvsBetaZpPanel1Total){

           HistoName="meResXvsBetaZpPanel1";
           TitleName="meResXvsBetaZpPanel1";
           ST="";

           totalMEs->meResXvsBetaZpPanel1 = bookME2D(ibooker, "TH2ResXvsBetaZpPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchResYvsBetaZpPanel1Total){

           HistoName="meResYvsBetaZpPanel1";
           TitleName="meResYvsBetaZpPanel1";
           ST="";

           totalMEs->meResYvsBetaZpPanel1 = bookME2D(ibooker, "TH2ResYvsBetaZpPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     } 

    //Zp-Panel1

    if (switchPullXvsAlphaZpPanel1Total){

           HistoName="mePullXvsAlphaZpPanel1";
           TitleName="mePullXvsAlphaZpPanel1";
           ST="";

           totalMEs->mePullXvsAlphaZpPanel1=bookME2D(ibooker, "TH2PullXvsAlphaZpPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullYvsAlphaZpPanel1Total){

           HistoName="mePullYvsAlphaZpPanel1";
           TitleName="mePullYvsAlphaZpPanel1";
           ST="";

           totalMEs->mePullYvsAlphaZpPanel1=bookME2D(ibooker, "TH2PullYvsAlphaZpPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullXvsBetaZpPanel1Total){

           HistoName="mePullXvsBetaZpPanel1";
           TitleName="mePullXvsBetaZpPanel1";
           ST="";

           totalMEs->mePullXvsBetaZpPanel1=bookME2D(ibooker, "TH2PullXvsBetaZpPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullYvsBetaZpPanel1Total){

           HistoName="mePullYvsBetaZpPanel1";
           TitleName="mePullYvsBetaZpPanel1";
           ST="";

           totalMEs->mePullYvsBetaZpPanel1=bookME2D(ibooker, "TH2PullYvsBetaZpPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullXvsPhiZpPanel1Total){

           HistoName="mePullXvsPhiZpPanel1";
           TitleName="mePullXvsPhiZpPanel1";
           ST="";

           totalMEs->mePullXvsPhiZpPanel1=bookME2D(ibooker, "TH2PullXvsPhiZpPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullYvsPhiZpPanel1Total){

           HistoName="mePullYvsPhiZpPanel1";
           TitleName="mePullYvsPhiZpPanel1";
           ST="";

           totalMEs->mePullYvsPhiZpPanel1=bookME2D(ibooker, "TH2PullYvsPhiZpPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullXvsEtaZpPanel1Total){

           HistoName="mePullXvsEtaZpPanel1";
           TitleName="mePullXvsEtaZpPanel1";
           ST="";

           totalMEs->mePullXvsEtaZpPanel1=bookME2D(ibooker, "TH2PullXvsEtaZpPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullYvsEtaZpPanel1Total){

           HistoName="mePullYvsEtaZpPanel1";
           TitleName="mePullYvsEtaZpPanel1";
           ST="";

           totalMEs->mePullYvsEtaZpPanel1=bookME2D(ibooker, "TH2PullYvsEtaZpPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
 
     if (switchWPullXvsAlphaZpPanel1Total){

           HistoName="meWPullXvsAlphaZpPanel1";
           TitleName="meWPullXvsAlphaZpPanel1";
           ST="";

           totalMEs->meWPullXvsAlphaZpPanel1=bookME2D(ibooker, "TH2WPullXvsAlphaZpPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchWPullYvsAlphaZpPanel1Total){

           HistoName="meWPullYvsAlphaZpPanel1";
           TitleName="meWPullYvsAlphaZpPanel1";
           ST="";

           totalMEs->meWPullYvsAlphaZpPanel1=bookME2D(ibooker, "TH2WPullYvsAlphaZpPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchWPullXvsBetaZpPanel1Total){

           HistoName="meWPullXvsBetaZpPanel1";
           TitleName="meWPullXvsBetaZpPanel1";
           ST="";

           totalMEs->meWPullXvsBetaZpPanel1=bookME2D(ibooker, "TH2WPullXvsBetaZpPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchWPullYvsBetaZpPanel1Total){

           HistoName="meWPullYvsBetaZpPanel1";
           TitleName="meWPullYvsBetaZpPanel1";
           ST="";

           totalMEs->meWPullYvsBetaZpPanel1=bookME2D(ibooker, "TH2WPullYvsBetaZpPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }

    if(switchPosxZmPanel2Total){

          HistoName="mePosxZmPanel2";
          TitleName="mePosxZmPanel2";

          totalMEs->mePosxZmPanel2 = bookME1D(ibooker, "TH1PosxZmPanel2", HistoName.c_str(), TitleName.c_str());

    }

    if(switchPosyZmPanel2Total){

          HistoName="mePosyZmPanel2";
          TitleName="mePosyZmPanel2";

          totalMEs->mePosyZmPanel2 = bookME1D(ibooker, "TH1PosyZmPanel2", HistoName.c_str(), TitleName.c_str());

    }

    if(switchErrxZmPanel2Total){

          HistoName="meErrxZmPanel2";
          TitleName="meErrxZmPanel2";

          totalMEs->meErrxZmPanel2 = bookME1D(ibooker, "TH1ErrxZmPanel2", HistoName.c_str(), TitleName.c_str());

    }

    if(switchErryZmPanel2Total){

          HistoName="meErryZmPanel2";
          TitleName="meErryZmPanel2";

          totalMEs->meErryZmPanel2 = bookME1D(ibooker, "TH1ErryZmPanel2", HistoName.c_str(), TitleName.c_str());

    }

    if(switchResxZmPanel2Total){

          HistoName="meResxZmPanel2";
          TitleName="meResxZmPanel2";

          totalMEs->meResxZmPanel2 = bookME1D(ibooker, "TH1ResxZmPanel2", HistoName.c_str(), TitleName.c_str());

    }

    if(switchResyZmPanel2Total){

          HistoName="meResyZmPanel2";
          TitleName="meResyZmPanel2";

          totalMEs->meResyZmPanel2 = bookME1D(ibooker, "TH1ResyZmPanel2", HistoName.c_str(), TitleName.c_str());

    }

    if(switchPullxZmPanel2Total){

          HistoName="mePullxZmPanel2";
          TitleName="mePullxZmPanel2";

          totalMEs->mePullxZmPanel2 = bookME1D(ibooker, "TH1PullxZmPanel2", HistoName.c_str(), TitleName.c_str());

    }

    if(switchPullyZmPanel2Total){

          HistoName="mePullyZmPanel2";
          TitleName="mePullyZmPanel2";

          totalMEs->mePullyZmPanel2 = bookME1D(ibooker, "TH1PullyZmPanel2", HistoName.c_str(), TitleName.c_str());

    }

    if(switchNpixZmPanel2Total){

          HistoName="meNpixZmPanel2";
          TitleName="meNpixZmPanel2";

          totalMEs->meNpixZmPanel2 = bookME1D(ibooker, "TH1NpixZmPanel2", HistoName.c_str(), TitleName.c_str());

    }

    if(switchNxpixZmPanel2Total){

          HistoName="meNxpixZmPanel2";
          TitleName="meNxpixZmPanel2";

          totalMEs->meNxpixZmPanel2 = bookME1D(ibooker, "TH1NxpixZmPanel2", HistoName.c_str(), TitleName.c_str());

    }

    if(switchNypixZmPanel2Total){

          HistoName="meNypixZmPanel2";
          TitleName="meNypixZmPanel2";

          totalMEs->meNypixZmPanel2 = bookME1D(ibooker, "TH1NypixZmPanel2", HistoName.c_str(), TitleName.c_str());

    }

    if(switchChargeZmPanel2Total){

          HistoName="meChargeZmPanel2";
          TitleName="meChargeZmPanel2";

          totalMEs->meChargeZmPanel2 = bookME1D(ibooker, "TH1ChargeZmPanel2", HistoName.c_str(), TitleName.c_str());

    }

    if (switchResXvsAlphaZmPanel2Total){

           HistoName="meResXvsAlphaZmPanel2";
           TitleName="meResXvsAlphaZmPanel2";
           ST="";

           totalMEs->meResXvsAlphaZmPanel2=bookME2D(ibooker, "TH2ResXvsAlphaZmPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }

    if (switchResYvsAlphaZmPanel2Total){

           HistoName="meResYvsAlphaZmPanel2";
           TitleName="meResYvsAlphaZmPanel2";
           ST="";

           totalMEs->meResYvsAlphaZmPanel2=bookME2D(ibooker, "TH2ResYvsAlphaZmPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }

  
     if (switchResXvsBetaZmPanel2Total){

           HistoName="meResXvsBetaZmPanel2";
           TitleName="meResXvsBetaZmPanel2";
           ST="";

           totalMEs->meResXvsBetaZmPanel2=bookME2D(ibooker, "TH2ResXvsBetaZmPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }

    if (switchResYvsBetaZmPanel2Total){

           HistoName="meResYvsBetaZmPanel2";
           TitleName="meResYvsBetaZmPanel2";
           ST="";

           totalMEs->meResYvsBetaZmPanel2=bookME2D(ibooker, "TH2ResYvsBetaZmPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }

     if (switchPullXvsAlphaZmPanel2Total){

           HistoName="mePullXvsAlphaZmPanel2";
           TitleName="mePullXvsAlphaZmPanel2";
           ST="";

           totalMEs->mePullXvsAlphaZmPanel2=bookME2D(ibooker, "TH2PullXvsAlphaZmPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }

     if (switchPullYvsAlphaZmPanel2Total){

           HistoName="mePullYvsAlphaZmPanel2";
           TitleName="mePullYvsAlphaZmPanel2";
           ST="";

           totalMEs->mePullYvsAlphaZmPanel2=bookME2D(ibooker, "TH2PullYvsAlphaZmPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }

     if (switchPullXvsBetaZmPanel2Total){

           HistoName="mePullXvsBetaZmPanel2";
           TitleName="mePullXvsBetaZmPanel2";
           ST="";

           totalMEs->mePullXvsBetaZmPanel2=bookME2D(ibooker, "TH2PullXvsBetaZmPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }

     if (switchPullYvsBetaZmPanel2Total){

           HistoName="mePullYvsBetaZmPanel2";
           TitleName="mePullYvsBetaZmPanel2";
           ST="";

           totalMEs->mePullYvsBetaZmPanel2=bookME2D(ibooker, "TH2PullYvsBetaZmPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }

      if (switchPullXvsPhiZmPanel2Total){

           HistoName="mePullXvsPhiZmPanel2";
           TitleName="mePullXvsPhiZmPanel2";
           ST="";

           totalMEs->mePullXvsPhiZmPanel2=bookME2D(ibooker, "TH2PullXvsPhiZmPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }

     if (switchPullYvsPhiZmPanel2Total){

           HistoName="mePullYvsPhiZmPanel2";
           TitleName="mePullYvsPhiZmPanel2";
           ST="";

           totalMEs->mePullYvsPhiZmPanel2=bookME2D(ibooker, "TH2PullYvsPhiZmPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }

      if (switchPullXvsEtaZmPanel2Total){

           HistoName="mePullXvsEtaZmPanel2";
           TitleName="mePullXvsEtaZmPanel2";
           ST="";

           totalMEs->mePullXvsEtaZmPanel2=bookME2D(ibooker, "TH2PullXvsEtaZmPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }

     if (switchPullYvsEtaZmPanel2Total){

           HistoName="mePullYvsEtaZmPanel2";
           TitleName="mePullYvsEtaZmPanel2";
           ST="";

           totalMEs->mePullYvsEtaZmPanel2=bookME2D(ibooker, "TH2PullYvsEtaZmPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }


    if(switchPosxZpPanel2Total){

          HistoName="mePosxZpPanel2";
          TitleName="mePosxZpPanel2";

          totalMEs->mePosxZpPanel2 = bookME1D(ibooker, "TH1PosxZpPanel2", HistoName.c_str(), TitleName.c_str());

    }

    if(switchPosyZpPanel2Total){

          HistoName="mePosyZpPanel2";
          TitleName="mePosyZpPanel2";

          totalMEs->mePosyZpPanel2 = bookME1D(ibooker, "TH1PosyZpPanel2", HistoName.c_str(), TitleName.c_str());

    }

    if(switchErrxZpPanel2Total){

          HistoName="meErrxZpPanel2";
          TitleName="meErrxZpPanel2";

          totalMEs->meErrxZpPanel2 = bookME1D(ibooker, "TH1ErrxZpPanel2", HistoName.c_str(), TitleName.c_str());

    }

    if(switchErryZpPanel2Total){

          HistoName="meErryZpPanel2";
          TitleName="meErryZpPanel2";

          totalMEs->meErryZpPanel2 = bookME1D(ibooker, "TH1ErryZpPanel2", HistoName.c_str(), TitleName.c_str());

    }

    if(switchResxZpPanel2Total){

          HistoName="meResxZpPanel2";
          TitleName="meResxZpPanel2";

          totalMEs->meResxZpPanel2 = bookME1D(ibooker, "TH1ResxZpPanel2", HistoName.c_str(), TitleName.c_str());

    }

    if(switchResyZpPanel2Total){

          HistoName="meResyZpPanel2";
          TitleName="meResyZpPanel2";

          totalMEs->meResyZpPanel2 = bookME1D(ibooker, "TH1ResyZpPanel2", HistoName.c_str(), TitleName.c_str());

    }

    if(switchPullxZpPanel2Total){

          HistoName="mePullxZpPanel2";
          TitleName="mePullxZpPanel2";

          totalMEs->mePullxZpPanel2 = bookME1D(ibooker, "TH1PullxZpPanel2", HistoName.c_str(), TitleName.c_str());

    }

    if(switchPullyZpPanel2Total){

          HistoName="mePullyZpPanel2";
          TitleName="mePullyZpPanel2";

          totalMEs->mePullyZpPanel2 = bookME1D(ibooker, "TH1PullyZpPanel2", HistoName.c_str(), TitleName.c_str());

    }

    if(switchNpixZpPanel2Total){

          HistoName="meNpixZpPanel2";
          TitleName="meNpixZpPanel2";

          totalMEs->meNpixZpPanel2 = bookME1D(ibooker, "TH1NpixZpPanel2", HistoName.c_str(), TitleName.c_str());

    }

    if(switchNxpixZpPanel2Total){

          HistoName="meNxpixZpPanel2";
          TitleName="meNxpixZpPanel2";

          totalMEs->meNxpixZpPanel2 = bookME1D(ibooker, "TH1NxpixZpPanel2", HistoName.c_str(), TitleName.c_str());

    }

    if(switchNypixZpPanel2Total){

          HistoName="meNypixZpPanel2";
          TitleName="meNypixZpPanel2";

          totalMEs->meNypixZpPanel2 = bookME1D(ibooker, "TH1NypixZpPanel2", HistoName.c_str(), TitleName.c_str());

    }

    if(switchChargeZpPanel2Total){

          HistoName="meChargeZpPanel2";
          TitleName="meChargeZpPanel2";

          totalMEs->meChargeZpPanel2 = bookME1D(ibooker, "TH1ChargeZpPanel2", HistoName.c_str(), TitleName.c_str());

    }


    if (switchWPullXvsAlphaZmPanel2Total){

           HistoName="meWPullXvsAlphaZmPanel2";
           TitleName="meWPullXvsAlphaZmPanel2";
           ST="";

           totalMEs->meWPullXvsAlphaZmPanel2=bookME2D(ibooker, "TH2WPullXvsAlphaZmPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }


     if (switchWPullYvsAlphaZmPanel2Total){

           HistoName="meWPullYvsAlphaZmPanel2";
           TitleName="meWPullYvsAlphaZmPanel2";
           ST="";

           totalMEs->meWPullYvsAlphaZmPanel2=bookME2D(ibooker, "TH2WPullYvsAlphaZmPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }

 
     if (switchWPullXvsBetaZmPanel2Total){

           HistoName="meWPullXvsBetaZmPanel2";
           TitleName="meWPullXvsBetaZmPanel2";
           ST="";

           totalMEs->meWPullXvsBetaZmPanel2=bookME2D(ibooker, "TH2WPullXvsBetaZmPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }


     if (switchWPullYvsBetaZmPanel2Total){

           HistoName="meWPullYvsBetaZmPanel2";
           TitleName="meWPullYvsBetaZmPanel2";
           ST="";

           totalMEs->meWPullYvsBetaZmPanel2=bookME2D(ibooker, "TH2WPullYvsBetaZmPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }


     if (switchResXvsAlphaZpPanel2Total){

           HistoName="meResXvsAlphaZpPanel2";
           TitleName="meResXvsAlphaZpPanel2";
           ST="";

           totalMEs->meResXvsAlphaZpPanel2=bookME2D(ibooker, "TH2ResXvsAlphaZpPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }

     if (switchResYvsAlphaZpPanel2Total){

           HistoName="meResYvsAlphaZpPanel2";
           TitleName="meResYvsAlphaZpPanel2";
           ST="";

           totalMEs->meResYvsAlphaZpPanel2=bookME2D(ibooker, "TH2ResYvsAlphaZpPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }

     if (switchResXvsBetaZpPanel2Total){

           HistoName="meResXvsBetaZpPanel2";
           TitleName="meResXvsBetaZpPanel2";
           ST="";

           totalMEs->meResXvsBetaZpPanel2=bookME2D(ibooker, "TH2ResXvsBetaZpPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }

     if (switchResYvsBetaZpPanel2Total){

           HistoName="meResYvsBetaZpPanel2";
           TitleName="meResYvsBetaZpPanel2";
           ST="";

           totalMEs->meResYvsBetaZpPanel2=bookME2D(ibooker, "TH2ResYvsBetaZpPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }



     if (switchPullXvsAlphaZpPanel2Total){

           HistoName="mePullXvsAlphaZpPanel2";
           TitleName="mePullXvsAlphaZpPanel2";
           ST="";

           totalMEs->mePullXvsAlphaZpPanel2=bookME2D(ibooker, "TH2PullXvsAlphaZpPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }

     if (switchPullYvsAlphaZpPanel2Total){

           HistoName="mePullYvsAlphaZpPanel2";
           TitleName="mePullYvsAlphaZpPanel2";
           ST="";

           totalMEs->mePullYvsAlphaZpPanel2=bookME2D(ibooker, "TH2PullYvsAlphaZpPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }

     if (switchPullXvsBetaZpPanel2Total){

           HistoName="mePullXvsBetaZpPanel2";
           TitleName="mePullXvsBetaZpPanel2";
           ST="";

           totalMEs->mePullXvsBetaZpPanel2=bookME2D(ibooker, "TH2PullXvsBetaZpPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }

     if (switchPullYvsBetaZpPanel2Total){

           HistoName="mePullYvsBetaZpPanel2";
           TitleName="mePullYvsBetaZpPanel2";
           ST="";

           totalMEs->mePullYvsBetaZpPanel2=bookME2D(ibooker, "TH2PullYvsBetaZpPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }

     if (switchPullXvsPhiZpPanel2Total){

           HistoName="mePullXvsPhiZpPanel2";
           TitleName="mePullXvsPhiZpPanel2";
           ST="";

           totalMEs->mePullXvsPhiZpPanel2=bookME2D(ibooker, "TH2PullXvsPhiZpPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }

     if (switchPullYvsPhiZpPanel2Total){

           HistoName="mePullYvsPhiZpPanel2";
           TitleName="mePullYvsPhiZpPanel2";
           ST="";

           totalMEs->mePullYvsPhiZpPanel2=bookME2D(ibooker, "TH2PullYvsPhiZpPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }

     if (switchPullXvsEtaZpPanel2Total){

           HistoName="mePullXvsEtaZpPanel2";
           TitleName="mePullXvsEtaZpPanel2";
           ST="";

           totalMEs->mePullXvsEtaZpPanel2=bookME2D(ibooker, "TH2PullXvsEtaZpPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }

     if (switchPullYvsEtaZpPanel2Total){

           HistoName="mePullYvsEtaZpPanel2";
           TitleName="mePullYvsEtaZpPanel2";
           ST="";

           totalMEs->mePullYvsEtaZpPanel2=bookME2D(ibooker, "TH2PullYvsEtaZpPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }

     if (switchPosxBarrel_all_hits){

           HistoName="mePosxBarrel_all_hits";
           TitleName="mePosxBarrel_all_hits";

           totalMEs->mePosxBarrel_all_hits = bookME1D(ibooker, "TH1PosxBarrel_all_hits", HistoName.c_str(), TitleName.c_str());

     }

     if (switchPosyBarrel_all_hits){

           HistoName="mePosyBarrel_all_hits";
           TitleName="mePosyBarrel_all_hits";

           totalMEs->mePosyBarrel_all_hits = bookME1D(ibooker, "TH1PosyBarrel_all_hits", HistoName.c_str(), TitleName.c_str());

     }

     if (switchPosxZmPanel1_all_hits){

           HistoName="mePosxZmPanel1_all_hits";
           TitleName="mePosxZmPanel1_all_hits";

           totalMEs->mePosxZmPanel1_all_hits = bookME1D(ibooker, "TH1PosxZmPanel1_all_hits", HistoName.c_str(), TitleName.c_str());

     }

     if (switchPosyZmPanel1_all_hits){

           HistoName="mePosyZmPanel1_all_hits";
           TitleName="mePosyZmPanel1_all_hits";

           totalMEs->mePosyZmPanel1_all_hits = bookME1D(ibooker, "TH1PosyZmPanel1_all_hits", HistoName.c_str(), TitleName.c_str());

     }

     if (switchPosxZmPanel2_all_hits){

           HistoName="mePosxZmPanel2_all_hits";
           TitleName="mePosxZmPanel2_all_hits";

           totalMEs->mePosxZmPanel2_all_hits = bookME1D(ibooker, "TH1PosxZmPanel2_all_hits", HistoName.c_str(), TitleName.c_str());

     }

     if (switchPosyZmPanel2_all_hits){

           HistoName="mePosyZmPanel2_all_hits";
           TitleName="mePosyZmPanel2_all_hits";

           totalMEs->mePosyZmPanel2_all_hits = bookME1D(ibooker, "TH1PosyZmPanel2_all_hits", HistoName.c_str(), TitleName.c_str());

     }
   
      if (switchPosxZpPanel1_all_hits){

           HistoName="mePosxZpPanel1_all_hits";
           TitleName="mePosxZpPanel1_all_hits";

           totalMEs->mePosxZpPanel1_all_hits = bookME1D(ibooker, "TH1PosxZpPanel1_all_hits", HistoName.c_str(), TitleName.c_str());

     }

     if (switchPosyZpPanel1_all_hits){

           HistoName="mePosyZpPanel1_all_hits";
           TitleName="mePosyZpPanel1_all_hits";

           totalMEs->mePosyZpPanel1_all_hits = bookME1D(ibooker, "TH1PosyZpPanel1_all_hits", HistoName.c_str(), TitleName.c_str());

     }

     if (switchPosxZpPanel2_all_hits){

           HistoName="mePosxZpPanel2_all_hits";
           TitleName="mePosxZpPanel2_all_hits";

           totalMEs->mePosxZpPanel2_all_hits = bookME1D(ibooker, "TH1PosxZpPanel2_all_hits", HistoName.c_str(), TitleName.c_str());

     }

     if (switchPosyZpPanel2_all_hits){

           HistoName="mePosyZpPanel2_all_hits";
           TitleName="mePosyZpPanel2_all_hits";

           totalMEs->mePosyZpPanel2_all_hits = bookME1D(ibooker, "TH1PosyZpPanel2_all_hits", HistoName.c_str(), TitleName.c_str());

     }

     if (switchTracksPerEvent){

           HistoName="meTracksPerEvent";
           TitleName="meTracksPerEvent";

           totalMEs->meTracksPerEvent = bookME1D(ibooker, "TH1TracksPerEvent", HistoName.c_str(), TitleName.c_str());

     }

     if (switchPixRecHitsPerTrack){

           HistoName="mePixRecHitsPerTrack";
           TitleName="mePixRecHitsPerTrack";

           totalMEs->mePixRecHitsPerTrack = bookME1D(ibooker, "TH1PixRecHitsPerTrack", HistoName.c_str(), TitleName.c_str());

     }


}

void SiPixelTrackingRecHitsValid::createLadderLayersMEs(DQMStore::IBooker &ibooker, std::string label){

        LadderLayersMEs *LadderLayerObject = new LadderLayersMEs();

        std::string HistoName;
        std::string TitleName;
        std::string ST;

        //Flipped
        if (switchResXvsAlphaBarrelFlippedLaddersLayer){

             HistoName="meResXvsAlphaBarrel"+label;
             TitleName="meResXvsAlphaBarrel"+label;
             ST="";

             LadderLayerObject->meResXvsAlphaBarrel=bookME2D(ibooker, "TH2ResXvsAlphaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

        }
    
        if (switchResYvsAlphaBarrelFlippedLaddersLayer){

             HistoName="meResYvsAlphaBarrel"+label;
             TitleName="meResYvsAlphaBarrel"+label;
             ST="";

             LadderLayerObject->meResYvsAlphaBarrel=bookME2D(ibooker, "TH2ResYvsAlphaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

        }

        if (switchResXvsBetaBarrelFlippedLaddersLayer){

             HistoName="meResXvsBetaBarrel"+label;
             TitleName="meResXvsBetaBarrel"+label;
             ST="";

             LadderLayerObject->meResXvsBetaBarrel=bookME2D(ibooker, "TH2ResXvsBetaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

        }

        if (switchResYvsBetaBarrelFlippedLaddersLayer){

             HistoName="meResYvsBetaBarrel"+label;
             TitleName="meResYvsBetaBarrel"+label;
             ST="";

             LadderLayerObject->meResYvsBetaBarrel=bookME2D(ibooker, "TH2ResYvsBetaBarrelFlippedLaddersLayer", HistoName.c_str(), TitleName.c_str(), ST.c_str());

        }

        //Non-Flipped

        if (switchResXvsAlphaBarrelNonFlippedLaddersLayer){

             HistoName="meResXvsAlphaBarrel"+label;
             TitleName="meResXvsAlphaBarrel"+label;
             ST="";

             LadderLayerObject->meResXvsAlphaBarrel=bookME2D(ibooker, "TH2ResXvsAlphaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

        }

        if (switchResYvsAlphaBarrelNonFlippedLaddersLayer){

             HistoName="meResYvsAlphaBarrel"+label;
             TitleName="meResYvsAlphaBarrel"+label;
             ST="";

             LadderLayerObject->meResYvsAlphaBarrel=bookME2D(ibooker, "TH2ResYvsAlphaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

        }

        if (switchResXvsBetaBarrelNonFlippedLaddersLayer){

             HistoName="meResXvsBetaBarrel"+label;
             TitleName="meResXvsBetaBarrel"+label;
             ST="";

             LadderLayerObject->meResXvsBetaBarrel=bookME2D(ibooker, "TH2ResXvsBetaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

        }

        if (switchResYvsBetaBarrelNonFlippedLaddersLayer){

             HistoName="meResYvsBetaBarrel"+label;
             TitleName="meResYvsBetaBarrel"+label;
             ST="";

             LadderLayerObject->meResYvsBetaBarrel=bookME2D(ibooker, "TH2ResYvsBetaBarrelNonFlippedLaddersLayer", HistoName.c_str(), TitleName.c_str(), ST.c_str());

        }


       	ladderLayersMEsMap[label] = LadderLayerObject;

}//END--------createLadderLayersMEs


void SiPixelTrackingRecHitsValid::createLayersMEs(DQMStore::IBooker &ibooker, std::string label){

      LayerMEs* LayerObject=new LayerMEs();

      std::string HistoName;
      std::string TitleName;

      if (switchResxBarrelLayer){

           HistoName="meResxBarrel"+label;
           TitleName="meResxBarrel"+label;

           LayerObject->meResxBarrel=bookME1D(ibooker, "TH1ResxBarrel",HistoName.c_str(), TitleName.c_str());

      }

      if (switchResyBarrelLayer){

           HistoName="meResyBarrel"+label;
           TitleName="meResyBarrel"+label;

           LayerObject->meResyBarrel=bookME1D(ibooker, "TH1ResyBarrel",HistoName.c_str(), TitleName.c_str());

      }

      if (switchPullxBarrelLayer){

           HistoName="mePullxBarrel"+label;
           TitleName="mePullxBarrel"+label;

           LayerObject->mePullxBarrel=bookME1D(ibooker, "TH1PullxBarrel",HistoName.c_str(), TitleName.c_str());

      }
      if (switchPullyBarrelLayer){

           HistoName="mePullyBarrel"+label;
           TitleName="mePullyBarrel"+label;

           LayerObject->mePullyBarrel=bookME1D(ibooker, "TH1PullyBarrel",HistoName.c_str(), TitleName.c_str());

      }

      layerMEsMap[label] = LayerObject;

}


void SiPixelTrackingRecHitsValid::createLayerModulesMEs(DQMStore::IBooker &ibooker, std::string label){

    LayerModulesMEs* LayerModulesObject=new LayerModulesMEs();
    
    std::string HistoName;
    std::string TitleName;
    std::string ST;
    
    if (switchPosxBarrelLayerModule){

           HistoName="mePosxBarrel"+label;
           TitleName="mePosxBarrel"+label;

           LayerModulesObject->mePosxBarrel=bookME1D(ibooker, "TH1PosxBarrel",HistoName.c_str(), TitleName.c_str());

    }  

    if (switchPosyBarrelLayerModule){

           HistoName="mePosyBarrel"+label;
           TitleName="mePosyBarrel"+label;

           LayerModulesObject->mePosxBarrel=bookME1D(ibooker, "TH1PosyBarrel",HistoName.c_str(), TitleName.c_str());

    }

    if (switchErrxBarrelLayerModule){

           HistoName="meErrxBarrel"+label;
           TitleName="meErrxBarrel"+label;

           LayerModulesObject->meErrxBarrel=bookME1D(ibooker, "TH1ErrxBarrel",HistoName.c_str(), TitleName.c_str());

    }

    if (switchErryBarrelLayerModule){

           HistoName="meErryBarrel"+label;
           TitleName="meErryBarrel"+label;

           LayerModulesObject->meErrxBarrel=bookME1D(ibooker, "TH1ErryBarrel",HistoName.c_str(), TitleName.c_str());

    }

     if (switchResxBarrelLayerModule){

           HistoName="meResxBarrel"+label;
           TitleName="meResxBarrel"+label;

           LayerModulesObject->meResxBarrel=bookME1D(ibooker, "TH1ResxBarrel",HistoName.c_str(), TitleName.c_str());

     }
   
     if (switchResyBarrelLayerModule){

           HistoName="meResyBarrel"+label;
           TitleName="meResyBarrel"+label;

           LayerModulesObject->meResyBarrel=bookME1D(ibooker, "TH1ResyBarrel",HistoName.c_str(), TitleName.c_str());

     }

     if (switchPullxBarrelLayerModule){

           HistoName="mePullxBarrel"+label;
           TitleName="mePullxBarrel"+label;

           LayerModulesObject->mePullxBarrel=bookME1D(ibooker, "TH1PullxBarrel",HistoName.c_str(), TitleName.c_str());

     }

     if (switchPullyBarrelLayerModule){

           HistoName="mePullyBarrel"+label;
           TitleName="mePullyBarrel"+label;

           LayerModulesObject->mePullyBarrel=bookME1D(ibooker, "TH1PullyBarrel",HistoName.c_str(), TitleName.c_str());

     }

      if (switchNpixBarrelLayerModule){

           HistoName="meNpixBarrel"+label;
           TitleName="meNpixBarrel"+label;

           LayerModulesObject->meNpixBarrel=bookME1D(ibooker, "TH1NpixBarrel",HistoName.c_str(), TitleName.c_str());

     }
   
     if (switchNxpixBarrelLayerModule){

           HistoName="meNxpixBarrel"+label;
           TitleName="meNxpixBarrel"+label;

           LayerModulesObject->meNxpixBarrel=bookME1D(ibooker, "TH1NxpixBarrel",HistoName.c_str(), TitleName.c_str());

     }

     if (switchNypixBarrelLayerModule){

           HistoName="meNypixBarrel"+label;
           TitleName="meNypixBarrel"+label;

           LayerModulesObject->meNypixBarrel=bookME1D(ibooker, "TH1NypixBarrel",HistoName.c_str(), TitleName.c_str());

     }

     if (switchChargeBarrelLayerModule){

           HistoName="meChargeBarrel"+label;
           TitleName="meChargeBarrel"+label;

           LayerModulesObject->meNypixBarrel=bookME1D(ibooker, "TH1ChargeBarrel",HistoName.c_str(), TitleName.c_str());

     }

     if (switchResXvsAlphaBarrelLayerModule){

           HistoName="meResXvsAlphaBarrel"+label;
           TitleName="meResXvsAlphaBarrel"+label;
           ST="";

           LayerModulesObject->meResXvsAlphaBarrel=bookME2D(ibooker, "TH2ResXvsAlphaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchResYvsAlphaBarrelLayerModule){

           HistoName="meResYvsAlphaBarrel"+label;
           TitleName="meResYvsAlphaBarrel"+label;
           ST="";

           LayerModulesObject->meResYvsAlphaBarrel=bookME2D(ibooker, "TH2ResYvsAlphaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }


     if (switchResXvsBetaBarrelLayerModule){

           HistoName="meResXvsBetaBarrel"+label;
           TitleName="meResXvsBetaBarrel"+label;
           ST="";

           LayerModulesObject->meResXvsBetaBarrel=bookME2D(ibooker, "TH2ResXvsBetaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchResYvsBetaBarrelLayerModule){

           HistoName="meResYvsBetaBarrel"+label;
           TitleName="meResYvsBetaBarrel"+label;
           ST="";

           LayerModulesObject->meResYvsBetaBarrel=bookME2D(ibooker, "TH2ResYvsBetaBarrelLayerModule", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }


     if (switchPullXvsAlphaBarrelLayerModule){

           HistoName="mePullXvsAlphaBarrel"+label;
           TitleName="mePullXvsAlphaBarrel"+label;
           ST="";

           LayerModulesObject->mePullXvsAlphaBarrel=bookME2D(ibooker, "TH2PullXvsAlphaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchPullYvsAlphaBarrelLayerModule){

           HistoName="mePullYvsAlphaBarrel"+label;
           TitleName="mePullYvsAlphaBarrel"+label;
           ST="";

           LayerModulesObject->mePullYvsAlphaBarrel=bookME2D(ibooker, "TH2PullYvsAlphaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchPullXvsBetaBarrelLayerModule){

           HistoName="mePullXvsBetaBarrel"+label;
           TitleName="mePullXvsBetaBarrel"+label;
           ST="";

           LayerModulesObject->mePullXvsBetaBarrel=bookME2D(ibooker, "TH2PullXvsBetaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchPullYvsBetaBarrelLayerModule){

           HistoName="mePullYvsBetaBarrel"+label;
           TitleName="mePullYvsBetaBarrel"+label;
           ST="";

           LayerModulesObject->mePullYvsBetaBarrel=bookME2D(ibooker, "TH2PullYvsBetaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

    
     if (switchPullXvsPhiBarrelLayerModule){

           HistoName="mePullXvsPhiBarrel"+label;
           TitleName="mePullXvsPhiBarrel"+label;
           ST="";

           LayerModulesObject->mePullXvsPhiBarrel=bookME2D(ibooker, "TH2PullXvsPhiBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchPullYvsPhiBarrelLayerModule){

           HistoName="mePullYvsPhiBarrel"+label;
           TitleName="mePullYvsPhiBarrel"+label;
           ST="";

           LayerModulesObject->mePullYvsPhiBarrel=bookME2D(ibooker, "TH2PullYvsPhiBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }


     if (switchPullXvsEtaBarrelLayerModule){

           HistoName="mePullXvsEtaBarrel"+label;
           TitleName="mePullXvsEtaBarrel"+label;
           ST="";

           LayerModulesObject->mePullXvsEtaBarrel=bookME2D(ibooker, "TH2PullXvsEtaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }

     if (switchPullYvsEtaBarrelLayerModule){

           HistoName="mePullYvsEtaBarrel"+label;
           TitleName="mePullYvsEtaBarrel"+label;
           ST="";

           LayerModulesObject->mePullYvsEtaBarrel=bookME2D(ibooker, "TH2PullYvsEtaBarrel", HistoName.c_str(), TitleName.c_str(), ST.c_str());

     }


     layerModulesMEsMap[label]=LayerModulesObject; 

}

void SiPixelTrackingRecHitsValid::createDiskPlaquettesMEs(DQMStore::IBooker &ibooker, std::string label){

    DiskPlaquettesMEs* DiskPlaquettesObject=new DiskPlaquettesMEs();
    
    std::string HistoName;
    std::string TitleName;
    std::string ST;

    //////////////////Zm
    if (switchPosxZmPanel1DiskPlaquette){

           HistoName="mePosxZmPanel1"+label;
           TitleName="mePosxZmPanel1"+label;

           DiskPlaquettesObject->mePosxZmPanel1=bookME1D(ibooker, "TH1PosxZmPanel1",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchPosyZmPanel1DiskPlaquette){

           HistoName="mePosyZmPanel1"+label;
           TitleName="mePosyZmPanel1"+label;

           DiskPlaquettesObject->mePosyZmPanel1=bookME1D(ibooker, "TH1PosyZmPanel1",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchErrxZmPanel1DiskPlaquette){

           HistoName="meErrxZmPanel1"+label;
           TitleName="meErrxZmPanel1"+label;

           DiskPlaquettesObject->meErrxZmPanel1=bookME1D(ibooker, "TH1ErrxZmPanel1",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchErryZmPanel1DiskPlaquette){

           HistoName="meErryZmPanel1"+label;
           TitleName="meErryZmPanel1"+label;

           DiskPlaquettesObject->meErryZmPanel1=bookME1D(ibooker, "TH1ErryZmPanel1",HistoName.c_str(), TitleName.c_str());

    }  
    if (switchResxZmPanel1DiskPlaquette){

           HistoName="meResxZmPanel1"+label;
           TitleName="meResxZmPanel1"+label;

           DiskPlaquettesObject->meResxZmPanel1=bookME1D(ibooker, "TH1ResxZmPanel1",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchResyZmPanel1DiskPlaquette){

           HistoName="meResyZmPanel1"+label;
           TitleName="meResyZmPanel1"+label;

           DiskPlaquettesObject->meResyZmPanel1=bookME1D(ibooker, "TH1ResyZmPanel1",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchPullxZmPanel1DiskPlaquette){

           HistoName="mePullxZmPanel1"+label;
           TitleName="mePullxZmPanel1"+label;

           DiskPlaquettesObject->mePullxZmPanel1=bookME1D(ibooker, "TH1PullxZmPanel1",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchPullyZmPanel1DiskPlaquette){

           HistoName="mePullyZmPanel1"+label;
           TitleName="mePullyZmPanel1"+label;

           DiskPlaquettesObject->mePullyZmPanel1=bookME1D(ibooker, "TH1PullyZmPanel1",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchNpixZmPanel1DiskPlaquette){

           HistoName="meNpixZmPanel1"+label;
           TitleName="meNpixZmPanel1"+label;

           DiskPlaquettesObject->meNpixZmPanel1=bookME1D(ibooker, "TH1NpixZmPanel1",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchNxpixZmPanel1DiskPlaquette){

           HistoName="meNxpixZmPanel1"+label;
           TitleName="meNxpixZmPanel1"+label;

           DiskPlaquettesObject->meNxpixZmPanel1=bookME1D(ibooker, "TH1NxpixZmPanel1",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchNypixZmPanel1DiskPlaquette){

           HistoName="meNypixZmPanel1"+label;
           TitleName="meNypixZmPanel1"+label;

           DiskPlaquettesObject->meNypixZmPanel1=bookME1D(ibooker, "TH1NypixZmPanel1",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchChargeZmPanel1DiskPlaquette){

           HistoName="meChargeZmPanel1"+label;
           TitleName="meChargeZmPanel1"+label;

           DiskPlaquettesObject->meChargeZmPanel1=bookME1D(ibooker, "TH1ChargeZmPanel1",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchResXvsAlphaZmPanel1DiskPlaquette){

           HistoName="meResXvsAlphaZmPanel1"+label;
           TitleName="meResXvsAlphaZmPanel1"+label;
           ST="";

           DiskPlaquettesObject->meResXvsAlphaZmPanel1=bookME2D(ibooker, "TH2ResXvsAlphaZmPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }

    if (switchResYvsAlphaZmPanel1DiskPlaquette){

           HistoName="meResYvsAlphaZmPanel1"+label;
           TitleName="meResYvsAlphaZmPanel1"+label;
           ST="";

           DiskPlaquettesObject->meResYvsAlphaZmPanel1=bookME2D(ibooker, "TH2ResYvsAlphaZmPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchResXvsBetaZmPanel1DiskPlaquette){

           HistoName="meResXvsBetaZmPanel1"+label;
           TitleName="meResXvsBetaZmPanel1"+label;
           ST="";

           DiskPlaquettesObject->meResXvsBetaZmPanel1=bookME2D(ibooker, "TH2ResXvsBetaZmPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchResYvsBetaZmPanel1DiskPlaquette){

           HistoName="meResYvsBetaZmPanel1"+label;
           TitleName="meResYvsBetaZmPanel1"+label;
           ST="";

           DiskPlaquettesObject->meResYvsBetaZmPanel1=bookME2D(ibooker, "TH2ResYvsBetaZmPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullXvsAlphaZmPanel1DiskPlaquette){

           HistoName="mePullXvsAlphaZmPanel1"+label;
           TitleName="mePullXvsAlphaZmPanel1"+label;
           ST="";

           DiskPlaquettesObject->mePullXvsAlphaZmPanel1=bookME2D(ibooker, "TH2PullXvsAlphaZmPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullYvsAlphaZmPanel1DiskPlaquette){

           HistoName="mePullYvsAlphaZmPanel1"+label;
           TitleName="mePullYvsAlphaZmPanel1"+label;
           ST="";

           DiskPlaquettesObject->mePullYvsAlphaZmPanel1=bookME2D(ibooker, "TH2PullYvsAlphaZmPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullXvsBetaZmPanel1DiskPlaquette){

           HistoName="mePullXvsBetaZmPanel1"+label;
           TitleName="mePullXvsBetaZmPanel1"+label;
           ST="";

           DiskPlaquettesObject->mePullXvsBetaZmPanel1=bookME2D(ibooker, "TH2PullXvsBetaZmPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullYvsBetaZmPanel1DiskPlaquette){

           HistoName="mePullYvsBetaZmPanel1"+label;
           TitleName="mePullYvsBetaZmPanel1"+label;
           ST="";

           DiskPlaquettesObject->mePullYvsBetaZmPanel1=bookME2D(ibooker, "TH2PullYvsBetaZmPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullXvsPhiZmPanel1DiskPlaquette){

           HistoName="mePullXvsPhiZmPanel1"+label;
           TitleName="mePullXvsPhiZmPanel1"+label;
           ST="";

           DiskPlaquettesObject->mePullXvsPhiZmPanel1=bookME2D(ibooker, "TH2PullXvsPhiZmPanel1DiskPlaquette", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullYvsPhiZmPanel1DiskPlaquette){

           HistoName="mePullYvsPhiZmPanel1"+label;
           TitleName="mePullYvsPhiZmPanel1"+label;
           ST="";

           DiskPlaquettesObject->mePullYvsPhiZmPanel1=bookME2D(ibooker, "TH2PullYvsPhiZmPanel1DiskPlaquette", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullXvsEtaZmPanel1DiskPlaquette){

           HistoName="mePullXvsEtaZmPanel1"+label;
           TitleName="mePullXvsEtaZmPanel1"+label;
           ST="";

           DiskPlaquettesObject->mePullXvsEtaZmPanel1=bookME2D(ibooker, "TH2PullXvsEtaZmPanel1DiskPlaquette", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullYvsEtaZmPanel1DiskPlaquette){

           HistoName="mePullYvsEtaZmPanel1"+label;
           TitleName="mePullYvsEtaZmPanel1"+label;
           ST="";

           DiskPlaquettesObject->mePullYvsEtaZmPanel1=bookME2D(ibooker, "TH2PullYvsEtaZmPanel1DiskPlaquette", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
     /////////////////END-Zm

    //////////////////Zp-Panel1
    if (switchPosxZpPanel1DiskPlaquette){

           HistoName="mePosxZpPanel1"+label;
           TitleName="mePosxZpPanel1"+label;

           DiskPlaquettesObject->mePosxZpPanel1=bookME1D(ibooker, "TH1PosxZpPanel1",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchPosyZpPanel1DiskPlaquette){

           HistoName="mePosyZpPanel1"+label;
           TitleName="mePosyZpPanel1"+label;

           DiskPlaquettesObject->mePosyZpPanel1=bookME1D(ibooker, "TH1PosyZpPanel1",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchErrxZpPanel1DiskPlaquette){

           HistoName="meErrxZpPanel1"+label;
           TitleName="meErrxZpPanel1"+label;

           DiskPlaquettesObject->meErrxZpPanel1=bookME1D(ibooker, "TH1ErrxZpPanel1",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchErryZpPanel1DiskPlaquette){

           HistoName="meErryZpPanel1"+label;
           TitleName="meErryZpPanel1"+label;

           DiskPlaquettesObject->meErryZpPanel1=bookME1D(ibooker, "TH1ErryZpPanel1",HistoName.c_str(), TitleName.c_str());

    }  
    if (switchResxZpPanel1DiskPlaquette){

           HistoName="meResxZpPanel1"+label;
           TitleName="meResxZpPanel1"+label;

           DiskPlaquettesObject->meResxZpPanel1=bookME1D(ibooker, "TH1ResxZpPanel1",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchResyZpPanel1DiskPlaquette){

           HistoName="meResyZpPanel1"+label;
           TitleName="meResyZpPanel1"+label;

           DiskPlaquettesObject->meResyZpPanel1=bookME1D(ibooker, "TH1ResyZpPanel1",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchPullxZpPanel1DiskPlaquette){

           HistoName="mePullxZpPanel1"+label;
           TitleName="mePullxZpPanel1"+label;

           DiskPlaquettesObject->mePullxZpPanel1=bookME1D(ibooker, "TH1PullxZpPanel1",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchPullyZpPanel1DiskPlaquette){

           HistoName="mePullyZpPanel1"+label;
           TitleName="mePullyZpPanel1"+label;

           DiskPlaquettesObject->mePullyZpPanel1=bookME1D(ibooker, "TH1PullyZpPanel1",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchNpixZpPanel1DiskPlaquette){

           HistoName="meNpixZpPanel1"+label;
           TitleName="meNpixZpPanel1"+label;

           DiskPlaquettesObject->meNpixZpPanel1=bookME1D(ibooker, "TH1NpixZpPanel1",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchNxpixZpPanel1DiskPlaquette){

           HistoName="meNxpixZpPanel1"+label;
           TitleName="meNxpixZpPanel1"+label;

           DiskPlaquettesObject->meNxpixZpPanel1=bookME1D(ibooker, "TH1NxpixZpPanel1",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchNypixZpPanel1DiskPlaquette){

           HistoName="meNypixZpPanel1"+label;
           TitleName="meNypixZpPanel1"+label;

           DiskPlaquettesObject->meNypixZpPanel1=bookME1D(ibooker, "TH1NypixZpPanel1",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchChargeZpPanel1DiskPlaquette){

           HistoName="meChargeZpPanel1"+label;
           TitleName="meChargeZpPanel1"+label;

           DiskPlaquettesObject->meChargeZpPanel1=bookME1D(ibooker, "TH1ChargeZpPanel1",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchResXvsAlphaZpPanel1DiskPlaquette){

           HistoName="meResXvsAlphaZpPanel1"+label;
           TitleName="meResXvsAlphaZpPanel1"+label;
           ST="";

           DiskPlaquettesObject->meResXvsAlphaZpPanel1=bookME2D(ibooker, "TH2ResXvsAlphaZpPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }

    if (switchResYvsAlphaZpPanel1DiskPlaquette){

           HistoName="meResYvsAlphaZpPanel1"+label;
           TitleName="meResYvsAlphaZpPanel1"+label;
           ST="";

           DiskPlaquettesObject->meResYvsAlphaZpPanel1=bookME2D(ibooker, "TH2ResYvsAlphaZpPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchResXvsBetaZpPanel1DiskPlaquette){

           HistoName="meResXvsBetaZpPanel1"+label;
           TitleName="meResXvsBetaZpPanel1"+label;
           ST="";

           DiskPlaquettesObject->meResXvsBetaZpPanel1=bookME2D(ibooker, "TH2ResXvsBetaZpPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchResYvsBetaZpPanel1DiskPlaquette){

           HistoName="meResYvsBetaZpPanel1"+label;
           TitleName="meResYvsBetaZpPanel1"+label;
           ST="";

           DiskPlaquettesObject->meResYvsBetaZpPanel1=bookME2D(ibooker, "TH2ResYvsBetaZpPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullXvsAlphaZpPanel1DiskPlaquette){

           HistoName="mePullXvsAlphaZpPanel1"+label;
           TitleName="mePullXvsAlphaZpPanel1"+label;
           ST="";

           DiskPlaquettesObject->mePullXvsAlphaZpPanel1=bookME2D(ibooker, "TH2PullXvsAlphaZpPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullYvsAlphaZpPanel1DiskPlaquette){

           HistoName="mePullYvsAlphaZpPanel1"+label;
           TitleName="mePullYvsAlphaZpPanel1"+label;
           ST="";

           DiskPlaquettesObject->mePullYvsAlphaZpPanel1=bookME2D(ibooker, "TH2PullYvsAlphaZpPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullXvsBetaZpPanel1DiskPlaquette){

           HistoName="mePullXvsBetaZpPanel1"+label;
           TitleName="mePullXvsBetaZpPanel1"+label;
           ST="";

           DiskPlaquettesObject->mePullXvsBetaZpPanel1=bookME2D(ibooker, "TH2PullXvsBetaZpPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullYvsBetaZpPanel1DiskPlaquette){

           HistoName="mePullYvsBetaZpPanel1"+label;
           TitleName="mePullYvsBetaZpPanel1"+label;
           ST="";

           DiskPlaquettesObject->mePullYvsBetaZpPanel1=bookME2D(ibooker, "TH2PullYvsBetaZpPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullXvsPhiZpPanel1DiskPlaquette){

           HistoName="mePullXvsPhiZpPanel1"+label;
           TitleName="mePullXvsPhiZpPanel1"+label;
           ST="";

           DiskPlaquettesObject->mePullXvsPhiZpPanel1=bookME2D(ibooker, "TH2PullXvsPhiZpPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullYvsPhiZpPanel1DiskPlaquette){

           HistoName="mePullYvsPhiZpPanel1"+label;
           TitleName="mePullYvsPhiZpPanel1"+label;
           ST="";

           DiskPlaquettesObject->mePullYvsPhiZpPanel1=bookME2D(ibooker, "TH2PullYvsPhiZpPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullXvsEtaZpPanel1DiskPlaquette){

           HistoName="mePullXvsEtaZpPanel1"+label;
           TitleName="mePullXvsEtaZpPanel1"+label;
           ST="";

           DiskPlaquettesObject->mePullXvsEtaZpPanel1=bookME2D(ibooker, "TH2PullXvsEtaZpPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullYvsEtaZpPanel1DiskPlaquette){

           HistoName="mePullYvsEtaZpPanel1"+label;
           TitleName="mePullYvsEtaZpPanel1"+label;
           ST="";

           DiskPlaquettesObject->mePullYvsEtaZpPanel1=bookME2D(ibooker, "TH2PullYvsEtaZpPanel1", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
     //END-Zp 
    
     //Zm-Panel2
    if (switchPosxZmPanel2DiskPlaquette){

           HistoName="mePosxZmPanel2"+label;
           TitleName="mePosxZmPanel2"+label;

           DiskPlaquettesObject->mePosxZmPanel2=bookME1D(ibooker, "TH1PosxZmPanel2",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchPosyZmPanel2DiskPlaquette){

           HistoName="mePosyZmPanel2"+label;
           TitleName="mePosyZmPanel2"+label;

           DiskPlaquettesObject->mePosyZmPanel2=bookME1D(ibooker, "TH1PosyZmPanel2",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchErrxZmPanel2DiskPlaquette){

           HistoName="meErrxZmPanel2"+label;
           TitleName="meErrxZmPanel2"+label;

           DiskPlaquettesObject->meErrxZmPanel2=bookME1D(ibooker, "TH1ErrxZmPanel2",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchErryZmPanel2DiskPlaquette){

           HistoName="meErryZmPanel2"+label;
           TitleName="meErryZmPanel2"+label;

           DiskPlaquettesObject->meErryZmPanel2=bookME1D(ibooker, "TH1ErryZmPanel2",HistoName.c_str(), TitleName.c_str());

    }  
    if (switchResxZmPanel2DiskPlaquette){

           HistoName="meResxZmPanel2"+label;
           TitleName="meResxZmPanel2"+label;

           DiskPlaquettesObject->meResxZmPanel2=bookME1D(ibooker, "TH1ResxZmPanel2",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchResyZmPanel2DiskPlaquette){

           HistoName="meResyZmPanel2"+label;
           TitleName="meResyZmPanel2"+label;

           DiskPlaquettesObject->meResyZmPanel2=bookME1D(ibooker, "TH1ResyZmPanel2",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchPullxZmPanel2DiskPlaquette){

           HistoName="mePullxZmPanel2"+label;
           TitleName="mePullxZmPanel2"+label;

           DiskPlaquettesObject->mePullxZmPanel2=bookME1D(ibooker, "TH1PullxZmPanel2",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchPullyZmPanel2DiskPlaquette){

           HistoName="mePullyZmPanel2"+label;
           TitleName="mePullyZmPanel2"+label;

           DiskPlaquettesObject->mePullyZmPanel2=bookME1D(ibooker, "TH1PullyZmPanel2",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchNpixZmPanel2DiskPlaquette){

           HistoName="meNpixZmPanel2"+label;
           TitleName="meNpixZmPanel2"+label;

           DiskPlaquettesObject->meNpixZmPanel2=bookME1D(ibooker, "TH1NpixZmPanel2",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchNxpixZmPanel2DiskPlaquette){

           HistoName="meNxpixZmPanel2"+label;
           TitleName="meNxpixZmPanel2"+label;

           DiskPlaquettesObject->meNxpixZmPanel2=bookME1D(ibooker, "TH1NxpixZmPanel2",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchNypixZmPanel2DiskPlaquette){

           HistoName="meNypixZmPanel2"+label;
           TitleName="meNypixZmPanel2"+label;

           DiskPlaquettesObject->meNypixZmPanel2=bookME1D(ibooker, "TH1NypixZmPanel2",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchChargeZmPanel2DiskPlaquette){

           HistoName="meChargeZmPanel2"+label;
           TitleName="meChargeZmPanel2"+label;

           DiskPlaquettesObject->meChargeZmPanel2=bookME1D(ibooker, "TH1ChargeZmPanel2",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchResXvsAlphaZmPanel2DiskPlaquette){

           HistoName="meResXvsAlphaZmPanel2"+label;
           TitleName="meResXvsAlphaZmPanel2"+label;
           ST="";

           DiskPlaquettesObject->meResXvsAlphaZmPanel2=bookME2D(ibooker, "TH2ResXvsAlphaZmPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }

    if (switchResYvsAlphaZmPanel2DiskPlaquette){

           HistoName="meResYvsAlphaZmPanel2"+label;
           TitleName="meResYvsAlphaZmPanel2"+label;
           ST="";

           DiskPlaquettesObject->meResYvsAlphaZmPanel2=bookME2D(ibooker, "TH2ResYvsAlphaZmPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchResXvsBetaZmPanel2DiskPlaquette){

           HistoName="meResXvsBetaZmPanel2"+label;
           TitleName="meResXvsBetaZmPanel2"+label;
           ST="";

           DiskPlaquettesObject->meResXvsBetaZmPanel2=bookME2D(ibooker, "TH2ResXvsBetaZmPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchResYvsBetaZmPanel2DiskPlaquette){

           HistoName="meResYvsBetaZmPanel2"+label;
           TitleName="meResYvsBetaZmPanel2"+label;
           ST="";

           DiskPlaquettesObject->meResYvsBetaZmPanel2=bookME2D(ibooker, "TH2ResYvsBetaZmPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullXvsAlphaZmPanel2DiskPlaquette){

           HistoName="mePullXvsAlphaZmPanel2"+label;
           TitleName="mePullXvsAlphaZmPanel2"+label;
           ST="";

           DiskPlaquettesObject->mePullXvsAlphaZmPanel2=bookME2D(ibooker, "TH2PullXvsAlphaZmPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullYvsAlphaZmPanel2DiskPlaquette){

           HistoName="mePullYvsAlphaZmPanel2"+label;
           TitleName="mePullYvsAlphaZmPanel2"+label;
           ST="";

           DiskPlaquettesObject->mePullYvsAlphaZmPanel2=bookME2D(ibooker, "TH2PullYvsAlphaZmPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullXvsBetaZmPanel2DiskPlaquette){

           HistoName="mePullXvsBetaZmPanel2"+label;
           TitleName="mePullXvsBetaZmPanel2"+label;
           ST="";

           DiskPlaquettesObject->mePullXvsBetaZmPanel2=bookME2D(ibooker, "TH2PullXvsBetaZmPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullYvsBetaZmPanel2DiskPlaquette){

           HistoName="mePullYvsBetaZmPanel2"+label;
           TitleName="mePullYvsBetaZmPanel2"+label;
           ST="";

           DiskPlaquettesObject->mePullYvsBetaZmPanel2=bookME2D(ibooker, "TH2PullYvsBetaZmPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullXvsPhiZmPanel2DiskPlaquette){

           HistoName="mePullXvsPhiZmPanel2"+label;
           TitleName="mePullXvsPhiZmPanel2"+label;
           ST="";

           DiskPlaquettesObject->mePullXvsPhiZmPanel2=bookME2D(ibooker, "TH2PullXvsPhiZmPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullYvsPhiZmPanel2DiskPlaquette){

           HistoName="mePullYvsPhiZmPanel2"+label;
           TitleName="mePullYvsPhiZmPanel2"+label;
           ST="";

           DiskPlaquettesObject->mePullYvsPhiZmPanel2=bookME2D(ibooker, "TH2PullYvsPhiZmPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullXvsEtaZmPanel2DiskPlaquette){

           HistoName="mePullXvsEtaZmPanel2"+label;
           TitleName="mePullXvsEtaZmPanel2"+label;
           ST="";

           DiskPlaquettesObject->mePullXvsEtaZmPanel2=bookME2D(ibooker, "TH2PullXvsEtaZmPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullYvsEtaZmPanel2DiskPlaquette){

           HistoName="mePullYvsEtaZmPanel2"+label;
           TitleName="mePullYvsEtaZmPanel2"+label;
           ST="";

           DiskPlaquettesObject->mePullYvsEtaZmPanel2=bookME2D(ibooker, "TH2PullYvsEtaZmPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }

     //END-Zm-Panel2
   

    //Zp_-Panel2
    if (switchPosxZpPanel2DiskPlaquette){

           HistoName="mePosxZpPanel2"+label;
           TitleName="mePosxZpPanel2"+label;

           DiskPlaquettesObject->mePosxZpPanel2=bookME1D(ibooker, "TH1PosxZpPanel2",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchPosyZpPanel2DiskPlaquette){

           HistoName="mePosyZpPanel2"+label;
           TitleName="mePosyZpPanel2"+label;

           DiskPlaquettesObject->mePosyZpPanel2=bookME1D(ibooker, "TH1PosyZpPanel2",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchErrxZpPanel2DiskPlaquette){

           HistoName="meErrxZpPanel2"+label;
           TitleName="meErrxZpPanel2"+label;

           DiskPlaquettesObject->meErrxZpPanel2=bookME1D(ibooker, "TH1ErrxZpPanel2",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchErryZpPanel2DiskPlaquette){

           HistoName="meErryZpPanel2"+label;
           TitleName="meErryZpPanel2"+label;

           DiskPlaquettesObject->meErryZpPanel2=bookME1D(ibooker, "TH1ErryZpPanel2",HistoName.c_str(), TitleName.c_str());

    }  
    if (switchResxZpPanel2DiskPlaquette){

           HistoName="meResxZpPanel2"+label;
           TitleName="meResxZpPanel2"+label;

           DiskPlaquettesObject->meResxZpPanel2=bookME1D(ibooker, "TH1ResxZpPanel2",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchResyZpPanel2DiskPlaquette){

           HistoName="meResyZpPanel2"+label;
           TitleName="meResyZpPanel2"+label;

           DiskPlaquettesObject->meResyZpPanel2=bookME1D(ibooker, "TH1ResyZpPanel2",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchPullxZpPanel2DiskPlaquette){

           HistoName="mePullxZpPanel2"+label;
           TitleName="mePullxZpPanel2"+label;

           DiskPlaquettesObject->mePullxZpPanel2=bookME1D(ibooker, "TH1PullxZpPanel2",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchPullyZpPanel2DiskPlaquette){

           HistoName="mePullyZpPanel2"+label;
           TitleName="mePullyZpPanel2"+label;

           DiskPlaquettesObject->mePullyZpPanel2=bookME1D(ibooker, "TH1PullyZpPanel2",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchNpixZpPanel2DiskPlaquette){

           HistoName="meNpixZpPanel2"+label;
           TitleName="meNpixZpPanel2"+label;

           DiskPlaquettesObject->meNpixZpPanel2=bookME1D(ibooker, "TH1NpixZpPanel2",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchNxpixZpPanel2DiskPlaquette){

           HistoName="meNxpixZpPanel2"+label;
           TitleName="meNxpixZpPanel2"+label;

           DiskPlaquettesObject->meNxpixZpPanel2=bookME1D(ibooker, "TH1NxpixZpPanel2",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchNypixZpPanel2DiskPlaquette){

           HistoName="meNypixZpPanel2"+label;
           TitleName="meNypixZpPanel2"+label;

           DiskPlaquettesObject->meNypixZpPanel2=bookME1D(ibooker, "TH1NypixZpPanel2",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchChargeZpPanel2DiskPlaquette){

           HistoName="meChargeZpPanel2"+label;
           TitleName="meChargeZpPanel2"+label;

           DiskPlaquettesObject->meChargeZpPanel2=bookME1D(ibooker, "TH1ChargeZpPanel2",HistoName.c_str(), TitleName.c_str());
    }  
    if (switchResXvsAlphaZpPanel2DiskPlaquette){

           HistoName="meResXvsAlphaZpPanel2"+label;
           TitleName="meResXvsAlphaZpPanel2"+label;
           ST="";

           DiskPlaquettesObject->meResXvsAlphaZpPanel2=bookME2D(ibooker, "TH2ResXvsAlphaZpPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }

    if (switchResYvsAlphaZpPanel2DiskPlaquette){

           HistoName="meResYvsAlphaZpPanel2"+label;
           TitleName="meResYvsAlphaZpPanel2"+label;
           ST="";

           DiskPlaquettesObject->meResYvsAlphaZpPanel2=bookME2D(ibooker, "TH2ResYvsAlphaZpPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchResXvsBetaZpPanel2DiskPlaquette){

           HistoName="meResXvsBetaZpPanel2"+label;
           TitleName="meResXvsBetaZpPanel2"+label;
           ST="";

           DiskPlaquettesObject->meResXvsBetaZpPanel2=bookME2D(ibooker, "TH2ResXvsBetaZpPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchResYvsBetaZpPanel1DiskPlaquette){

           HistoName="meResYvsBetaZpPanel1"+label;
           TitleName="meResYvsBetaZpPanel1"+label;
           ST="";

           DiskPlaquettesObject->meResYvsBetaZpPanel1=bookME2D(ibooker, "TH2ResYvsBetaZpPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullXvsAlphaZpPanel2DiskPlaquette){

           HistoName="mePullXvsAlphaZpPanel2"+label;
           TitleName="mePullXvsAlphaZpPanel2"+label;
           ST="";

           DiskPlaquettesObject->mePullXvsAlphaZpPanel2=bookME2D(ibooker, "TH2PullXvsAlphaZpPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullYvsAlphaZpPanel2DiskPlaquette){

           HistoName="mePullYvsAlphaZpPanel2"+label;
           TitleName="mePullYvsAlphaZpPanel2"+label;
           ST="";

           DiskPlaquettesObject->mePullYvsAlphaZpPanel2=bookME2D(ibooker, "TH2PullYvsAlphaZpPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullXvsBetaZpPanel2DiskPlaquette){

           HistoName="mePullXvsBetaZpPanel2"+label;
           TitleName="mePullXvsBetaZpPanel2"+label;
           ST="";

           DiskPlaquettesObject->mePullXvsBetaZpPanel2=bookME2D(ibooker, "TH2PullXvsBetaZpPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullYvsBetaZpPanel2DiskPlaquette){

           HistoName="mePullYvsBetaZpPanel2"+label;
           TitleName="mePullYvsBetaZpPanel2"+label;
           ST="";

           DiskPlaquettesObject->mePullYvsBetaZpPanel2=bookME2D(ibooker, "TH2PullYvsBetaZpPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullXvsPhiZpPanel2DiskPlaquette){

           HistoName="mePullXvsPhiZpPanel2"+label;
           TitleName="mePullXvsPhiZpPanel2"+label;
           ST="";

           DiskPlaquettesObject->mePullXvsPhiZpPanel2=bookME2D(ibooker, "TH2PullXvsPhiZpPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullYvsPhiZpPanel2DiskPlaquette){

           HistoName="mePullYvsPhiZpPanel2"+label;
           TitleName="mePullYvsPhiZpPanel2"+label;
           ST="";

           DiskPlaquettesObject->mePullYvsPhiZpPanel2=bookME2D(ibooker, "TH2PullYvsPhiZpPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullXvsEtaZpPanel2DiskPlaquette){

           HistoName="mePullXvsEtaZpPanel2"+label;
           TitleName="mePullXvsEtaZpPanel2"+label;
           ST="";

           DiskPlaquettesObject->mePullXvsEtaZpPanel2=bookME2D(ibooker, "TH2PullXvsEtaZpPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
    if (switchPullYvsEtaZpPanel2DiskPlaquette){

           HistoName="mePullYvsEtaZpPanel2"+label;
           TitleName="mePullYvsEtaZpPanel2"+label;
           ST="";

           DiskPlaquettesObject->mePullYvsEtaZpPanel2=bookME2D(ibooker, "TH2PullYvsEtaZpPanel2", HistoName.c_str(), TitleName.c_str(), ST.c_str());
     }
     

     diskPlaquettesMEsMap[label]=DiskPlaquettesObject; 

}


void SiPixelTrackingRecHitsValid::bookHistograms(DQMStore::IBooker & ibooker,const edm::Run& run, const edm::EventSetup& es){
  
  // Book histograms
  dbe_ = edm::Service<DQMStore>().operator->();


  std::cout<<"Damir..."<<std::endl;
  
  edm::ESHandle<TrackerGeometry> pDD;
  es.get<TrackerDigiGeometryRecord>().get( pDD );

  PixelBarrelName *pbn;
  PixelEndcapName *pecn;

  std::stringstream labelss;
  std::string pathName;
 
  SiPixelFolderOrganizer theSiPixelFolder;


  ibooker.setCurrentFolder("Tracking/TrackingRecHits/Pixel/Histograms_per_ring-layer_or_disk-plaquette");

  for(TrackerGeometry::DetContainer::const_iterator it = pDD->dets().begin(); it != pDD->dets().end(); it++){ //LOOP DET-UNITS TO SET GEOMETRY EXTENSIONS FOR HISTOGRAMS

           DetId detId = (*it)->geographicalId();
           uint32_t id = detId();
           
       
           if(detId.subdetId() == static_cast<int>(PixelSubdetector::PixelBarrel)) {

                                 pbn = new PixelBarrelName(DetId(id));

                                 if (
                                         (pbn->shell() == 1) // mO
                                      && ((pbn->ladderName() == 1) || (pbn->ladderName() == 2))
                                      && (pbn->moduleName() == 1)
                                    ) {
                                         labelss.str("");
                                         if ((pbn->layerName()%2)^(pbn->ladderName()%2))
                                            labelss << "Flipped" << "Ladders" << "Layer_" << pbn->layerName();
                                         else
                                            labelss << "NonFlipper" << "Ladders" << "Layer_" << pbn->layerName();

                                         createLadderLayersMEs(ibooker,labelss.str());
                                         std::cout<<"Label="<<labelss.str()<<std::endl;
                                      }
                                    

                               	if (
                                      ((pbn->shell() == 1) || (pbn->shell() == 3)) // mO or pO
                                   && (pbn->ladderName() == 1)
                                   ) {
                                        labelss.str("");
                                        if (pbn->shell() == 1) // is mO
                                                      labelss << "Layer"<< "Module_" << pbn->layerName()<<"_"<<pbn->moduleName();
                                        else
                                                      labelss << "Layer"<< "Module_" << pbn->layerName()<<"_"<<(pbn->moduleName()+4);

                                        createLayerModulesMEs(ibooker, labelss.str());
                                        std::cout<<"Label="<<labelss.str()<<std::endl;
                                     }

                                 
                                if (
                                    ((pbn->shell() == 1)) // mO or pO
                                 && (pbn->ladderName() == 1)
                                 && (pbn->moduleName() == 1)
                                   ) {
                                           labelss.str("");
                                           labelss << "Layer_" << pbn->layerName();
                                           std::cout<<"Label="<<labelss.str()<<std::endl;
                                           createLayersMEs(ibooker, labelss.str());
                                     }

            }//if---PixelBarrel

            if(detId.subdetId() == static_cast<int>(PixelSubdetector::PixelEndcap)){

                	theSiPixelFolder.getModuleFolder(id, pathName, false); // isUpgrade = false
                        std::cout<<pathName<<std::endl;

                        pecn = new PixelEndcapName(DetId(id));

                    	if (
                              (pecn->halfCylinder() == 1) // mO
                           && (pecn->bladeName() == 1)
                           ) {
                                labelss.str("");

                                if ((pecn->pannelName() == 1) && (pecn->diskName()==1))
                                                   labelss << "Disk" << "Plaq_"<<pecn->diskName() << "_" << pecn->plaquetteName();
                                if ((pecn->pannelName() == 2) && (pecn->diskName()==2))
                                                   labelss << "Disk" << "Plaq_"<<pecn->diskName() << "_" << pecn->plaquetteName();

                                createDiskPlaquettesMEs(ibooker, labelss.str());
                                std::cout<<labelss.str()<<std::endl;
                             
                             }

            }

  }//END--LOOPING DET-UNITS
  
  ibooker.setCurrentFolder("Tracking/TrackingRecHits/Pixel/Histograms_all");

  createTotalMEs(ibooker, "");//no geometry extensions---outside det-units for loop


}

// Virtual destructor needed.
SiPixelTrackingRecHitsValid::~SiPixelTrackingRecHitsValid() 
{  
  //save local root file only in standalone mode
  if ( runStandalone && outputFile_.size() != 0 && dbe_ ) dbe_->save(outputFile_);
}  

// Functions that gets called by framework every event
void SiPixelTrackingRecHitsValid::analyze(const edm::Event& e, const edm::EventSetup& es)
{
  //Retrieve tracker topology from geometry
  edm::ESHandle<TrackerTopology> tTopo;
  es.get<IdealGeometryRecord>().get(tTopo);

  std::stringstream labelss;

  run = e.id().run();
  evt = e.id().event();

  //  if ( evt%1000 == 0 ) 
    //cout << "evt = " << evt << endl;
  
  float math_pi = 3.14159265;
  float radtodeg = 180.0 / math_pi;
    
  DetId detId;

  LocalPoint position;
  LocalError error;
  float mindist = 999999.9;

  std::vector<PSimHit> matched;
  TrackerHitAssociator associate(e,conf_);

  edm::ESHandle<TrackerGeometry> pDD;
  es.get<TrackerDigiGeometryRecord> ().get (pDD);
  const TrackerGeometry* tracker = &(* pDD);

  if ( !MTCCtrack_ )
    {
      // --------------------------------------- all hits -----------------------------------------------------------
      //--- Fetch Pixel RecHits
      edm::Handle<SiPixelRecHitCollection> recHitColl;
      e.getByToken( siPixelRecHitCollectionToken_, recHitColl );
      
      //cout <<" ----- Found " 
      //   << const_cast<SiPixelRecHitCollection*>(recHitColl.product())->size()
      //   << " Pixel RecHits" << std::endl;
  
      //-----Iterate over detunits
      for (TrackerGeometry::DetContainer::const_iterator it = pDD->dets().begin(); it != pDD->dets().end(); it++) 
	{
	  DetId detId = ((*it)->geographicalId());
	 
	  unsigned int subid = detId.subdetId();
	  if ( !((subid==1) || (subid==2)) ) 
	    continue; // end subid if

          SiPixelRecHitCollection::const_iterator match = recHitColl->find(detId);
          if (match == recHitColl->end()) continue;

          SiPixelRecHitCollection::DetSet pixelrechitRange = *match;
          SiPixelRecHitCollection::DetSet::const_iterator pixelrechitRangeIteratorBegin = pixelrechitRange.begin();
          SiPixelRecHitCollection::DetSet::const_iterator pixelrechitRangeIteratorEnd = pixelrechitRange.end();
          SiPixelRecHitCollection::DetSet::const_iterator pixeliter = pixelrechitRangeIteratorBegin;
	  std::vector<PSimHit> matched;
	  
	  //----Loop over rechits for this detId
	  for ( ; pixeliter != pixelrechitRangeIteratorEnd; ++pixeliter) 
	    {
	      LocalPoint lp = pixeliter->localPosition();
	      float rechitx = lp.x();
	      float rechity = lp.y();
	     
	      detId = (*it)->geographicalId();
	      subdetId = (int)detId.subdetId();
	      if ( (int)detId.subdetId() == (int)PixelSubdetector::PixelBarrel ) 
		{
		  fillME(totalMEs->mePosxBarrel_all_hits, rechitx );
		  fillME(totalMEs->mePosyBarrel_all_hits, rechity );
		}
	      else if ( (int)detId.subdetId() == (int)PixelSubdetector::PixelEndcap )
		{
		  
		  side  = tTopo->pxfSide(detId);
		  disk  = tTopo->pxfDisk(detId);
		  blade = tTopo->pxfBlade(detId);
		  panel = tTopo->pxfPanel(detId);
		  plaq  = tTopo->pxfModule(detId); // also known as plaquette
		  
		  if ( side==1 ) 
		    {
		      if ( panel==1 )
			{
			  fillME(totalMEs->mePosxZmPanel1_all_hits, rechitx );
			  fillME(totalMEs->mePosyZmPanel1_all_hits, rechity );
			}
		      else if ( panel==2 )
			{
			  fillME(totalMEs->mePosxZmPanel2_all_hits, rechitx );
			  fillME(totalMEs->mePosyZmPanel2_all_hits, rechity );
			}
		      else edm::LogWarning("SiPixelTrackingRecHitsValid") << "..............................................Wrong panel number !"; 
		    } // if ( side==1 ) 
		  else if ( side==2 )
		    {
		      if ( panel==1 )
			{
			  fillME(totalMEs->mePosxZpPanel1_all_hits, rechitx );
			  fillME(totalMEs->mePosyZpPanel1_all_hits, rechity );
			}
		       else if ( panel==2 )
			 {
			   fillME(totalMEs->mePosxZpPanel2_all_hits, rechitx );
			   fillME(totalMEs->mePosyZpPanel2_all_hits, rechity );
			 }
		       else  edm::LogWarning("SiPixelTrackingRecHitsValid")<< "..............................................Wrong panel number !";
		    } //else if ( side==2 )
		  else edm::LogWarning("SiPixelTrackingRecHitsValid") << ".......................................................Wrong side !" ;
		  
		} // else if ( detId.subdetId()==PixelSubdetector::PixelEndcap )
	      else edm::LogWarning("SiPixelTrackingRecHitsValid") << "Pixel rechit collection but we are not in the pixel detector" << (int)detId.subdetId() ;
	      
	    }
	}
      // ------------------------------------------------ all hits ---------------------------------------------------------------
       
      // Get tracks
      edm::Handle<reco::TrackCollection> trackCollection;
      e.getByToken( recoTrackCollectionToken_, trackCollection );
      const reco::TrackCollection *tracks = trackCollection.product();
      reco::TrackCollection::const_iterator tciter;

      int n_tracks = (int)tracks->size(); // number of tracks in this event
      fillME(totalMEs->meTracksPerEvent, n_tracks );

      if ( tracks->size() > 0 )
	{
	  // Loop on tracks
	  for ( tciter=tracks->begin(); tciter!=tracks->end(); tciter++)
	    {
	      phi = tciter->momentum().phi() / math_pi*180.0;
	      eta = tciter->momentum().eta();
	      
	      int n_hits = 0;
	      // First loop on hits: find matched hits
	      for ( trackingRecHit_iterator it = tciter->recHitsBegin(); it != tciter->recHitsEnd(); it++) 
		{
		  const TrackingRecHit &thit = **it;
		  // Is it a matched hit?
		  const SiPixelRecHit* matchedhit = dynamic_cast<const SiPixelRecHit*>(&thit);
		  
		  if ( matchedhit ) 
		    {
		      ++n_hits;
		      
		      layer  = -9999; 
		      ladder = -9999; 
		      mod    = -9999; 
		      side   = -9999;  
		      disk   = -9999;  
		      blade  = -9999; 
		      panel  = -9999; 
		      plaq   = -9999; 

		      rechitx = -9999.9;
		      rechity = -9999.9;
		      rechitz = -9999.9;
		      rechiterrx = -9999.9;
		      rechiterry = -9999.9;		      
		      rechitresx = -9999.9;
		      rechitresy = -9999.9;
		      rechitpullx = -9999.9;
		      rechitpully = -9999.9;
		      
		      npix = -9999;
		      nxpix = -9999;
		      nypix = -9999;
		      charge = -9999.9;
		      
		      alpha = -9999.9;
		      beta  = -9999.9;

		      half = -9999;
		      flipped = -9999;
		 
		      nsimhit = -9999;
		         
		      simhitx = -9999.9;
		      simhity = -9999.9;

		      position = (*it)->localPosition();
		      error = (*it)->localPositionError();

		      rechitx = position.x();
		      rechity = position.y();
		      rechitz = position.z();
		      rechiterrx = sqrt(error.xx());
		      rechiterry = sqrt(error.yy());

		      npix = (*matchedhit).cluster()->size();
		      nxpix = (*matchedhit).cluster()->sizeX();
		      nypix = (*matchedhit).cluster()->sizeY();
		      charge = (*matchedhit).cluster()->charge();

		      //Association of the rechit to the simhit
		      matched.clear();
		      matched = associate.associateHit(*matchedhit);

		      nsimhit = (int)matched.size();

		      if ( !matched.empty() ) 
			{
			  mindist = 999999.9;
			  float distx, disty, dist;
			  bool found_hit_from_generated_particle = false;
			  
			  int n_assoc_muon = 0;

			  std::vector<PSimHit>::const_iterator closestit = matched.begin();
			  for (std::vector<PSimHit>::const_iterator m=matched.begin(); m<matched.end(); m++)
			    {
			      if ( checkType_ )
				{
				  int pid = (*m).particleType();
				  if ( abs(pid) != genType_ )
				    continue;
				  
				}
			      
			      float simhitx = 0.5 * ( (*m).entryPoint().x() + (*m).exitPoint().x() );
			      float simhity = 0.5 * ( (*m).entryPoint().y() + (*m).exitPoint().y() );
			      
			      distx = fabs(rechitx - simhitx);
			      disty = fabs(rechity - simhity);
			      dist = sqrt( distx*distx + disty*disty );
	
			      if ( dist < mindist )
				{
				  n_assoc_muon++;

				  mindist = dist;
				  closestit = m;
				  found_hit_from_generated_particle = true;
				}
			    } // for (std::vector<PSimHit>::const_iterator m=matched.begin(); m<matched.end(); m++)
			  
			  // This recHit does not have any simHit with the same particleType as the particles generated
			  // Ignore it as most probably come from delta rays.
			  if ( checkType_ && !found_hit_from_generated_particle )
			    continue; 
			  
			  if ( n_assoc_muon > 1 )
			    {
			      edm::LogWarning("SiPixelTrackingRecHitsValid") << " ----- This is not good: n_assoc_muon = " << n_assoc_muon ;
			      edm::LogWarning("SiPixelTrackingRecHitsValid") << "evt = " << evt ;
			    }

			  pidhit = (*closestit).particleType();

			  simhitx = 0.5*( (*closestit).entryPoint().x() + (*closestit).exitPoint().x() );
			  simhity = 0.5*( (*closestit).entryPoint().y() + (*closestit).exitPoint().y() );
			  
			  rechitresx = rechitx - simhitx;
			  rechitresy = rechity - simhity;
			  rechitpullx = ( rechitx - simhitx ) / sqrt(error.xx());
			  rechitpully = ( rechity - simhity ) / sqrt(error.yy());

			  float simhitpx = (*closestit).momentumAtEntry().x();
			  float simhitpy = (*closestit).momentumAtEntry().y();
			  float simhitpz = (*closestit).momentumAtEntry().z();
			  			  
			  //beta  = atan2(simhitpz, simhitpy) * radtodeg;
			  //alpha = atan2(simhitpz, simhitpx) * radtodeg;
			  
			  beta  = fabs(atan2(simhitpz, simhitpy)) * radtodeg;
			  alpha = fabs(atan2(simhitpz, simhitpx)) * radtodeg;
		
			  detId = (*it)->geographicalId();

			  subdetId = (int)detId.subdetId();

			  if ( (int)detId.subdetId() == (int)PixelSubdetector::PixelBarrel ) 
			    {
			      fillME(totalMEs->mePosxBarrel, rechitx );
			      fillME(totalMEs->mePosyBarrel, rechity );
			      fillME(totalMEs->meErrxBarrel, rechiterrx );			      
			      fillME(totalMEs->meErryBarrel, rechiterry );
			      fillME(totalMEs->meResxBarrel, rechitresx );
			      fillME(totalMEs->meResyBarrel, rechitresy );
			      fillME(totalMEs->mePullxBarrel, rechitpullx );
			      fillME(totalMEs->mePullyBarrel, rechitpully );
			      fillME(totalMEs->meNpixBarrel, npix );
			      fillME(totalMEs->meNxpixBarrel, nxpix );
			      fillME(totalMEs->meNypixBarrel, nypix );
			      fillME(totalMEs->meChargeBarrel, charge );
			      fillME(totalMEs->meResXvsAlphaBarrel, alpha, fabs(rechitresx) );
			      fillME(totalMEs->meResYvsAlphaBarrel, alpha, fabs(rechitresy) );
			      fillME(totalMEs->meResXvsBetaBarrel, beta, fabs(rechitresx) );
			      fillME(totalMEs->meResYvsBetaBarrel, beta, fabs(rechitresy) );
			      fillME(totalMEs->mePullXvsAlphaBarrel, alpha, rechitpullx );
			      fillME(totalMEs->mePullYvsAlphaBarrel, alpha, rechitpully );
			      fillME(totalMEs->mePullXvsBetaBarrel, beta, rechitpullx );
			      fillME(totalMEs->mePullYvsBetaBarrel, beta, rechitpully );
			      fillME(totalMEs->mePullXvsPhiBarrel, phi, rechitpullx );
			      fillME(totalMEs->mePullYvsPhiBarrel, phi, rechitpully );
			      fillME(totalMEs->mePullXvsEtaBarrel, eta, rechitpullx );
			      fillME(totalMEs->mePullYvsEtaBarrel, eta, rechitpully );
			      
			      const PixelGeomDetUnit * theGeomDet 
				= dynamic_cast<const PixelGeomDetUnit*> ( tracker->idToDet(detId) );
			      //const PixelTopology * topol = (&(theGeomDet->specificTopology()));
			      
			      int tmp_nrows = theGeomDet->specificTopology().nrows();
			      
			      if ( tmp_nrows == 80 ) 
				{
				  fillME(totalMEs->mePosxBarrelHalfModule, rechitx );
				  fillME(totalMEs->mePosyBarrelHalfModule, rechity );

				  half = 1;
				}
			      else if ( tmp_nrows == 160 ) 
				{
				  fillME( totalMEs->mePosxBarrelFullModule, rechitx );
				  fillME( totalMEs->mePosyBarrelFullModule, rechity );
				  half = 0;
				}
			      else 
				edm::LogWarning("SiPixelTrackingRecHitsValid") << "-------------------------------------------------- Wrong module size !!!";

			      float tmp1 = theGeomDet->surface().toGlobal(Local3DPoint(0.,0.,0.)).perp();
			      float tmp2 = theGeomDet->surface().toGlobal(Local3DPoint(0.,0.,1.)).perp();
			      
			      if ( tmp2<tmp1 ) 
				{ // flipped
				  fillME(totalMEs->mePosxBarrelFlippedLadders, rechitx );
				  fillME(totalMEs->mePosyBarrelFlippedLadders, rechity );
				  flipped = 1;
				
				  fillME(totalMEs->meResXvsAlphaBarrelFlippedLadders, alpha, fabs(rechitresx) );
				  fillME(totalMEs->meResYvsAlphaBarrelFlippedLadders, alpha, fabs(rechitresy) );
				  fillME(totalMEs->meResXvsBetaBarrelFlippedLadders, beta, fabs(rechitresx) );
				  fillME(totalMEs->meResYvsBetaBarrelFlippedLadders, beta, fabs(rechitresy) );
				  fillME(totalMEs->mePullXvsAlphaBarrelFlippedLadders, alpha, rechitpullx );
				  fillME(totalMEs->mePullYvsAlphaBarrelFlippedLadders, alpha, rechitpully );
				  fillME(totalMEs->mePullXvsBetaBarrelFlippedLadders, beta, rechitpullx );
				  fillME(totalMEs->mePullYvsBetaBarrelFlippedLadders, beta, rechitpully );
				  fillME(totalMEs->mePullXvsPhiBarrelFlippedLadders, phi, rechitpullx );
				  fillME(totalMEs->mePullYvsPhiBarrelFlippedLadders, phi, rechitpully );
				  fillME(totalMEs->mePullXvsEtaBarrelFlippedLadders, eta, rechitpullx );
				  fillME(totalMEs->mePullYvsEtaBarrelFlippedLadders, eta, rechitpully );
				
				  fillME(totalMEs->meWPullXvsAlphaBarrelFlippedLadders, alpha, fabs(rechitpullx) );
				  fillME(totalMEs->meWPullYvsAlphaBarrelFlippedLadders, alpha, fabs(rechitpully) );
				  fillME(totalMEs->meWPullXvsBetaBarrelFlippedLadders, beta, fabs(rechitpullx) );
				  fillME(totalMEs->meWPullYvsBetaBarrelFlippedLadders, beta, fabs(rechitpully) );
				}
			      else 
				{ // not flipped
				  fillME(totalMEs->mePosxBarrelNonFlippedLadders, rechitx );
				  fillME(totalMEs->mePosyBarrelNonFlippedLadders, rechity );
				  flipped = 0;
				
				  fillME(totalMEs->meResXvsAlphaBarrelNonFlippedLadders, alpha, fabs(rechitresx) );
				  fillME(totalMEs->meResYvsAlphaBarrelNonFlippedLadders, alpha, fabs(rechitresy) );
				  fillME(totalMEs->meResXvsBetaBarrelNonFlippedLadders, beta, fabs(rechitresx) );
				  fillME(totalMEs->meResYvsBetaBarrelNonFlippedLadders, beta, fabs(rechitresy) );
				  fillME(totalMEs->mePullXvsAlphaBarrelNonFlippedLadders, alpha, rechitpullx );
				  fillME(totalMEs->mePullYvsAlphaBarrelNonFlippedLadders, alpha, rechitpully );
				  fillME(totalMEs->mePullXvsBetaBarrelNonFlippedLadders, beta, rechitpullx );
				  fillME(totalMEs->mePullYvsBetaBarrelNonFlippedLadders, beta, rechitpully );
				  fillME(totalMEs->mePullXvsPhiBarrelNonFlippedLadders, phi, rechitpullx );
				  fillME(totalMEs->mePullYvsPhiBarrelNonFlippedLadders, phi, rechitpully );
				  fillME(totalMEs->mePullXvsEtaBarrelNonFlippedLadders, eta, rechitpullx );
				  fillME(totalMEs->mePullYvsEtaBarrelNonFlippedLadders, eta, rechitpully );

				  fillME(totalMEs->meWPullXvsAlphaBarrelNonFlippedLadders, alpha, fabs(rechitpullx) );
				  fillME(totalMEs->meWPullYvsAlphaBarrelNonFlippedLadders, alpha, fabs(rechitpully) );
				  fillME(totalMEs->meWPullXvsBetaBarrelNonFlippedLadders, beta, fabs(rechitpullx) );
				  fillME(totalMEs->meWPullYvsBetaBarrelNonFlippedLadders, beta, fabs(rechitpully) );
				}
			          
			      
			      layer  = tTopo->pxbLayer(detId);   // Layer: 1,2,3.
			      ladder = tTopo->pxbLadder(detId);  // Ladder: 1-20, 32, 44. 
			      mod   = tTopo->pxbModule(detId);  // Mod: 1-8.
			      

                              labelss.str("");
                              labelss << "Layer"<< "Module_" << layer<<"_"<<mod;
                            
                              fillME(layerModulesMEsMap[labelss.str()]->mePosxBarrel, rechitx );
                              fillME(layerModulesMEsMap[labelss.str()]->mePosyBarrel, rechity );
                              fillME(layerModulesMEsMap[labelss.str()]->meErrxBarrel, rechiterrx );
                              fillME(layerModulesMEsMap[labelss.str()]->meErryBarrel, rechiterry );
                              fillME(layerModulesMEsMap[labelss.str()]->meResxBarrel, rechitresx );
                              fillME(layerModulesMEsMap[labelss.str()]->meResyBarrel, rechitresy );
                              fillME(layerModulesMEsMap[labelss.str()]->mePullxBarrel, rechitpullx );
                              fillME(layerModulesMEsMap[labelss.str()]->mePullyBarrel, rechitpully );
                              fillME(layerModulesMEsMap[labelss.str()]->meNpixBarrel, npix );
                              fillME(layerModulesMEsMap[labelss.str()]->meNxpixBarrel, nxpix );
                              fillME(layerModulesMEsMap[labelss.str()]->meNypixBarrel, nypix );
                              fillME(layerModulesMEsMap[labelss.str()]->meChargeBarrel, charge );
                              fillME(layerModulesMEsMap[labelss.str()]->meResXvsAlphaBarrel, alpha, fabs(rechitresx) );
                              fillME(layerModulesMEsMap[labelss.str()]->meResYvsAlphaBarrel, alpha, fabs(rechitresy) );
                              fillME(layerModulesMEsMap[labelss.str()]->meResXvsBetaBarrel, beta, fabs(rechitresx) );
                              fillME(layerModulesMEsMap[labelss.str()]->meResYvsBetaBarrel, beta, fabs(rechitresy) );
                              fillME(layerModulesMEsMap[labelss.str()]->mePullXvsAlphaBarrel, alpha, rechitpullx );
                              fillME(layerModulesMEsMap[labelss.str()]->mePullYvsAlphaBarrel, alpha, rechitpully );
                              fillME(layerModulesMEsMap[labelss.str()]->mePullXvsBetaBarrel, beta, rechitpullx );
                              fillME(layerModulesMEsMap[labelss.str()]->mePullYvsBetaBarrel, beta, rechitpully );
                              fillME(layerModulesMEsMap[labelss.str()]->mePullXvsPhiBarrel, phi, rechitpullx );
                              fillME(layerModulesMEsMap[labelss.str()]->mePullYvsPhiBarrel, phi, rechitpully );
                              fillME(layerModulesMEsMap[labelss.str()]->mePullXvsEtaBarrel, eta, rechitpullx );
                              fillME(layerModulesMEsMap[labelss.str()]->mePullYvsEtaBarrel, eta, rechitpully );

			    
                              labelss.str("");
                              labelss << "Layer_"<< layer;
                            
                              fillME(layerMEsMap[labelss.str()]->meResxBarrel, rechitresx );
                              fillME(layerMEsMap[labelss.str()]->meResyBarrel, rechitresy );
                              fillME(layerMEsMap[labelss.str()]->mePullxBarrel, rechitpullx );
                              fillME(layerMEsMap[labelss.str()]->mePullyBarrel, rechitpully );


			      if ( tmp2<tmp1 ) 
				{ // flipped

                               	  labelss.str("");
                                  labelss << "Flipped" << "Ladders" << "Layer_" << layer;                                   

                                  fillME(ladderLayersMEsMap[labelss.str()]->meResXvsAlphaBarrel, alpha, fabs(rechitresx) );
                                  fillME(ladderLayersMEsMap[labelss.str()]->meResYvsAlphaBarrel, alpha, fabs(rechitresy) );
				  fillME(ladderLayersMEsMap[labelss.str()]->meResXvsBetaBarrel, beta, fabs(rechitresx) );//-----to add
				  fillME(ladderLayersMEsMap[labelss.str()]->meResYvsBetaBarrel, beta, fabs(rechitresy) );//-----to add
				}
			      else
				{ // not flipped

                               	  labelss.str("");
                                  labelss << "NonFlipped" << "Ladders" << "Layer_" << layer;                                   

                                  fillME(ladderLayersMEsMap[labelss.str()]->meResXvsAlphaBarrel, alpha, fabs(rechitresx) );
                                  fillME(ladderLayersMEsMap[labelss.str()]->meResYvsAlphaBarrel, alpha, fabs(rechitresy) );
				  fillME(ladderLayersMEsMap[labelss.str()]->meResXvsBetaBarrel, beta, fabs(rechitresx) );//-----to add
				  fillME(ladderLayersMEsMap[labelss.str()]->meResYvsBetaBarrel, beta, fabs(rechitresy) );//-----to add
				}
			      
			    }
			  else if ( (int)detId.subdetId() == (int)PixelSubdetector::PixelEndcap )
			    {
			      std::cout<<"Inside of endcap"<<std::endl;
			      side  = tTopo->pxfSide(detId);
			      disk  = tTopo->pxfDisk(detId);
			      blade = tTopo->pxfBlade(detId);
			      panel = tTopo->pxfPanel(detId);
			      plaq  = tTopo->pxfModule(detId); // also known as plaquette

			      if ( side==1 ) 
				{
				  if ( panel==1 )
				    {
				      fillME(totalMEs->mePosxZmPanel1, rechitx );
				      fillME(totalMEs->mePosyZmPanel1, rechity );
				      fillME(totalMEs->meErrxZmPanel1, rechiterrx );
				      fillME(totalMEs->meErryZmPanel1, rechiterry );
				      fillME(totalMEs->meResxZmPanel1, rechitresx );
				      fillME(totalMEs->meResyZmPanel1, rechitresy );
				      fillME(totalMEs->mePullxZmPanel1, rechitpullx );
				      fillME(totalMEs->mePullyZmPanel1, rechitpully );
				      fillME(totalMEs->meNpixZmPanel1, npix );
				      fillME(totalMEs->meNxpixZmPanel1, nxpix );
				      fillME(totalMEs->meNypixZmPanel1, nypix );
				      fillME(totalMEs->meChargeZmPanel1, charge );
				      fillME(totalMEs->meResXvsAlphaZmPanel1, alpha, fabs(rechitresx) );
				      fillME(totalMEs->meResYvsAlphaZmPanel1, alpha, fabs(rechitresy) );
				      fillME(totalMEs->meResXvsBetaZmPanel1, beta, fabs(rechitresx) );
				      fillME(totalMEs->meResYvsBetaZmPanel1, beta, fabs(rechitresy) );
				      fillME(totalMEs->mePullXvsAlphaZmPanel1, alpha, rechitpullx );
				      fillME(totalMEs->mePullYvsAlphaZmPanel1, alpha, rechitpully );
				      fillME(totalMEs->mePullXvsBetaZmPanel1, beta, rechitpullx );
				      fillME(totalMEs->mePullYvsBetaZmPanel1, beta, rechitpully );
				      fillME(totalMEs->mePullXvsPhiZmPanel1, phi, rechitpullx );
				      fillME(totalMEs->mePullYvsPhiZmPanel1, phi, rechitpully );
				      fillME(totalMEs->mePullXvsEtaZmPanel1, eta, rechitpullx );
				      fillME(totalMEs->mePullYvsEtaZmPanel1, eta, rechitpully );
				      
				      fillME(totalMEs->meWPullXvsAlphaZmPanel1, alpha, fabs(rechitpullx) );
				      fillME(totalMEs->meWPullYvsAlphaZmPanel1, alpha, fabs(rechitpully) );
				      fillME(totalMEs->meWPullXvsBetaZmPanel1, beta, fabs(rechitpullx) );
				      fillME(totalMEs->meWPullYvsBetaZmPanel1, beta, fabs(rechitpully) );
                               	      labelss.str("");
                                      labelss << "Disk" << "Plaq_"<<disk<< "_" <<plaq;
      
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePosxZmPanel1, rechitx);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePosyZmPanel1, rechity);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meErrxZmPanel1, rechiterrx);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meErryZmPanel1, rechiterry);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meResxZmPanel1, rechitresx);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meResyZmPanel1, rechitresy);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullxZmPanel1, rechitpullx);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullyZmPanel1, rechitpully);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meNpixZmPanel1, npix);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meNxpixZmPanel1, nxpix);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meNypixZmPanel1, nypix);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meChargeZmPanel1, charge);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meResXvsAlphaZmPanel1, alpha, fabs(rechitresx));
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meResYvsAlphaZmPanel1, alpha, fabs(rechitresy));
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meResXvsBetaZmPanel1, beta, fabs(rechitresx));
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meResYvsBetaZmPanel1, beta, fabs(rechitresy));
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullXvsAlphaZmPanel1, alpha, rechitpullx);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullYvsAlphaZmPanel1, alpha, rechitpully);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullXvsBetaZmPanel1, beta, rechitpullx);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullYvsBetaZmPanel1, beta, rechitpully);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullXvsPhiZmPanel1, phi, rechitpullx);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullYvsPhiZmPanel1, phi, rechitpully);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullXvsEtaZmPanel1, eta, rechitpullx);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullYvsEtaZmPanel1, eta, rechitpully);

				      
				    }
				  else if ( panel==2 )
				    {
				      fillME(totalMEs->mePosxZmPanel2, rechitx );
				      fillME(totalMEs->mePosyZmPanel2, rechity );
				      fillME(totalMEs->meErrxZmPanel2, rechiterrx );
				      fillME(totalMEs->meErryZmPanel2, rechiterry );
				      fillME(totalMEs->meResxZmPanel2, rechitresx );
				      fillME(totalMEs->meResyZmPanel2, rechitresy );
				      fillME(totalMEs->mePullxZmPanel2, rechitpullx );
				      fillME(totalMEs->mePullyZmPanel2, rechitpully );
				      fillME(totalMEs->meNpixZmPanel2, npix );
				      fillME(totalMEs->meNxpixZmPanel2, nxpix );
				      fillME(totalMEs->meNypixZmPanel2, nypix );
				      fillME(totalMEs->meChargeZmPanel2, charge );

				      fillME(totalMEs->meResXvsAlphaZmPanel2, alpha, fabs(rechitresx) );
				      fillME(totalMEs->meResYvsAlphaZmPanel2, alpha, fabs(rechitresy) );
				      fillME(totalMEs->meResXvsBetaZmPanel2, beta, fabs(rechitresx) );
				      fillME(totalMEs->meResYvsBetaZmPanel2, beta, fabs(rechitresy) );
				      fillME(totalMEs->mePullXvsAlphaZmPanel2, alpha, rechitpullx );
				      fillME(totalMEs->mePullYvsAlphaZmPanel2, alpha, rechitpully );
				      fillME(totalMEs->mePullXvsBetaZmPanel2, beta, rechitpullx );
				      fillME(totalMEs->mePullYvsBetaZmPanel2, beta, rechitpully );
				      fillME(totalMEs->mePullXvsPhiZmPanel2, phi, rechitpullx );
				      fillME(totalMEs->mePullYvsPhiZmPanel2, phi, rechitpully );
				      fillME(totalMEs->mePullXvsEtaZmPanel2, eta, rechitpullx );
				      fillME(totalMEs->mePullYvsEtaZmPanel2, eta, rechitpully );

				      fillME(totalMEs->meWPullXvsAlphaZmPanel2, alpha, fabs(rechitpullx) );
				      fillME(totalMEs->meWPullYvsAlphaZmPanel2, alpha, fabs(rechitpully) );
				      fillME(totalMEs->meWPullXvsBetaZmPanel2, beta, fabs(rechitpullx) );
				      fillME(totalMEs->meWPullYvsBetaZmPanel2, beta, fabs(rechitpully) );



                                      labelss.str("");
                                      labelss << "Disk" << "Plaq_"<<disk<< "_" <<plaq;

                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePosxZmPanel2, rechitx);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePosyZmPanel2, rechity);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meErrxZmPanel2, rechiterrx);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meErryZmPanel2, rechiterry);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meResxZmPanel2, rechitresx);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meResyZmPanel2, rechitresy);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullxZmPanel2, rechitpullx);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullyZmPanel2, rechitpully);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meNpixZmPanel2, npix);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meNxpixZmPanel2, nxpix);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meNypixZmPanel2, nypix);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meChargeZmPanel2, charge);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meResXvsAlphaZmPanel2, alpha, fabs(rechitresx));
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meResYvsAlphaZmPanel2, alpha, fabs(rechitresy));
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meResXvsBetaZmPanel2, beta, fabs(rechitresx));
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meResYvsBetaZmPanel2, beta, fabs(rechitresy));
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullXvsAlphaZmPanel2, alpha, rechitpullx);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullYvsAlphaZmPanel2, alpha, rechitpully);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullXvsBetaZmPanel2, beta, rechitpullx);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullYvsBetaZmPanel2, beta, rechitpully);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullXvsPhiZmPanel2, phi, rechitpullx);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullYvsPhiZmPanel2, phi, rechitpully);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullXvsEtaZmPanel2, eta, rechitpullx);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullYvsEtaZmPanel2, eta, rechitpully);


				    }
				  else edm::LogWarning("SiPixelTrackingRecHitsValid") << "..............................................Wrong panel number !"; 
				} // if ( side==1 ) 
			      else if ( side==2 )
				{
				  if ( panel==1 )
				    {
				      fillME(totalMEs->mePosxZpPanel1, rechitx );
				      fillME(totalMEs->mePosyZpPanel1, rechity );
				      fillME(totalMEs->meErrxZpPanel1, rechiterrx );
				      fillME(totalMEs->meErryZpPanel1, rechiterry );
				      fillME(totalMEs->meResxZpPanel1, rechitresx );
				      fillME(totalMEs->meResyZpPanel1, rechitresy );
				      fillME(totalMEs->mePullxZpPanel1, rechitpullx );
				      fillME(totalMEs->mePullyZpPanel1, rechitpully );
				      fillME(totalMEs->meNpixZpPanel1, npix );
				      fillME(totalMEs->meNxpixZpPanel1, nxpix );
				      fillME(totalMEs->meNypixZpPanel1, nypix );
				      fillME(totalMEs->meChargeZpPanel1, charge );
				      fillME(totalMEs->meResXvsAlphaZpPanel1, alpha, fabs(rechitresx) );
				      fillME(totalMEs->meResYvsAlphaZpPanel1, alpha, fabs(rechitresy) );
				      fillME(totalMEs->meResXvsBetaZpPanel1, beta, fabs(rechitresx) );
				      fillME(totalMEs->meResYvsBetaZpPanel1, beta, fabs(rechitresy) );
				      fillME(totalMEs->mePullXvsAlphaZpPanel1, alpha, rechitpullx );
				      fillME(totalMEs->mePullYvsAlphaZpPanel1, alpha, rechitpully );
				      fillME(totalMEs->mePullXvsBetaZpPanel1, beta, rechitpullx );
				      fillME(totalMEs->mePullYvsBetaZpPanel1, beta, rechitpully );
				      fillME(totalMEs->mePullXvsPhiZpPanel1, phi, rechitpullx );
				      fillME(totalMEs->mePullYvsPhiZpPanel1, phi, rechitpully );
				      fillME(totalMEs->mePullXvsEtaZpPanel1, eta, rechitpullx );
				      fillME(totalMEs->mePullYvsEtaZpPanel1, eta, rechitpully );

				      fillME(totalMEs->meWPullXvsAlphaZpPanel1, alpha, fabs(rechitpullx) );
				      fillME(totalMEs->meWPullYvsAlphaZpPanel1, alpha, fabs(rechitpully) );
				      fillME(totalMEs->meWPullXvsBetaZpPanel1, beta, fabs(rechitpullx) );
				      fillME(totalMEs->meWPullYvsBetaZpPanel1, beta, fabs(rechitpully) );

                               	      labelss.str("");
                                      labelss << "Disk" << "Plaq_"<<disk<< "_" <<plaq;

                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePosxZpPanel1, rechitx);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePosyZpPanel1, rechity);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meErrxZpPanel1, rechiterrx);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meErryZpPanel1, rechiterry);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meResxZpPanel1, rechitresx);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meResyZpPanel1, rechitresy);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullxZpPanel1, rechitpullx);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullyZpPanel1, rechitpully);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meNpixZpPanel1, npix);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meNxpixZpPanel1, nxpix);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meNypixZpPanel1, nypix);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meChargeZpPanel1, charge); 
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meResXvsAlphaZpPanel1, alpha, fabs(rechitresx));
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meResYvsAlphaZpPanel1, alpha, fabs(rechitresy));
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meResXvsBetaZpPanel1, beta, fabs(rechitresx));
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meResYvsBetaZpPanel1, beta, fabs(rechitresy));
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullXvsAlphaZpPanel1, alpha, rechitpullx);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullYvsAlphaZpPanel1, alpha, rechitpully);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullXvsBetaZpPanel1, beta, rechitpullx);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullYvsBetaZpPanel1, beta, rechitpully);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullXvsPhiZpPanel1, phi, rechitpullx);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullYvsPhiZpPanel1, phi, rechitpully);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullXvsEtaZpPanel1, eta, rechitpullx);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullYvsEtaZpPanel1, eta, rechitpully);

				      
				    }
				  else if ( panel==2 )
				    {

				      fillME(totalMEs->mePosxZpPanel2, rechitx );
				      fillME(totalMEs->mePosyZpPanel2, rechity );
				      fillME(totalMEs->meErrxZpPanel2, rechiterrx );
				      fillME(totalMEs->meErryZpPanel2, rechiterry );
				      fillME(totalMEs->meResxZpPanel2, rechitresx );
				      fillME(totalMEs->meResyZpPanel2, rechitresy );
				      fillME(totalMEs->mePullxZpPanel2, rechitpullx );
				      fillME(totalMEs->mePullyZpPanel2, rechitpully );
				      fillME(totalMEs->meNpixZpPanel2, npix );
				      fillME(totalMEs->meNxpixZpPanel2, nxpix );
				      fillME(totalMEs->meNypixZpPanel2, nypix );
				      fillME(totalMEs->meChargeZpPanel2, charge );


				      fillME(totalMEs->meResXvsAlphaZpPanel2, alpha, fabs(rechitresx) );
				      fillME(totalMEs->meResYvsAlphaZpPanel2, alpha, fabs(rechitresy) );
				      fillME(totalMEs->meResXvsBetaZpPanel2, beta, fabs(rechitresx) );
				      fillME(totalMEs->meResYvsBetaZpPanel2, beta, fabs(rechitresy) );
				      fillME(totalMEs->mePullXvsAlphaZpPanel2, alpha, rechitpullx );
				      fillME(totalMEs->mePullYvsAlphaZpPanel2, alpha, rechitpully );
				      fillME(totalMEs->mePullXvsBetaZpPanel2, beta, rechitpullx );
				      fillME(totalMEs->mePullYvsBetaZpPanel2, beta, rechitpully );
				      fillME(totalMEs->mePullXvsPhiZpPanel2, phi, rechitpullx );
				      fillME(totalMEs->mePullYvsPhiZpPanel2, phi, rechitpully );
				      fillME(totalMEs->mePullXvsEtaZpPanel2, eta, rechitpullx );
				      fillME(totalMEs->mePullYvsEtaZpPanel2, eta, rechitpully );

				      fillME(totalMEs->meWPullXvsAlphaZpPanel2, alpha, fabs(rechitpullx) );
				      fillME(totalMEs->meWPullYvsAlphaZpPanel2, alpha, fabs(rechitpully) );
				      fillME(totalMEs->meWPullXvsBetaZpPanel2, beta, fabs(rechitpullx) );
				      fillME(totalMEs->meWPullYvsBetaZpPanel2, beta, fabs(rechitpully) );

                                      labelss.str("");
                                      labelss << "Disk" << "Plaq_"<<disk<< "_" <<plaq;

                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePosxZpPanel2, rechitx);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePosyZpPanel2, rechity);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meErrxZpPanel2, rechiterrx);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meErryZpPanel2, rechiterry);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meResxZpPanel2, rechitresx);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meResyZpPanel2, rechitresy);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullxZpPanel2, rechitpullx);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullyZpPanel2, rechitpully);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meNpixZpPanel2, npix);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meNxpixZpPanel2, nxpix);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meNypixZpPanel2, nypix);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meChargeZpPanel2, charge);

                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meResXvsAlphaZpPanel2, alpha, fabs(rechitresx));
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meResYvsAlphaZpPanel2, alpha, fabs(rechitresy));
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meResXvsBetaZpPanel2, beta, fabs(rechitresx));
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->meResYvsBetaZpPanel2, beta, fabs(rechitresy));
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullXvsAlphaZpPanel2, alpha, rechitpullx);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullYvsAlphaZpPanel2, alpha, rechitpully);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullXvsBetaZpPanel2, beta, rechitpullx);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullYvsBetaZpPanel2, beta, rechitpully);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullXvsPhiZpPanel2, phi, rechitpullx);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullYvsPhiZpPanel2, phi, rechitpully);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullXvsEtaZpPanel2, eta, rechitpullx);
                                      fillME(diskPlaquettesMEsMap[labelss.str()]->mePullYvsEtaZpPanel2, eta, rechitpully);


				    }
				  else edm::LogWarning("SiPixelTrackingRecHitsValid") << "..............................................Wrong panel number !"; 
				} //else if ( side==2 )
			      else edm::LogWarning("SiPixelTrackingRecHitsValid") << ".......................................................Wrong side !" ;
			      
			    } // else if ( detId.subdetId()==PixelSubdetector::PixelEndcap )
			  else edm::LogWarning("SiPixelTrackingRecHitsValid") << "Pixel rechit but we are not in the pixel detector" << (int)detId.subdetId() ;
			  
			  if(debugNtuple_.size()!=0)t_->Fill();

			} // if ( !matched.empty() )
		      //else
		      //cout << "---------------- RecHit with no associated SimHit !!! -------------------------- " << endl;
		      
		    } // matchedhit.
		  
		} // end of loop on hits
	      
	      fillME(totalMEs->mePixRecHitsPerTrack, n_hits );
	      //cout << "n_hits = " << n_hits << endl;
	      
	    } //end of loop on track 
	  
	} // tracks > 0.
      
    } //end of MTCCTrack

}
