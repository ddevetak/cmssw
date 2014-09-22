
// File: SiPixelTrackingRecHitsValid.hh
// // Authors:  Xingtao Huang (Puerto Rico Univ.)
//              Gavril Giurgiu (JHU)
// Creation Date:  Oct. 2006.
//
//--------------------------------------------

#ifndef Validation_RecoTrack_SiPixelTrackingRecHitsValid_h
#define Validation_RecoTrack_SiPixelTrackingRecHitsValid_h

//DQM services for histogram
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "TrackingTools/TrackFitters/interface/KFTrajectoryFitter.h"
#include "TrackingTools/TrackFitters/interface/KFTrajectorySmoother.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h" 
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
//#include "Validation/RecoTrack/interface/TrackLocalAngle.h"
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TProfile.h>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//--- for SimHit association
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"  
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h" 
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h" 
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h" 
#include "Geometry/TrackerGeometryBuilder/interface/GluedGeomDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"

#include <string>

class TTree;
class TFile;

class SiPixelTrackingRecHitsValid : public thread_unsafe::DQMEDAnalyzer
{
 public:
  
  explicit SiPixelTrackingRecHitsValid(const edm::ParameterSet& conf);
  
  virtual ~SiPixelTrackingRecHitsValid();

  virtual void analyze(const edm::Event& e, const edm::EventSetup& c);
  void bookHistograms(DQMStore::IBooker & ibooker,const edm::Run& run, const edm::EventSetup& es);
  virtual void beginJob();
  virtual void endJob();

  //xt std::pair<LocalPoint,LocalVector> projectHit( const PSimHit& hit, const StripGeomDetUnit* stripDet,const BoundPlane& plane);
  std::pair<LocalPoint,LocalVector> projectHit( const PSimHit& hit, const PixelGeomDetUnit* pixelDet,const BoundPlane& plane);

  MonitorElement* bookME2D(DQMStore::IBooker &ibooker, const char* ParameterSetLabel, const char* HistoName, const char* HistoTitle, const char* ST);
  MonitorElement* bookME1D(DQMStore::IBooker &ibooker, const char* ParameterSetLabel, const char* HistoName, const char* HistoTitle);  

  inline void fillME(MonitorElement* ME,float value1,float value2){if (ME!=0) ME->Fill(value1,value2);}//need this function for using the switchon parameter
  inline void fillME(MonitorElement* ME,float value1){if (ME!=0) ME->Fill(value1);}

  void createLadderLayersMEs(DQMStore::IBooker &ibooker, std::string label);
  void createLayerModulesMEs(DQMStore::IBooker &ibooker, std::string label);
  void createDiskPlaquettesMEs(DQMStore::IBooker &ibooker, std::string label);
  void createLayersMEs(DQMStore::IBooker &ibooker, std::string label);
  void createTotalMEs(DQMStore::IBooker &ibooker, std::string label);

 private:

  edm::ParameterSet conf_;
  edm::ParameterSet Parameters;
  //TrackLocalAngle *anglefinder_;
  DQMStore* dbe_;
  bool runStandalone;
  std::string outputFile_;
  std::string debugNtuple_;
  std::string builderName_;
  edm::EDGetTokenT<SiPixelRecHitCollection> siPixelRecHitCollectionToken_;
  edm::EDGetTokenT<reco::TrackCollection> recoTrackCollectionToken_;
  bool MTCCtrack_;

  bool checkType_; // do we check that the simHit associated with recHit is of the expected particle type ?
  int genType_; // the type of particle that the simHit associated with recHits should be

  //SWITCH BOOL PARAMETERS---Layers
   
  bool switchResxBarrelLayer;
  bool switchResyBarrelLayer;
  bool switchPullxBarrelLayer;
  bool switchPullyBarrelLayer;


  //SWITCH BOOL PARAMETERS--FlippedLadders
  bool switchResXvsAlphaBarrelFlippedLaddersLayer;
  bool switchResYvsAlphaBarrelFlippedLaddersLayer;
  bool switchResXvsBetaBarrelFlippedLaddersLayer;
  bool switchResYvsBetaBarrelFlippedLaddersLayer;

  //SWITCH BOOL PARAMETERS--NonFlippedLadders
  bool switchResXvsAlphaBarrelNonFlippedLaddersLayer;
  bool switchResYvsAlphaBarrelNonFlippedLaddersLayer;
  bool switchResXvsBetaBarrelNonFlippedLaddersLayer;
  bool switchResYvsBetaBarrelNonFlippedLaddersLayer;

  //SWITCH BOOL PARAMETERS--BarrelLayerModule
  bool switchPosxBarrelLayerModule;
  bool switchPosyBarrelLayerModule;
  bool switchErrxBarrelLayerModule;
  bool switchErryBarrelLayerModule;
  bool switchResxBarrelLayerModule;
  bool switchResyBarrelLayerModule;
  bool switchPullxBarrelLayerModule;
  bool switchPullyBarrelLayerModule;
  bool switchNpixBarrelLayerModule;
  bool switchNxpixBarrelLayerModule;
  bool switchNypixBarrelLayerModule;
  bool switchChargeBarrelLayerModule;
  bool switchResXvsAlphaBarrelLayerModule;
  bool switchResYvsAlphaBarrelLayerModule;
  bool switchResXvsBetaBarrelLayerModule;
  bool switchResYvsBetaBarrelLayerModule;
  bool switchPullXvsAlphaBarrelLayerModule;
  bool switchPullYvsAlphaBarrelLayerModule;
  bool switchPullXvsBetaBarrelLayerModule;
  bool switchPullYvsBetaBarrelLayerModule;
  bool switchPullXvsPhiBarrelLayerModule;
  bool switchPullYvsPhiBarrelLayerModule;
  bool switchPullXvsEtaBarrelLayerModule;
  bool switchPullYvsEtaBarrelLayerModule;

  // SWITCH BOOL PARAMETERS--TOTAL
  bool switchPosxBarrelTotal;
  bool switchPosyBarrelTotal;
  bool switchErrxBarrelTotal;
  bool switchErryBarrelTotal;
  bool switchResxBarrelTotal;
  bool switchResyBarrelTotal;
  bool switchPullxBarrelTotal;
  bool switchPullyBarrelTotal;
  bool switchNpixBarrelTotal;
  bool switchNxpixBarrelTotal;
  bool switchNypixBarrelTotal;
  bool switchChargeBarrelTotal;
  bool switchResXvsAlphaBarrelTotal;
  bool switchResYvsAlphaBarrelTotal;
  bool switchResXvsBetaBarrelTotal;
  bool switchResYvsBetaBarrelTotal;
  bool switchPullXvsAlphaBarrelTotal;
  bool switchPullYvsAlphaBarrelTotal;
  bool switchPullXvsBetaBarrelTotal;
  bool switchPullYvsBetaBarrelTotal;
  bool switchPullXvsPhiBarrelTotal;
  bool switchPullYvsPhiBarrelTotal;
  bool switchPullXvsEtaBarrelTotal;
  bool switchPullYvsEtaBarrelTotal;

  bool switchPosxBarrelHalfModuleTotal;
  bool switchPosxBarrelFullModuleTotal;
  bool switchPosxBarrelFlippedLaddersTotal;
  bool switchPosxBarrelNonFlippedLaddersTotal;
  bool switchPosyBarrelHalfModuleTotal;
  bool switchPosyBarrelFullModuleTotal;
  bool switchPosyBarrelFlippedLaddersTotal;
  bool switchPosyBarrelNonFlippedLaddersTotal;
  
  bool switchResXvsAlphaBarrelFlippedLaddersTotal;
  bool switchResYvsAlphaBarrelFlippedLaddersTotal;
  bool switchResXvsBetaBarrelFlippedLaddersTotal;
  bool switchResYvsBetaBarrelFlippedLaddersTotal;

  bool switchPullXvsAlphaBarrelFlippedLaddersTotal;
  bool switchPullYvsAlphaBarrelFlippedLaddersTotal;
  bool switchPullXvsBetaBarrelFlippedLaddersTotal;
  bool switchPullYvsBetaBarrelFlippedLaddersTotal;
  bool switchPullXvsPhiBarrelFlippedLaddersTotal;
  bool switchPullYvsPhiBarrelFlippedLaddersTotal;
  bool switchPullXvsEtaBarrelFlippedLaddersTotal;
  bool switchPullYvsEtaBarrelFlippedLaddersTotal;

  bool switchWPullXvsAlphaBarrelFlippedLaddersTotal;
  bool switchWPullYvsAlphaBarrelFlippedLaddersTotal;
  bool switchWPullXvsBetaBarrelFlippedLaddersTotal;
  bool switchWPullYvsBetaBarrelFlippedLaddersTotal;
  bool switchResXvsAlphaBarrelNonFlippedLaddersTotal;
  bool switchResYvsAlphaBarrelNonFlippedLaddersTotal;
  bool switchResXvsBetaBarrelNonFlippedLaddersTotal;
  bool switchResYvsBetaBarrelNonFlippedLaddersTotal;

  bool switchPullXvsAlphaBarrelNonFlippedLaddersTotal;
  bool switchPullYvsAlphaBarrelNonFlippedLaddersTotal;
  bool switchPullXvsBetaBarrelNonFlippedLaddersTotal;
  bool switchPullYvsBetaBarrelNonFlippedLaddersTotal;
  bool switchPullXvsPhiBarrelNonFlippedLaddersTotal;
  bool switchPullYvsPhiBarrelNonFlippedLaddersTotal;
  bool switchPullXvsEtaBarrelNonFlippedLaddersTotal;
  bool switchPullYvsEtaBarrelNonFlippedLaddersTotal;

  bool switchWPullXvsAlphaBarrelNonFlippedLaddersTotal;
  bool switchWPullYvsAlphaBarrelNonFlippedLaddersTotal;
  bool switchWPullXvsBetaBarrelNonFlippedLaddersTotal;
  bool switchWPullYvsBetaBarrelNonFlippedLaddersTotal;

  bool switchPosxZmPanel1Total;
  bool switchPosyZmPanel1Total;
  bool switchErrxZmPanel1Total;
  bool switchErryZmPanel1Total;
  bool switchResxZmPanel1Total;
  bool switchResyZmPanel1Total;
  bool switchPullxZmPanel1Total;
  bool switchPullyZmPanel1Total;
  bool switchNpixZmPanel1Total;
  bool switchNxpixZmPanel1Total;
  bool switchNypixZmPanel1Total;
  bool switchChargeZmPanel1Total;
  bool switchResXvsAlphaZmPanel1Total;
  bool switchResYvsAlphaZmPanel1Total;
  bool switchResXvsBetaZmPanel1Total;
  bool switchResYvsBetaZmPanel1Total;

  bool switchPullXvsAlphaZmPanel1Total;
  bool switchPullYvsAlphaZmPanel1Total;
  bool switchPullXvsBetaZmPanel1Total;
  bool switchPullYvsBetaZmPanel1Total;
  bool switchPullXvsPhiZmPanel1Total;
  bool switchPullYvsPhiZmPanel1Total;
  bool switchPullXvsEtaZmPanel1Total;
  bool switchPullYvsEtaZmPanel1Total;

  bool switchWPullXvsAlphaZmPanel1Total;
  bool switchWPullYvsAlphaZmPanel1Total;
  bool switchWPullXvsBetaZmPanel1Total;
  bool switchWPullYvsBetaZmPanel1Total;

  bool switchPosxZpPanel1Total;
  bool switchPosyZpPanel1Total;
  bool switchErrxZpPanel1Total;
  bool switchErryZpPanel1Total;
  bool switchResxZpPanel1Total;
  bool switchResyZpPanel1Total;
  bool switchPullxZpPanel1Total;
  bool switchPullyZpPanel1Total;
  bool switchNpixZpPanel1Total;
  bool switchNxpixZpPanel1Total;
  bool switchNypixZpPanel1Total;
  bool switchChargeZpPanel1Total;
  bool switchResXvsAlphaZpPanel1Total;
  bool switchResYvsAlphaZpPanel1Total;
  bool switchResXvsBetaZpPanel1Total;
  bool switchResYvsBetaZpPanel1Total;

  bool switchPullXvsAlphaZpPanel1Total;
  bool switchPullYvsAlphaZpPanel1Total;
  bool switchPullXvsBetaZpPanel1Total;
  bool switchPullYvsBetaZpPanel1Total;
  bool switchPullXvsPhiZpPanel1Total;
  bool switchPullYvsPhiZpPanel1Total;
  bool switchPullXvsEtaZpPanel1Total;
  bool switchPullYvsEtaZpPanel1Total;

  bool switchWPullXvsAlphaZpPanel1Total;
  bool switchWPullYvsAlphaZpPanel1Total;
  bool switchWPullXvsBetaZpPanel1Total;
  bool switchWPullYvsBetaZpPanel1Total;
   
  bool switchPosxZmPanel2Total;
  bool switchPosyZmPanel2Total;
  bool switchErrxZmPanel2Total;
  bool switchErryZmPanel2Total;
  bool switchResxZmPanel2Total;
  bool switchResyZmPanel2Total;
  bool switchPullxZmPanel2Total;
  bool switchPullyZmPanel2Total;
  bool switchNpixZmPanel2Total;
  bool switchNxpixZmPanel2Total;
  bool switchNypixZmPanel2Total;
  bool switchChargeZmPanel2Total;

  bool switchResXvsAlphaZmPanel2Total;
  bool switchResYvsAlphaZmPanel2Total;
  bool switchResXvsBetaZmPanel2Total;
  bool switchResYvsBetaZmPanel2Total;
  bool switchPullXvsAlphaZmPanel2Total;
  bool switchPullYvsAlphaZmPanel2Total;
  bool switchPullXvsBetaZmPanel2Total;
  bool switchPullYvsBetaZmPanel2Total;
  bool switchPullXvsPhiZmPanel2Total;
  bool switchPullYvsPhiZmPanel2Total;
  bool switchPullXvsEtaZmPanel2Total;
  bool switchPullYvsEtaZmPanel2Total;

  bool switchPosxZpPanel2Total;
  bool switchPosyZpPanel2Total;
  bool switchErrxZpPanel2Total;
  bool switchErryZpPanel2Total;
  bool switchResxZpPanel2Total;
  bool switchResyZpPanel2Total;
  bool switchPullxZpPanel2Total;
  bool switchPullyZpPanel2Total;
  bool switchNpixZpPanel2Total;
  bool switchNxpixZpPanel2Total;
  bool switchNypixZpPanel2Total;
  bool switchChargeZpPanel2Total;

  bool switchWPullXvsAlphaZmPanel2Total;
  bool switchWPullYvsAlphaZmPanel2Total;
  bool switchWPullXvsBetaZmPanel2Total;
  bool switchWPullYvsBetaZmPanel2Total;

  bool switchResXvsAlphaZpPanel2Total;
  bool switchResYvsAlphaZpPanel2Total;
  bool switchResXvsBetaZpPanel2Total;
  bool switchResYvsBetaZpPanel2Total;
  bool switchPullXvsAlphaZpPanel2Total;
  bool switchPullYvsAlphaZpPanel2Total;
  bool switchPullXvsBetaZpPanel2Total;
  bool switchPullYvsBetaZpPanel2Total;
  bool switchPullXvsPhiZpPanel2Total;
  bool switchPullYvsPhiZpPanel2Total;
  bool switchPullXvsEtaZpPanel2Total;
  bool switchPullYvsEtaZpPanel2Total;


  bool switchWPullXvsAlphaZpPanel2Total;
  bool switchWPullYvsAlphaZpPanel2Total;
  bool switchWPullXvsBetaZpPanel2Total;
  bool switchWPullYvsBetaZpPanel2Total;

  bool switchPosxBarrel_all_hits;
  bool switchPosyBarrel_all_hits;

  bool switchPosxZmPanel1_all_hits;
  bool switchPosyZmPanel1_all_hits;
  bool switchPosxZmPanel2_all_hits;
  bool switchPosyZmPanel2_all_hits;
  
  bool switchPosxZpPanel1_all_hits;
  bool switchPosyZpPanel1_all_hits;
  bool switchPosxZpPanel2_all_hits;
  bool switchPosyZpPanel2_all_hits;

  bool switchTracksPerEvent;
  bool switchPixRecHitsPerTrack;
  



///////////////////////////////////////////////////////

  //SWITCH BOOL PARAMETERS---DiskPlaquette
  /////Zm-Panel1
  bool switchPosxZmPanel1DiskPlaquette;
  bool switchPosyZmPanel1DiskPlaquette;
  bool switchErrxZmPanel1DiskPlaquette;
  bool switchErryZmPanel1DiskPlaquette;
  bool switchResxZmPanel1DiskPlaquette;
  bool switchResyZmPanel1DiskPlaquette;
  bool switchPullxZmPanel1DiskPlaquette;
  bool switchPullyZmPanel1DiskPlaquette;
  bool switchNpixZmPanel1DiskPlaquette;
  bool switchNxpixZmPanel1DiskPlaquette;
  bool switchNypixZmPanel1DiskPlaquette;
  bool switchChargeZmPanel1DiskPlaquette;
  bool switchResXvsAlphaZmPanel1DiskPlaquette;
  bool switchResYvsAlphaZmPanel1DiskPlaquette;
  bool switchResXvsBetaZmPanel1DiskPlaquette;
  bool switchResYvsBetaZmPanel1DiskPlaquette;
  bool switchPullXvsAlphaZmPanel1DiskPlaquette;
  bool switchPullYvsAlphaZmPanel1DiskPlaquette;
  bool switchPullXvsBetaZmPanel1DiskPlaquette;
  bool switchPullYvsBetaZmPanel1DiskPlaquette;
  bool switchPullXvsPhiZmPanel1DiskPlaquette;
  bool switchPullYvsPhiZmPanel1DiskPlaquette;
  bool switchPullXvsEtaZmPanel1DiskPlaquette;
  bool switchPullYvsEtaZmPanel1DiskPlaquette;
  //////Zp-Panel1
  bool switchPosxZpPanel1DiskPlaquette;
  bool switchPosyZpPanel1DiskPlaquette;
  bool switchErrxZpPanel1DiskPlaquette;
  bool switchErryZpPanel1DiskPlaquette;
  bool switchResxZpPanel1DiskPlaquette;
  bool switchResyZpPanel1DiskPlaquette;
  bool switchPullxZpPanel1DiskPlaquette;
  bool switchPullyZpPanel1DiskPlaquette;
  bool switchNpixZpPanel1DiskPlaquette;
  bool switchNxpixZpPanel1DiskPlaquette;
  bool switchNypixZpPanel1DiskPlaquette;
  bool switchChargeZpPanel1DiskPlaquette;
  bool switchResXvsAlphaZpPanel1DiskPlaquette;
  bool switchResYvsAlphaZpPanel1DiskPlaquette;
  bool switchResXvsBetaZpPanel1DiskPlaquette;
  bool switchResYvsBetaZpPanel1DiskPlaquette;
  bool switchPullXvsAlphaZpPanel1DiskPlaquette;
  bool switchPullYvsAlphaZpPanel1DiskPlaquette;
  bool switchPullXvsBetaZpPanel1DiskPlaquette;
  bool switchPullYvsBetaZpPanel1DiskPlaquette;
  bool switchPullXvsPhiZpPanel1DiskPlaquette;
  bool switchPullYvsPhiZpPanel1DiskPlaquette;
  bool switchPullXvsEtaZpPanel1DiskPlaquette;
  bool switchPullYvsEtaZpPanel1DiskPlaquette;
  //Zm-Panel2
  bool switchPosxZmPanel2DiskPlaquette;
  bool switchPosyZmPanel2DiskPlaquette;
  bool switchErrxZmPanel2DiskPlaquette;
  bool switchErryZmPanel2DiskPlaquette;
  bool switchResxZmPanel2DiskPlaquette;
  bool switchResyZmPanel2DiskPlaquette;
  bool switchPullxZmPanel2DiskPlaquette;
  bool switchPullyZmPanel2DiskPlaquette;
  bool switchNpixZmPanel2DiskPlaquette;
  bool switchNxpixZmPanel2DiskPlaquette;
  bool switchNypixZmPanel2DiskPlaquette;
  bool switchChargeZmPanel2DiskPlaquette;
  bool switchResXvsAlphaZmPanel2DiskPlaquette;
  bool switchResYvsAlphaZmPanel2DiskPlaquette;
  bool switchResXvsBetaZmPanel2DiskPlaquette;
  bool switchResYvsBetaZmPanel2DiskPlaquette;
  bool switchPullXvsAlphaZmPanel2DiskPlaquette;
  bool switchPullYvsAlphaZmPanel2DiskPlaquette;
  bool switchPullXvsBetaZmPanel2DiskPlaquette;
  bool switchPullYvsBetaZmPanel2DiskPlaquette;
  bool switchPullXvsPhiZmPanel2DiskPlaquette;
  bool switchPullYvsPhiZmPanel2DiskPlaquette;
  bool switchPullXvsEtaZmPanel2DiskPlaquette;
  bool switchPullYvsEtaZmPanel2DiskPlaquette;
  //////Zp-Panel2
  bool switchPosxZpPanel2DiskPlaquette;
  bool switchPosyZpPanel2DiskPlaquette;
  bool switchErrxZpPanel2DiskPlaquette;
  bool switchErryZpPanel2DiskPlaquette;
  bool switchResxZpPanel2DiskPlaquette;
  bool switchResyZpPanel2DiskPlaquette;
  bool switchPullxZpPanel2DiskPlaquette;
  bool switchPullyZpPanel2DiskPlaquette;
  bool switchNpixZpPanel2DiskPlaquette;
  bool switchNxpixZpPanel2DiskPlaquette;
  bool switchNypixZpPanel2DiskPlaquette;
  bool switchChargeZpPanel2DiskPlaquette;
  bool switchResXvsAlphaZpPanel2DiskPlaquette;
  bool switchResYvsAlphaZpPanel2DiskPlaquette;
  bool switchResXvsBetaZpPanel2DiskPlaquette;
  bool switchResYvsBetaZpPanel2DiskPlaquette;
  bool switchPullXvsAlphaZpPanel2DiskPlaquette;
  bool switchPullYvsAlphaZpPanel2DiskPlaquette;
  bool switchPullXvsBetaZpPanel2DiskPlaquette;
  bool switchPullYvsBetaZpPanel2DiskPlaquette;
  bool switchPullXvsPhiZpPanel2DiskPlaquette;
  bool switchPullYvsPhiZpPanel2DiskPlaquette;
  bool switchPullXvsEtaZpPanel2DiskPlaquette;
  bool switchPullYvsEtaZpPanel2DiskPlaquette;

  


  // Pixel barrel detector has 3 layers and 8 modules; book histograms for each module = (layer, ring) pair
 
  struct LayerMEs{

      MonitorElement* meResxBarrel;
      MonitorElement* meResyBarrel;
      MonitorElement* mePullxBarrel;
      MonitorElement* mePullyBarrel;
   
      LayerMEs() { meResxBarrel=0; meResyBarrel=0; mePullxBarrel=0; mePullyBarrel=0;};

      

  };

  struct LadderLayersMEs{

    MonitorElement* meResXvsAlphaBarrel;
    MonitorElement* meResYvsAlphaBarrel;
    MonitorElement* meResXvsBetaBarrel;
    MonitorElement* meResYvsBetaBarrel;

    LadderLayersMEs() {meResXvsAlphaBarrel = 0; meResYvsAlphaBarrel = 0; meResXvsBetaBarrel = 0; meResYvsBetaBarrel = 0;};

  };

  struct LayerModulesMEs{
    
    MonitorElement* mePosxBarrel;
    MonitorElement* mePosyBarrel;
    MonitorElement* meErrxBarrel;
    MonitorElement* meErryBarrel;
    MonitorElement* meResxBarrel;
    MonitorElement* meResyBarrel;
    MonitorElement* mePullxBarrel;
    MonitorElement* mePullyBarrel;
    MonitorElement* meNpixBarrel;
    MonitorElement* meNxpixBarrel;
    MonitorElement* meNypixBarrel;
    MonitorElement* meChargeBarrel;
    MonitorElement* meResXvsAlphaBarrel;
    MonitorElement* meResYvsAlphaBarrel;
    MonitorElement* meResXvsBetaBarrel;
    MonitorElement* meResYvsBetaBarrel;
    MonitorElement* mePullXvsAlphaBarrel;
    MonitorElement* mePullYvsAlphaBarrel;
    MonitorElement* mePullXvsBetaBarrel;
    MonitorElement* mePullYvsBetaBarrel;
    MonitorElement* mePullXvsPhiBarrel;
    MonitorElement* mePullYvsPhiBarrel;
    MonitorElement* mePullXvsEtaBarrel;
    MonitorElement* mePullYvsEtaBarrel;


    LayerModulesMEs() { 

                mePosxBarrel=0;  
                mePosyBarrel=0; 
                meErrxBarrel=0; 
                meErryBarrel=0; 
                meResxBarrel=0; 
                meResyBarrel=0; 
                mePullxBarrel=0; 
                mePullyBarrel=0; 
                meNpixBarrel=0; 
                meNxpixBarrel=0; 
                meChargeBarrel=0; 
                meResXvsAlphaBarrel=0;
                meResYvsAlphaBarrel=0;
                meResXvsBetaBarrel=0;
                meResYvsBetaBarrel=0;
                mePullXvsAlphaBarrel=0;
                mePullYvsAlphaBarrel=0;
                mePullXvsBetaBarrel=0;
                mePullYvsBetaBarrel=0;
                mePullXvsPhiBarrel=0;
                mePullYvsPhiBarrel=0;
                mePullXvsEtaBarrel=0;
                mePullYvsEtaBarrel=0;

    };

  };


  struct DiskPlaquettesMEs {
   //Zm-Panel1
   MonitorElement*mePosxZmPanel1;
   MonitorElement*mePosyZmPanel1;
   MonitorElement*meErrxZmPanel1;
   MonitorElement*meErryZmPanel1;
   MonitorElement*meResxZmPanel1;
   MonitorElement*meResyZmPanel1;
   MonitorElement*mePullxZmPanel1;
   MonitorElement*mePullyZmPanel1;
   MonitorElement*meNpixZmPanel1;
   MonitorElement*meNxpixZmPanel1;
   MonitorElement*meNypixZmPanel1;
   MonitorElement*meChargeZmPanel1;
   MonitorElement*meResXvsAlphaZmPanel1;
   MonitorElement*meResYvsAlphaZmPanel1;
   MonitorElement*meResXvsBetaZmPanel1;
   MonitorElement*meResYvsBetaZmPanel1;
   MonitorElement*mePullXvsAlphaZmPanel1;
   MonitorElement*mePullYvsAlphaZmPanel1;
   MonitorElement*mePullXvsBetaZmPanel1;
   MonitorElement*mePullYvsBetaZmPanel1;
   MonitorElement*mePullXvsPhiZmPanel1;
   MonitorElement*mePullYvsPhiZmPanel1;
   MonitorElement*mePullXvsEtaZmPanel1;
   MonitorElement*mePullYvsEtaZmPanel1;

   //Zp-Panel1
   MonitorElement*mePosxZpPanel1;
   MonitorElement*mePosyZpPanel1;
   MonitorElement*meErrxZpPanel1;
   MonitorElement*meErryZpPanel1;
   MonitorElement*meResxZpPanel1;
   MonitorElement*meResyZpPanel1;
   MonitorElement*mePullxZpPanel1;
   MonitorElement*mePullyZpPanel1;
   MonitorElement*meNpixZpPanel1;
   MonitorElement*meNxpixZpPanel1;
   MonitorElement*meNypixZpPanel1;
   MonitorElement*meChargeZpPanel1;
   MonitorElement*meResXvsAlphaZpPanel1;
   MonitorElement*meResYvsAlphaZpPanel1;
   MonitorElement*meResXvsBetaZpPanel1;
   MonitorElement*meResYvsBetaZpPanel1;
   MonitorElement*mePullXvsAlphaZpPanel1;
   MonitorElement*mePullYvsAlphaZpPanel1;
   MonitorElement*mePullXvsBetaZpPanel1;
   MonitorElement*mePullYvsBetaZpPanel1;
   MonitorElement*mePullXvsPhiZpPanel1;
   MonitorElement*mePullYvsPhiZpPanel1;
   MonitorElement*mePullXvsEtaZpPanel1;
   MonitorElement*mePullYvsEtaZpPanel1;

   //Zm-Panel2
   MonitorElement*mePosxZmPanel2;
   MonitorElement*mePosyZmPanel2;
   MonitorElement*meErrxZmPanel2;
   MonitorElement*meErryZmPanel2;
   MonitorElement*meResxZmPanel2;
   MonitorElement*meResyZmPanel2;
   MonitorElement*mePullxZmPanel2;
   MonitorElement*mePullyZmPanel2;
   MonitorElement*meNpixZmPanel2;
   MonitorElement*meNxpixZmPanel2;
   MonitorElement*meNypixZmPanel2;
   MonitorElement*meChargeZmPanel2;
   MonitorElement*meResXvsAlphaZmPanel2;
   MonitorElement*meResYvsAlphaZmPanel2;
   MonitorElement*meResXvsBetaZmPanel2;
   MonitorElement*meResYvsBetaZmPanel2;
   MonitorElement*mePullXvsAlphaZmPanel2;
   MonitorElement*mePullYvsAlphaZmPanel2;
   MonitorElement*mePullXvsBetaZmPanel2;
   MonitorElement*mePullYvsBetaZmPanel2;
   MonitorElement*mePullXvsPhiZmPanel2;
   MonitorElement*mePullYvsPhiZmPanel2;
   MonitorElement*mePullXvsEtaZmPanel2;
   MonitorElement*mePullYvsEtaZmPanel2;

   //Zp-Panel2
   MonitorElement*mePosxZpPanel2;
   MonitorElement*mePosyZpPanel2;
   MonitorElement*meErrxZpPanel2;
   MonitorElement*meErryZpPanel2;
   MonitorElement*meResxZpPanel2;
   MonitorElement*meResyZpPanel2;
   MonitorElement*mePullxZpPanel2;
   MonitorElement*mePullyZpPanel2;
   MonitorElement*meNpixZpPanel2;
   MonitorElement*meNxpixZpPanel2;
   MonitorElement*meNypixZpPanel2;
   MonitorElement*meChargeZpPanel2;

   MonitorElement*meResXvsAlphaZpPanel2;
   MonitorElement*meResYvsAlphaZpPanel2;
   MonitorElement*meResXvsBetaZpPanel2;
   MonitorElement*meResYvsBetaZpPanel2;
   MonitorElement*mePullXvsAlphaZpPanel2;
   MonitorElement*mePullYvsAlphaZpPanel2;
   MonitorElement*mePullXvsBetaZpPanel2;
   MonitorElement*mePullYvsBetaZpPanel2;
   MonitorElement*mePullXvsPhiZpPanel2;
   MonitorElement*mePullYvsPhiZpPanel2;
   MonitorElement*mePullXvsEtaZpPanel2;
   MonitorElement*mePullYvsEtaZpPanel2;
  

  DiskPlaquettesMEs(){
             
               //Zm-Panel1
               mePosxZmPanel1=0;
               mePosyZmPanel1=0;
               meErrxZmPanel1=0;
               meErryZmPanel1=0;
               meResxZmPanel1=0;
               meResyZmPanel1=0;
               mePullxZmPanel1=0;
               mePullyZmPanel1=0;
               meNpixZmPanel1=0;
               meNxpixZmPanel1=0;
               meNypixZmPanel1=0;
               meChargeZmPanel1=0;
               meResXvsAlphaZmPanel1=0;
               meResYvsAlphaZmPanel1=0;
               meResXvsBetaZmPanel1=0;
               mePullXvsAlphaZmPanel1=0;
               mePullYvsAlphaZmPanel1=0;
               mePullXvsBetaZmPanel1=0;
               mePullYvsBetaZmPanel1=0;
               mePullXvsPhiZmPanel1=0;
               mePullYvsPhiZmPanel1=0;
               mePullXvsEtaZmPanel1=0;
               mePullYvsEtaZmPanel1=0;

               //Zp-Panel1
               mePosxZpPanel1=0;
               mePosyZpPanel1=0;
               meErrxZpPanel1=0;
               meErryZpPanel1=0;
               meResxZpPanel1=0;
               meResyZpPanel1=0;
               mePullxZpPanel1=0;
               mePullyZpPanel1=0;
               meNpixZpPanel1=0;
               meNxpixZpPanel1=0;
               meNypixZpPanel1=0;
               meChargeZpPanel1=0;
               meResXvsAlphaZpPanel1=0;
               meResYvsAlphaZpPanel1=0;
               meResXvsBetaZpPanel1=0;
               mePullXvsAlphaZpPanel1=0;
               mePullYvsAlphaZpPanel1=0;
               mePullXvsBetaZpPanel1=0;
               mePullYvsBetaZpPanel1=0;
               mePullXvsPhiZpPanel1=0;
               mePullYvsPhiZpPanel1=0;
               mePullXvsEtaZpPanel1=0;
               mePullYvsEtaZpPanel1=0;

              //Zm-Panel2
               mePosxZmPanel2=0;
               mePosyZmPanel2=0;
               meErrxZmPanel2=0;
               meErryZmPanel2=0;
               meResxZmPanel2=0;
               meResyZmPanel2=0;
               mePullxZmPanel2=0;
               mePullyZmPanel2=0;
               meNpixZmPanel2=0;
               meNxpixZmPanel2=0;
               meNypixZmPanel2=0;
               meChargeZmPanel2=0;
               meResXvsAlphaZmPanel2=0;
               meResYvsAlphaZmPanel2=0;
               meResXvsBetaZmPanel2=0;
               mePullXvsAlphaZmPanel2=0;
               mePullYvsAlphaZmPanel2=0;
               mePullXvsBetaZmPanel2=0;
               mePullYvsBetaZmPanel2=0;
               mePullXvsPhiZmPanel2=0;
               mePullYvsPhiZmPanel2=0;
               mePullXvsEtaZmPanel2=0;
               mePullYvsEtaZmPanel2=0;

                 //Zp-Panel2
               mePosxZpPanel2=0;
               mePosyZpPanel2=0;
               meErrxZpPanel2=0;
               meErryZpPanel2=0;
               meResxZpPanel2=0;
               meResyZpPanel2=0;
               mePullxZpPanel2=0;
               mePullyZpPanel2=0;
               meNpixZpPanel2=0;
               meNxpixZpPanel2=0;
               meNypixZpPanel2=0;
               meChargeZpPanel2=0;
               meResXvsAlphaZpPanel2=0;
               meResYvsAlphaZpPanel2=0;
               meResXvsBetaZpPanel2=0;
               mePullXvsAlphaZpPanel2=0;
               mePullYvsAlphaZpPanel2=0;
               mePullXvsBetaZpPanel2=0;
               mePullYvsBetaZpPanel2=0;
               mePullXvsPhiZpPanel2=0;
               mePullYvsPhiZpPanel2=0;
               mePullXvsEtaZpPanel2=0;
               mePullYvsEtaZpPanel2=0;


  };


  };

  struct TotalMEs{

    MonitorElement* mePosxBarrel;
    MonitorElement* mePosyBarrel;
    MonitorElement* meErrxBarrel;
    MonitorElement* meErryBarrel;
    MonitorElement* meResxBarrel;
    MonitorElement* meResyBarrel;
    MonitorElement* mePullxBarrel;
    MonitorElement* mePullyBarrel;
    MonitorElement* meNpixBarrel;
    MonitorElement* meNxpixBarrel;
    MonitorElement* meNypixBarrel;
    MonitorElement* meChargeBarrel;
    MonitorElement* meResXvsAlphaBarrel;
    MonitorElement* meResYvsAlphaBarrel;
    MonitorElement* meResXvsBetaBarrel;
    MonitorElement* meResYvsBetaBarrel;
    MonitorElement* mePullXvsAlphaBarrel;
    MonitorElement* mePullYvsAlphaBarrel;
    MonitorElement* mePullXvsBetaBarrel;
    MonitorElement* mePullYvsBetaBarrel;
    MonitorElement* mePullXvsPhiBarrel;
    MonitorElement* mePullYvsPhiBarrel;
    MonitorElement* mePullXvsEtaBarrel;
    MonitorElement* mePullYvsEtaBarrel;
    
    MonitorElement* mePosxBarrelHalfModule;
    MonitorElement* mePosxBarrelFullModule;
    MonitorElement* mePosxBarrelFlippedLadders;
    MonitorElement* mePosxBarrelNonFlippedLadders;
    MonitorElement* mePosyBarrelHalfModule;
    MonitorElement* mePosyBarrelFullModule;
    MonitorElement* mePosyBarrelFlippedLadders;
    MonitorElement* mePosyBarrelNonFlippedLadders;

    MonitorElement* meResXvsAlphaBarrelFlippedLadders;
    MonitorElement* meResYvsAlphaBarrelFlippedLadders;
    MonitorElement* meResXvsBetaBarrelFlippedLadders;
    MonitorElement* meResYvsBetaBarrelFlippedLadders;
 
    MonitorElement* mePullXvsAlphaBarrelFlippedLadders;
    MonitorElement* mePullYvsAlphaBarrelFlippedLadders;
    MonitorElement* mePullXvsBetaBarrelFlippedLadders;
    MonitorElement* mePullYvsBetaBarrelFlippedLadders;
    MonitorElement* mePullXvsPhiBarrelFlippedLadders;
    MonitorElement* mePullYvsPhiBarrelFlippedLadders;
    MonitorElement* mePullXvsEtaBarrelFlippedLadders;
    MonitorElement* mePullYvsEtaBarrelFlippedLadders;
   
    MonitorElement* meWPullXvsAlphaBarrelFlippedLadders;
    MonitorElement* meWPullYvsAlphaBarrelFlippedLadders;
    MonitorElement* meWPullXvsBetaBarrelFlippedLadders;
    MonitorElement* meWPullYvsBetaBarrelFlippedLadders;
    MonitorElement* meResXvsAlphaBarrelNonFlippedLadders;
    MonitorElement* meResYvsAlphaBarrelNonFlippedLadders;
    MonitorElement* meResXvsBetaBarrelNonFlippedLadders;
    MonitorElement* meResYvsBetaBarrelNonFlippedLadders;
   
    MonitorElement* mePullXvsAlphaBarrelNonFlippedLadders;
    MonitorElement* mePullYvsAlphaBarrelNonFlippedLadders;
    MonitorElement* mePullXvsBetaBarrelNonFlippedLadders;
    MonitorElement* mePullYvsBetaBarrelNonFlippedLadders;
    MonitorElement* mePullXvsPhiBarrelNonFlippedLadders;
    MonitorElement* mePullYvsPhiBarrelNonFlippedLadders;
    MonitorElement* mePullXvsEtaBarrelNonFlippedLadders;
    MonitorElement* mePullYvsEtaBarrelNonFlippedLadders;

    MonitorElement* meWPullXvsAlphaBarrelNonFlippedLadders;
    MonitorElement* meWPullYvsAlphaBarrelNonFlippedLadders;
    MonitorElement* meWPullXvsBetaBarrelNonFlippedLadders;
    MonitorElement* meWPullYvsBetaBarrelNonFlippedLadders;
   
    //Zm-Panel1
    MonitorElement* mePosxZmPanel1;
    MonitorElement* mePosyZmPanel1;
    MonitorElement* meErrxZmPanel1;
    MonitorElement* meErryZmPanel1;
    MonitorElement* meResxZmPanel1;
    MonitorElement* meResyZmPanel1;
    MonitorElement* mePullxZmPanel1;
    MonitorElement* mePullyZmPanel1;
    MonitorElement* meNpixZmPanel1;
    MonitorElement* meNxpixZmPanel1;
    MonitorElement* meNypixZmPanel1;
    MonitorElement* meChargeZmPanel1;
    MonitorElement* meResXvsAlphaZmPanel1;
    MonitorElement* meResYvsAlphaZmPanel1;
    MonitorElement* meResXvsBetaZmPanel1;
    MonitorElement* meResYvsBetaZmPanel1;

    MonitorElement* mePullXvsAlphaZmPanel1;
    MonitorElement* mePullYvsAlphaZmPanel1;
    MonitorElement* mePullXvsBetaZmPanel1;
    MonitorElement* mePullYvsBetaZmPanel1;
    MonitorElement* mePullXvsPhiZmPanel1;
    MonitorElement* mePullYvsPhiZmPanel1;
    MonitorElement* mePullXvsEtaZmPanel1;
    MonitorElement* mePullYvsEtaZmPanel1;

    MonitorElement* meWPullXvsAlphaZmPanel1;
    MonitorElement* meWPullYvsAlphaZmPanel1;
    MonitorElement* meWPullXvsBetaZmPanel1;
    MonitorElement* meWPullYvsBetaZmPanel1;

   //Zp-Panel1
    MonitorElement* mePosxZpPanel1;
    MonitorElement* mePosyZpPanel1;
    MonitorElement* meErrxZpPanel1;
    MonitorElement* meErryZpPanel1;
    MonitorElement* meResxZpPanel1;
    MonitorElement* meResyZpPanel1;
    MonitorElement* mePullxZpPanel1;
    MonitorElement* mePullyZpPanel1;
    MonitorElement* meNpixZpPanel1;
    MonitorElement* meNxpixZpPanel1;
    MonitorElement* meNypixZpPanel1;
    MonitorElement* meChargeZpPanel1;
    MonitorElement* meResXvsAlphaZpPanel1;
    MonitorElement* meResYvsAlphaZpPanel1;
    MonitorElement* meResXvsBetaZpPanel1;
    MonitorElement* meResYvsBetaZpPanel1;

    MonitorElement* mePullXvsAlphaZpPanel1;
    MonitorElement* mePullYvsAlphaZpPanel1;
    MonitorElement* mePullXvsBetaZpPanel1;
    MonitorElement* mePullYvsBetaZpPanel1;
    MonitorElement* mePullXvsPhiZpPanel1;
    MonitorElement* mePullYvsPhiZpPanel1;
    MonitorElement* mePullXvsEtaZpPanel1;
    MonitorElement* mePullYvsEtaZpPanel1;

    MonitorElement* meWPullXvsAlphaZpPanel1;
    MonitorElement* meWPullYvsAlphaZpPanel1;
    MonitorElement* meWPullXvsBetaZpPanel1;
    MonitorElement* meWPullYvsBetaZpPanel1;

    //Zm-Panel2

    MonitorElement* mePosxZmPanel2;
    MonitorElement* mePosyZmPanel2;
    MonitorElement* meErrxZmPanel2;
    MonitorElement* meErryZmPanel2;
    MonitorElement* meResxZmPanel2;
    MonitorElement* meResyZmPanel2;
    MonitorElement* mePullxZmPanel2;
    MonitorElement* mePullyZmPanel2;
    MonitorElement* meNpixZmPanel2;
    MonitorElement* meNxpixZmPanel2;
    MonitorElement* meNypixZmPanel2;
    MonitorElement* meChargeZmPanel2;

    MonitorElement* meResXvsAlphaZmPanel2;
    MonitorElement* meResYvsAlphaZmPanel2;
    MonitorElement* meResXvsBetaZmPanel2;
    MonitorElement* meResYvsBetaZmPanel2;
    MonitorElement* mePullXvsAlphaZmPanel2;
    MonitorElement* mePullYvsAlphaZmPanel2;
    MonitorElement* mePullXvsBetaZmPanel2;
    MonitorElement* mePullYvsBetaZmPanel2;
    MonitorElement* mePullXvsPhiZmPanel2;
    MonitorElement* mePullYvsPhiZmPanel2;
    MonitorElement* mePullXvsEtaZmPanel2;
    MonitorElement* mePullYvsEtaZmPanel2;

    MonitorElement* mePosxZpPanel2;
    MonitorElement* mePosyZpPanel2;
    MonitorElement* meErrxZpPanel2;
    MonitorElement* meErryZpPanel2;
    MonitorElement* meResxZpPanel2;
    MonitorElement* meResyZpPanel2;
    MonitorElement* mePullxZpPanel2;
    MonitorElement* mePullyZpPanel2;
    MonitorElement* meNpixZpPanel2;
    MonitorElement* meNxpixZpPanel2;
    MonitorElement* meNypixZpPanel2;
    MonitorElement* meChargeZpPanel2;

    MonitorElement* meWPullXvsAlphaZmPanel2;
    MonitorElement* meWPullYvsAlphaZmPanel2;
    MonitorElement* meWPullXvsBetaZmPanel2;
    MonitorElement* meWPullYvsBetaZmPanel2;

    //Zp-Panel2
    
    MonitorElement*meResXvsAlphaZpPanel2;
    MonitorElement*meResYvsAlphaZpPanel2;
    MonitorElement*meResXvsBetaZpPanel2;
    MonitorElement*meResYvsBetaZpPanel2;
    MonitorElement*mePullXvsAlphaZpPanel2;
    MonitorElement*mePullYvsAlphaZpPanel2;
    MonitorElement*mePullXvsBetaZpPanel2;
    MonitorElement*mePullYvsBetaZpPanel2;
    MonitorElement*mePullXvsPhiZpPanel2;
    MonitorElement*mePullYvsPhiZpPanel2;
    MonitorElement*mePullXvsEtaZpPanel2;
    MonitorElement*mePullYvsEtaZpPanel2;

    MonitorElement* meWPullXvsAlphaZpPanel2;
    MonitorElement* meWPullYvsAlphaZpPanel2;
    MonitorElement* meWPullXvsBetaZpPanel2;
    MonitorElement* meWPullYvsBetaZpPanel2;

    MonitorElement* mePosxBarrel_all_hits;
    MonitorElement* mePosyBarrel_all_hits;

    MonitorElement* mePosxZmPanel1_all_hits;
    MonitorElement* mePosyZmPanel1_all_hits;
    MonitorElement* mePosxZmPanel2_all_hits;
    MonitorElement* mePosyZmPanel2_all_hits;

    MonitorElement* mePosxZpPanel1_all_hits;
    MonitorElement* mePosyZpPanel1_all_hits;
    MonitorElement* mePosxZpPanel2_all_hits;
    MonitorElement* mePosyZpPanel2_all_hits;

    MonitorElement* meTracksPerEvent;
    MonitorElement* mePixRecHitsPerTrack;


  TotalMEs(){

          mePosxBarrel=0;
          mePosyBarrel=0;
          meErrxBarrel=0;
          meErryBarrel=0;
          meResxBarrel=0;
          meResyBarrel=0;
          mePullxBarrel=0;
          mePullyBarrel=0;
          meNpixBarrel=0;
          meNxpixBarrel=0;
          meChargeBarrel=0;
          meResXvsAlphaBarrel=0;
          meResYvsAlphaBarrel=0;
          meResXvsBetaBarrel=0;
          meResYvsBetaBarrel=0;
          mePullXvsAlphaBarrel=0;
          mePullYvsAlphaBarrel=0;
          mePullXvsBetaBarrel=0;
          mePullYvsBetaBarrel=0;
          mePullXvsPhiBarrel=0;
          mePullYvsPhiBarrel=0;
          mePullXvsEtaBarrel=0;
          mePullYvsEtaBarrel=0;
  
          mePosxBarrelHalfModule=0;
          mePosxBarrelFullModule=0;
          mePosxBarrelFlippedLadders=0;
          mePosxBarrelNonFlippedLadders=0;
          mePosyBarrelHalfModule=0;
          mePosyBarrelFullModule=0;
          mePosyBarrelFlippedLadders=0;
          mePosyBarrelNonFlippedLadders=0;
 
          meResXvsAlphaBarrelFlippedLadders=0;
          meResYvsAlphaBarrelFlippedLadders=0;
          meResXvsBetaBarrelFlippedLadders=0;
          meResYvsBetaBarrelFlippedLadders=0;
   
          mePullXvsAlphaBarrelFlippedLadders=0;
          mePullYvsAlphaBarrelFlippedLadders=0;
          mePullXvsBetaBarrelFlippedLadders=0;
          mePullYvsBetaBarrelFlippedLadders=0;
          mePullXvsPhiBarrelFlippedLadders=0;
          mePullYvsPhiBarrelFlippedLadders=0;
          mePullXvsEtaBarrelFlippedLadders=0;
          mePullYvsEtaBarrelFlippedLadders=0;

          meWPullXvsAlphaBarrelFlippedLadders=0;
          meWPullYvsAlphaBarrelFlippedLadders=0;
          meWPullXvsBetaBarrelFlippedLadders=0;
          meWPullYvsBetaBarrelFlippedLadders=0;
          meResXvsAlphaBarrelNonFlippedLadders=0;
          meResYvsAlphaBarrelNonFlippedLadders=0;
          meResXvsBetaBarrelNonFlippedLadders=0;
          meResYvsBetaBarrelNonFlippedLadders=0;

          mePullXvsAlphaBarrelNonFlippedLadders=0;
          mePullYvsAlphaBarrelNonFlippedLadders=0;
          mePullXvsBetaBarrelNonFlippedLadders=0;
          mePullYvsBetaBarrelNonFlippedLadders=0;
          mePullXvsPhiBarrelNonFlippedLadders=0;
          mePullYvsPhiBarrelNonFlippedLadders=0;
          mePullXvsEtaBarrelNonFlippedLadders=0;
          mePullYvsEtaBarrelNonFlippedLadders=0;

          meWPullXvsAlphaBarrelNonFlippedLadders=0;
          meWPullYvsAlphaBarrelNonFlippedLadders=0;
          meWPullXvsBetaBarrelNonFlippedLadders=0;
          meWPullYvsBetaBarrelNonFlippedLadders=0;

          mePosxZmPanel1=0;
          mePosyZmPanel1=0;
          meErrxZmPanel1=0;
          meErryZmPanel1=0;
          meResxZmPanel1=0;
          meResyZmPanel1=0;
          mePullxZmPanel1=0;
          mePullyZmPanel1=0;
          meNpixZmPanel1=0;
          meNxpixZmPanel1=0;
          meNypixZmPanel1=0;
          meChargeZmPanel1=0;
          meResXvsAlphaZmPanel1=0;
          meResYvsAlphaZmPanel1=0;
          meResXvsBetaZmPanel1=0;
          meResYvsBetaZmPanel1=0;

          mePullXvsAlphaZmPanel1=0;
          mePullYvsAlphaZmPanel1=0;
          mePullXvsBetaZmPanel1=0;
          mePullYvsBetaZmPanel1=0;
          mePullXvsPhiZmPanel1=0;
          mePullYvsPhiZmPanel1=0;
          mePullXvsEtaZmPanel1=0;
          mePullYvsEtaZmPanel1=0;

          meWPullXvsAlphaZmPanel1=0;
          meWPullYvsAlphaZmPanel1=0;
          meWPullXvsBetaZmPanel1=0;
          meWPullYvsBetaZmPanel1=0;

          mePosxZpPanel1=0;
          mePosyZpPanel1=0;
          meErrxZpPanel1=0;
          meErryZpPanel1=0;
          meResxZpPanel1=0;
          meResyZpPanel1=0;
          mePullxZpPanel1=0;
          mePullyZpPanel1=0;
          meNpixZpPanel1=0;
          meNxpixZpPanel1=0;
          meNypixZpPanel1=0;
          meChargeZpPanel1=0;
          meResXvsAlphaZpPanel1=0;
          meResYvsAlphaZpPanel1=0;
          meResXvsBetaZpPanel1=0;
          meResYvsBetaZpPanel1=0;

          mePullXvsAlphaZpPanel1=0;
          mePullYvsAlphaZpPanel1=0;
          mePullXvsBetaZpPanel1=0;
          mePullYvsBetaZpPanel1=0;
          mePullXvsPhiZpPanel1=0;
          mePullYvsPhiZpPanel1=0;
          mePullXvsEtaZpPanel1=0;
          mePullYvsEtaZpPanel1=0;

          meWPullXvsAlphaZpPanel1=0;
          meWPullYvsAlphaZpPanel1=0;
          meWPullXvsBetaZpPanel1=0;
          meWPullYvsBetaZpPanel1=0;

          mePosxZmPanel2=0;
          mePosyZmPanel2=0;
          meErrxZmPanel2=0;
          meErryZmPanel2=0;
          meResxZmPanel2=0;
          meResyZmPanel2=0;
          mePullxZmPanel2=0;
          mePullyZmPanel2=0;
          meNpixZmPanel2=0;
          meNxpixZmPanel2=0;
          meNypixZmPanel2=0;
          meChargeZmPanel2=0;
    
          meResXvsAlphaZmPanel2=0;
          meResYvsAlphaZmPanel2=0;
          meResXvsBetaZmPanel2=0;
          meResYvsBetaZmPanel2=0;
          mePullXvsAlphaZmPanel2=0;
          mePullYvsAlphaZmPanel2=0;
          mePullXvsBetaZmPanel2=0;
          mePullYvsBetaZmPanel2=0;
          mePullXvsPhiZmPanel2=0;
          mePullYvsPhiZmPanel2=0;
          mePullXvsEtaZmPanel2=0;
          mePullYvsEtaZmPanel2=0;

          mePosxZpPanel2=0;
          mePosyZpPanel2=0;
          meErrxZpPanel2=0;
          meErryZpPanel2=0;
          meResxZpPanel2=0;
          meResyZpPanel2=0;
          mePullxZpPanel2=0;
          mePullyZpPanel2=0;
          meNpixZpPanel2=0;
          meNxpixZpPanel2=0;
          meNypixZpPanel2=0;
          meChargeZpPanel2=0;

          meWPullXvsAlphaZmPanel2=0;
          meWPullYvsAlphaZmPanel2=0;
          meWPullXvsBetaZmPanel2=0;
          meWPullYvsBetaZmPanel2=0;

          meResXvsAlphaZpPanel2=0;
          meResYvsAlphaZpPanel2=0;
          meResXvsBetaZpPanel2=0;
          mePullXvsAlphaZpPanel2=0;
          mePullYvsAlphaZpPanel2=0;
          mePullXvsBetaZpPanel2=0;
          mePullYvsBetaZpPanel2=0;
          mePullXvsPhiZpPanel2=0;
          mePullYvsPhiZpPanel2=0;
          mePullXvsEtaZpPanel2=0;
          mePullYvsEtaZpPanel2=0;

          meWPullXvsAlphaZpPanel2=0;
          meWPullYvsAlphaZpPanel2=0;
          meWPullXvsBetaZpPanel2=0;
          meWPullYvsBetaZpPanel2=0;

          mePosxBarrel_all_hits=0;
          mePosyBarrel_all_hits=0;

          mePosxZmPanel1_all_hits=0;
          mePosyZmPanel1_all_hits=0;
          mePosxZmPanel2_all_hits=0;
          mePosyZmPanel2_all_hits=0;

          mePosxZpPanel1_all_hits=0;
          mePosyZpPanel1_all_hits=0;
          mePosxZpPanel2_all_hits=0;
          mePosyZpPanel2_all_hits=0;

          meTracksPerEvent=0;
          mePixRecHitsPerTrack=0;



  };

  };

  std::map<std::string, LadderLayersMEs*> ladderLayersMEsMap;
  std::map<std::string, LayerModulesMEs*> layerModulesMEsMap;
  std::map<std::string, DiskPlaquettesMEs*> diskPlaquettesMEsMap;
  std::map<std::string, LayerMEs*> layerMEsMap;
  TotalMEs *totalMEs = new TotalMEs();



  // variables that go in the output tree
  float rechitx; // x position of hit 
  float rechity; // y position of hit
  float rechitz; // z position of hit
  float rechiterrx; // x position error of hit (error not squared)
  float rechiterry; // y position error of hit (error not squared)

  float rechitresx; // difference between reconstructed hit x position and 'true' x position
  float rechitresy; // difference between reconstructed hit y position and 'true' y position
  float rechitpullx; // x residual divideded by error
  float rechitpully; // y residual divideded by error

  int npix; // number of pixel in the cluster
  int nxpix; // size of cluster (number of pixels) along x direction
  int nypix; // size of cluster (number of pixels) along y direction
  float charge; // total charge in cluster

  float alpha; // track angle in the xz plane of the module local coordinate system  
  float beta;  // track angle in the yz plane of the module local coordinate system  

  float phi;   // polar track angle
  float eta;   // pseudo-rapidity (function of theta, the azimuthal angle)

  int subdetId;
  int layer;  
  int ladder; 
  int mod;    
  int side;    
  int disk;   
  int blade;  
  int panel;  
  int plaq;   
  int half; // half = 1 if the barrel module is half size and 0 if it is full size (only defined for barrel) 
  int flipped; // flipped = 1 if the module is flipped and 0 if non-flipped (only defined for barrel) 

  int nsimhit; // number of simhits associated with a rechit
  int pidhit; // PID of the particle that produced the simHit associated with the recHit

  float simhitx; // true x position of hit 
  float simhity; // true y position of hit

  int evt;
  int run;

  TFile * tfile_;
  TTree * t_;
};

#endif
