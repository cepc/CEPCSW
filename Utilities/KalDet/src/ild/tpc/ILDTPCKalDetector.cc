
#include "ILDTPCKalDetector.h"
#include "ILDCylinderMeasLayer.h"
#include "ILDCylinderHit.h"

#include "TMath.h"
#include "TTUBE.h"

#include "MaterialDataBase.h"

#include <sstream>

#include "DetInterface/IGeomSvc.h"
#include "DD4hep/Detector.h"
#include "DDRec/DetectorData.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "DD4hep/DD4hepUnits.h"

#include "gear/GEAR.h"
#include "gear/BField.h"
#include "gear/TPCParameters.h"
#include "gear/PadRowLayout2D.h"
#include "gearimpl/Util.h"

#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>

// #include "streamlog/streamlog.h"


ILDTPCKalDetector::ILDTPCKalDetector( const gear::GearMgr& gearMgr, IGeomSvc* geoSvc ) : 
TVKalDetector(250) // SJA:FIXME initial size, 250 looks reasonable for ILD, though this would be better stored as a const somewhere
{
  Double_t bz;
  Int_t    nlayers;
  Double_t lhalf, rstep, rmin, rtub, outerr, inthick, outthick;
  // streamlog_out(DEBUG1) << "ILDTPCKalDetector building TPC detector using GEAR " << std::endl ;
  if(geoSvc){
    dd4hep::DetElement world = geoSvc->getDD4HepGeo();
    dd4hep::DetElement tpc;
    const std::map<std::string, dd4hep::DetElement>& subs = world.children();
    for(std::map<std::string, dd4hep::DetElement>::const_iterator it=subs.begin();it!=subs.end();it++){
      if(it->first!="TPC") continue;
      tpc = it->second;
    }
    dd4hep::rec::FixedPadSizeTPCData* tpcData = nullptr;
    try{
      tpcData = tpc.extension<dd4hep::rec::FixedPadSizeTPCData>();
    }
    catch(std::runtime_error& e){
      std::cout << e.what() << " " << tpcData << std::endl;
      throw GaudiException(e.what(), "FATAL", StatusCode::FAILURE);
    }

    const dd4hep::Direction& field = geoSvc->lcdd()->field().magneticField(dd4hep::Position(0,0,0));
    bz = field.z()/dd4hep::tesla;
    nlayers = tpcData->maxRow;
    lhalf = tpcData->driftLength*CLHEP::cm;
    rstep = tpcData->padHeight*CLHEP::cm;
    rmin = tpcData->rMinReadout*CLHEP::cm + 0.5*rstep;
    rtub = tpcData->rMin*CLHEP::cm;
    outerr = tpcData->rMax*CLHEP::cm;
    inthick = tpcData->innerWallThickness*CLHEP::cm;
    outthick = tpcData->outerWallThickness*CLHEP::cm;
    //std::cout << "TPC: " << nlayers << " " << lhalf << " " << rstep << " " << rmin << " " << rtub << " " << outerr << " " << inthick << " " << outthick << std::endl;
  }
  else{
    const gear::TPCParameters& tpcParams = gearMgr.getTPCParameters();
  
    const gear::PadRowLayout2D& pL = tpcParams.getPadLayout() ; 

    // std::cout << "ILDTPCKalDetector - got padlayout with nLayers = " <<  pL.getNRows()  <<  std::endl ;
    
    bz = gearMgr.getBField().at( gear::Vector3D( 0.,0.,0.)  ).z() ;
  
    nlayers   =  pL.getNRows() ;   // n rows
    lhalf     =  tpcParams.getMaxDriftLength() ;  // half length
    rstep     =  pL.getRowHeight(0) ;     // step length of radius
    
    // assuming that this is the radius of the first measurment layer ....
    rmin      =  tpcParams.getPlaneExtent()[0]   + rstep/2. ;   // minimum radius
    
    // std::cout << tpcParams << std::endl ;
    
    rtub      = tpcParams.getDoubleVal("tpcInnerRadius")  ; // inner r of support tube
    outerr    = tpcParams.getDoubleVal("tpcOuterRadius")  ; // outer radius of TPC
    
    inthick   =  tpcParams.getDoubleVal("tpcInnerWallThickness")  ;   // thickness of inner shell
    outthick  =  tpcParams.getDoubleVal("tpcOuterWallThickness")  ;   // thickness of outer shell
    //std::cout << "TPC: " << nlayers << " " << lhalf << " " << rstep << " " << rmin << " " << rtub << " " << outerr << " " << inthick << " " << outthick << std::endl;
  }    

  MaterialDataBase::Instance().registerForService(gearMgr, geoSvc);
  TMaterial & air          = *MaterialDataBase::Instance().getMaterial("air");
  TMaterial & tpcgas       = *MaterialDataBase::Instance().getMaterial("tpcgas");
  //  TMaterial & aluminium    = *MaterialDataBase::Instance().getMaterial("aluminium");
  TMaterial & tpcinnerfieldcage = *MaterialDataBase::Instance().getMaterial("tpcinnerfieldcage");
  TMaterial & tpcouterfieldcage = *MaterialDataBase::Instance().getMaterial("tpcouterfieldcage");
  
  Bool_t active = true;
  Bool_t dummy  = false;
  
  std::string name = "TPC";
  
  const double x0 = 0.0;
  const double y0 = 0.0;
  const double z0 = 0.0;

  
  // add inner field cage
  Add( new ILDCylinderMeasLayer(air, tpcinnerfieldcage , rtub, lhalf, x0, y0, z0, bz, dummy,-1,"TPCInnerFCInr" ) );
  // streamlog_out( DEBUG0 )   << " *** adding " << name << " Measurement layer using CellID: [ inner field cage ] at R = " << rtub
  // << " X0_in = " << air.GetRadLength() << "  X0_out = " <<  tpcinnerfieldcage.GetRadLength()    
  // << std::endl ;  
  
  Add( new ILDCylinderMeasLayer(tpcinnerfieldcage , tpcgas, rtub+inthick, lhalf, x0, y0, z0, bz, dummy,-1,"TPCInnerFCOtr" ) );
  // streamlog_out( DEBUG0 )   << " *** adding " << name << " Measurement layer using CellID: [ inner field cage ] at R = " << rtub+inthick
  // << " X0_in = " << tpcinnerfieldcage.GetRadLength() << "  X0_out = " <<  tpcgas.GetRadLength()    
  // << std::endl ;  
  
  
  // streamlog_out( DEBUG0 )   << " *** Inner Field Cage =  " << int( (inthick/(tpcinnerfieldcage.GetRadLength()*10.0) /*cm*/ )*1000) / 10.0  << "% of a radiation length " << std::endl ;  
  
  // create measurement layers
  Double_t r = rmin;
  
  UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ; 
  
  for (Int_t layer = 0; layer < nlayers; layer++) {
    
    encoder.reset() ;  // reset to 0
    
    encoder[lcio::ILDCellID0::subdet] = lcio::ILDDetID::TPC ;
    encoder[lcio::ILDCellID0::layer] = layer ;
    
    int CellID = encoder.lowWord() ;
    
    ILDCylinderMeasLayer* tpcL =  new ILDCylinderMeasLayer(tpcgas, tpcgas, r, lhalf, x0, y0, z0, bz, active, CellID, "TPCMeasLayer") ;
    
    Add( tpcL ) ;  
    
    int nth_layers(10) ;
    
    if( layer % nth_layers == 0 ){
      
      // streamlog_out( DEBUG0 )   << " *** for TPC Gas printing only every " << nth_layers << "th layer"  << std::endl ; 
      // streamlog_out( DEBUG0 )   << " *** adding " << name << " Measurement layer using CellID: [" << CellID <<  "] at R = " << r
      // << " X0_in = " << tpcgas.GetRadLength() << "  X0_out = " <<  tpcgas.GetRadLength()    
      // << std::endl ;  
    }
    
    r += rstep;

  }
  
  // add outer field cage
  Add( new ILDCylinderMeasLayer(tpcgas, tpcouterfieldcage, outerr-outthick, lhalf, x0, y0, z0, bz, dummy,-1,"TPCOuterFCInr") ) ;

  // streamlog_out( DEBUG0 )   << " *** adding " << name << " Measurement layer using CellID: [ outer field cage ] at R = " << outerr-outthick
  // << " X0_in = " << tpcgas.GetRadLength() << "  X0_out = " <<  tpcouterfieldcage.GetRadLength()    
  // << std::endl ;  
  
  Add( new ILDCylinderMeasLayer(tpcouterfieldcage, air, outerr, lhalf, x0, y0, z0, bz, dummy,-1,"TPCOuterFCOtr") ) ;

  // streamlog_out( DEBUG0 )   << " *** adding " << name << " Measurement layer using CellID: [ outer field cage ] at R = " << outerr
  // << " X0_in = " << tpcouterfieldcage.GetRadLength() << "  X0_out = " <<  air.GetRadLength()    
  // << std::endl ;  
  
  // streamlog_out( DEBUG0 )   << " *** Outer Field Cage =  " << int( (outthick/(tpcouterfieldcage.GetRadLength()*10.0) /*cm*/ )*1000) / 10.0  << "% of a radiation length " << std::endl ; 
  
  
  SetOwner();
}


