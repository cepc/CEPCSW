
#include "kaldet/ILDSITCylinderKalDetector.h"
#include "kaldet/ILDCylinderMeasLayer.h"
#include "kaldet/ILDCylinderHit.h"

#include "kaldet/MaterialDataBase.h"

#include <sstream>
#include <cstdlib>

#include "gear/GEAR.h"
#include "gear/BField.h"
#include "gearimpl/Util.h"

#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>

// #include "streamlog/streamlog.h"


ILDSITCylinderKalDetector::ILDSITCylinderKalDetector( const gear::GearMgr& gearMgr ) : 
TVKalDetector(100)
{
  
  // streamlog_out(DEBUG1) << "ILDSITCylinderKalDetector building SIT Simple Cylinder Based Detector using GEAR " << std::endl ;

 ;
  
  this->setupGearGeom(gearMgr);
  
  MaterialDataBase::Instance().registerForService(gearMgr);
  TMaterial & air      = *MaterialDataBase::Instance().getMaterial("air");
  TMaterial & silicon  = *MaterialDataBase::Instance().getMaterial("silicon");
  TMaterial & carbon   = *MaterialDataBase::Instance().getMaterial("carbon");
  
  Bool_t active = true;
  Bool_t dummy  = false;
  
  std::string name = "SIT";

  UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ; 
  
  for ( unsigned int ilayer = 0 ; ilayer<_nLayers ; ++ilayer) {

    encoder[lcio::ILDCellID0::subdet] = lcio::ILDDetID::SIT ;
    encoder[lcio::ILDCellID0::layer]  = ilayer ;
    encoder[lcio::ILDCellID0::module] = 0 ;
    encoder[lcio::ILDCellID0::sensor] = 0 ;
    
    int CellID = encoder.lowWord() ;
    
    double rad   = _SITgeo[ilayer].radius - 0.5 * _SITgeo[ilayer].senThickness;
    double lhalf = _SITgeo[ilayer].half_length;
    
    const double x0 = 0.0;
    const double y0 = 0.0;
    const double z0 = 0.0;
    
    Add( new ILDCylinderMeasLayer(air, silicon , rad, lhalf, x0, y0, z0, _bZ, dummy,-1,"SITSenFront" ) );

    // streamlog_out( DEBUG0 )   << " *** adding " << name << " Front face Measurement layer at R = " << rad << " and half length = " << lhalf 
    // << " X0_in = " << air.GetRadLength() << "  X0_out = " <<  silicon.GetRadLength()    
    // << std::endl ;  
    
    rad += 0.5 * _SITgeo[ilayer].senThickness;
    
    Add( new ILDCylinderMeasLayer(silicon , silicon, rad, lhalf, x0, y0, z0, _bZ, active,CellID,"SITMeaslayer" ) );
    
    // streamlog_out( DEBUG0 )   << " *** adding " << name << " Measurement layer using CellID: " << CellID << " at R = " << rad << " and half length = " << lhalf 
    // << " X0_in = " << silicon.GetRadLength() << "  X0_out = " <<  silicon.GetRadLength()    
    // << std::endl ;  

    rad += 0.5 * _SITgeo[ilayer].senThickness;

    Add( new ILDCylinderMeasLayer(silicon , carbon, rad, lhalf, x0, y0, z0, _bZ, dummy, -1,"SITSenSupInterFace" ) );
    
    // streamlog_out( DEBUG0 )   << " *** adding " << name << " Interface between senstive and support Measurement layer at R = " << rad << " and half length = " << lhalf 
    // << " X0_in = " << silicon.GetRadLength() << "  X0_out = " <<  carbon.GetRadLength()    
    // << std::endl ;  

    rad += _SITgeo[ilayer].supThickness;
    
    Add( new ILDCylinderMeasLayer(carbon , air, rad, lhalf, x0, y0, z0, _bZ, dummy, -1,"SITSupRear" ) );
    
    // streamlog_out( DEBUG0 )   << " *** adding " << name << " Rear face Measurement layer using at R = " << rad << " and half length = " << lhalf 
    // << " X0_in = " << carbon.GetRadLength() << "  X0_out = " <<  air.GetRadLength()    
    // << std::endl ;  

    
    
  }
  
  SetOwner();
}


void ILDSITCylinderKalDetector::setupGearGeom( const gear::GearMgr& gearMgr ){
  
  
  const gear::GearParameters& pSIT = gearMgr.getGearParameters("SIT");
  
  _bZ = gearMgr.getBField().at( gear::Vector3D( 0.,0.,0.)  ).z() ;
  
  double SIT_si  =  pSIT.getDoubleVal("SITLayerThickness" )  ;
  double SIT_sp  =  pSIT.getDoubleVal("SITSupportLayerThickness" )  ;
  const std::vector<double>& SIT_r   =  pSIT.getDoubleVals("SITLayerRadius" )  ;
  const std::vector<double>& SIT_hl  =  pSIT.getDoubleVals("SITSupportLayerHalfLength" )  ;
  
  _nLayers = SIT_r.size() ; 
  _SITgeo.resize(_nLayers);
  
  
  if (_nLayers != SIT_r.size() || _nLayers != SIT_hl.size()) {

    // streamlog_out( ERROR ) << "ILDSITCylinderKalDetector miss-match between DoubleVec and nlayers exit(1) called from file " << __FILE__ << " line " << __LINE__  << std::endl ;
    exit(1);

  }

  
  for(unsigned int layer=0; layer< _nLayers; ++layer){

    _SITgeo[layer].senThickness =  SIT_si;
    _SITgeo[layer].supThickness =  SIT_sp;
    _SITgeo[layer].radius       =  SIT_r[layer];
    _SITgeo[layer].half_length  =  SIT_hl[layer];


  }
  
  
}

