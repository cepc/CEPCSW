
#include "ILDFTDDiscBasedKalDetector.h"

#include "kaldet/MaterialDataBase.h"

#include <sstream>

#include "gear/GEAR.h"
#include "gear/BField.h"
#include "gearimpl/Util.h"
#include "gear/FTDLayerLayout.h"

#include "kaldet/ILDDiscMeasLayer.h"

#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>

// #include "streamlog/streamlog.h"

#include "TVector3.h"

ILDFTDDiscBasedKalDetector::ILDFTDDiscBasedKalDetector( const gear::GearMgr& gearMgr ) : 
TVKalDetector(30), _nDisks(0) // SJA:FIXME initial size, 300 looks reasonable for ILD, though this would be better stored as a const somewhere
{
  
  // streamlog_out(DEBUG1) << "ILDFTDDiscBasedKalDetector building Simple Disc Based FTD detector using GEAR " << std::endl ;
  
  MaterialDataBase::Instance().registerForService(gearMgr);
  setupGearGeom( gearMgr ) ; 
  
  TMaterial & air       = *MaterialDataBase::Instance().getMaterial("air") ;
  TMaterial & silicon   = *MaterialDataBase::Instance().getMaterial("silicon") ;
  TMaterial & carbon    = *MaterialDataBase::Instance().getMaterial("carbon") ;
  
  
  Bool_t active = true ;
  Bool_t dummy  = false ;
  
  std::string name = "FTD" ;
  
  double eps1 = 1.0e-04 ; // disk  
  double eps3 = 1.0e-06 ; // layer in disk 
  double eps4 = 1.0e-08 ; // forward or backwards

  
  UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ; 
  
  for (int idisk = 0; idisk < _nDisks; ++idisk) {
    
    
    double rOuter = _FTDgeo[idisk].rOuter ;
    double rInner = _FTDgeo[idisk].rInner ;
    double senThickness = _FTDgeo[idisk].senThickness ;
    double supThickness = _FTDgeo[idisk].supThickness ;
    double zPos = _FTDgeo[idisk].zPos;
    
    encoder[lcio::ILDCellID0::subdet] = lcio::ILDDetID::FTD ;
    encoder[lcio::ILDCellID0::side] = 1 ;
    encoder[lcio::ILDCellID0::layer]  = idisk ;
    encoder[lcio::ILDCellID0::module] = 0 ;
    encoder[lcio::ILDCellID0::sensor] = 0 ;
    
    int CellID_FWD = encoder.lowWord() ;
    
    encoder[lcio::ILDCellID0::side] = -1 ;
    
    int CellID_BWD = encoder.lowWord() ;
    
    // note the z position given in gear is actually the mid point (z) of the sensitive i.e. the z of the measurement plane
    TVector3 sen_front_face_centre_fwd( 0.0, 0.0, zPos - senThickness*0.5); // for +z  
    
    TVector3 measurement_plane_centre_fwd( sen_front_face_centre_fwd.X(), 
                                          sen_front_face_centre_fwd.Y(), 
                                          sen_front_face_centre_fwd.Z() + senThickness*0.5 ); 
    
    TVector3 sen_rear_face_centre_fwd( sen_front_face_centre_fwd.X(), 
                                      sen_front_face_centre_fwd.Y(), 
                                      sen_front_face_centre_fwd.Z() + senThickness ); 
    
    TVector3 sup_rear_face_centre_fwd( sen_rear_face_centre_fwd.X(), 
                                      sen_rear_face_centre_fwd.Y(), 
                                      sen_rear_face_centre_fwd.Z() + supThickness ); 
    
    TVector3 normal_fwd(sen_front_face_centre_fwd) ;    
    normal_fwd.SetMag(1.0) ;    
    
    
    double sort_policy = rOuter + eps1 * idisk ;
    
    
    sort_policy += eps3 * 1;
    Add(new ILDDiscMeasLayer( air, silicon, sen_front_face_centre_fwd, normal_fwd, _bZ, sort_policy, rInner, rOuter, dummy ) );
    
    sort_policy += eps3 * 2;
    Add(new ILDDiscMeasLayer( silicon, silicon, measurement_plane_centre_fwd, normal_fwd, _bZ, sort_policy, rInner, rOuter, active, CellID_FWD ) );

    sort_policy += eps3 * 3;
    Add(new ILDDiscMeasLayer( silicon, carbon, sen_rear_face_centre_fwd, normal_fwd, _bZ, sort_policy, rInner, rOuter, dummy ) );

    sort_policy += eps3 * 4;
    Add(new ILDDiscMeasLayer( carbon, air, sup_rear_face_centre_fwd, normal_fwd, _bZ, sort_policy, rInner, rOuter, dummy ) );
    
    
    // note the z position given in gear is actually the mid point (z) of the sensitive i.e. the z of the measurement plane
    TVector3 sen_front_face_centre_bwd( 0.0, 0.0, -zPos + senThickness*0.5 );         // for -z  
    
    TVector3 measurement_plane_centre_bwd( sen_front_face_centre_bwd.X(), 
                                          sen_front_face_centre_bwd.Y(), 
                                          sen_front_face_centre_bwd.Z() - senThickness*0.5 ); 
    
    TVector3 sen_rear_face_centre_bwd( sen_front_face_centre_bwd.X(), 
                                      sen_front_face_centre_bwd.Y(), 
                                      sen_front_face_centre_bwd.Z() - senThickness ); 
    
    TVector3 sup_rear_face_centre_bwd( sen_rear_face_centre_bwd.X(), 
                                      sen_rear_face_centre_bwd.Y(), 
                                      sen_rear_face_centre_bwd.Z() - supThickness ); 
    
    TVector3 normal_bwd(sen_front_face_centre_bwd) ;
    normal_bwd.SetMag(1.0) ;
    
    
    sort_policy += eps4 ; // for backward 
    
    sort_policy += eps3 * 1;
    Add(new ILDDiscMeasLayer( air, silicon, sen_front_face_centre_bwd, normal_bwd, _bZ, sort_policy, rInner, rOuter, dummy ) );
    
    sort_policy += eps3 * 2;
    Add(new ILDDiscMeasLayer( silicon, silicon, measurement_plane_centre_bwd, normal_bwd, _bZ, sort_policy, rInner, rOuter, active, CellID_BWD ) );
    
    sort_policy += eps3 * 3;
    Add(new ILDDiscMeasLayer( silicon, carbon, sen_rear_face_centre_bwd, normal_bwd, _bZ, sort_policy, rInner, rOuter, dummy ) );
    
    sort_policy += eps3 * 4;
    Add(new ILDDiscMeasLayer( carbon, air, sup_rear_face_centre_bwd, normal_bwd, _bZ, sort_policy, rInner, rOuter, dummy ) );
    
    
  }
  
  SetOwner();
  
}



void ILDFTDDiscBasedKalDetector::setupGearGeom( const gear::GearMgr& gearMgr ){
  
  
  const gear::GearParameters& pFTD = gearMgr.getGearParameters("FTD");
  
  _bZ = gearMgr.getBField().at( gear::Vector3D( 0.,0.,0.)  ).z() ;
  
  const std::vector<double>& FTD_si  =  pFTD.getDoubleVals("FTDDiskSiThickness" )  ;
  const std::vector<double>& FTD_sp  =  pFTD.getDoubleVals("FTDDiskSupportThickness" )  ;
  const std::vector<double>& FTD_ri  =  pFTD.getDoubleVals("FTDInnerRadius" )  ;
  const std::vector<double>& FTD_ro  =  pFTD.getDoubleVals("FTDOuterRadius" )  ;
  const std::vector<double>& FTD_z   =  pFTD.getDoubleVals("FTDZCoordinate" )  ;
  
  _nDisks = FTD_si.size() ; 
  _FTDgeo.resize(_nDisks);
  
  for(int disk=0; disk< _nDisks; ++disk){
    
    _FTDgeo[disk].rInner = FTD_ri[disk];
    _FTDgeo[disk].rOuter = FTD_ro[disk];
    _FTDgeo[disk].senThickness =  FTD_si[disk];
    _FTDgeo[disk].supThickness =  FTD_sp[disk];
    _FTDgeo[disk].zPos = FTD_z[disk];
    
    
  }
  
  
}
