
#include "ILDVXDKalDetector.h"

#include "MaterialDataBase.h"

#include "ILDParallelPlanarMeasLayer.h"
#include "ILDCylinderMeasLayer.h"
#include "ILDDiscMeasLayer.h"

#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>

#include "DetInterface/IGeomSvc.h"
#include "DD4hep/Detector.h"
#include "DDRec/DetectorData.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "DD4hep/DD4hepUnits.h"

#include <gear/GEAR.h>
#include "gear/BField.h"
#include <gearimpl/ZPlanarParametersImpl.h>
#include <gear/VXDParameters.h>
#include <gear/VXDLayerLayout.h>
#include "gearimpl/Util.h"
#include "DetInterface/IGeomSvc.h"

#include "TMath.h"

#include "math.h"
#include <sstream>

// #include "streamlog/streamlog.h"

ILDVXDKalDetector::ILDVXDKalDetector( const gear::GearMgr& gearMgr, IGeomSvc* geoSvc )
: TVKalDetector(300) // SJA:FIXME initial size, 300 looks reasonable for ILD, though this would be better stored as a const somewhere
{
  
  // streamlog_out(DEBUG1) << "ILDVXDKalDetector building VXD detector using GEAR " << std::endl ;
  
  MaterialDataBase::Instance().registerForService(gearMgr, geoSvc);
  
  TMaterial & air       = *MaterialDataBase::Instance().getMaterial("air");
  TMaterial & silicon   = *MaterialDataBase::Instance().getMaterial("silicon");
  TMaterial & carbon    = *MaterialDataBase::Instance().getMaterial("VXDSupportMaterial");
  TMaterial & beryllium = *MaterialDataBase::Instance().getMaterial("beryllium");

  // needed for cryostat
  
  TMaterial & aluminium = *MaterialDataBase::Instance().getMaterial("aluminium");
  
  _vxd_Cryostat.exists = false;

  if(geoSvc){
    this->setupGearGeom(geoSvc) ;
  }
  else{
    this->setupGearGeom(gearMgr) ;
  }

  //--The Ladder structure (realistic ladder)--
  int nLadders;
  
  Bool_t active = true;
  Bool_t dummy  = false;
  
  static const double eps = 1e-6; 
  
  UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ; 
  
  for (int layer=0; layer<_nLayers; ++layer) {
    
    nLadders = _VXDgeo[layer].nLadders ;
    
    double phi0 = _VXDgeo[layer].phi0 ;
    
    double ladder_distance = _VXDgeo[layer].supRMin ;
    double ladder_thickness = _VXDgeo[layer].supThickness ;
    
    double sensitive_distance = _VXDgeo[layer].senRMin ;
    double sensitive_thickness = _VXDgeo[layer].senThickness ;
    
    double width = _VXDgeo[layer].width ;
    double length = _VXDgeo[layer].length;
    double offset = _VXDgeo[layer].offset;
    
    double pos_xi_nonoverlap_width = (2.0 * (( width / 2.0 ) - fabs(offset))); 
    
    double currPhi;
    double dphi = _VXDgeo[layer].dphi ;
    
    static const double z_offset = 0.0; // all VXD planes are centred at zero 
    
    for (int ladder=0; ladder<nLadders; ++ladder) {
      
      currPhi = phi0 + (dphi * ladder);
      
      encoder.reset() ;  // reset to 0
      
      encoder[lcio::ILDCellID0::subdet] = lcio::ILDDetID::VXD ;
      encoder[lcio::ILDCellID0::side] = 0 ;
      encoder[lcio::ILDCellID0::layer]  = layer ;
      encoder[lcio::ILDCellID0::module] = ladder ;
      encoder[lcio::ILDCellID0::sensor] = 0 ;
      
      int CellID = encoder.lowWord() ;
      
      // even layers have the senstive side facing the IP
      if(layer%2 == 0 ){ // overlap section of ladder0 is defined after the last ladder,
        
        
        double sen_front_sorting_policy         = sensitive_distance  + (4 * ladder+0) * eps ;
        double measurement_plane_sorting_policy = sensitive_distance  + (4 * ladder+1) * eps ;
        double sen_back_sorting_policy          = sensitive_distance  + (4 * ladder+2) * eps ;
        double sup_back_sorting_policy          = sensitive_distance  + (4 * ladder+3) * eps ;
        
        
        if(ladder==0){   // bacause overlap section of ladder0 is further outer than the last ladder.
          
          // streamlog_out(DEBUG0) << "ILDVXDKalDetector add surface with CellID = "
          // << CellID
          // << std::endl ;
                            
          // non overlapping region
          // air - sensitive boundary
          Add(new ILDParallelPlanarMeasLayer(air, silicon, sensitive_distance, currPhi, _bZ, sen_front_sorting_policy, pos_xi_nonoverlap_width, length, 0.0, z_offset, offset, dummy,-1,"VXDSenFront_non_overlap_even" )) ;
          
          // measurement plane defined as the middle of the sensitive volume  - unless "relative_position_of_measurement_surface" parameter given in GEAR
          Add(new ILDParallelPlanarMeasLayer(silicon, silicon, sensitive_distance+sensitive_thickness*_relative_position_of_measurement_surface, currPhi, _bZ, measurement_plane_sorting_policy, pos_xi_nonoverlap_width, length, 0.0, z_offset, offset, active, CellID, "VXDMeasLayer_non_overlap_even" )) ;          
          
          // sensitive - support boundary 
          Add(new ILDParallelPlanarMeasLayer(silicon, carbon, sensitive_distance+sensitive_thickness, currPhi, _bZ, sen_back_sorting_policy, pos_xi_nonoverlap_width, length, 0.0, z_offset, offset, dummy,-1,"VXDSenSuppportIntf_non_overlap_even" )) ; 
          
          // support - air boundary
          Add(new ILDParallelPlanarMeasLayer(carbon, air, ladder_distance+ladder_thickness, currPhi, _bZ, sup_back_sorting_policy, pos_xi_nonoverlap_width, length, 0.0, z_offset, offset, dummy,-1,"VXDSupRear_non_overlap_even" )) ;           
          
          
          // overlapping region
          double overlap_region_width  = width - pos_xi_nonoverlap_width ;
          double overlap_region_offset = -(overlap_region_width/2.0) - (pos_xi_nonoverlap_width)/2.0 ;

          
          // overlap sorting policy uses nLadders as the overlapping "ladder" is the order i.e. there will now be nLadders+1 
          double overlap_front_sorting_policy                = sensitive_distance + (4* nLadders+0) * eps;
          double overlap_measurement_plane_sorting_policy    = sensitive_distance + (4* nLadders+1) * eps;
          double overlap_back_sorting_policy                 = sensitive_distance + (4* nLadders+2) * eps;
          double overlap_sup_back_sorting_policy             = sensitive_distance + (4* nLadders+3) * eps;
          
          // streamlog_out(DEBUG0) << "ILDVXDKalDetector add surface with CellID = "
          // << CellID
          // << std::endl ;
          
          // air - sensitive boundary
          Add(new ILDParallelPlanarMeasLayer(air, silicon, sensitive_distance, currPhi, _bZ, overlap_front_sorting_policy, overlap_region_width, length, overlap_region_offset, z_offset, offset, dummy,-1,"VXDSenFront_overlap_even")) ;
          
          // measurement plane defined as the middle of the sensitive volume  - unless "relative_position_of_measurement_surface" parameter given in GEAR
          Add(new ILDParallelPlanarMeasLayer(silicon, silicon, sensitive_distance+sensitive_thickness*_relative_position_of_measurement_surface, currPhi, _bZ, overlap_measurement_plane_sorting_policy, overlap_region_width, length, overlap_region_offset, z_offset, offset, active, CellID, "VXDMeasLayer_overlap_even" )) ;
          
          
          // sensitive - support boundary 
          Add(new ILDParallelPlanarMeasLayer(silicon, carbon, sensitive_distance+sensitive_thickness, currPhi, _bZ, overlap_back_sorting_policy, overlap_region_width, length, overlap_region_offset, z_offset, offset, dummy,-1,"VXDSenSuppportIntf_overlap_even")) ; 
          
          // support - air boundary
          Add(new ILDParallelPlanarMeasLayer(carbon, air, ladder_distance+ladder_thickness, currPhi, _bZ, overlap_sup_back_sorting_policy, overlap_region_width, length, overlap_region_offset, z_offset, offset, dummy,-1,"VXDSupRear_overlap_even")) ; 
          
        }
        else{
          
          // streamlog_out(DEBUG0) << "ILDVXDKalDetector (ILDParallelPlanarMeasLayer) add surface with CellID = "
          // << CellID
          // << std::endl ;                                        
          
          
          // air - sensitive boundary
          Add(new ILDParallelPlanarMeasLayer(air, silicon, sensitive_distance, currPhi, _bZ, sen_front_sorting_policy, width, length, offset, z_offset, offset, dummy,-1, "VXDSenFront_even")) ;
                    

          // measurement plane defined as the middle of the sensitive volume  - unless "relative_position_of_measurement_surface" parameter given in GEAR - even layers face outwards ! 
          Add(new ILDParallelPlanarMeasLayer(silicon, silicon, sensitive_distance+sensitive_thickness*( 1.-_relative_position_of_measurement_surface ), currPhi, _bZ, measurement_plane_sorting_policy, width, length, offset, z_offset, offset, active, CellID, "VXDMeaslayer_even" )) ;
          
          // sensitive - support boundary 
          Add(new ILDParallelPlanarMeasLayer(silicon, carbon, sensitive_distance+sensitive_thickness, currPhi, _bZ, sen_back_sorting_policy, width, length, offset, z_offset, offset, dummy,-1,"VXDSenSuppportIntf_even" )) ; 
          
          // support - air boundary
          Add(new ILDParallelPlanarMeasLayer(carbon, air, ladder_distance+ladder_thickness, currPhi, _bZ, sup_back_sorting_policy, width, length, offset, z_offset, offset, dummy,-1,"VXDSupRear_even" )) ; 
          
        }        
      }
      else{ // counting from 0, odd numbered layers are placed with the support closer to the IP than the sensitive
                                
        
        
        double sup_forward_sorting_policy        = ladder_distance + (4 * ladder+0) * eps ;
        double sup_back_sorting_policy           = ladder_distance + (4 * ladder+1) * eps ;
        double measurement_plane_sorting_policy  = ladder_distance + (4 * ladder+2) * eps ;
        double sen_back_sorting_policy           = ladder_distance + (4 * ladder+3) * eps ;
        
        // streamlog_out(DEBUG0) << "ILDVXDKalDetector (ILDPlanarMeasLayer) add surface with CellID = "
        // << CellID
        // << std::endl ;
                
        // air - support boundary
        Add(new ILDParallelPlanarMeasLayer(air, carbon, ladder_distance, currPhi, _bZ, sup_forward_sorting_policy, width, length, offset, z_offset, offset, dummy,-1,"VXDSupFront_odd" )) ; 
        
        // support - sensitive boundary 
        Add(new ILDParallelPlanarMeasLayer(carbon, silicon, (ladder_distance+ladder_thickness), currPhi, _bZ, sup_back_sorting_policy, width, length, offset, z_offset, offset, dummy,-1,"VXDSenSuppportIntf_odd")) ; 
        
        // measurement plane defined as the middle of the sensitive volume
        Add(new ILDParallelPlanarMeasLayer(silicon, silicon, (sensitive_distance+sensitive_thickness*0.5), currPhi, _bZ, measurement_plane_sorting_policy, width, length, offset, z_offset, offset, active, CellID, "VXDMeaslayer_odd")) ; 
        
        // sensitive air - sensitive boundary
        Add(new ILDParallelPlanarMeasLayer(silicon, air, (sensitive_distance+sensitive_thickness), currPhi, _bZ, sen_back_sorting_policy, width, length, offset, z_offset, offset, dummy,-1,"VXDSenRear_odd")) ;
        
        
      }
    }
  }
  
  if (_vxd_Cryostat.exists) {
    // build Cryostat according to mokka driver vxd04.cc
    
    // beryllium shell
    
    double rtub  = _vxd_Cryostat.shellInnerR;
    double halfz = _vxd_Cryostat.shelllHalfZ;

    const double x0 = 0.0;
    const double y0 = 0.0;
    const double z0 = 0.0;

    
    // beryllium cylinder inner wall
    Add( new ILDCylinderMeasLayer(air, beryllium , rtub, halfz, x0, y0, z0, _bZ, dummy,-1,"VXDShellInnerWall" ) );
    
    // streamlog_out( DEBUG0 )   << " *** adding " << "VXDShellInnerWall" << " Measurement layer using CellID: [ VXDShellInnerWall ] at R = " << rtub
    // << " X0_in = " << air.GetRadLength() << "  X0_out = " <<  beryllium.GetRadLength()
    // << std::endl;
    

    rtub  += _vxd_Cryostat.shellThickness;
    
    // beryllium cylinder outer wall
    Add( new ILDCylinderMeasLayer(beryllium, air , rtub, halfz, x0, y0, z0, _bZ, dummy,-1,"VXDShellOuterWall" ) );
    
    // streamlog_out( DEBUG0 )   << " *** adding " << "VXDShellOuterWall" << " Measurement layer using CellID: [ VXDShellOuterWall ] at R = " << rtub
    // << " X0_in = " << beryllium.GetRadLength() << "  X0_out = " <<  air.GetRadLength()
    // << std::endl;

    
    rtub  = _vxd_Cryostat.alRadius;
    halfz = _vxd_Cryostat.alHalfZ;

    
    // aluminum cylinder inner wall
    Add( new ILDCylinderMeasLayer(air, aluminium , _vxd_Cryostat.alRadius, halfz, x0, y0, z0, _bZ, dummy,-1,"VXDCryoAlInnerWall" ) );


    // streamlog_out( DEBUG0 )   << " *** adding " << "VXDCryoAlInnerWall" << " Measurement layer using CellID: [ VXDCryoAlInnerWall ] at R = " << rtub
    // << " X0_in = " << air.GetRadLength() << "  X0_out = " <<  aluminium.GetRadLength()
    // << std::endl;
    
    rtub  += 1.1 * _vxd_Cryostat.alThickness; // SJA:FIXME: increase the thickness as we don't have the information on the foam in the GEAR file.

    // aluminum cylinder outer wall
    Add( new ILDCylinderMeasLayer(aluminium, air , rtub, halfz, x0, y0, z0, _bZ, dummy,-1,"VXDCryoAlOuterWall" ) );

    
    // streamlog_out( DEBUG0 )   << " *** adding " << "VXDCryoAlOuterWall" << " Measurement layer using CellID: [ VXDCryoAlOuterWall ] at R = " << rtub
    // << " X0_in = " << aluminium.GetRadLength() << "  X0_out = " <<  air.GetRadLength()
    // << std::endl;
    

    // aluminum endcaps
   
    
    const double z_front_face = _vxd_Cryostat.alZEndCap;
    const double z_rear_face = z_front_face + _vxd_Cryostat.alThickness;
 
    double rOuter = _vxd_Cryostat.alRadius - 0.1 ; // make sure we don't collide with the aluminium cryostat cylinder
    double rInner = _vxd_Cryostat.alInnerR ;

    
    double eps_face = 1.0e-06 ; // offset for forwards and backwards
    double eps_side = 1.0e-08 ; // offset for forwards and backwards
    
    double sort_policy = rOuter ;

    
    TVector3 xc_fwd(0.0, 0.0, z_front_face) ;
    TVector3 normal_fwd(xc_fwd) ;
    normal_fwd.SetMag(1.0) ;
    
    // streamlog_out(DEBUG) << "VXDCryoAleEndCap created disk at " << xc_fwd.z() << " sort_policy = " << sort_policy << std::endl;
    
    Add(new ILDDiscMeasLayer( air, aluminium, xc_fwd, normal_fwd, _bZ, sort_policy,
                             rInner, rOuter, dummy,-1, "VXDCryoAleEndCapDiscFrontFwd" ) );
    
    sort_policy += eps_face;
    xc_fwd.SetZ(z_rear_face);
    
    Add(new ILDDiscMeasLayer( aluminium, air, xc_fwd, normal_fwd, _bZ, sort_policy,
                             rInner, rOuter, dummy,-1, "VXDCryoAleEndCapDiscBackFwd" ) );

    
    TVector3 xc_bwd(0.0, 0.0, -z_front_face) ;
    TVector3 normal_bwd(xc_bwd) ;
    normal_bwd.SetMag(1.0) ;
    
    // offset needed for rear disks
    sort_policy = rOuter + eps_side ;
    
    // streamlog_out(DEBUG) << "VXDCryoAleEndCap created disk at " <<  xc_bwd.z() << " sort_policy = " << sort_policy << std::endl;

    Add(new ILDDiscMeasLayer( air, aluminium, xc_bwd, normal_bwd, _bZ, sort_policy,
                             rInner, rOuter, dummy,-1, "VXDCryoAleEndCapDiscFrontBwd" ) );
    
    sort_policy += eps_face;
    xc_fwd.SetZ(-z_rear_face);
    
    Add(new ILDDiscMeasLayer( aluminium, air, xc_bwd, normal_bwd, _bZ, sort_policy,
                             rInner, rOuter, dummy,-1, "VXDCryoAleEndCapDiscFrontBwd" ) );

    
    
  }
  
  SetOwner();                   
}

void ILDVXDKalDetector::setupGearGeom( const gear::GearMgr& gearMgr ){
  const gear::VXDParameters& pVXDDetMain = gearMgr.getVXDParameters();
  const gear::VXDLayerLayout& pVXDLayerLayout = pVXDDetMain.getVXDLayerLayout();

  _bZ = gearMgr.getBField().at( gear::Vector3D( 0.,0.,0.)  ).z() ;

  _nLayers = pVXDLayerLayout.getNLayers();
  _VXDgeo.resize(_nLayers);

  for( int layer=0; layer < _nLayers; ++layer){
    _VXDgeo[layer].nLadders = pVXDLayerLayout.getNLadders(layer);
    _VXDgeo[layer].phi0 = pVXDLayerLayout.getPhi0(layer);
    _VXDgeo[layer].dphi = 2*M_PI / _VXDgeo[layer].nLadders;
    _VXDgeo[layer].senRMin = pVXDLayerLayout.getSensitiveDistance(layer);
    _VXDgeo[layer].supRMin = pVXDLayerLayout.getLadderDistance(layer);
    _VXDgeo[layer].length = pVXDLayerLayout.getSensitiveLength(layer) * 2.0 ; // note: gear for historical reasons uses the halflength
    _VXDgeo[layer].width = pVXDLayerLayout.getSensitiveWidth(layer);
    _VXDgeo[layer].offset = pVXDLayerLayout.getSensitiveOffset(layer);
    _VXDgeo[layer].senThickness = pVXDLayerLayout.getSensitiveThickness(layer);
    _VXDgeo[layer].supThickness = pVXDLayerLayout.getLadderThickness(layer);
    //std::cout << layer << ": " << _VXDgeo[layer].nLadders << " " << _VXDgeo[layer].phi0 << " " << _VXDgeo[layer].dphi << " " << _VXDgeo[layer].senRMin
    //          << " " << _VXDgeo[layer].supRMin << " " << _VXDgeo[layer].length << " " << _VXDgeo[layer].width << " " << _VXDgeo[layer].offset
    //          << " " << _VXDgeo[layer].senThickness << " " << _VXDgeo[layer].supThickness << std::endl;
  }
  
  _relative_position_of_measurement_surface = 0.5 ;

  try {
    _relative_position_of_measurement_surface =  pVXDDetMain.getDoubleVal( "relative_position_of_measurement_surface"  );
  }
  catch (gear::UnknownParameterException& e) {}
  try {
    const gear::GearParameters& pVXDInfra = gearMgr.getGearParameters("VXDInfra");                                                                                                  
    _vxd_Cryostat.alRadius    = pVXDInfra.getDoubleVal( "CryostatAlRadius"  );
    _vxd_Cryostat.alThickness = pVXDInfra.getDoubleVal( "CryostatAlThickness"  );
    _vxd_Cryostat.alInnerR    = pVXDInfra.getDoubleVal( "CryostatAlInnerR"  );
    _vxd_Cryostat.alZEndCap   = pVXDInfra.getDoubleVal( "CryostatAlZEndCap"  );
    _vxd_Cryostat.alHalfZ     = pVXDInfra.getDoubleVal( "CryostatAlHalfZ"  );

    _vxd_Cryostat.shellInnerR    = pVXDDetMain.getShellInnerRadius();
    _vxd_Cryostat.shellThickness = pVXDDetMain.getShellOuterRadius() - _vxd_Cryostat.shellInnerR;
    _vxd_Cryostat.shelllHalfZ    = pVXDDetMain.getShellHalfLength();

    _vxd_Cryostat.exists = true;
    //std::cout << "VXDInfra: " << _vxd_Cryostat.alRadius << " " << _vxd_Cryostat.alThickness << " " << _vxd_Cryostat.alInnerR << " " << _vxd_Cryostat.alZEndCap << " " 
    //	      << _vxd_Cryostat.alHalfZ << " " << _vxd_Cryostat.shellInnerR << " " << _vxd_Cryostat.shellThickness << " " << _vxd_Cryostat.shelllHalfZ << std::endl;
  }    
  catch (gear::UnknownParameterException& e) {
    std::cout << e.what() << std::endl ;
    _vxd_Cryostat.exists = false;

  }
}


void ILDVXDKalDetector::setupGearGeom( IGeomSvc* geoSvc){
  /*
  dd4hep::DetElement world = geoSvc->getDD4HepGeo();
  dd4hep::DetElement vxd;
  const std::map<std::string, dd4hep::DetElement>& subs = world.children();
  for(std::map<std::string, dd4hep::DetElement>::const_iterator it=subs.begin();it!=subs.end();it++){
    if(it->first!="VXD") continue;
    vxd = it->second;
  }
  dd4hep::rec::ZPlanarData* vxdData = nullptr;
  try{
    vxdData = vxd.extension<dd4hep::rec::ZPlanarData>();
  }
  catch(std::runtime_error& e){
    std::cout << e.what() << " " << vxdData << std::endl;
    throw GaudiException(e.what(), "FATAL", StatusCode::FAILURE);
  }
  */
  const dd4hep::Direction& field = geoSvc->lcdd()->field().magneticField(dd4hep::Position(0,0,0));
  _bZ = field.z()/dd4hep::tesla;
  
  const gear::ZPlanarParametersImpl* pVXDDetMain = geoSvc->getVXDParameters();
  const gear::VXDLayerLayout& pVXDLayerLayout = pVXDDetMain->getVXDLayerLayout();
  
  _nLayers = pVXDLayerLayout.getNLayers(); 
  _VXDgeo.resize(_nLayers);
  
  //SJA:FIXME: for now the support is taken as the same size the sensitive
  //           if this is not done then the exposed areas of the support would leave a carbon - air boundary,
  //           which if traversed in the reverse direction to the next boundary then the track will be propagated through carbon
  //           for a significant distance 
  
  for( int layer=0; layer < _nLayers; ++layer){
    _VXDgeo[layer].nLadders = pVXDLayerLayout.getNLadders(layer); 
    _VXDgeo[layer].phi0 = pVXDLayerLayout.getPhi0(layer); 
    _VXDgeo[layer].dphi = 2*M_PI / _VXDgeo[layer].nLadders; 
    _VXDgeo[layer].senRMin = pVXDLayerLayout.getSensitiveDistance(layer); 
    _VXDgeo[layer].supRMin = pVXDLayerLayout.getLadderDistance(layer); 
    _VXDgeo[layer].length = pVXDLayerLayout.getSensitiveLength(layer) * 2.0 ; // note: gear for historical reasons uses the halflength 
    _VXDgeo[layer].width = pVXDLayerLayout.getSensitiveWidth(layer); 
    _VXDgeo[layer].offset = pVXDLayerLayout.getSensitiveOffset(layer); 
    _VXDgeo[layer].senThickness = pVXDLayerLayout.getSensitiveThickness(layer); 
    _VXDgeo[layer].supThickness = pVXDLayerLayout.getLadderThickness(layer); 
    //std::cout << layer << ": " << _VXDgeo[layer].nLadders << " " << _VXDgeo[layer].phi0 << " " << _VXDgeo[layer].dphi << " " << _VXDgeo[layer].senRMin 
    //	      << " " << _VXDgeo[layer].supRMin << " " << _VXDgeo[layer].length << " " << _VXDgeo[layer].width << " " << _VXDgeo[layer].offset
    //	      << " " << _VXDgeo[layer].senThickness << " " << _VXDgeo[layer].supThickness << std::endl; 
  }
  // by default, we put the measurement surface in the middle of the sensitive
  // layer, this can optionally be changed, e.g. in the case of the FPCCD where the 
  // epitaxial layer is 15 mu thick (in a 50 mu wafer)
  _relative_position_of_measurement_surface = 0.5 ;
  // Cryostat
  try {
    _vxd_Cryostat.alRadius    = geoSvc->getDetParameter("VXDInfra","CryostatAlRadius");
    _vxd_Cryostat.alThickness = geoSvc->getDetParameter("VXDInfra","CryostatAlThickness");
    _vxd_Cryostat.alInnerR    = geoSvc->getDetParameter("VXDInfra","CryostatAlInnerR");
    _vxd_Cryostat.alZEndCap   = geoSvc->getDetParameter("VXDInfra","CryostatAlZEndCap");
    _vxd_Cryostat.alHalfZ     = geoSvc->getDetParameter("VXDInfra","CryostatAlHalfZ");

    _vxd_Cryostat.shellInnerR    = pVXDDetMain->getShellInnerRadius();
    _vxd_Cryostat.shellThickness = pVXDDetMain->getShellOuterRadius() - _vxd_Cryostat.shellInnerR;    
    _vxd_Cryostat.shelllHalfZ    = pVXDDetMain->getShellHalfLength();
    
    _vxd_Cryostat.exists = true;
    //std::cout << "VXDInfra: " << _vxd_Cryostat.alRadius << " " << _vxd_Cryostat.alThickness << " " << _vxd_Cryostat.alInnerR << " " << _vxd_Cryostat.alZEndCap << " "
    //          << _vxd_Cryostat.alHalfZ << " " << _vxd_Cryostat.shellInnerR << " " << _vxd_Cryostat.shellThickness << " " << _vxd_Cryostat.shelllHalfZ << std::endl;
  }
  catch (std::runtime_error& e) {
    std::cout << e.what() << std::endl ;
    _vxd_Cryostat.exists = false;
  
  }
  
}
