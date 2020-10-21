
#include "kaldet/ILDSETKalDetector.h"

#include "kaldet/MaterialDataBase.h"

#include "kaldet/ILDParallelPlanarStripMeasLayer.h"


#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>

#include "DetInterface/IGeomSvc.h"
#include "DD4hep/Detector.h"
#include "DDRec/DetectorData.h"

#include <gear/GEAR.h>
#include "gear/BField.h"
#include <gear/ZPlanarParameters.h>
#include <gear/ZPlanarLayerLayout.h>
#include "gearimpl/Util.h"

#include "CLHEP/Units/SystemOfUnits.h"
#include "DD4hep/DD4hepUnits.h"
#include "TMath.h"

#include "math.h"
#include <sstream>

// #include "streamlog/streamlog.h"

ILDSETKalDetector::ILDSETKalDetector( const gear::GearMgr& gearMgr, IGeomSvc* geoSvc )
: TVKalDetector(300) // SJA:FIXME initial size, 300 looks reasonable for ILD, though this would be better stored as a const somewhere
{
  
  
  // streamlog_out(DEBUG4) << "ILDSETKalDetector building SET detector using GEAR " << std::endl ;
  
  MaterialDataBase::Instance().registerForService(gearMgr);
  
  TMaterial & air       = *MaterialDataBase::Instance().getMaterial("air");
  TMaterial & silicon   = *MaterialDataBase::Instance().getMaterial("silicon");
  TMaterial & carbon    = *MaterialDataBase::Instance().getMaterial("carbon");
  
  if(geoSvc){
    setupGearGeom( geoSvc ) ;
  }                                                                                                                                                                                    
  else{
    setupGearGeom( gearMgr );
  }

  if (_isStripDetector) {
    // streamlog_out(DEBUG4) << "\t\t building SET detector as STRIP Detector." << std::endl ;
  } else {
    // streamlog_out(DEBUG4) << "\t\t building SET detector as PIXEL Detector." << std::endl ;
  }

  
  //--The Ladder structure (realistic ladder)--
  int nLadders;
  
//  Bool_t active = true;
  Bool_t dummy  = false;
  
  static const double eps_layer  = 1e-6; 
  static const double eps_sensor = 1e-8; 
  
  UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ; 
  
  for (int layer=0; layer<_nLayers; ++layer) {
    
    double offset = _SETgeo[layer].offset ;
    
    if( offset!=0.0 ) {
      std::cerr << "SET can not have offsets using the current implementation: exit(1)" << std::endl ; 
      exit(1) ;
    }
    
    nLadders = _SETgeo[layer].nLadders ;
    
    const double phi0 = _SETgeo[layer].phi0 ;
    
    const double ladder_distance = _SETgeo[layer].supRMin ;
    const double ladder_thickness = _SETgeo[layer].supThickness ;
    
    const double sensitive_distance = _SETgeo[layer].senRMin ;
    const double sensitive_thickness = _SETgeo[layer].senThickness ;
    
    const double width = _SETgeo[layer].width ;
    const double length = _SETgeo[layer].length;
    
    double currPhi;
    const double dphi = _SETgeo[layer].dphi ;
    
    const double stripAngle = pow(-1,layer) *_SETgeo[layer].stripAngle;
    
    const int nsensors = _SETgeo[layer].nSensorsPerLadder;
    
    const double sensor_length = _SETgeo[layer].sensorLength;
    
    
    for (int ladder=0; ladder<nLadders; ++ladder) {
      
      
      currPhi = phi0 + (dphi * ladder);
      
      encoder.reset() ;  // reset to 0
      
      encoder[lcio::ILDCellID0::subdet] = lcio::ILDDetID::SET ;
      encoder[lcio::ILDCellID0::side] = 0 ;
      encoder[lcio::ILDCellID0::layer]  = layer ;
      encoder[lcio::ILDCellID0::module] = ladder ;
      
      double z_centre_support = 0.0;
      
      // check if the sensitive is inside or outside for the support 
      if( sensitive_distance < ladder_distance  ) {
        
        double sen_front_sorting_policy         = sensitive_distance  + (4 * ladder+0) * eps_layer ;
        double sen_back_sorting_policy          = sensitive_distance  + (4 * ladder+2) * eps_layer ;
        double sup_back_sorting_policy          = ladder_distance     + (4 * ladder+3) * eps_layer ;
        
        // air - sensitive boundary
        Add(new ILDParallelPlanarMeasLayer(air, silicon, sensitive_distance, currPhi, _bZ, sen_front_sorting_policy, width, length, offset,z_centre_support, offset, dummy,-1,"SETSenFront")) ;
        
        for (int isensor=0; isensor<nsensors; ++isensor) {

          encoder[lcio::ILDCellID0::sensor] = isensor ;          
          int CellID = encoder.lowWord() ;
          
          double measurement_plane_sorting_policy = sensitive_distance  + (4 * ladder+1) * eps_layer + eps_sensor * isensor ;
          
          double z_centre_sensor = -0.5*length + (0.5*sensor_length) + (isensor*sensor_length) ;
          
          
          if (_isStripDetector) {
            
            // measurement plane defined as the middle of the sensitive volume
            Add(new ILDParallelPlanarStripMeasLayer(silicon, silicon, sensitive_distance+sensitive_thickness*0.5, currPhi, _bZ, measurement_plane_sorting_policy, width, sensor_length, offset, z_centre_sensor, offset, stripAngle, CellID, "SETMeaslayer" )) ;
            
          } else {
            // measurement plane defined as the middle of the sensitive volume
            Add(new ILDParallelPlanarMeasLayer(silicon, silicon, sensitive_distance+sensitive_thickness*0.5, currPhi, _bZ, measurement_plane_sorting_policy, width, sensor_length, offset, z_centre_sensor, offset, true, CellID, "SETMeaslayer" )) ;
          }

          
          // streamlog_out(DEBUG0) << "ILDSETKalDetector add surface with CellID = "
          // << CellID
          // << std::endl ;
          
        }
        
       
        // sensitive - support boundary 
        Add(new ILDParallelPlanarMeasLayer(silicon, carbon, sensitive_distance+sensitive_thickness, currPhi, _bZ, sen_back_sorting_policy, width, length, offset,z_centre_support, offset, dummy,-1,"SETSenSupportIntf" )) ; 
        
        // support - air boundary
        Add(new ILDParallelPlanarMeasLayer(carbon, air, ladder_distance+ladder_thickness, currPhi, _bZ, sup_back_sorting_policy, width, length, offset,z_centre_support, offset, dummy,-1,"SETSupRear" )) ; 
      }
      else {
        
        double sup_front_sorting_policy         = ladder_distance     + (4 * ladder+0) * eps_layer ;
        double sen_front_sorting_policy         = sensitive_distance  + (4 * ladder+1) * eps_layer ;
        double sen_back_sorting_policy          = sensitive_distance  + (4 * ladder+3) * eps_layer ;
        
        // air - support boundary
        Add(new ILDParallelPlanarMeasLayer(air, carbon, ladder_distance, currPhi, _bZ, sup_front_sorting_policy, width, length, offset,z_centre_support, offset, dummy,-1,"SETSupFront")) ;
        
        // support boundary - sensitive
        Add(new ILDParallelPlanarMeasLayer(carbon, silicon, sensitive_distance, currPhi, _bZ, sen_front_sorting_policy, width, length, offset,z_centre_support, offset, dummy,-1,"SETSenSupportIntf" )) ; 
        
        
        for (int isensor=0; isensor<nsensors; ++isensor) {

          encoder[lcio::ILDCellID0::sensor] = isensor ;          
          int CellID = encoder.lowWord() ;
          
          double measurement_plane_sorting_policy = sensitive_distance  + (4 * ladder+2) * eps_layer + eps_sensor * isensor ;

          double z_centre_sensor = -0.5*length + (0.5*sensor_length) + (isensor*sensor_length) ;

          
          if (_isStripDetector) {
            
            // measurement plane defined as the middle of the sensitive volume
            Add(new ILDParallelPlanarStripMeasLayer(silicon, silicon, sensitive_distance+sensitive_thickness*0.5, currPhi, _bZ, measurement_plane_sorting_policy, width, sensor_length, offset, z_centre_sensor, offset, stripAngle, CellID, "SETMeaslayer" )) ;
            
          } else {
            // measurement plane defined as the middle of the sensitive volume
            Add(new ILDParallelPlanarMeasLayer(silicon, silicon, sensitive_distance+sensitive_thickness*0.5, currPhi, _bZ, measurement_plane_sorting_policy, width, sensor_length, offset, z_centre_sensor, offset, true, CellID, "SETMeaslayer" )) ;
          }

          
          // streamlog_out(DEBUG0) << "ILDSETKalDetector add surface with CellID = "
          // << CellID
          // << std::endl ;

          
        }

        
        // support - air boundary
        Add(new ILDParallelPlanarMeasLayer(silicon, air, sensitive_distance+sensitive_thickness, currPhi, _bZ, sen_back_sorting_policy, width, length, offset,z_centre_support, offset, dummy,-1,"SETSenRear" )) ;  
      }
      
      
    }    
    
  }
  
  SetOwner();                   
}



void ILDSETKalDetector::setupGearGeom( const gear::GearMgr& gearMgr ){
  
  const gear::ZPlanarParameters& pSETDetMain = gearMgr.getSETParameters();
  const gear::ZPlanarLayerLayout& pSETLayerLayout = pSETDetMain.getZPlanarLayerLayout();
  
  _bZ = gearMgr.getBField().at( gear::Vector3D( 0.,0.,0.)  ).z() ;
  
  _nLayers = pSETLayerLayout.getNLayers(); 
  _SETgeo.resize(_nLayers);
  
  bool n_sensors_per_ladder_present = true;
  
  try {

    std::vector<int> v = pSETDetMain.getIntVals("n_sensors_per_ladder");

  } catch (gear::UnknownParameterException& e) {

    n_sensors_per_ladder_present = false;

  }

  double strip_angle_deg = 0.0;
  
  try {
    
    strip_angle_deg = pSETDetMain.getDoubleVal("strip_angle_deg");
    _isStripDetector = true;
    
  } catch (gear::UnknownParameterException& e) {
    
    _isStripDetector = false;
    
  }
  
  
  //SJA:FIXME: for now the support is taken as the same size the sensitive
  //           if this is not done then the exposed areas of the support would leave a carbon - air boundary,
  //           which if traversed in the reverse direction to the next boundary then the track would be propagated through carbon
  //           for a significant distance 
  //std::cout << "=============SET strip angle: " << strip_angle_deg << "==============" << std::endl;
  for( int layer=0; layer < _nLayers; ++layer){
    
      
    _SETgeo[layer].nLadders = pSETLayerLayout.getNLadders(layer); 
    _SETgeo[layer].phi0 = pSETLayerLayout.getPhi0(layer); 
    _SETgeo[layer].dphi = 2*M_PI / _SETgeo[layer].nLadders; 
    _SETgeo[layer].senRMin = pSETLayerLayout.getSensitiveDistance(layer); 
    _SETgeo[layer].supRMin = pSETLayerLayout.getLadderDistance(layer); 
    _SETgeo[layer].length = pSETLayerLayout.getSensitiveLength(layer)*2.0; // note: gear for historical reasons uses the halflength 
    _SETgeo[layer].width = pSETLayerLayout.getSensitiveWidth(layer); 
    _SETgeo[layer].offset = pSETLayerLayout.getSensitiveOffset(layer); 
    _SETgeo[layer].senThickness = pSETLayerLayout.getSensitiveThickness(layer); 
    _SETgeo[layer].supThickness = pSETLayerLayout.getLadderThickness(layer); 

    if (n_sensors_per_ladder_present) {
      _SETgeo[layer].nSensorsPerLadder =   pSETDetMain.getIntVals("n_sensors_per_ladder")[layer];
    }
    else{
      _SETgeo[layer].nSensorsPerLadder = 1 ;
    }
    
    _SETgeo[layer].sensorLength = _SETgeo[layer].length / _SETgeo[layer].nSensorsPerLadder;
    
    
    if (_isStripDetector) {
      _SETgeo[layer].stripAngle = strip_angle_deg * M_PI/180 ;
    } else {
      _SETgeo[layer].stripAngle = 0.0 ;
    }
    //std::cout << _SETgeo[layer].nLadders << " " << _SETgeo[layer].phi0 << " "<< _SETgeo[layer].dphi << " " << _SETgeo[layer].senRMin << " " << _SETgeo[layer].supRMin << " "
    //          << _SETgeo[layer].length << " " << _SETgeo[layer].width << " " << _SETgeo[layer].offset << " " << _SETgeo[layer].senThickness << " " << _SETgeo[layer].supThickness << " "
    //          << _SETgeo[layer].nSensorsPerLadder << " " << _SETgeo[layer].sensorLength << std::endl;
    // streamlog_out(DEBUG0) << " layer  = " << layer << std::endl;
    // streamlog_out(DEBUG0) << " nSensorsPerLadder  = " << _SETgeo[layer].nSensorsPerLadder << std::endl;
    // streamlog_out(DEBUG0) << " sensorLength  = " << _SETgeo[layer].sensorLength << std::endl;
    // streamlog_out(DEBUG0) << " stripAngle  = " << _SETgeo[layer].stripAngle << std::endl;
    // streamlog_out(DEBUG0) << " _isStripDetector  = " << _isStripDetector << std::endl;

  }
}

void ILDSETKalDetector::setupGearGeom( IGeomSvc* geoSvc ){
  dd4hep::DetElement world = geoSvc->getDD4HepGeo();
  dd4hep::DetElement set;
  const std::map<std::string, dd4hep::DetElement>& subs = world.children();
  for(std::map<std::string, dd4hep::DetElement>::const_iterator it=subs.begin();it!=subs.end();it++){
    if(it->first!="SET") continue;
    set = it->second;
  }
  dd4hep::rec::ZPlanarData* setData = nullptr;
  try{
    setData = set.extension<dd4hep::rec::ZPlanarData>();
  }
  catch(std::runtime_error& e){
    std::cout << e.what() << " " << setData << std::endl;
    throw GaudiException(e.what(), "FATAL", StatusCode::FAILURE);
  }

  const dd4hep::Direction& field = geoSvc->lcdd()->field().magneticField(dd4hep::Position(0,0,0));
  _bZ = field.z()/dd4hep::tesla;

  std::vector<dd4hep::rec::ZPlanarData::LayerLayout>& setlayers = setData->layers;
  _nLayers = setlayers.size();
  _SETgeo.resize(_nLayers);

  double strip_angle_deg = setData->angleStrip/CLHEP::degree;
  //std::cout << "=============SET strip angle: " << strip_angle_deg << "==============" << std::endl;
  for( int layer=0; layer < _nLayers; ++layer){
    dd4hep::rec::ZPlanarData::LayerLayout& pSETLayerLayout = setlayers[layer];

    _SETgeo[layer].nLadders = pSETLayerLayout.ladderNumber;
    _SETgeo[layer].phi0 = pSETLayerLayout.phi0;
    _SETgeo[layer].dphi = 2*M_PI / _SETgeo[layer].nLadders;
    _SETgeo[layer].senRMin = pSETLayerLayout.distanceSensitive*CLHEP::cm;
    _SETgeo[layer].supRMin = pSETLayerLayout.distanceSupport*CLHEP::cm;
    _SETgeo[layer].length = pSETLayerLayout.zHalfSensitive*2.0*CLHEP::cm; // note: gear for historical reasons uses the halflength
    _SETgeo[layer].width = pSETLayerLayout.widthSensitive*CLHEP::cm;
    _SETgeo[layer].offset = pSETLayerLayout.offsetSensitive*CLHEP::cm;
    _SETgeo[layer].senThickness = pSETLayerLayout.thicknessSensitive*CLHEP::cm;
    _SETgeo[layer].supThickness = pSETLayerLayout.thicknessSupport*CLHEP::cm;
    _SETgeo[layer].nSensorsPerLadder = pSETLayerLayout.sensorsPerLadder;
    _SETgeo[layer].sensorLength = _SETgeo[layer].length / _SETgeo[layer].nSensorsPerLadder;

    if (_isStripDetector) {
      _SETgeo[layer].stripAngle = strip_angle_deg * M_PI/180 ;
    } else {
      _SETgeo[layer].stripAngle = 0.0 ;
    }
    //std::cout << _SETgeo[layer].nLadders << " " << _SETgeo[layer].phi0 << " "<< _SETgeo[layer].dphi << " " << _SETgeo[layer].senRMin << " " << _SETgeo[layer].supRMin << " "
    //        << _SETgeo[layer].length << " " << _SETgeo[layer].width << " " << _SETgeo[layer].offset << " " << _SETgeo[layer].senThickness << " " << _SETgeo[layer].supThickness << " "
    //        << _SETgeo[layer].nSensorsPerLadder << " " << _SETgeo[layer].sensorLength << " " << pSETLayerLayout.lengthSensor << std::endl;
  }
}
