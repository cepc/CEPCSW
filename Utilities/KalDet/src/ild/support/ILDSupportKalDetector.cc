
#include "ILDSupportKalDetector.h"
#include "ILDCylinderMeasLayer.h"
#include "ILDConeMeasLayer.h"
#include "ILDPolygonBarrelMeasLayer.h"
#include "ILDDiscMeasLayer.h"

#include "TMath.h"
#include "TTUBE.h"

#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>

#include "MaterialDataBase.h"

#include <sstream>
#include <cmath>

#include "DD4hep/Detector.h"
#include "DDRec/DetectorData.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "DD4hep/DD4hepUnits.h"

#include "gear/GEAR.h"
#include "gear/BField.h"
#include "gearimpl/Util.h"
#include "gear/CalorimeterParameters.h"
// #include "streamlog/streamlog.h"
#include "DetInterface/IGeomSvc.h"

ILDSupportKalDetector::ILDSupportKalDetector( const gear::GearMgr& gearMgr, IGeomSvc* geoSvc ) : 
TVKalDetector(10) 
{
  Double_t bz;
  std::vector<double> z, rInner, rOuter;
  // streamlog_out(DEBUG1) << "ILDSupportKalDetector building beampipe using GEAR " << std::endl ;
  if(geoSvc){
    const dd4hep::rec::ConicalSupportData* pBeamPipeData = geoSvc->getBeamPipeData();
    const std::vector<dd4hep::rec::ConicalSupportData::Section>& sections = pBeamPipeData->sections;
    const dd4hep::Direction& field = geoSvc->lcdd()->field().magneticField(dd4hep::Position(0,0,0));
    bz = field.z()/dd4hep::tesla;
    for(int i=0;i<sections.size();i++){
      z.push_back(sections[i].zPos*CLHEP::cm  );
      rInner.push_back(sections[i].rInner*CLHEP::cm );
      rOuter.push_back(sections[i].rOuter*CLHEP::cm );
    }
  }
  else{
    const gear::GearParameters& pBeamPipe = gearMgr.getGearParameters("BeamPipe");
    bz = gearMgr.getBField().at( gear::Vector3D( 0.,0.,0.)  ).z() ;
    
    // the beampipe is supposed to be a chain of cut cones (cut, that means the spike is cut off, also called a cone frustum).
    // as they are connected, RStart[i] == REnd[i-1]. With this all we need are the z values and the radii at the place.
    z = pBeamPipe.getDoubleVals("Z");
    rInner = pBeamPipe.getDoubleVals("RInner"); //inner radius of the cone
    rOuter = pBeamPipe.getDoubleVals("ROuter"); //outer radius of the cone
  }
  //for(int i=0;i<z.size();i++){
  //std::cout << z[i] << " " << rInner[i] << " " << rOuter[i] << std::endl;
  //}
  MaterialDataBase::Instance().registerForService(gearMgr, geoSvc);
  TMaterial & beam      = *MaterialDataBase::Instance().getMaterial("beam");
  TMaterial & air       = *MaterialDataBase::Instance().getMaterial("air");
  //  TMaterial & aluminium = *MaterialDataBase::Instance().getMaterial("aluminium");
  TMaterial & beryllium = *MaterialDataBase::Instance().getMaterial("beryllium");
  
  Bool_t dummy  = false;
    
  double eps_side    = 1.0e-8;
  double eps_section = 1.0e-6;
  
  // first check the sizes of the beam-pipe values r_in r_out and z

  if ( ( z.size() != rInner.size() ) || ( z.size() != rOuter.size() ) ) {
    // streamlog_out(ERROR) << "ILDSupportKalDetector::ILDSupportKalDetector miss-match in number of double values for beam-pipe: Z = " << z.size() << " RInner = " << rInner.size() << " ROuter = " << rOuter.size()  << " exit(1) called from "<< __FILE__ << "   line " << __LINE__ << std::endl;
    exit(1);
  }
  
  double epsilon = 1.e-6;
  
  double rInnerStart_last = -DBL_MAX;
  double rInnerEnd_last   = -DBL_MAX;
  double rOuterStart_last = -DBL_MAX;
  double rOuterEnd_last   = -DBL_MAX;
  
  // add beam pipe cones
  for( unsigned i=0; i<z.size()-1; i++){
    
    double zStart = z[i];
    double zEnd = z[i+1];
    double rInnerStart = rInner[i];
    double rInnerEnd = rInner[i+1];
    double rOuterStart = rOuter[i];
    double rOuterEnd = rOuter[i+1];
    
    
    //SJA:FIXME: HERE WE NEED TO MAKE SURE THAT THE BEAM PIPE IS NEVER DECREASING IN RADII OTHERWISE THE SORTING GOES WRONG
    if ( rInnerStart < rInnerStart_last ) rInnerStart = rInnerStart_last; rInnerStart_last = rInnerStart;
    if ( rInnerEnd   < rInnerEnd_last   ) rInnerEnd   = rInnerEnd_last  ; rInnerEnd_last   = rInnerEnd;
    if ( rOuterStart < rOuterStart_last ) rOuterStart = rOuterStart_last; rOuterStart_last = rOuterStart;
    if ( rOuterEnd   < rOuterEnd_last   ) rOuterEnd   = rOuterEnd_last;   rOuterEnd_last   = rOuterEnd;
    
    
    std::stringstream sname;

    const double x0 = 0.0;
    const double y0 = 0.0;
    const double z0 = 0.0;

       
    // central part has to be treated separately as it should not be composed of two halves but a single section.
    if (i==0) {
      // check that z values start a zero
      if (fabs(zStart) > epsilon) {
        // streamlog_out(ERROR) << "ILDSupportKalDetector::ILDSupportKalDetector first z value for beam-pipe must equal zero: Z = " << zStart << " exit(1) called from "<< __FILE__ << "   line " << __LINE__ << std::endl;
        exit(1);
      }
      
      // just make a minimal beam pipe (tube) no cone for now as the sorting policy does not works as expected
      
      double zhalf = zEnd;      

           
      sname << "BeamPipeInnerWallCentralSection_" << i;
      _surface_names.push_back(sname.str());
      
      
      Add( new ILDCylinderMeasLayer(beam, beryllium , rInnerStart , zhalf, x0, y0, z0, bz, dummy,-1, _surface_names.back().c_str() ));
      
      // streamlog_out( DEBUG0 )   << " *** adding " << _surface_names.back() << " Measurement layer using CellID: [ beampipe ] at R = " << rInnerStart
      // << " zHalf = " << zhalf << " X0_in = " << beam.GetRadLength() << "  X0_out = " <<  beryllium.GetRadLength()
      // << std::endl ;
      
      sname.str(std::string());
      sname << "BeamPipeOuterWallCentralSection_" << i;
      _surface_names.push_back(sname.str());
      
      Add( new ILDCylinderMeasLayer(beryllium, air , rOuterStart , zhalf, x0, y0, z0, bz, dummy,-1, _surface_names.back().c_str() ));
      
      // streamlog_out( DEBUG0 )   << " *** adding " << _surface_names.back() << " Measurement layer using CellID: [ beampipe ] at R = " << rOuterStart
      // << " zHalf = " << zhalf << " X0_in = " << beryllium.GetRadLength() << "  X0_out = " <<  air.GetRadLength()
      // << std::endl ;

      
    
    } else {
      
      
      if( fabs( zEnd-zStart ) > epsilon ){
        
                
        double SortingPolicy = rInnerStart+i*eps_section;
         
        if (fabs(rInnerEnd-rInnerStart) > epsilon) { // use cone 

          sname.str(std::string());
          sname << "BeamPipeConeInnerFwd" << i;
          _surface_names.push_back(sname.str());
          
          Add( new ILDConeMeasLayer(beam, beryllium ,  zStart, rInnerStart,  zEnd, rInnerEnd, bz, SortingPolicy,          dummy,-1, _surface_names.back().c_str() ) );
          
          sname.str(std::string());
          sname << "BeamPipeConeInnerBwd" << i;
          _surface_names.push_back(sname.str());
          
          Add( new ILDConeMeasLayer(beam, beryllium , -zStart, rInnerStart, -zEnd, rInnerEnd, bz, SortingPolicy+eps_side, dummy,-1, _surface_names.back().c_str() ) );

        } else { // use cylinder
          
          sname.str(std::string());
          sname << "BeamPipeCylindeInnerFwd" << i;
          _surface_names.push_back(sname.str());
          
          const double zhalf   = fabs(zEnd-zStart);
          
          const double zoffset = zStart + zhalf;
          
          Add( new ILDCylinderMeasLayer(beam, beryllium , rInnerStart+i*eps_section          , zhalf, x0, y0,  zoffset, bz, dummy,-1, _surface_names.back().c_str() ));
          
          sname.str(std::string());
          sname << "BeamPipeCylindeInnerBwd" << i;
          _surface_names.push_back(sname.str());
          
          Add( new ILDCylinderMeasLayer(beam, beryllium , rInnerStart+i*eps_section+eps_side , zhalf, x0, y0, -zoffset, bz, dummy,-1, _surface_names.back().c_str() ));
          
          
        }
        
        
        // streamlog_out( DEBUG0 )   << " *** adding inner " << _surface_names[_surface_names.size()-2] << " and " << _surface_names.back() << " Measurement layer using CellID: [ beampipe ] at"
        // << " z1 = +-" << zStart
        // << " z2 = +-" << zEnd
        // << " r1 = " << rInnerStart
        // << " r2 = " << rInnerEnd
        // << " X0_in = " << beam.GetRadLength() << "  X0_out = " <<  beryllium.GetRadLength()
        // << std::endl ;
        
        SortingPolicy = rOuterStart+i*eps_section;

        if (fabs(rOuterEnd-rOuterStart) > epsilon) { // use cone
        
          sname.str(std::string());
          sname << "BeamPipeConeOuterFwd" << i;
          _surface_names.push_back(sname.str());
          
          Add( new ILDConeMeasLayer(beryllium , air ,  zStart, rOuterStart,  zEnd, rOuterEnd, bz, SortingPolicy,          dummy,-1, _surface_names.back().c_str() ) );
          
          sname.str(std::string());
          sname << "BeamPipeConeOuterBwd" << i;
          _surface_names.push_back(sname.str());
          
          Add( new ILDConeMeasLayer(beryllium , air , -zStart, rOuterStart, -zEnd, rOuterEnd, bz, SortingPolicy+eps_side, dummy,-1, _surface_names.back().c_str() ) );
                    
          
        } else { // use cylinder
          
          sname.str(std::string());
          sname << "BeamPipeCylindeOuterFwd" << i;
          _surface_names.push_back(sname.str());
        
          const double zhalf   = fabs(zEnd-zStart);
          
          const double zoffset = zStart + zhalf;
          
          Add( new ILDCylinderMeasLayer(beryllium, air , rOuterStart+i*eps_section          , zhalf, x0, y0,  zoffset, bz, dummy,-1, _surface_names.back().c_str() ));
          
          sname.str(std::string());
          sname << "BeamPipeCylindeOuterBwd" << i;
          _surface_names.push_back(sname.str());
          
          Add( new ILDCylinderMeasLayer(beryllium, air , rOuterStart+i*eps_section+eps_side , zhalf, x0, y0, -zoffset, bz, dummy,-1, _surface_names.back().c_str() ));
                  
        }

        
        // streamlog_out( DEBUG0 )   << " *** adding outer " << _surface_names[_surface_names.size()-2] << " and " << _surface_names.back() << " Measurement layer using CellID: [ beampipe ] at"
        // << " z1 = +-" << zStart
        // << " z2 = +-" << zEnd
        // << " r1 = " << rOuterStart
        // << " r2 = " << rOuterEnd
        // << " X0_in = " << beryllium.GetRadLength() << "  X0_out = " <<  air.GetRadLength()
        // << std::endl ;
        
      }
      
    }
    
  }

  
  // add vacuum layer 1mm inside the beam pipe to assist propagation to the IP
  // therefore make a cylinder that is 1mm smaller than the lowest RInner value of the cones
  // and make it so long that it gets as long as the beamtube
  const double rvacuum = *min_element(rInner.begin(), rInner.end()) - 1.0; 
  const double zHalfVacuum = z.back();

  const double x0 = 0.0;
  const double y0 = 0.0;
  const double z0 = 0.0;

  
  _ipLayer = new ILDCylinderMeasLayer(beam, beam , rvacuum , zHalfVacuum, x0, y0, z0, bz, dummy,-1,"IPLayer" );
  
  Add( _ipLayer );
  // streamlog_out( DEBUG0 )   << " *** adding " << "IPLayer" << " Measurement layer using CellID: [ beampipe ] at R = " << rvacuum
  // << " zHalf = " << zHalfVacuum << " X0_in = " << beam.GetRadLength() << "  X0_out = " <<  beam.GetRadLength()    
  // << std::endl ;  
  
  
  // add calo bounding inner surface as ILDPolygonBarrelMeasLayer and two planes at the inner z of the calo endcaps with ILDDiscMeasLayer

  // SJA:FIXME: OK for now we will just use a set of planes as there is not yet an implementation to propagate to a multilayer
  
  const gear::CalorimeterParameters& ecalB = gearMgr.getEcalBarrelParameters();
  const gear::CalorimeterParameters& ecalE = gearMgr.getEcalEndcapParameters();
  
  if (ecalB.getSymmetryOrder()!=8) {
    // streamlog_out(ERROR) << "ILDSupportKalDetector::ILDSupportKalDetector ECal barrel is not eightfold symmetry: exit(1) called from " << __FILE__ << "   line " << __LINE__ << std::endl; 
    exit(1);

  }
  
  double phi0  = ecalB.getPhi0();
  double r_min_ecal_bar = ecalB.getExtent()[0];
  double z_max_ecal_bar = ecalB.getExtent()[3];

  UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ; 
  encoder.reset() ;  // reset to 0
  
  encoder[lcio::ILDCellID0::subdet] = lcio::ILDDetID::ECAL ;
  encoder[lcio::ILDCellID0::side] = lcio::ILDDetID::barrel;
  encoder[lcio::ILDCellID0::layer]  = 0 ;
  
  std::vector<int> module_ids;
  
  for (int i=0; i<8; ++i) {

    encoder[lcio::ILDCellID0::module] = i;
    module_ids.push_back(encoder.lowWord());

    double segment_dphi = 2.0*M_PI / 8; 
    double phi = i*segment_dphi+phi0;
    double width = 2.0*(r_min_ecal_bar*tan(segment_dphi*0.5));
    double length = 2.0*z_max_ecal_bar;
    Add ( new ILDParallelPlanarMeasLayer(air,air,r_min_ecal_bar,phi,bz,r_min_ecal_bar+i*1.0e-06,width,length,0.0,0.0,0.0,true,encoder.lowWord(),"ECalBarrelFace"));    
    
  }
  
  //  Add ( new ILDPolygonBarrelMeasLayer(air,air,bz,r_min_ecal_bar,r_min_ecal_bar,z_max_ecal_bar*0.5,8,0.0,phi0,module_ids,"ECalBarrelFace"));  
  
  // streamlog_out( DEBUG0 )   << " *** adding ECalBarrelFace Measurement layer at Rmin = " << r_min_ecal_bar << std::endl ;
  
  double r_max_ecal_ecap = ecalE.getExtent()[1];
  double z_min_ecal_ecap = ecalE.getExtent()[2];

  encoder[lcio::ILDCellID0::module] = 0;
  
  encoder[lcio::ILDCellID0::side] = lcio::ILDDetID::fwd;
  
  TVector3 front_face_centre_fwd( 0.0, 0.0, z_min_ecal_ecap); // for +z  
  TVector3 front_face_normal_fwd(front_face_centre_fwd);
  front_face_normal_fwd.SetMag(1.0);
  
  Add (new ILDDiscMeasLayer(air, air, front_face_centre_fwd,front_face_normal_fwd, bz, z_min_ecal_ecap, 0., r_max_ecal_ecap/cos(M_PI/8.0), true, encoder.lowWord(),"ECalEndcapFace+Z") );
  
  // streamlog_out( DEBUG0 )   << " *** adding ECalEndcapFace+Z Measurement layer at Zmin = " << front_face_centre_fwd.z() << " and Rmax = " << r_max_ecal_ecap/cos(M_PI/8.0) << std::endl ;
  
  encoder[lcio::ILDCellID0::side] = lcio::ILDDetID::bwd;
  TVector3 front_face_centre_bwd( -front_face_centre_fwd ); // for -z  
  TVector3 front_face_normal_bwd(front_face_centre_bwd);
  front_face_normal_bwd.SetMag(1.0);
  
  Add (new ILDDiscMeasLayer(air, air, front_face_centre_bwd,front_face_normal_bwd, bz, z_min_ecal_ecap+1.0e-06, 0., r_max_ecal_ecap/cos(M_PI/8.0), true, encoder.lowWord(),"ECalEndcapFace-Z"));

  // streamlog_out( DEBUG0 )   << " *** adding ECalEndcapFace-Z Measurement layer at Zmin = " << front_face_centre_bwd.z() << " and Rmax = " << r_max_ecal_ecap/cos(M_PI/8.0) << std::endl ;
  
  SetOwner();
}

