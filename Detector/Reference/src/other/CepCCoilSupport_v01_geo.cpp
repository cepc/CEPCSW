//====================================================================
// CepC BeamPipe models in DD4hep 
//--------------------------------------------------------------------
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetType.h"
#include "XML/Utilities.h"
//#include "XMLHandlerDB.h"
#include <cmath>
#include <map>
#include <string>
#include "TGeoMatrix.h"

using dd4hep::Assembly;
using dd4hep::BUILD_ENVELOPE;
using dd4hep::DetElement;
using dd4hep::Detector;
using dd4hep::SensitiveDetector;

//Coil support
static dd4hep::Ref_t create_detector(Detector& theDet,
				     xml_h element,
				     SensitiveDetector /*sens*/) {

  std::cout << "This is the Coil support:"  << std::endl;

  //Access to the XML File
  xml_det_t x_support = element;
  const std::string name = x_support.nameStr();

  DetElement support(  name, x_support.id()  ) ;

  // --- create an envelope volume and position it into the world ---------------------
  dd4hep::Volume envelope = dd4hep::xml::createPlacedEnvelope( theDet,  element , support ) ;
  
  dd4hep::xml::setDetectorTypeFlag( element, support ) ;

  if( theDet.buildType() == BUILD_ENVELOPE ) return support ;
    
  dd4hep::Material holeMaterial    = theDet.material("G4_AIR");

  int ip = 0;
  for(xml_coll_t pi( x_support,Unicode("plate")); pi; ++pi, ++ip) {
    xml_comp_t x_p(pi);

    const double rmin       = x_p.attr< double > (_Unicode(cellRmin));
    const double thickness  = x_p.attr< double > (_Unicode(cellthickness));
    const double height     = x_p.attr< double > (_Unicode(cellheight));
    const int numx          = x_p.attr< int > (_Unicode(numx));
    const int numy          = x_p.attr< int > (_Unicode(numy));
        
    const std::string volName      = "Coil_" + x_p.nameStr();
    dd4hep::Material holderMaterial    = theDet.material(x_p.materialStr());
    
    xml_comp_t x_pos( x_p.child( _U(position) ) ); 
    const double x = x_pos.x();
    const double y = x_pos.y();
    const double z = x_pos.z();
    
    xml_comp_t x_rot( x_p.child( _U(rotation) ) );
    const double rotx = x_rot.x();
    const double roty = x_rot.y();
    const double rotz = x_rot.z();

    std::cout << "plate: "
	      << std::setw(8) << rmin      /dd4hep::mm
	      << std::setw(8) << thickness	     /dd4hep::mm
	      << std::setw(8) << height /dd4hep::mm
	      << std::setw(8) << numx   /dd4hep::mm
	      << std::setw(8) << numy        /dd4hep::mm
	      << std::setw(15) << x_p.materialStr()
	      << std::setw(35) << volName
	      << std::setw(8) << x           /dd4hep::mm
	      << std::setw(8) << y           /dd4hep::mm
	      << std::setw(8) << z           /dd4hep::mm
	      << std::setw(8) << rotx           /dd4hep::degree
	      << std::setw(8) << roty           /dd4hep::degree
	      << std::setw(8) << rotz           /dd4hep::degree
	      << std::endl;    

    double xhalf = numx*(rmin+thickness);
    double yhalf = ((int)(numy/2)+1)*(rmin+thickness)/cos(30*dd4hep::degree)+(int)(numy/2)*(rmin+thickness)*tan(30*dd4hep::degree);
    dd4hep::Box plate(xhalf, yhalf, 0.5*height);
    dd4hep::Volume vol_plate(volName+"_air", plate, holeMaterial);
    dd4hep::PlacedVolume pv_plate = envelope.placeVolume(vol_plate, dd4hep::Transform3D(dd4hep::RotationZYX(rotz,roty,rotx), dd4hep::Position(x,y,z)));
    vol_plate.setVisAttributes(theDet, "SeeThrough");
    std::string detName =  dd4hep::_toString( ip , "plate%02d") ;
    DetElement plateDE( support, detName, x_support.id() );
    plateDE.setPlacement( pv_plate ) ;

    dd4hep::PolyhedraRegular poly(6,rmin,rmin+thickness,0.5*height);
    dd4hep::Volume vol_holder(volName+"holder", poly, holderMaterial);
    vol_holder.setVisAttributes(theDet, "HolderVis");

    double yc = -yhalf;
    for(int j=0;j<numy;j++){
      yc += (rmin+thickness)/cos(30*dd4hep::degree);
      if(j%2==0){
	for(int i=0;i<numx;i++){
	  double xc = -xhalf + (2*i+1)*(rmin+thickness);
	  vol_plate.placeVolume(vol_holder, dd4hep::Transform3D(dd4hep::RotationZ(30*dd4hep::degree),dd4hep::Position(xc,yc,0)));
	}
      }
      else{
	for(int i=0;i<numx-1;i++){
	  double xc = -xhalf + (2*i+2)*(rmin+thickness);
          vol_plate.placeVolume(vol_holder, dd4hep::Transform3D(dd4hep::RotationZ(30*dd4hep::degree),dd4hep::Position(xc,yc,0)));
	}
      }
      yc += (rmin+thickness)*tan(30*dd4hep::degree);
    }
    
  }  
  return support;
}
DECLARE_DETELEMENT(DD4hep_CepCCoilSupport_v01, create_detector)
