// $Id: $
//==========================================================================
//  Detector description implementation for CEPC
//--------------------------------------------------------------------------
// All rights reserved.
//
// For the licensing terms see $DD4hepINSTALL/LICENSE.
// For the list of contributors see $DD4hepINSTALL/doc/CREDITS.
// 
// Author     : FU Chengdong
//==========================================================================
// Specialized generic detector constructor for skew planar
//              _______________
//              \      |      /
//               \     |     /
//                \    |    /
//                 \   |   /
//                  \__|__/
//                  center                       
//
// z       :      z position at maximum radius of center
// dz      :      z gap between odd and even modules
// gap     :      rphi gap between neighbouring modules
// rmin    :      minimum radius at center
// rmax    :      maximum radius at center
// phi0    :      start phi for #0 module
// skew    :      skew angle
// nmodules:      module number at phi direction
// is_pixel:      pixel tag for reconstruction
//==========================================================================
#include <DD4hep/Detector.h>
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Utilities.h"
#include "DDRec/Surface.h"
#include "DDRec/DetectorData.h"
#include "Math/AxisAngle.h"

#include <map>

using namespace std;

using dd4hep::DetElement;
using dd4hep::Detector;
using dd4hep::Material;
using dd4hep::PlacedVolume;
using dd4hep::SensitiveDetector;
using dd4hep::Volume;
using dd4hep::Trap;
using dd4hep::_toString;
using dd4hep::Position;
using dd4hep::Transform3D;
using dd4hep::Rotation3D;
using dd4hep::RotationZYX;

static dd4hep::Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens){
  typedef vector<PlacedVolume> Placements;
  xml_det_t   x_det     = e;
  Material    air       = description.air();
  int         det_id    = x_det.id();
  string      name      = x_det.nameStr();
  bool        reflect   = x_det.reflect(false);
  DetElement  tracker(name, det_id);

  Volume envelope = dd4hep::xml::createPlacedEnvelope(description, e, tracker);
  dd4hep::xml::setDetectorTypeFlag(e, tracker) ;
  if(description.buildType()==dd4hep::BUILD_ENVELOPE) return tracker;
  envelope.setVisAttributes(description.visAttributes("SeeThrough"));

  sens.setType("tracker");
  std::cout << " ** building SiTrackerEndcapRing_v01 ..." << std::endl ;

  dd4hep::rec::ZDiskPetalsData*  zDiskPetalsData = new dd4hep::rec::ZDiskPetalsData;

  PlacedVolume pv;

  for(xml_coll_t li(x_det,_U(layer)); li; ++li)  {
    xml_comp_t  x_layer(li);
    int layer_id    = x_layer.id();

    double zstart   = x_layer.z();
    double dz       = x_layer.dz();
    double rmin     = x_layer.inner_r();
    double rmax     = x_layer.outer_r();
    double phi0     = x_layer.phi0(0);
    double skew     = x_layer.skew(0);
    double gap      = x_layer.gap();
    bool   is_pixel = x_layer.attr<bool>(_Unicode(is_pixel)); 
    int    nmodules = x_layer.nmodules();

    double dphi       = 2*M_PI/nmodules;
    double phi        = phi0;
    double x1         = rmin*tan(0.5*dphi) - 0.5*gap;
    double x2         = rmax*tan(0.5*dphi) - 0.5*gap;
    double half_width = (rmax-rmin)/2/cos(skew);
    double rpos       = (rmax+rmin)/2;
    double zpos       = zstart + half_width*sin(skew);

    double layer_thickness        = 0.;
    double sensitive_thickness[2] = {0.,0.};
    double support_thickness      = 0.;
    int    nsensor                = 0;
    for(xml_coll_t ci(x_layer,_U(component)); ci; ++ci){
      xml_comp_t c     = ci;
      double thickness = c.thickness();
      layer_thickness  += thickness;
      if(c.isSensitive()){
	sensitive_thickness[nsensor] = thickness;
	nsensor++;
      }
      else      support_thickness   += thickness;
    }
    double half_thickness = layer_thickness/2;
    double zshift_support = dz + 0.5*(sensitive_thickness[0] + sensitive_thickness[1] + support_thickness);

    //std::cout << " ****Layer: " << layer_id << " nmodules = " << nmodules << " z = " << zpos << " zshift = " << zshift_support
    //	      << " rmin = " << rmin << " rmax = " << rmax << " x1 = " << x1 << " x2 = " << x2 << " thickness = " << layer_thickness << std::endl;    
    
    Trap   moduleSolid(half_thickness, 0, 0, half_width, x1, x2, 0, half_width, x1, x2, 0);
    Volume moduleVol(_toString(layer_id,"layer%d"), moduleSolid, air);
    moduleVol.setVisAttributes(description.visAttributes(x_layer.visStr()));
    
    Placements sensitives;
    int    sensor_id = 1;
    int    c_id      = 0;
    double c_pos     = -half_thickness;
    for(xml_coll_t ci(x_layer,_U(component)); ci; ++ci, c_id++){
      xml_comp_t c       = ci;
      double     c_thickness = c.thickness();
      Material   c_mat   = description.material(c.materialStr());
      string     c_name  = _toString(c_id,"component%d");
      Trap       c_solid(c_thickness/2, 0, 0, half_width, x1, x2, 0, half_width, x1, x2, 0);
      Volume     c_vol(c_name, c_solid, c_mat);
      c_vol.setVisAttributes(description.visAttributes(c.visStr()));
      pv = moduleVol.placeVolume(c_vol,Position(0,0,c_pos+c_thickness/2));
      if(c.isSensitive()){
        tracker.check(sensor_id > 2," fromCompact: "+c_name+" Max of 2 modules allowed!");
        pv.addPhysVolID("sensor", sensor_id);
        c_vol.setSensitiveDetector(sens);
	sensitives.push_back(pv);
        ++sensor_id;
      }
      c_pos += c_thickness;
    }
    
    for(int module_id=0; module_id<nmodules; module_id++){
      
      string m_base = _toString(layer_id,"layer%d") + _toString(module_id,"_module%d");
      
      double x = rpos*std::cos(phi);
      double y = rpos*std::sin(phi);
      double zshift_layer   = zshift_support - 0.5*sensitive_thickness[0] + 0.5*sensitive_thickness[1];
      //std::cout << "****** module_id = " << module_id << " phi = " << phi << " x = " << x << " y = " << y
      //	<< " zsup = " << zpos+zshift_support << " zsens = " << zpos+zshift_support-0.5*sensitive_thickness[0]-0.5*support_thickness << std::endl;
      DetElement module(tracker, m_base+"_pos", det_id);
      Rotation3D rot = Rotation3D(ROOT::Math::AxisAngle(dd4hep::PositionPolar(1,M_PI/2,-M_PI/2+phi),-skew))*Rotation3D(RotationZYX(-M_PI/2+phi,0,0));
      pv = envelope.placeVolume(moduleVol, Transform3D(rot, Position(x,y,zpos+zshift_layer)));
      pv.addPhysVolID("side",1).addPhysVolID("layer", layer_id).addPhysVolID("module",module_id);
      module.setPlacement(pv);
      for(size_t ic=0; ic<sensitives.size(); ++ic)  {
	PlacedVolume sens_pv = sensitives[ic];
	DetElement comp_elt(module, sens_pv.volume().name(), module_id);
	comp_elt.setPlacement(sens_pv);
      }

      if(reflect){
	Rotation3D rotRef = Rotation3D(ROOT::Math::AxisAngle(dd4hep::PositionPolar(1,M_PI/2,-M_PI/2+phi),skew))*Rotation3D(RotationZYX(M_PI/2-phi,M_PI,0));
	pv = envelope.placeVolume(moduleVol, Transform3D(rotRef, Position(x,y,-zpos-zshift_layer)));
	pv.addPhysVolID("side",-1).addPhysVolID("layer",layer_id).addPhysVolID("module",module_id);
	DetElement r_module(tracker, m_base+"_neg", det_id);
	r_module.setPlacement(pv);
	for(size_t ic=0; ic<sensitives.size(); ++ic)  {
	  PlacedVolume sens_pv = sensitives[ic];
	  DetElement comp_elt(r_module, sens_pv.volume().name(), module_id);
	  comp_elt.setPlacement(sens_pv);
	}
      }
      zshift_support = -zshift_support;
      phi   +=  dphi;
    }

    dd4hep::rec::ZDiskPetalsData::LayerLayout thisLayer;
    thisLayer.typeFlags[ dd4hep::rec::ZDiskPetalsData::SensorType::DoubleSided ] = bool(sensor_id>2);
    thisLayer.typeFlags[ dd4hep::rec::ZDiskPetalsData::SensorType::Pixel ]       = is_pixel;
    thisLayer.petalHalfAngle      = dphi/2;
    thisLayer.alphaPetal          = skew;
    thisLayer.zPosition           = zpos;
    thisLayer.petalNumber         = nmodules;
    thisLayer.sensorsPerPetal     = sensor_id-1;
    thisLayer.phi0                = phi0;
    thisLayer.zOffsetSupport      = fabs(zshift_support);
    thisLayer.distanceSupport     = rmin;
    thisLayer.thicknessSupport    = support_thickness;
    thisLayer.widthInnerSupport   = 2.*x1;
    thisLayer.widthOuterSupport   = 2.*x2;
    thisLayer.lengthSupport       = 2.*half_width;
    thisLayer.zOffsetSensitive    = fabs(zshift_support) - 0.5*support_thickness - 0.5*sensitive_thickness[0];
    thisLayer.distanceSensitive   = rmin;
    thisLayer.thicknessSensitive  = sensitive_thickness[0];
    thisLayer.widthInnerSensitive = 2.*x1;
    thisLayer.widthOuterSensitive = 2.*x2;
    thisLayer.lengthSensitive     = 2.*half_width;

    zDiskPetalsData->layers.push_back(thisLayer);
  }

  dd4hep::xml::Component recPar = x_det.child(_Unicode(reconstruction));
  const double strip_width  = recPar.attr< double >(_Unicode(strip_width));
  const double strip_length = recPar.attr< double >(_Unicode(strip_length));
  const double strip_pitch  = recPar.attr< double >(_Unicode(strip_pitch));
  const double strip_angle  = recPar.attr< double >(_Unicode(strip_angle));

  zDiskPetalsData->widthStrip  = strip_width;
  zDiskPetalsData->lengthStrip = strip_length;
  zDiskPetalsData->pitchStrip  = strip_pitch;
  zDiskPetalsData->angleStrip  = strip_angle;

  std::cout << (*zDiskPetalsData) << std::endl;
  tracker.addExtension<dd4hep::rec::ZDiskPetalsData>( zDiskPetalsData );
  
  if ( x_det.hasAttr(_U(combineHits)) ) {
    tracker.setCombineHits(x_det.attr<bool>(_U(combineHits)),sens);
  }
  
  return tracker;
}

DECLARE_DETELEMENT(SiTrackerSkewRing_v01, create_detector)
