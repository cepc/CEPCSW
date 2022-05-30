//====================================================================
//  Detector description implementation of the Drift Chamber with only straight
//--------------------------------------------------------------------
//
//  Author: Chengdong FU
//
//====================================================================

#include <DD4hep/Detector.h>
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Utilities.h"
#include "DDRec/Surface.h"
#include "DDRec/DetectorData.h"

#include <map>

using namespace std;

using dd4hep::DetElement;
using dd4hep::Detector;
using dd4hep::Material;
using dd4hep::PlacedVolume;
using dd4hep::SensitiveDetector;
using dd4hep::Volume;
using dd4hep::Tube;
using dd4hep::_toString;
using dd4hep::Position;
using dd4hep::rec::FixedPadSizeTPCData;
using dd4hep::rec::Vector3D;
using dd4hep::rec::VolCylinder;
using dd4hep::rec::SurfaceType;
using dd4hep::rec::volSurfaceList;
using dd4hep::rec::VolPlane;

static dd4hep::Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens) {
  xml_det_t   x_det     = e;
  Material    air       = description.air();
  int         det_id    = x_det.id();
  string      name      = x_det.nameStr();
  DetElement  tracker(name, det_id);

  Volume envelope = dd4hep::xml::createPlacedEnvelope(description, e, tracker);
  dd4hep::xml::setDetectorTypeFlag(e, tracker) ;
  if(description.buildType()==dd4hep::BUILD_ENVELOPE) return tracker;
  envelope.setVisAttributes(description.visAttributes("SeeThrough"));

  sens.setType("tracker");
  std::cout << " ** building StraightDriftChamber_v01 ..." << std::endl ;
  
  xml_coll_t c_shell(x_det,_Unicode(shell));
  xml_comp_t x_shell = c_shell;

  xml_coll_t c_inner(x_shell,_U(inner));
  xml_comp_t x_inner = c_inner;
  double thickness_inner = x_inner.thickness();
  Material mat_inner(description.material(x_inner.materialStr())); 

  xml_coll_t c_outer(x_shell,_U(outer));
  xml_comp_t x_outer = c_outer;
  double thickness_outer = x_outer.thickness();
  Material mat_outer(description.material(x_outer.materialStr()));

  xml_coll_t c_side(x_shell,_U(side));
  xml_comp_t x_side = c_side;
  double thickness_side = x_side.thickness();
  Material mat_side(description.material(x_side.materialStr()));

  xml_coll_t c_chamber(x_det,_U(chamber));
  xml_comp_t x_chamber = c_chamber;
  double rmin_chamber = x_chamber.rmin();
  double rmax_chamber = x_chamber.rmax();
  double halflength_chamber = x_chamber.zhalf();
  Material mat_chamber(description.material(x_chamber.materialStr()));

  xml_coll_t c_signal(x_det,_Unicode(signal_wire));
  xml_comp_t x_signal = c_signal;
  double rmin_signal = x_signal.rmin();
  double rmax_signal = x_signal.rmax();
  Material mat_signal = (description.material(x_signal.materialStr()));
  Tube solid_signal(rmin_signal,rmax_signal,halflength_chamber);
  Volume volume_signal(name+"_signal",solid_signal,mat_signal);
  for(xml_coll_t c(x_signal, _U(tubs)); c; ++c) {
    xml_comp_t x_tub = c;
    string name_tub = x_tub.nameStr();
    double rmin_tub = x_tub.rmin();
    double rmax_tub = x_tub.rmax();
    Material mat_tub(description.material(x_tub.materialStr()));
    Tube solid_tub(rmin_tub,rmax_tub,halflength_chamber);
    Volume volume_tub(name+"_signal"+name_tub,solid_tub,mat_tub);
    volume_tub.setVisAttributes(description.visAttributes(x_tub.visStr()));
    PlacedVolume phy_tub = volume_signal.placeVolume(volume_tub);
  }

  xml_coll_t c_field(x_det,_Unicode(field_wire));
  xml_comp_t x_field = c_field;
  double rmin_field = x_field.rmin();
  double rmax_field = x_field.rmax();
  Material mat_field = (description.material(x_field.materialStr()));
  Tube solid_field(rmin_field,rmax_field,halflength_chamber);
  Volume volume_field(name+"_field",solid_field,mat_field);
  for(xml_coll_t c(x_field, _U(tubs)); c; ++c) {
    xml_comp_t x_tub = c;
    string name_tub = x_tub.nameStr();
    double rmin_tub = x_tub.rmin();
    double rmax_tub = x_tub.rmax();
    Material mat_tub(description.material(x_tub.materialStr()));
    Tube solid_tub(rmin_tub,rmax_tub,halflength_chamber);
    Volume volume_tub(name+"_field"+name_tub,solid_tub,mat_tub);
    volume_tub.setVisAttributes(description.visAttributes(x_tub.visStr()));
    PlacedVolume phy_tub = volume_field.placeVolume(volume_tub);
  }
		    
  Tube solid_chamber(rmin_chamber,rmax_chamber,halflength_chamber);
  Volume volume_chamber(name+"_chamber",solid_chamber,mat_chamber);
  volume_chamber.setVisAttributes(description.visAttributes(x_chamber.visStr()));
  PlacedVolume phy_chamber = envelope.placeVolume(volume_chamber);
  if(x_det.hasAttr(_U(id))){
    phy_chamber.addPhysVolID("system",x_det.id());
  }
  DetElement det_chamber(tracker, name+"_chamber", 0);
  det_chamber.setPlacement(phy_chamber);

  Tube solid_inner(rmin_chamber-thickness_inner, rmin_chamber, halflength_chamber);
  Volume volume_inner(name+"_inner_wall",solid_inner,mat_inner);
  volume_inner.setVisAttributes(description.visAttributes(x_inner.visStr()));
  PlacedVolume phy_inner = envelope.placeVolume(volume_inner);
  DetElement det_inner(tracker, name+"_inner_wall", 0);
  det_inner.setPlacement(phy_inner);
  Vector3D ocyl_inner(rmin_chamber-0.5*thickness_inner, 0., 0.);
  VolCylinder surfI(volume_inner, SurfaceType(SurfaceType::Helper), 0.5*thickness_inner, 0.5*thickness_inner, ocyl_inner);
  volSurfaceList(tracker)->push_back(surfI);

  Tube solid_outer(rmax_chamber, rmax_chamber+thickness_outer, halflength_chamber);
  Volume volume_outer(name+"_outer_wall",solid_outer,mat_outer);
  volume_outer.setVisAttributes(description.visAttributes(x_outer.visStr()));
  PlacedVolume phy_outer = envelope.placeVolume(volume_outer);
  DetElement det_outer(tracker, name+"_outer_wall", 0);
  det_outer.setPlacement(phy_outer);
  Vector3D ocyl_outer(rmax_chamber+0.5*thickness_inner, 0., 0.);
  VolCylinder surfO(volume_outer, SurfaceType(SurfaceType::Helper), 0.5*thickness_outer, 0.5*thickness_outer, ocyl_outer);
  volSurfaceList(tracker)->push_back(surfO);

  Tube solid_side(rmin_chamber-thickness_inner, rmax_chamber+thickness_inner, 0.5*thickness_side);
  Volume volume_side(name+"_side_wall",solid_side,mat_side);
  volume_side.setVisAttributes(description.visAttributes(x_side.visStr()));
  PlacedVolume phy_plus  = envelope.placeVolume(volume_side,Position(0,0,halflength_chamber+0.5*thickness_side));
  PlacedVolume phy_minus = envelope.placeVolume(volume_side,Position(0,0,-halflength_chamber-0.5*thickness_side));
  DetElement det_plus(tracker, name+"_plusside_wall", 0);
  det_plus.setPlacement(phy_plus);
  DetElement det_minus(tracker, name+"_minusside_wall", 1);
  det_minus.setPlacement(phy_minus);
  Vector3D u(0., 1., 0.);
  Vector3D v(1., 0., 0.);
  Vector3D n(0., 0., 1.);
  double mid_r = 0.5*(rmin_chamber-thickness_inner+rmax_chamber+thickness_inner);
  Vector3D o(0., mid_r, 0.);
  VolPlane surfS(volume_side, SurfaceType(SurfaceType::Helper), 0.5*thickness_side, 0.5*thickness_side, u, v, n, o);
  volSurfaceList(det_plus)->push_back(surfS);
  volSurfaceList(det_minus)->push_back(surfS);

  int chamber_id=0, max_layer=0;
  double rmin_sensitive=10000, rmax_sensitive=0, cell_height=0;
  for(xml_coll_t c(x_chamber, _U(layer)); c; ++c) {
    xml_comp_t x_layer = c;
    string name_layer = x_layer.nameStr();
    double rmin_layer = x_layer.rmin();
    double rmax_layer = x_layer.rmax();
    Material mat_layer(description.material(x_layer.materialStr()));
    Tube solid_sub(rmin_layer, rmax_layer, halflength_chamber);
    Volume volume_sub(name+name_layer,solid_sub,mat_layer);
    volume_sub.setVisAttributes(description.visAttributes("SeeThrough"));
    PlacedVolume phy_sub = volume_chamber.placeVolume(volume_sub);
    DetElement det_sub(det_chamber, name+name_layer, x_layer.id());
    det_sub.setPlacement(phy_sub);

    int nlayer = x_layer.number();
    double height = (rmax_layer-rmin_layer)/nlayer;
    double radius = rmin_layer;
    for(int layer_id=0; layer_id<nlayer; layer_id++){
      Tube solid_layer(radius, radius+height, halflength_chamber);
      Volume volume_layer(name+name_layer+_toString(layer_id, "_%d"),solid_layer,mat_layer);
      volume_layer.setVisAttributes(description.visAttributes(x_layer.visStr()));
      
      PlacedVolume phy = volume_sub.placeVolume(volume_layer);
      if(x_layer.isSensitive()){
	volume_layer.setSensitiveDetector(sens);
	volume_layer.setLimitSet(description,x_det.limitsStr());
	phy.addPhysVolID("chamber", chamber_id).addPhysVolID("layer", layer_id);
	if(radius<rmin_sensitive){
	  rmin_sensitive = radius;
	  // TODO: more than one chamber, with different height; now only minimum radius chamber include
	  rmax_sensitive = radius + nlayer*height;
	  max_layer = nlayer;
	  cell_height = height;
	}
	double radius_center = radius+0.5*height;
	double radius_edge   = radius+rmax_field;
	int ncell = floor(dd4hep::twopi*radius_center/height);
	double dphi = dd4hep::twopi / ncell;
	double offset = 0.;
	if(layer_id%2!=0) offset = 0.5*dphi;
	//std::cout << "debug: " << layer_id << " rmid=" << radius_center << " redge=" << radius_edge
	//	  << " ncell=" << ncell << " dphi=" << dphi << " offset=" << offset << std::endl;
	for(int icell=0;icell<ncell;icell++){
	  double phi = (icell+0.5)*dphi + offset;
	  volume_layer.placeVolume(volume_signal, Position(radius_center*cos(phi),radius_center*sin(phi),0));
	  volume_layer.placeVolume(volume_field, Position(radius_center*cos(phi+0.5*dphi),radius_center*sin(phi+0.5*dphi),0));
	  volume_layer.placeVolume(volume_field, Position(radius_edge*cos(phi+0.5*dphi),radius_edge*sin(phi+0.5*dphi),0));
	  volume_layer.placeVolume(volume_field, Position(radius_edge*cos(phi+0.25*dphi),radius_edge*sin(phi+0.25*dphi),0));
	  volume_layer.placeVolume(volume_field, Position(radius_edge*cos(phi),radius_edge*sin(phi),0));
	  volume_layer.placeVolume(volume_field, Position(radius_edge*cos(phi-0.25*dphi),radius_edge*sin(phi-0.25*dphi),0));
	  
	  if(layer_id==nlayer){
	    double radius_max = radius+height-rmax_field;
	    volume_layer.placeVolume(volume_field, Position(radius_max*cos(phi+0.5*dphi),radius_max*sin(phi+0.5*dphi),0));
	    volume_layer.placeVolume(volume_field, Position(radius_max*cos(phi+0.25*dphi),radius_max*sin(phi+0.25*dphi),0));
	    volume_layer.placeVolume(volume_field, Position(radius_max*cos(phi),radius_max*sin(phi),0));
	    volume_layer.placeVolume(volume_field, Position(radius_max*cos(phi-0.25*dphi),radius_max*sin(phi-0.25*dphi),0));
	  }
	}
      }
      DetElement det_layer(det_sub, name+name_layer+_toString(layer_id, "_%d"), layer_id);
      det_layer.setPlacement(phy);
      Vector3D ol(radius+0.5*height, 0., 0.);
      SurfaceType type = x_layer.isSensitive()?SurfaceType(SurfaceType::Sensitive, SurfaceType::Invisible):SurfaceType(SurfaceType::Helper);
      VolCylinder surf(volume_layer, type, 0.5*height, 0.5*height, ol);
      volSurfaceList(det_layer)->push_back(surf);

      radius += height;
    }
    if(x_layer.isSensitive()) chamber_id++;
  }

  FixedPadSizeTPCData* dcData = new FixedPadSizeTPCData;
  dcData->zHalf = halflength_chamber+thickness_side;
  dcData->rMin = rmin_chamber-thickness_inner;
  dcData->rMax = rmax_chamber+thickness_outer;
  dcData->innerWallThickness = thickness_inner;
  dcData->outerWallThickness = thickness_outer;
  dcData->rMinReadout = rmin_sensitive;
  dcData->rMaxReadout = rmax_sensitive;
  dcData->maxRow = max_layer;
  dcData->padHeight = cell_height;
  dcData->padWidth = cell_height;
  dcData->driftLength = halflength_chamber;
  dcData->zMinReadout = thickness_side/2.0;
  tracker.addExtension<FixedPadSizeTPCData>(dcData);
  
  return tracker;
}

DECLARE_DETELEMENT(StraightDriftChamber_v01, create_detector)
