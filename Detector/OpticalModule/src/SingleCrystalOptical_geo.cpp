#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Detector.h"

#include <iostream>
#include <map>

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

static Ref_t create_detector(Detector &description, xml_h e, SensitiveDetector sens)  {
  xml_det_t x_det     = e;
  string det_name     = x_det.nameStr();

  xml_det_t x_crystal = x_det.child(_Unicode(crystal));
  xml_det_t x_wrap    = x_det.child(_Unicode(wrap));
  xml_det_t x_pmt     = x_det.child(_Unicode(PM));
  xml_det_t x_cath    = x_pmt.child(_Unicode(cathode));
  xml_det_t x_shell   = x_pmt.child(_Unicode(shell));

  DetElement sdet(det_name,x_det.id());

  sens.setType("calorimeter");
  Assembly vol_hall("Hall");

  std::string type = x_pmt.attr< std::string >(_Unicode(type));
  
  double gap = 0.5*x_det.attr< double >(_Unicode(gap));
  double xhalf_crystal = 0.5*x_crystal.dx();
  double yhalf_crystal = 0.5*x_crystal.dy();
  double lhalf_crystal = 0.5*x_crystal.length();
  gap += 0.5*x_crystal.attr< double >(_Unicode(gap));

  bool active_cath  = x_cath.attr<bool>(_Unicode(active));
  bool active_shell = x_shell.attr<bool>(_Unicode(active)); 

  double t_PM   = x_pmt.thickness();
  double t_cath  = x_cath.thickness();
  double t_shell = x_shell.thickness();
  double t_glass = t_PM-t_shell;
  double t_wrap = 0;
  for(xml_coll_t li( x_wrap ,Unicode("layer")); li; ++li) {
    xml_comp_t x_layer(li);
    t_wrap += x_layer.thickness();
  }
  double depth_PM = x_pmt.depth();
  if(t_glass-t_cath<depth_PM+t_wrap) throw std::runtime_error("SingleCrystalOptical_geo.cpp: photocathode inside crystal, invalid!");
  double l_shell = t_PM-depth_PM-t_wrap;
  
  double x_PM = x_pmt.x();
  double y_PM = x_pmt.y();
  double z_PM = lhalf_crystal-depth_PM+0.5*t_glass;
  dd4hep::Position pos_right(x_PM,y_PM,z_PM);
  dd4hep::Position pos_left(x_PM,y_PM,-z_PM);

  dd4hep::Material material_crystal = description.material(x_crystal.attr<std::string>(_Unicode(material)));
  dd4hep::Material material_air = description.air();
  dd4hep::Material material_PM    = description.material(x_pmt.attr<std::string>(_Unicode(material)));
  dd4hep::Material material_cath = description.material(x_cath.attr<std::string>(_Unicode(material)));
  dd4hep::Material material_shell = description.material(x_shell.attr<std::string>(_Unicode(material)));

  Solid solid_PM;
  Volume vol_PM, vol_cath, vol_glass;
  std::vector<PlacedVolume> pvs;
  if(type=="Square"){
    Box box_shell(0.5*x_pmt.dx(), 0.5*x_pmt.dy(), 0.5*l_shell);
    Box box_glass(0.5*x_pmt.dx()-t_shell, 0.5*x_pmt.dy()-t_shell, 0.5*t_glass);
    solid_PM = UnionSolid(box_glass, box_shell, dd4hep::Position(0,0,0.5*(depth_PM+t_wrap)+0.5*t_shell));
    
    vol_glass = Volume("PM_Glass", box_glass, material_PM);
    vol_PM = Volume("PM_Shell", solid_PM, material_shell);

    Box box_cath(0.5*x_cath.dx(), 0.5*x_cath.dy(), 0.5*t_cath);
    vol_cath = Volume("Photocathode", box_cath, material_cath);
  }
  else if(type=="Circular"){
    Tube tube_shell(0, x_pmt.dr(), 0.5*l_shell, 0, 360*dd4hep::degree);
    Tube tube_glass(0, x_pmt.dr()-t_shell, 0.5*t_glass, 0, 360*dd4hep::degree);
    solid_PM = UnionSolid(tube_glass, tube_shell, dd4hep::Position(0,0,0.5*(depth_PM+t_wrap)+0.5*t_shell));
    
    vol_glass = Volume("PM_Glass", tube_glass, material_PM);
    vol_PM = Volume("PM_Shell", solid_PM, material_shell);
    
    Tube tube_cath(0, x_cath.dr(), 0.5*t_cath, 0, 360*dd4hep::degree);
    vol_cath = Volume("Photocathode", tube_cath, material_cath);
  }
  else{
    throw std::runtime_error("SingleCrystalOptical_geo.cpp: invalid PM type");
  }
  PlacedVolume pv_glass = vol_PM.placeVolume(vol_glass, dd4hep::Position(0,0,0));
  PlacedVolume pv_cath  = vol_glass.placeVolume(vol_cath, dd4hep::Position(0,0,0.5*t_glass-0.5*t_cath));

  vol_PM.setVisAttributes(description, x_pmt.visStr());
  vol_cath.setVisAttributes(description, x_cath.visStr());

  if(active_cath) vol_cath.setSensitiveDetector(sens);
  else            vol_glass.setSensitiveDetector(sens);
  
  Box box_crystal(xhalf_crystal,yhalf_crystal,lhalf_crystal);
  SubtractionSolid tmp_crystal(box_crystal, solid_PM, pos_right);
  SubtractionSolid solid_crystal(tmp_crystal, solid_PM, dd4hep::Transform3D(dd4hep::RotationY(180*dd4hep::degree),pos_left));
  Volume vol_crystal("Crystal", solid_crystal, material_crystal);
  vol_crystal.setVisAttributes(description, x_crystal.visStr());
  
  double thickness = 0;
  std::vector<Volume> vols;
  vols.push_back(vol_crystal);
  for(xml_coll_t li( x_wrap ,Unicode("layer")); li; ++li) {
    xml_comp_t x_layer(li);
    thickness += x_layer.thickness();
    
    Box box_wrap(xhalf_crystal+thickness, yhalf_crystal+thickness, lhalf_crystal+thickness);
    SubtractionSolid tmp_wrap(box_wrap, solid_PM, pos_right);
    SubtractionSolid solid_wrap(tmp_wrap, solid_PM, dd4hep::Transform3D(dd4hep::RotationY(180*dd4hep::degree),pos_left));
    std::string name = std::string("Wrap_")+x_layer.materialStr(); 
    Volume vol_wrap(name, solid_wrap, description.material(x_layer.attr<std::string>(_Unicode(material))));
    PlacedVolume pv_into = vol_wrap.placeVolume(vols.back(), dd4hep::Position(0,0,0));
    
    vol_wrap.setVisAttributes(description, x_wrap.visStr());
    vols.push_back(vol_wrap);
    pvs.push_back(pv_into);
  }
  Box box_module(xhalf_crystal+t_wrap+gap, yhalf_crystal+t_wrap+gap, lhalf_crystal+t_wrap+l_shell+gap);
  Volume vol_module("Module", box_module, material_air);
  PlacedVolume pv_wrap = vol_module.placeVolume(vols.back());
  vol_module.setVisAttributes(description, x_det.visStr());

  PlacedVolume pv_PM_right = vol_module.placeVolume(vol_PM, pos_right);
  pv_PM_right.addPhysVolID("pmt",1);
  PlacedVolume pv_PM_left = vol_module.placeVolume(vol_PM, dd4hep::Transform3D(dd4hep::RotationY(180*dd4hep::degree),pos_left));
  pv_PM_left.addPhysVolID("pmt",2);
    
  PlacedVolume pv_module = vol_hall.placeVolume(vol_module, dd4hep::Position(x_det.x(),x_det.y(),x_det.z()));
  pv_module.addPhysVolID("crystal",1);

  OpticalSurfaceManager surfMgr = description.surfaceManager();
  OpticalSurface surf_crystal   = surfMgr.opticalSurface("/world/"+det_name+"#CrystalSurface");
  BorderSurface  border_crystal2wrap = BorderSurface(description, sdet, "Crystal2Wrap", surf_crystal, pvs[0], pvs[1]);
  border_crystal2wrap.isValid();
  
  if(active_cath){
    OpticalSurface surf_cath      = surfMgr.opticalSurface("/world/"+det_name+"#CathodeSurface");
    SkinSurface    skin_cath = SkinSurface(description, sdet, "PMCathode", surf_cath, vol_cath);
    skin_cath.isValid();
  }
  else{
    OpticalSurface surf_PM       = surfMgr.opticalSurface("/world/"+det_name+"#PMSurface");
    BorderSurface  border_crystal2PM  = BorderSurface(description, sdet, "Crystal2PM", surf_PM, pvs[0], pv_glass);
    border_crystal2PM.isValid();
  }
  if(active_shell){
    OpticalSurface surf_shell     = surfMgr.opticalSurface("/world/"+det_name+"#ShellSurface");
    BorderSurface  border_shell_r = BorderSurface(description, sdet, "PM2Shell_r", surf_shell, pv_glass, pv_PM_right);
    BorderSurface  border_shell_l = BorderSurface(description, sdet, "PM2Shell_l", surf_shell, pv_glass, pv_PM_left);
    border_shell_r.isValid();
    border_shell_l.isValid();
  }

  PlacedVolume pv_hall   = description.pickMotherVolume(sdet).placeVolume(vol_hall);
  pv_hall.addPhysVolID("system",x_det.id());
  sdet.setPlacement(pv_hall);

  return sdet;
}
DECLARE_DETELEMENT(DD4hep_SingleCrystalOptical,create_detector)
