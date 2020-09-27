//====================================================================
//  Detector description implementation of the Drift Chamber
//--------------------------------------------------------------------
//
//  Author: Tao Lin
//          Mengyao Liu
//
//====================================================================

#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"
#include "DDSegmentation/Segmentation.h"
#include "DetSegmentation/GridDriftChamber.h"

using namespace dd4hep;
using namespace dd4hep::detail;
using namespace dd4hep::rec ;

#define MYDEBUG(x) std::cout << __FILE__ << ":" << __LINE__ << ": " << x << std::endl;
#define MYDEBUGVAL(x) std::cout << __FILE__ << ":" << __LINE__ << ": " << #x << ": " << x << std::endl;

static dd4hep::Ref_t create_detector(dd4hep::Detector& theDetector,
        xml_h e,
        dd4hep::SensitiveDetector sens) {
    // =======================================================================
    // Parameter Definition
    // =======================================================================

    xml_det_t x_det = e;

    std::string det_name = x_det.nameStr();
    std::string det_type = x_det.typeStr();

    // - global
    double chamber_radius_min = theDetector.constant<double>("SDT_radius_min");
    double chamber_radius_max = theDetector.constant<double>("SDT_radius_max");
    double chamber_length     = theDetector.constant<double>("SDT_length");

    // - inner chamber
    double inner_chamber_radius_min = theDetector.constant<double>("SDT_inner_chamber_radius_min");
    double inner_chamber_radius_max = theDetector.constant<double>("SDT_inner_chamber_radius_max");
    double inner_chamber_length     = theDetector.constant<double>("SDT_inner_chamber_length");

    // - outer chamber
    double outer_chamber_radius_min = theDetector.constant<double>("SDT_outer_chamber_radius_min");
    double outer_chamber_radius_max = theDetector.constant<double>("SDT_outer_chamber_radius_max");
    double outer_chamber_length     = theDetector.constant<double>("SDT_outer_chamber_length");

    // - layer
    int inner_chamber_layer_number = theDetector.constant<int>("SDT_inner_chamber_layer_number");
    int outer_chamber_layer_number = theDetector.constant<int>("SDT_outer_chamber_layer_number");
    double chamber_layer_width  = theDetector.constant<double>("SDT_chamber_layer_width");

    // =======================================================================
    // Detector Construction
    // =======================================================================

    dd4hep::DetElement sdet(det_name, x_det.id());

    dd4hep::Volume motherVol = theDetector.pickMotherVolume(sdet);

    dd4hep::Material det_mat(theDetector.material("Air"));

    // - global
    dd4hep::Tube det_solid(chamber_radius_min, chamber_radius_max, chamber_length*0.5);
    dd4hep::Volume det_vol(det_name+"_vol", det_solid, det_mat);

    // - inner
    dd4hep::Tube det_inner_chamber_solid(inner_chamber_radius_min, inner_chamber_radius_max, inner_chamber_length*0.5);
    dd4hep::Volume det_inner_chamber_vol(det_name+"_inner_chamber_vol", det_inner_chamber_solid, det_mat);

    // - outer
    dd4hep::Tube det_outer_chamber_solid(outer_chamber_radius_min, outer_chamber_radius_max, outer_chamber_length*0.5);
    dd4hep::Volume det_outer_chamber_vol(det_name+"_outer_chamber_vol", det_outer_chamber_solid, det_mat);


    // - layer
    for(int layer_id=0; layer_id<(inner_chamber_layer_number+outer_chamber_layer_number-1);layer_id++) {
      double rmin,rmax;
      std::string layer_name;
      dd4hep::Volume* current_vol_ptr = nullptr;

      if(layer_id<inner_chamber_layer_number){
        current_vol_ptr = &det_inner_chamber_vol;
        rmin = inner_chamber_radius_min+(layer_id*chamber_layer_width);
        rmax = rmin+chamber_layer_width;
        layer_name = det_name+"_inner_chamber_vol"+_toString(layer_id,"_layer%d");
      }else{
        current_vol_ptr = &det_outer_chamber_vol;
        rmin = outer_chamber_radius_min+((layer_id-inner_chamber_layer_number)*chamber_layer_width);
        rmax = rmin+chamber_layer_width;
        layer_name = det_name+"_outer_chamber_vol"+_toString(layer_id,"_layer%d");
      }
      /// Construction of drift chamber layers
      dd4hep::Tube layer_solid(rmin,rmax,chamber_length*0.5);
      dd4hep::Volume layer_vol(layer_name,layer_solid,det_mat);
      dd4hep::Transform3D transform_layer(dd4hep::Rotation3D(),dd4hep::Position(0.,0.,0.));
      dd4hep::PlacedVolume layer_phy = (*current_vol_ptr).placeVolume(layer_vol, transform_layer);
      layer_phy.addPhysVolID("layer",layer_id);

      /// Set drift chamber layers to sensitive detector
      dd4hep::SensitiveDetector sd = sens;
      layer_vol.setSensitiveDetector(sens);
      sd.setType("tracker");
   }

    // - place in det
    // inner
     dd4hep::Transform3D transform_inner_chamber(dd4hep::Rotation3D(),
            dd4hep::Position(0,0,0));
     dd4hep::PlacedVolume det_inner_chamber_phy = det_vol.placeVolume(det_inner_chamber_vol,
            transform_inner_chamber);
     det_inner_chamber_phy.addPhysVolID("chamber", 0);

    // outer
    dd4hep::Transform3D transform_outer_chamber(dd4hep::Rotation3D(),
            dd4hep::Position(0,0,0));
    dd4hep::PlacedVolume det_outer_chamber_phy = det_vol.placeVolume(det_outer_chamber_vol,
            transform_inner_chamber);
    det_outer_chamber_phy.addPhysVolID("chamber", 1);

    // - place in world
    dd4hep::Transform3D transform(dd4hep::Rotation3D(),
            dd4hep::Position(0,0,0));
    dd4hep::PlacedVolume phv = motherVol.placeVolume(det_vol,transform);


    if ( x_det.hasAttr(_U(id)) )  {
        phv.addPhysVolID("system",x_det.id());
    }

    sdet.setPlacement(phv);

    MYDEBUG("Build Detector Drift Chamber successfully.");
    return sdet;

}

DECLARE_DETELEMENT(DriftChamber, create_detector);
