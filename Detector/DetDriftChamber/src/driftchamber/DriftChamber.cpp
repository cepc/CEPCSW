//====================================================================
//  Detector description implementation of the Drift Chamber
//--------------------------------------------------------------------
//
//  Author: Tao Lin
//
//====================================================================

#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"
#include "DDSegmentation/Segmentation.h"

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

    // =======================================================================
    // Detector Construction
    // =======================================================================

    dd4hep::DetElement sdet(det_name, x_det.id());

    dd4hep::Volume motherVol = theDetector.pickMotherVolume(sdet);

    dd4hep::Material det_mat(theDetector.material("Air"));
    
    // - global 
    dd4hep::Tube det_solid(chamber_radius_min, chamber_radius_max, chamber_length*0.5);
    dd4hep::Volume det_vol(det_name+"_vol", det_solid, det_mat);

    if ( x_det.isSensitive() )   {
        dd4hep::SensitiveDetector sd = sens;
        det_vol.setSensitiveDetector(sens);
        sd.setType("tracker");
    }

    // - inner
    dd4hep::Tube det_inner_chamber_solid(inner_chamber_radius_min, inner_chamber_radius_max, inner_chamber_length*0.5);
    dd4hep::Volume det_inner_chamber_vol(det_name+"_inner_chamber_vol", det_inner_chamber_solid, det_mat);

    dd4hep::Transform3D transform_inner_chamber(dd4hep::Rotation3D(),
                                                dd4hep::Position(0,0,0));
    dd4hep::PlacedVolume det_inner_chamber_phy = det_vol.placeVolume(det_inner_chamber_vol,
                                                                     transform_inner_chamber);

    // - outer
    dd4hep::Tube det_outer_chamber_solid(outer_chamber_radius_min, outer_chamber_radius_max, outer_chamber_length*0.5);
    dd4hep::Volume det_outer_chamber_vol(det_name+"_outer_chamber_vol", det_outer_chamber_solid, det_mat);

    dd4hep::Transform3D transform_outer_chamber(dd4hep::Rotation3D(),
                                                dd4hep::Position(0,0,0));
    dd4hep::PlacedVolume det_outer_chamber_phy = det_vol.placeVolume(det_outer_chamber_vol,
                                                                     transform_outer_chamber);



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
