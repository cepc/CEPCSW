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
    xml_det_t x_det = e;

    std::string det_name = x_det.nameStr();
    std::string det_type = x_det.typeStr();

    dd4hep::DetElement sdet(det_name, x_det.id());

    dd4hep::Volume motherVol = theDetector.pickMotherVolume(sdet);

    dd4hep::Material det_mat(theDetector.material("Air"));
    
    dd4hep::Tube det_solid(23.5, 171.6, 222.5);
    dd4hep::Volume det_vol(det_name+"_vol", det_solid, det_mat);

    dd4hep::Transform3D transform(dd4hep::Rotation3D(),
                                  dd4hep::Position(0,0,0));
    dd4hep::PlacedVolume phv = motherVol.placeVolume(det_vol,transform);

    if ( x_det.isSensitive() )   {
        dd4hep::SensitiveDetector sd = sens;
        det_vol.setSensitiveDetector(sens);
        sd.setType("tracker");
    }

    if ( x_det.hasAttr(_U(id)) )  {
        phv.addPhysVolID("system",x_det.id());
    }

    sdet.setPlacement(phv);

    MYDEBUG("Build Detector Drift Chamber successfully.");
    return sdet;

}

DECLARE_DETELEMENT(DriftChamber, create_detector);
