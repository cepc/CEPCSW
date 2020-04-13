//====================================================================
//  Detector description implementation for Chunxiu Liu's EcalMatrix
//--------------------------------------------------------------------
//
//  Author     : Tao Lin
//               Examples from lcgeo
//                   lcgeo/detector/calorimeter/
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
    MYDEBUGVAL(det_name);
    MYDEBUGVAL(det_type);
    xml_dim_t   pos    (x_det.child(_U(position)));


    dd4hep::DetElement sdet(det_name, x_det.id());

    dd4hep::Volume motherVol = theDetector.pickMotherVolume(sdet);

    dd4hep::Material det_mat(theDetector.material("G4_BGO"));
    dd4hep::Volume det_vol(det_name+"_vol", dd4hep::Box(60, 60, 60), det_mat);

    dd4hep::Transform3D transform(dd4hep::Rotation3D(),
                                  dd4hep::Position(pos.x(),pos.y(),pos.z()));
    dd4hep::PlacedVolume phv = motherVol.placeVolume(det_vol,transform);

    if ( x_det.isSensitive() )   {
        dd4hep::SensitiveDetector sd = sens;
        det_vol.setSensitiveDetector(sens);
        sd.setType("calorimeter");
    }
    if ( x_det.hasAttr(_U(id)) )  {
        phv.addPhysVolID("system",x_det.id());
    }
    sdet.setPlacement(phv);
    
    MYDEBUG("create_detector DONE. ");
    return sdet;
}

DECLARE_DETELEMENT(EcalMatrix, create_detector)
