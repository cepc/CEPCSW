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

using dd4hep::rec::LayeredCalorimeterData;
static dd4hep::Ref_t create_detector(dd4hep::Detector& theDetector,
                                     xml_h e,
                                     dd4hep::SensitiveDetector sens) {

    xml_det_t x_det = e;

    std::string det_name = x_det.nameStr();
    std::string det_type = x_det.typeStr();
    MYDEBUGVAL(det_name);
    MYDEBUGVAL(det_type);
    xml_dim_t   pos    (x_det.child(_U(position)));
    xml_dim_t   dim    (x_det.child(_U(dimensions)));


    dd4hep::DetElement sdet(det_name, x_det.id());

    dd4hep::Volume motherVol = theDetector.pickMotherVolume(sdet);

    dd4hep::Material det_mat(theDetector.material("G4_BGO"));
    dd4hep::Volume det_vol(det_name+"_vol", dd4hep::Box(dim.dx(), dim.dy(), dim.dz()), det_mat);

    dd4hep::Transform3D transform(dd4hep::Rotation3D(),
                                  dd4hep::Position(pos.x(),pos.y(),pos.z()));
    dd4hep::PlacedVolume phv = motherVol.placeVolume(det_vol,transform);
    //Create caloData object to extend driver with data required for reconstruction
    LayeredCalorimeterData* caloData = new LayeredCalorimeterData ;
    caloData->layoutType = LayeredCalorimeterData::BarrelLayout ;
    caloData->inner_symmetry = 8 ;//nsides
    caloData->outer_symmetry = 8; //nsides;
    caloData->inner_phi0 = 0.;
    caloData->outer_phi0 = 0.; 
    caloData->gap0 = 0.; //FIXME
    caloData->gap1 = 0.; //FIXME
    caloData->gap2 = 0.; //FIXME
    // extent of the calorimeter in the r-z-plane [ rmin, rmax, zmin, zmax ] in mm.
    caloData->extent[0] = 10*( pos.x()-dim.dx() );// from cm to mm
    caloData->extent[1] = 10*( pos.x()+dim.dx() );// from cm to mm
    caloData->extent[2] = 0 ;
    caloData->extent[3] = 10*( pos.z()+dim.dz() );// from cm to mm
    std::cout<<"rmin="<<caloData->extent[0]<<",rmax="<<caloData->extent[1]<<",zmin="<<caloData->extent[2]<<",zmax="<<caloData->extent[3]<<std::endl;
    dd4hep::Readout readout = sens.readout();
    dd4hep::Segmentation seg = readout.segmentation();

    std::cout << "TAO: "
              << " field description: " << seg.segmentation()->fieldDescription()
              // << " grid size x/y/z: "
              // << seg.segmentation()->parameter("grid_size_x")->toString()
              // << seg.segmentation()->parameter("grid_size_x")->toString()
              // << seg.segmentation()->parameter("grid_size_x")->toString()
              << std::endl;


    double cellSize_x = ::atof( seg.segmentation()->parameter("grid_size_x")->value().c_str() ) * 10;// from cm to mm
    double cellSize_y = ::atof( seg.segmentation()->parameter("grid_size_y")->value().c_str() ) * 10;// from cm to mm
    double cellSize_z = ::atof( seg.segmentation()->parameter("grid_size_z")->value().c_str() ) * 10;// from cm to mm
    int n_layer = int(2*dim.dx()*10/cellSize_x) ; // here the calorimeter is placed in barrel, so x direaction is layer direction
    std::cout<<"cellx="<<cellSize_x<<",celly="<<cellSize_y<<",cellz="<<cellSize_z<<",dx="<<dim.dx()*10<<"mm,n_layer="<<n_layer<<std::endl;
    for(int i=1 ; i <= n_layer; i++)
    {
        LayeredCalorimeterData::Layer caloLayer ;
        caloLayer.distance = caloData->extent[0] + (i-0.5)*cellSize_x; //NEED TO START FROM ORIGIN, to mm
        caloLayer.sensitive_thickness = cellSize_x ;
        caloLayer.inner_thickness = cellSize_x ;
        caloLayer.outer_thickness = cellSize_x ;
        caloLayer.absorberThickness = cellSize_x;
        caloLayer.cellSize0 = cellSize_y;
        caloLayer.cellSize1 = cellSize_z;
        caloData->layers.push_back(caloLayer); 
    }

    if ( x_det.isSensitive() )   {
        dd4hep::SensitiveDetector sd = sens;
        det_vol.setSensitiveDetector(sens);
        sd.setType("calorimeter");
    }
    if ( x_det.hasAttr(_U(id)) )  {
        phv.addPhysVolID("system",x_det.id());
    }
    sdet.setPlacement(phv);
    sdet.addExtension< LayeredCalorimeterData >(caloData) ; 
    MYDEBUG("create_detector DONE. ");
    return sdet;
}

DECLARE_DETELEMENT(EcalMatrix, create_detector)
