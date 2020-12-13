//====================================================================
//  Detector description implementation of the Drift Chamber
//--------------------------------------------------------------------
//
//  Author: Tao Lin
//          Mengyao Liu
//
//====================================================================

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include "XML/XMLElements.h"
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
    // ------- Lambda functions ---- //
    auto delta_a_func = [](auto x, auto y) { return 0.5 * ( x + y ); };

    // =======================================================================
    // Parameter Definition
    // =======================================================================

    xml_det_t x_det = e;

    std::string det_name = x_det.nameStr();
    std::string det_type = x_det.typeStr();

    dd4hep::SensitiveDetector sd = sens;

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

    double epsilon = theDetector.constant<double>("Epsilon");

    // =======================================================================
    // Detector Construction
    // =======================================================================

    dd4hep::DetElement sdet(det_name, x_det.id());

    Volume envelope = dd4hep::xml::createPlacedEnvelope( theDetector,  e , sdet ) ;

    dd4hep::xml::setDetectorTypeFlag( e, sdet ) ;

    if( theDetector.buildType() == BUILD_ENVELOPE ) return sdet ;


    dd4hep::Material det_mat(theDetector.material("Air"));

    // - global
    Assembly det_vol( det_name+"_assembly" ) ;

    // - inner
    dd4hep::Tube det_inner_chamber_solid(inner_chamber_radius_min, inner_chamber_radius_max, inner_chamber_length*0.5);
    dd4hep::Volume det_inner_chamber_vol(det_name+"_inner_chamber_vol", det_inner_chamber_solid, det_mat);

    // - outer
    dd4hep::Tube det_outer_chamber_solid(outer_chamber_radius_min, outer_chamber_radius_max, outer_chamber_length*0.5);
    dd4hep::Volume det_outer_chamber_vol(det_name+"_outer_chamber_vol", det_outer_chamber_solid, det_mat);

    //Initialize the segmentation
    dd4hep::Readout readout = sd.readout();
    dd4hep::Segmentation geomseg = readout.segmentation();
    dd4hep::Segmentation* _geoSeg = &geomseg;

    auto DCHseg = dynamic_cast<dd4hep::DDSegmentation::GridDriftChamber*>(_geoSeg->segmentation());

    // - layer
      for(int layer_id = 0; layer_id < (inner_chamber_layer_number+outer_chamber_layer_number); layer_id++) {
        double rmin,rmax,offset;
        std::string layer_name;
        dd4hep::Volume* current_vol_ptr = nullptr;

        if( layer_id < inner_chamber_layer_number ) {
           current_vol_ptr = &det_inner_chamber_vol;
           rmin = inner_chamber_radius_min+(layer_id*chamber_layer_width);
           rmax = rmin+chamber_layer_width;
           layer_name = det_name+"_inner_chamber_vol"+_toString(layer_id,"_layer%d");
        }
        else {
           current_vol_ptr = &det_outer_chamber_vol;
           rmin = outer_chamber_radius_min+((layer_id-inner_chamber_layer_number)*chamber_layer_width);
           rmax = rmin+chamber_layer_width;
           layer_name = det_name+"_outer_chamber_vol"+_toString(layer_id,"_layer%d");
        }

        //Construction of drift chamber layers
        double rmid = delta_a_func(rmin,rmax);
        double ilayer_cir = 2 * M_PI * rmid;
        double ncell = ilayer_cir / chamber_layer_width;
        int ncell_layer = floor(ncell);
        int numWire = ncell_layer;
        double layer_Phi = 2*M_PI / ncell_layer;

        if(layer_id %2 ==0){ offset = 0.; }
        else { offset = 0.5 * layer_Phi; }

        DCHseg->setGeomParams(layer_id, layer_Phi, rmid, epsilon, offset);
        DCHseg->setWiresInLayer(layer_id, numWire);

        dd4hep::Tube layer_solid(rmin,rmax,chamber_length*0.5);
        dd4hep::Volume layer_vol(layer_name,layer_solid,det_mat);
        layer_vol.setAttributes(theDetector,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());

        // - wire vol
//        if(layer_id==0) {
        for(int icell=0; icell< numWire; icell++) {
            double wire_phi = (icell+0.5)*layer_Phi + offset;
            for(xml_coll_t c(x_det,_U(module)); c; ++c) {
                xml_comp_t x_module = c;
                if(x_module.id()==0) {
                  double module_rmin = x_module.rmin();
                  double module_rmax = x_module.rmax();
                  std::string module_name = x_module.nameStr();
                  dd4hep::Tube module_solid(module_rmin,module_rmax,chamber_length*0.5);
                  std::string Module_name = layer_name + module_name;
                  dd4hep::Volume module_vol(Module_name,module_solid,det_mat);
                  for(xml_coll_t l(x_module,_U(tubs)); l; ++l) {
                      xml_comp_t x_tube =l;
                      double tube_rmin = x_tube.rmin();
                      double tube_rmax = x_tube.rmax();
                      std::string tube_name = x_tube.nameStr();
                      std::string wire_name= Module_name + _toString(icell,"_%d") + tube_name;
                      dd4hep::Material tube_mat = theDetector.material(x_tube.materialStr());
                      dd4hep::Tube wire_solid(tube_rmin,tube_rmax,chamber_length*0.5);
                      dd4hep::Volume wire_vol(wire_name,wire_solid,tube_mat);
                      dd4hep::Transform3D transform_wire(dd4hep::Rotation3D(),dd4hep::Position(0.,0.,0.));
                      dd4hep::PlacedVolume wire_phy = module_vol.placeVolume(wire_vol,transform_wire);
                  }
                dd4hep::Transform3D transform_module(dd4hep::Rotation3D(),dd4hep::Position(rmid*std::cos(wire_phi),rmid*std::sin(wire_phi),0.));
                dd4hep::PlacedVolume module_phy = layer_vol.placeVolume(module_vol,transform_module);
              }
           }
         }
//       }
        dd4hep::Transform3D transform_layer(dd4hep::Rotation3D(),dd4hep::Position(0.,0.,0.));
        dd4hep::PlacedVolume layer_phy = (*current_vol_ptr).placeVolume(layer_vol, transform_layer);
        layer_phy.addPhysVolID("layer",layer_id);

        //Set drift chamber layers to sensitive detector
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
    dd4hep::PlacedVolume phv = envelope.placeVolume(det_vol,transform); 

    if ( x_det.hasAttr(_U(id)) )  {
        phv.addPhysVolID("system",x_det.id());
    }
    DetElement assDE( sdet , det_name+"_assembly" , x_det.id() )  ;
    assDE.setPlacement(phv);

    MYDEBUG("Build Detector Drift Chamber successfully.");
    return sdet;

}

DECLARE_DETELEMENT(DriftChamber, create_detector);
