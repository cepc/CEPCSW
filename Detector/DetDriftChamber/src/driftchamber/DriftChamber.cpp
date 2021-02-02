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

    // - Control the number of drift chambersr
    int inner_chamber_enabled = theDetector.constant<int>("DC_inner_chamber_enabled");
    int outer_chamber_enabled = theDetector.constant<int>("DC_outer_chamber_enabled");

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

    // - wall
    double inner_chamber_inner_wall_rmin = theDetector.constant<double>("SDT_inner_chamber_inner_wall_radius_min");
    double inner_chamber_inner_wall_rmax = theDetector.constant<double>("SDT_inner_chamber_inner_wall_radius_max");
    double inner_chamber_outer_wall_rmin = theDetector.constant<double>("SDT_inner_chamber_outer_wall_radius_min");
    double inner_chamber_outer_wall_rmax = theDetector.constant<double>("SDT_inner_chamber_outer_wall_radius_max");
    double outer_chamber_outer_wall_rmin = theDetector.constant<double>("SDT_outer_chamber_outer_wall_radius_min");
    double outer_chamber_outer_wall_rmax = theDetector.constant<double>("SDT_outer_chamber_outer_wall_radius_max");
    double outer_chamber_inner_wall_rmin = theDetector.constant<double>("SDT_outer_chamber_inner_wall_radius_min");
    double outer_chamber_inner_wall_rmax = theDetector.constant<double>("SDT_outer_chamber_inner_wall_radius_max");

    dd4hep::Material wall_mat(theDetector.material("CarbonFiber"));

    double wall_rmin[4] = {inner_chamber_inner_wall_rmin,inner_chamber_outer_wall_rmin,outer_chamber_inner_wall_rmin,outer_chamber_outer_wall_rmin};
    double wall_rmax[4] = {inner_chamber_inner_wall_rmax,inner_chamber_outer_wall_rmax,outer_chamber_inner_wall_rmax,outer_chamber_outer_wall_rmax};

   // - wire
    dd4hep::Volume module_vol;
    dd4hep::Volume Module_vol;
    for(xml_coll_t c(x_det,_U(module)); c; ++c) {
      xml_comp_t x_module = c;
      double  module_rmin = x_module.rmin();
      double  module_rmax = x_module.rmax();
      std::string module_name = x_module.nameStr();
      dd4hep::Tube module_solid(module_rmin,module_rmax,chamber_length*0.5);
      if(x_module.id()==0) {
         module_vol = dd4hep::Volume(module_name,module_solid,det_mat);
         module_vol.setVisAttributes(theDetector.visAttributes(x_module.visStr()));
      } else {
         Module_vol = dd4hep::Volume(module_name,module_solid,det_mat);
         Module_vol.setVisAttributes(theDetector.visAttributes(x_module.visStr()));
      }

      for(xml_coll_t l(x_module,_U(tubs)); l; ++l) {
         xml_comp_t x_tube =l;
         double tube_rmin = x_tube.rmin();
         double tube_rmax = x_tube.rmax();
         std::string tube_name = x_tube.nameStr();
         std::string wire_name= module_name + tube_name;
         dd4hep::Material tube_mat = theDetector.material(x_tube.materialStr());
         dd4hep::Tube wire_solid(tube_rmin,tube_rmax,chamber_length*0.5);
         dd4hep::Volume wire_vol(wire_name,wire_solid,tube_mat);
         dd4hep::Transform3D transform_wire(dd4hep::Rotation3D(),dd4hep::Position(0.,0.,0.));
         dd4hep::PlacedVolume wire_phy;
         if(x_module.id()==0) {
            wire_phy = module_vol.placeVolume(wire_vol,transform_wire);
         } else {
            wire_phy = Module_vol.placeVolume(wire_vol,transform_wire);
         }
      }
  }

    // End cap
    double Endcap_rmin = theDetector.constant<double>("DC_Endcap_rmin");
    double Endcap_rmax = theDetector.constant<double>("DC_Endcap_rmax");
    double Endcap_z = theDetector.constant<double>("DC_Endcap_dz");
    dd4hep::Tube det_Endcap_solid(Endcap_rmin,Endcap_rmax,0.5*Endcap_z);
    dd4hep::Volume det_Endcap_vol(det_name+"Endcap",det_Endcap_solid,det_mat);
    det_Endcap_vol.setVisAttributes(theDetector,"YellowVis");

    //Initialize the segmentation
    dd4hep::Readout readout = sd.readout();
    dd4hep::Segmentation geomseg = readout.segmentation();
    dd4hep::Segmentation* _geoSeg = &geomseg;

    auto DCHseg = dynamic_cast<dd4hep::DDSegmentation::GridDriftChamber*>(_geoSeg->segmentation());

    // - layer
    for(int layer_id = 0; layer_id < (inner_chamber_layer_number+outer_chamber_layer_number); layer_id++) {
        double rmin,rmax,offset=0;
        std::string layer_name;
        dd4hep::Volume* current_vol_ptr = nullptr;
        dd4hep::Material layer_mat(theDetector.material("GasHe_90Isob_10"));
        if(inner_chamber_enabled && (layer_id < inner_chamber_layer_number)) {
           current_vol_ptr = &det_inner_chamber_vol;
           rmin = inner_chamber_radius_min+(layer_id*chamber_layer_width);
           rmax = rmin+chamber_layer_width;
           layer_name = det_name+"_inner_chamber_vol"+_toString(layer_id,"_layer%d");
        } else if(outer_chamber_enabled && (layer_id > (inner_chamber_layer_number-1))) {
           current_vol_ptr = &det_outer_chamber_vol;
           rmin = outer_chamber_radius_min+((layer_id-inner_chamber_layer_number)*chamber_layer_width);
           rmax = rmin+chamber_layer_width;
           layer_name = det_name+"_outer_chamber_vol"+_toString(layer_id,"_layer%d");
        } else continue;

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
        dd4hep::Volume layer_vol(layer_name,layer_solid,layer_mat);
        layer_vol.setAttributes(theDetector,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());

        // - wire vol
        //phi <-------------------> -phi
        //    |   F8    F7   F6   F5|  Only on the outermost layer.
        //    |                     |
        //    |         S         F4|
        //    |                     |
        //    |   F0    F1   F2   F3|
        //    -----------------------
//     if(layer_id == -1) {
        for(int icell=0; icell< numWire; icell++) {
            double wire_phi = (icell+0.5)*layer_Phi + offset;
            // - signal wire
            dd4hep::Transform3D transform_module(dd4hep::Rotation3D(),dd4hep::Position(rmid*std::cos(wire_phi),rmid*std::sin(wire_phi),0.));
            dd4hep::PlacedVolume module_phy = layer_vol.placeVolume(module_vol,transform_module);
           // - Field wire
            dd4hep::PlacedVolume Module_phy;
            double radius[9] = {rmid-chamber_layer_width*0.5,rmid-chamber_layer_width*0.5,rmid-chamber_layer_width*0.5,rmid-chamber_layer_width*0.5,rmid,rmid+chamber_layer_width*0.5,rmid+chamber_layer_width*0.5,rmid+chamber_layer_width*0.5,rmid+chamber_layer_width*0.5};
            double phi[9] = {wire_phi+layer_Phi*0.25,wire_phi,wire_phi-layer_Phi*0.25,wire_phi-layer_Phi*0.5,wire_phi-layer_Phi*0.5,wire_phi-layer_Phi*0.5,wire_phi-layer_Phi*0.25,wire_phi,wire_phi+layer_Phi*0.25};
            int num = 5;
            if(layer_id==(inner_chamber_layer_number-1)||layer_id==(outer_chamber_layer_number-1)) { num = 9; }
            for(int i=0; i<num ; i++) {
                dd4hep::Position tr3D = Position(radius[i]*std::cos(phi[i]),radius[i]*std::sin(phi[i]),0.);
                dd4hep::Transform3D transform_Module(dd4hep::Rotation3D(),tr3D);
                Module_phy = layer_vol.placeVolume(Module_vol,transform_Module);
            }
        }
//  }

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
    if(inner_chamber_enabled) {
         dd4hep::PlacedVolume det_inner_chamber_phy = det_vol.placeVolume(det_inner_chamber_vol,
                 transform_inner_chamber);

         det_inner_chamber_phy.addPhysVolID("chamber", 0);
    }
    // outer
    dd4hep::Transform3D transform_outer_chamber(dd4hep::Rotation3D(),
            dd4hep::Position(0,0,0));
    if(outer_chamber_enabled) {
       dd4hep::PlacedVolume det_outer_chamber_phy = det_vol.placeVolume(det_outer_chamber_vol,
               transform_inner_chamber);

       det_outer_chamber_phy.addPhysVolID("chamber", 1);
    }
    // - place in world
    dd4hep::Transform3D transform(dd4hep::Rotation3D(),
            dd4hep::Position(0,0,0));
    dd4hep::PlacedVolume phv = envelope.placeVolume(det_vol,transform);
    // - place wall
//    dd4hep::PlacedVolume wall_phy;
    for(int i=0; i<4; i++) {
       dd4hep::Tube wall_solid(wall_rmin[i],wall_rmax[i],chamber_length*0.5);
       dd4hep::Volume wall_vol(det_name+"_wall_vol",wall_solid,wall_mat);
       wall_vol.setVisAttributes(theDetector,"VisibleGreen");
       dd4hep::PlacedVolume wall_phy = envelope.placeVolume(wall_vol,transform);
    }

    // - place Endcap
    double endcap_pos[2] = {chamber_length*0.5+Endcap_z*0.5,-chamber_length*0.5-Endcap_z*0.5};
    dd4hep::PlacedVolume endcap_phy;
    for(int i=0; i<2; i++) {
        dd4hep::Transform3D Endcap_transform(dd4hep::Rotation3D(),dd4hep::Position(0,0,endcap_pos[i]));
        endcap_phy = envelope.placeVolume(det_Endcap_vol,Endcap_transform);
   }

    if ( x_det.hasAttr(_U(id)) )  {
        phv.addPhysVolID("system",x_det.id());
    }
    DetElement assDE( sdet , det_name+"_assembly" , x_det.id() )  ;
    assDE.setPlacement(phv);

    MYDEBUG("Build Detector Drift Chamber successfully.");
    return sdet;

}

DECLARE_DETELEMENT(DriftChamber, create_detector);
