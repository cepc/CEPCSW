#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include "XML/XMLElements.h"
#include "DDRec/DetectorData.h"
#include "DDSegmentation/Segmentation.h"
#include "DetSegmentation/GridDriftChamber.h"

#include <iostream>
#include <chrono>

using namespace dd4hep;
using namespace dd4hep::detail;
using namespace dd4hep::rec ;

#define MYDEBUG(x) std::cout << __FILE__ << ":" << __LINE__ << ": " << x << std::endl;
#define MYDEBUGVAL(x) std::cout << __FILE__ << ":" << __LINE__ << ": " << #x << ": " << x << std::endl;

static dd4hep::Ref_t create_detector(dd4hep::Detector& theDetector,
        xml_h e,
        dd4hep::SensitiveDetector sens) {
    // ------- Lambda functions ---- //
    auto delta_b_func = [](auto x, auto y) { return 2 * std::sqrt((x + y) * (x - y)); };
    auto epsilon_func = [](auto delta_a, auto L) { return std::atan(delta_a / L); };

    // =======================================================================
    // Parameter Definition
    // =======================================================================

    auto start = std::chrono::steady_clock::now();

    xml_det_t x_det = e;

    xml_coll_t c(x_det,_U(chamber));
    xml_comp_t x_chamber = c;

    xml_coll_t cc(x_det,_U(side));
    xml_comp_t x_side = cc;

    std::string det_name = x_det.nameStr();
    std::string det_type = x_det.typeStr();

    dd4hep::SensitiveDetector sd = sens;

    // - global
    double chamber_half_length = theDetector.constant<double>("Gas_half_length");
    double chamber_length = chamber_half_length*2;

    double alpha = theDetector.constant<double>("Alpha");
    double Pi = 0;
    if(0!=alpha) Pi = 0.5*M_PI;

    double safe_distance = theDetector.constant<double>("DC_safe_distance");

    double DC_inner_wall_thickness = theDetector.constant<double>("DC_inner_wall_thickness");
    double DC_outer_wall_thickness = theDetector.constant<double>("DC_outer_wall_thickness");

    // - chamber gas
    double DC_rend = theDetector.constant<double>("DC_rend");

    double chamber_radius_min = theDetector.constant<double>("Gas_radius_min");
    double chamber_radius_max = DC_rend-safe_distance-DC_outer_wall_thickness;
    double layer_radius_max = chamber_radius_max*std::cos(alpha);
    int chamberID = x_chamber.id();

    // - layer
    int DC_layer_number = theDetector.constant<int>("DC_layer_number");

    double layer_width = (layer_radius_max-chamber_radius_min)/DC_layer_number;

    double cell_width = theDetector.constant<double>("DC_cell_width");

    int DC_construct_wire = theDetector.constant<int>("DC_construct_wire");

    // =======================================================================
    // Detector Construction
    // =======================================================================

    dd4hep::DetElement sdet(det_name, x_det.id());

    Volume envelope = dd4hep::xml::createPlacedEnvelope( theDetector,  e , sdet ) ;

    dd4hep::xml::setDetectorTypeFlag( e, sdet ) ;

    if( theDetector.buildType() == BUILD_ENVELOPE ) return sdet ;

    dd4hep::Material det_mat(theDetector.material(x_det.materialStr()));
    dd4hep::Material chamber_mat = theDetector.material(x_chamber.materialStr());

    /// - global
    Assembly det_vol( det_name+"_assembly" ) ;

    /// - chamber volume
    dd4hep::Tube det_chamber_solid(chamber_radius_min, chamber_radius_max, chamber_half_length);
    dd4hep::Volume det_chamber_vol(det_name+"_chamber_vol", det_chamber_solid, chamber_mat);

    /// - wall
    double chamber_inner_wall_rmin = theDetector.constant<double>("DC_inner_wall_radius_min");
    double chamber_inner_wall_rmax = theDetector.constant<double>("DC_inner_wall_radius_max");
    double chamber_outer_wall_rmin = chamber_radius_max+safe_distance;
    double chamber_outer_wall_rmax = chamber_outer_wall_rmin+DC_outer_wall_thickness;

    dd4hep::Material wall_mat(theDetector.material(x_side.materialStr()));

    double wall_rmin[2] = {chamber_inner_wall_rmin, chamber_outer_wall_rmin};
    double wall_rmax[2] = {chamber_inner_wall_rmax, chamber_outer_wall_rmax};

    /// - construct wires
    dd4hep::Volume signalWireVolume;
    dd4hep::Volume fieldWireVolume;
    for(xml_coll_t dcModule(x_det,_U(module)); dcModule; ++dcModule) {
        xml_comp_t x_module = dcModule;
        std::string module_name = x_module.nameStr();
        dd4hep::Tube module_solid(x_module.rmin(),x_module.rmax(),chamber_half_length);
        if(0==x_module.id()) {
            signalWireVolume = dd4hep::Volume(module_name,module_solid,det_mat);
            signalWireVolume.setVisAttributes(theDetector.visAttributes(x_module.visStr()));
        } else {
            fieldWireVolume = dd4hep::Volume(module_name,module_solid,det_mat);
            fieldWireVolume.setVisAttributes(theDetector.visAttributes(x_module.visStr()));
        }

        /// construct wire tubes
        for(xml_coll_t l(x_module,_U(tubs)); l; ++l) {
            xml_comp_t x_tube =l;
            std::string wire_name= module_name + x_tube.nameStr();
            dd4hep::Material tube_mat = theDetector.material(x_tube.materialStr());
            dd4hep::Tube wire_solid(x_tube.rmin(),x_tube.rmax(),chamber_half_length);
            dd4hep::Volume wire_vol(wire_name,wire_solid,tube_mat);
            dd4hep::Transform3D transform_wire(dd4hep::Rotation3D(),dd4hep::Position(0.,0.,0.));
            if(0==x_module.id()) {
                signalWireVolume.placeVolume(wire_vol,transform_wire);
            } else {
                fieldWireVolume.placeVolume(wire_vol,transform_wire);
            }
        }//end of construct tubes
    }//end of construct wire

    /// construct End cap
    double endcap_rmin = theDetector.constant<double>("DC_Endcap_rmin");
    double endcap_rmax = theDetector.constant<double>("DC_Endcap_rmax");
    double endcap_z = theDetector.constant<double>("DC_Endcap_dz");
    dd4hep::Tube det_Endcap_solid(endcap_rmin,endcap_rmax,0.5*endcap_z);
    dd4hep::Volume det_Endcap_vol(det_name+"Endcap",det_Endcap_solid,det_mat);
    det_Endcap_vol.setVisAttributes(theDetector,"YellowVis");

    ///Initialize the segmentation
    auto segmentationDC =
        dynamic_cast<dd4hep::DDSegmentation::GridDriftChamber*>
        (sd.readout().segmentation().segmentation());

    /// - construct layers
    int chamber_id = 0;
    for(int layer_id = 0; layer_id < DC_layer_number; layer_id++) {
        dd4hep::Volume* current_vol_ptr = nullptr;
        double layer_rmin = chamber_radius_min+(layer_id*layer_width);
        double layer_rmax = layer_rmin+layer_width;
        double rmid_zZero = (layer_rmin+layer_rmax)/2.;  // z=0
        double rmid_zEnd = rmid_zZero/std::cos(alpha/2);  //  z=endcap
        int nCell = floor((2. * M_PI * rmid_zZero) / layer_width);
        int nWire = nCell;
        if(!DC_construct_wire) nWire =0;
        double cell_phi = 2*M_PI / nCell;
        double offset=0;//phi offset of first cell in each layer
        double sign_eps = 1;// setero angle sign
        if(0==(layer_id%2)) {
            offset = 0.;
            sign_eps = -1;
        } else {
            offset = 0.5 * cell_phi;
        }
        double epsilon = sign_eps*std::atan(rmid_zZero * std::tan(alpha / 2.0) / chamber_half_length);

        segmentationDC->setGeomParams(chamber_id, layer_id, cell_phi, rmid_zEnd , epsilon, offset);
        segmentationDC->setWiresInLayer(chamber_id, layer_id, nCell);

        double r_in_test = rmid_zZero*std::cos(alpha / 2.0);

        double r_in0 = rmid_zZero - layer_width / 2.0;//
        double r_in = r_in0 / std::cos(alpha / 2.0);//layer min at z=half length

        double r_out0 = rmid_zZero + layer_width / 2.0;
        double r_out = r_out0 / std::cos(alpha / 2.0);

        double delta_a_in = delta_b_func(r_in, r_in0);
        double delta_a_out = delta_b_func(r_out, r_out0);

        double eps_in = epsilon_func(delta_a_in, chamber_length );
        double eps_out = epsilon_func(delta_a_out, chamber_length );

        /// create hyper layer volume
        dd4hep::Hyperboloid layer_vol_solid(r_in0, eps_in, r_out0, eps_out, chamber_half_length);
        dd4hep::Volume layer_vol(det_name+"_layer_vol",layer_vol_solid,det_mat);
        current_vol_ptr = &layer_vol;

        if ( x_det.isSensitive() )   {
            layer_vol.setRegion(theDetector,x_det.regionStr());
            layer_vol.setLimitSet(theDetector,x_det.limitsStr());
            layer_vol.setSensitiveDetector(sens);
            sd.setType("tracker");
        }

        // - create wire volume
        //phi <---------------> -phi
        //    |    F4     F3|  Only on the outermost layer.
        //    |             |
        //    |    S      F2|
        //    |             |
        //    |    F0     F1|
        //--------------------
        // loop over cells
        for(int icell=0; icell<nWire; icell++) {

            double wire_phi = icell*cell_phi + offset;

            // - signal wire
            dd4hep::RotationZ rz(wire_phi+Pi);
            dd4hep::RotationY ry(epsilon);
            dd4hep::Position tr3D = Position(rmid_zZero*std::cos(wire_phi),rmid_zZero*std::sin(wire_phi),0.);
            dd4hep::Transform3D transform_signal_wire(rz*ry,tr3D);

            (*current_vol_ptr).placeVolume(signalWireVolume,transform_signal_wire);
            // - Field wire
            double radius[5] = {rmid_zZero-layer_width*0.5+safe_distance,rmid_zZero-layer_width*0.5+safe_distance,rmid_zZero,rmid_zZero+layer_width*0.5-safe_distance,rmid_zZero+layer_width*0.5-safe_distance};
            double phi[5] = {wire_phi,wire_phi-cell_phi*0.5,wire_phi-cell_phi*0.5,wire_phi-cell_phi*0.5,wire_phi};
            int num = 3;
            if(layer_id==(DC_layer_number-1)) {
                num = 5;
            }
            for(int i=0; i<num ; i++) {
                dd4hep::RotationZ rz_field(phi[i]+Pi);
                dd4hep::RotationY ry_field(epsilon);
                dd4hep::Position tr3D_field = Position(radius[i]*std::cos(phi[i]),radius[i]*std::sin(phi[i]),0.);

                dd4hep::Transform3D transform_field_wire(rz_field*ry_field,tr3D_field);

                (*current_vol_ptr).placeVolume(fieldWireVolume,transform_field_wire);
            }
        }//end of loop over cell
        dd4hep::Transform3D transform_layer(dd4hep::Rotation3D(),
                dd4hep::Position(0,0,0));
        dd4hep::PlacedVolume layer_phy = det_chamber_vol.placeVolume(layer_vol,transform_layer);
        layer_phy.addPhysVolID("layer", layer_id);
    }//end of loop over layers


    // - place in det
    // - chamber
    dd4hep::Transform3D transform_chamber(dd4hep::Rotation3D(),
            dd4hep::Position(0,0,0));
    dd4hep::PlacedVolume det_chamber_phy = det_vol.placeVolume(det_chamber_vol,
            transform_chamber);

    det_chamber_phy.addPhysVolID("chamber", chamberID);

    // - place in world
    dd4hep::Transform3D transform(dd4hep::Rotation3D(),
            dd4hep::Position(0,0,0));
    dd4hep::PlacedVolume phv = envelope.placeVolume(det_vol,transform);
    // - place wall
    dd4hep::PlacedVolume wall_phy;
    for(int i=0; i<2; i++) {
        dd4hep::Tube wall_solid(wall_rmin[i],wall_rmax[i],chamber_half_length);
        dd4hep::Volume wall_vol(det_name+"_wall_vol",wall_solid,wall_mat);
        wall_vol.setVisAttributes(theDetector,"VisibleGreen");
        wall_phy = envelope.placeVolume(wall_vol,transform);
    }

    // - place Endcap
    double endcap_pos[2] = {chamber_half_length+endcap_z*0.5,-chamber_half_length-endcap_z*0.5};
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
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n"; 
    return sdet;

}

DECLARE_DETELEMENT(DriftChamber, create_detector)
