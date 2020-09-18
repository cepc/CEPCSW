#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetType.h"

#include "XML/Utilities.h"
#include "XML/XMLElements.h"
#include "XMLHandlerDB.h"

static dd4hep::Ref_t create_element(dd4hep::Detector& theDetector,
                                    dd4hep::xml::Handle_t e, dd4hep::SensitiveDetector sens)
{
    // typedef dd4hep::xml::DetElement      xml_det_t;
    xml_det_t x_det( e );

    std::string det_name = x_det.nameStr();

    dd4hep::DetElement sdet( det_name, x_det.id() );

    // since createPlacedEnvelope function calls following, it is turned off.
    //dd4hep::Volume mother_vol = theDetector.pickMotherVolume( sdet );


    // --- create an envelope volume and position it into the world :
    // --- this function call "addPhysVolID("system", sdet.id()) inside of it.  
    dd4hep::Volume envelope_vol = dd4hep::xml::createPlacedEnvelope( theDetector, e, sdet );

    // --- Set detector type flag
    dd4hep::xml::setDetectorTypeFlag( e, sdet );
    
    // if only the envelope part is necessary,,,
    if( theDetector.buildType() == dd4hep::BUILD_ENVELOPE ){ return sdet; }
    
    
    // Set Sensitive-Detector (SD) type.
    // Following is just a type name, but it is used later ( i.e. geosvc in Gaidi framework )
    // to match its SD
    sens.setType("tracker");     // default type !

    // Read parameter values from the xml file
    xml_comp_t layer_params( x_det.child(_U(layer)) );

    int n_layer = layer_params.attr<int>( _Unicode(nLayer) );
    int n_cell = layer_params.attr<int>( _Unicode(nCell) );
    double cell_size = layer_params.attr<double>( _Unicode(CellSize) );
    double det_half_length = layer_params.attr<double>( _Unicode(HalfLength)  );
    double min_r0 = layer_params.attr<double>( _Unicode(r0) );
    
    // Temporally define a PI : Perhaps, it is defined in other place
    //const double PI=3.14159265358979323846;   ---> M_PI is defined in the dd4hep headers
    
    double theta = 2*M_PI/(double)(n_cell);

    //dd4hep::Material foamSpacerMaterial = theDetector.material( foam_spacer_material);

    dd4hep::PlacedVolume pv;

    
    // Define & Place each volume to a mother volume  
    for( int layer_id = 0; layer_id < n_layer; layer_id++)
    {
        double r = min_r0 + layer_id*cell_size;        
        //double width = ( r - cell_size/2.0 ) * std::tan( theta/2.0 );  // original
        double width = 2*( r - cell_size/2.0 ) * std::tan( theta/2.0 );    // bug fix

	// 1. Box shape
        //dd4hep::Box box( width/2.0,  det_half_length, cell_size/2.0 );

	// 2. Trapezoid shape
	double width_lowZ = 2*( r - cell_size/2.0 ) * std::tan( theta/2.0 );
	double width_highZ = 2*( r + cell_size/2.0 ) * std::tan( theta/2.0 );
	dd4hep::Trapezoid box( width_lowZ/2.0,  width_highZ/2.0, det_half_length, det_half_length, cell_size/2.0 );

        dd4hep::Volume box_vol( dd4hep::xml::_toString(layer_id, "DCH_box_layer%i"), box, theDetector.material("TDR_gas") );


        // use visualization attributes defined in <display> 
        box_vol.setVisAttributes( theDetector, "RedVis" );

        // Set this volume as "sensitive"
        box_vol.setSensitiveDetector( sens );

        // Rotate around Z (== beam pipe) axis 
        for( int n_rotation = 0;  n_rotation < n_cell; n_rotation++)
        {
            // Add each segment to the DCH volume
            double rot_angle = n_rotation*theta;
            pv = envelope_vol.placeVolume(
                box_vol,
                dd4hep::Transform3D(
                    dd4hep::RotationZYX( 0.0, rot_angle, M_PI/2.0 ),                    
                    //dd4hep::Position( r * std::sin(rot_angle), -r * std::cos(rot_angle), -det_half_length )
                    dd4hep::Position( r * std::sin(rot_angle), -r * std::cos(rot_angle), 0.0 )                       
                    )
                );
            pv.addPhysVolID("layer", layer_id).addPhysVolID( "module", n_rotation );
        }
    }

    // Set the DCH physical volume on the mother volume
    // ----> following is done by createPlacedEnvelope function, therefore, no need to repeat.
    
    //dd4hep::PlacedVolume phv = mother_vol.placeVolume(envelope_vol);
    //phv.addPhysVolID( "system", x_det.id() );
    //sdet.setPlacement( phv );
    
    return  sdet;
}

DECLARE_DETELEMENT( DCH, create_element )


