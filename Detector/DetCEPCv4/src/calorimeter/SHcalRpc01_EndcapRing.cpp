//====================================================================
//  SHcalRpc01 - Implementation from ILCSoft's Mokka version                              
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetType.h"

#include "DDRec/Surface.h"
#include "DDRec/DetectorData.h"
#include "XML/Utilities.h"

using namespace std;

using dd4hep::Ref_t;
using dd4hep::BUILD_ENVELOPE;
using dd4hep::DetElement;
using dd4hep::Detector;
using dd4hep::SensitiveDetector;
using dd4hep::Segmentation;
using dd4hep::Readout;
using dd4hep::Material;
using dd4hep::Volume;
using dd4hep::PlacedVolume;
using dd4hep::Position;
using dd4hep::RotationZYX;
using dd4hep::Transform3D;
using dd4hep::Box;
using dd4hep::Tube;
using dd4hep::PolyhedraRegular;
using dd4hep::SubtractionSolid;
using dd4hep::IntersectionSolid;
using dd4hep::_toString;
using dd4hep::pi;
using dd4hep::rec::LayeredCalorimeterData;

/** Construction of SHcalRpc01 detector, ported from Mokka driver SHcalRpc01.cc
 *
 *  Mokka History:
 * - first implementation from ILCSoft
 * - http://cepcgit.ihep.ac.cn/cepcsoft/MokkaC
 */
static Ref_t create_detector(Detector& theDetector, xml_h element, SensitiveDetector sens)  {
  cout << "------------------------------" << endl;
  cout << "creating SHcalRpc01_EndcapRing" << endl;
  cout << "------------------------------" << endl;

  xml_det_t    x_det = element;
  string       name  = x_det.nameStr();
  
  int          det_id    = x_det.id();
  DetElement   det(name, det_id) ;
  
  xml_comp_t    x_staves          = x_det.staves();
  string   Hcal_radiator_material = x_staves.materialStr();
  Material      stavesMaterial    = theDetector.material(Hcal_radiator_material);
  Material      air               = theDetector.air();

  Volume envelope = dd4hep::xml::createPlacedEnvelope( theDetector,  element , det ) ;

  dd4hep::xml::setDetectorTypeFlag( element, det ) ;

  if( theDetector.buildType() == BUILD_ENVELOPE ) return det ;

  sens.setType("calorimeter");

  DetElement module(det,"module0",det_id);
  DetElement layer(module, "stave_layer", det_id);
  DetElement slice(layer, "slice", det_id);

  Readout readout = sens.readout();
  Segmentation seg = readout.segmentation();
  
  std::vector<double> cellSizeVector = seg.segmentation()->cellDimensions(0);
  double cell_sizeX      = cellSizeVector[0];
  double cell_sizeY      = cellSizeVector[1];

  double Hcal_outer_radius            = theDetector.constant<double>("Hcal_outer_radius");
  int    Hcal_endcap_outer_symmetry   = theDetector.constant<int>("Hcal_endcap_outer_symmetry");
  double Hcal_stave_gaps              = theDetector.constant<double>("Hcal_stave_gaps");
  int    Hcal_nlayers                 = theDetector.constant<int>("Hcal_endcap_nlayers");
  double Hcal_start_z                 = theDetector.constant<double>("Hcal_endcap_zmin");
  double Hcal_lateral_plate_thickness = theDetector.constant<double>("Hcal_lateral_structure_thickness");
  double Ecal_endcap_zmin             = theDetector.constant<double>("Ecal_endcap_zmin");
  double Hcal_endcap_ecal_gap         = theDetector.constant<double>("Hcal_endcap_ecal_gap");
  double Ecal_endcap_outer_radius     = theDetector.constant<double>("EcalEndcap_outer_radius");
  double Hcal_radial_ring_inner_gap   = theDetector.constant<double>("Hcal_radial_ring_inner_gap");

  xml_coll_t c(x_det,_U(layer));
  xml_comp_t x_layer = c;

  double Hcal_radiator_thickness = 0;
  double layerThickness = 0.0;
  for(xml_coll_t k(x_layer,_U(slice)); k; ++k)  {
    xml_comp_t x_slice = k;
    layerThickness += x_slice.thickness();
    if(x_slice.materialStr()==Hcal_radiator_material) Hcal_radiator_thickness = x_slice.thickness();
  }
  cout << " layer_thickness (from slices) = " << layerThickness << " and radiator_thickness = " << Hcal_radiator_thickness << endl;
  double Hcal_chamber_thickness = layerThickness - Hcal_radiator_thickness;

  int numSide = Hcal_endcap_outer_symmetry;
  double Hcal_endcap_rmax = Hcal_outer_radius * cos(pi/numSide);

  LayeredCalorimeterData* caloData = new LayeredCalorimeterData ;
  caloData->layoutType = LayeredCalorimeterData::EndcapLayout ;
  caloData->inner_symmetry = numSide;
  caloData->outer_symmetry = numSide;
  caloData->phi0           = 0;
  
  double start_z        = Ecal_endcap_zmin;
  double SpaceForLayers = Hcal_start_z - Hcal_endcap_ecal_gap - Ecal_endcap_zmin - 2*Hcal_lateral_plate_thickness;
  int MaxNumberOfLayers = (int)(SpaceForLayers / (Hcal_chamber_thickness + Hcal_radiator_thickness));
  double stop_z         = start_z + MaxNumberOfLayers*layerThickness + 2*Hcal_lateral_plate_thickness;

  double pRMax = Hcal_endcap_rmax;
  double pDz   = (stop_z - start_z)/2.;
  double pRMin = Ecal_endcap_outer_radius + Hcal_radial_ring_inner_gap;

  cout << "Rings will have " << MaxNumberOfLayers << " layers." << endl;
  cout << " Z: " << start_z << " -> " << stop_z << endl;
  cout << " R: " << pRMin   << " -> " << pRMax  << endl;

  caloData->extent[0] = pRMin;
  caloData->extent[1] = pRMax;
  caloData->extent[2] = start_z;
  caloData->extent[3] = stop_z;
  
  PolyhedraRegular EndCapRingSolidPoly(numSide, -pi/numSide, pRMin, pRMax, 2*pDz);
  Volume EndCapRingLogical(name+"_radiator", EndCapRingSolidPoly, stavesMaterial);
  EndCapRingLogical.setAttributes(theDetector,x_staves.regionStr(),x_staves.limitsStr(),x_staves.visStr());

  int number_of_chambers = Hcal_nlayers;
  if(MaxNumberOfLayers < number_of_chambers) number_of_chambers = MaxNumberOfLayers;
  
  double rInner = pRMin + Hcal_lateral_plate_thickness;
  double rOuter = pRMax - Hcal_lateral_plate_thickness;
  
  PolyhedraRegular EndCapRingChamberPoly(numSide, -pi/numSide, rInner, rOuter, Hcal_chamber_thickness);
  Box IntersectionStaveBox(rOuter/2., rOuter/2., Hcal_chamber_thickness/2);
  Position IntersectPos(rOuter/2. + Hcal_stave_gaps/2., rOuter/2. + Hcal_stave_gaps/2., 0.);
  IntersectionSolid  EndCapRingStaveSolid(EndCapRingChamberPoly, IntersectionStaveBox, IntersectPos);
  Volume EndCapRingChamberLogical(name+"_chamber", EndCapRingStaveSolid, air);
  EndCapRingChamberLogical.setAttributes(theDetector,x_layer.regionStr(),x_layer.limitsStr(),x_layer.visStr());
  
  double nRadiationLengthsInside=0.;
  double nInteractionLengthsInside=0.;
  double inner_thickness=0;
  double sensitive_thickness=0;
  double nRadiationLengths=0.;
  double nInteractionLengths=0.;
  double thickness_sum=0;

  double slice_pos_z = -Hcal_chamber_thickness/2.;
  int slice_number = 0;
  for(xml_coll_t k(x_layer,_U(slice)); k; ++k)  {
    xml_comp_t x_slice = k;
    string   slice_name      = name + _toString(slice_number,"_slice%d");
    double   slice_thickness = x_slice.thickness();
    Material slice_material  = theDetector.material(x_slice.materialStr());
    cout<<"  Layer_slice:  " <<  slice_name << " slice_thickness:  " << slice_thickness<< endl;

    nRadiationLengths   += slice_thickness/(2.*slice_material.radLength());
    nInteractionLengths += slice_thickness/(2.*slice_material.intLength());
    thickness_sum       += slice_thickness/2;
    if(x_slice.materialStr()==Hcal_radiator_material) continue;

    slice_pos_z += slice_thickness/2.;

    PolyhedraRegular slicePoly(numSide, -pi/numSide, rInner, rOuter, slice_thickness);
    IntersectionSolid sliceStaveSolid(slicePoly, IntersectionStaveBox, IntersectPos);
    Volume sliceVol(name + _toString(slice_number,"_slice%d"), sliceStaveSolid, slice_material);
    sliceVol.setAttributes(theDetector,x_slice.regionStr(),x_slice.limitsStr(),x_slice.visStr());
    if(x_slice.isSensitive()){
      sliceVol.setSensitiveDetector(sens);
      nRadiationLengthsInside   = nRadiationLengths;
      nInteractionLengthsInside = nInteractionLengths;
      inner_thickness           = thickness_sum;
      sensitive_thickness       = slice_thickness;
      
      nRadiationLengths=0.;
      nInteractionLengths=0.;
      thickness_sum = 0.;
    }

    nRadiationLengths   += slice_thickness/(2.*slice_material.radLength());
    nInteractionLengths += slice_thickness/(2.*slice_material.intLength());
    thickness_sum       += slice_thickness/2;
    // slice PlacedVolume
    PlacedVolume slice_phv = EndCapRingChamberLogical.placeVolume(sliceVol,Position(0,0,slice_pos_z));
    //DetElement slice(layer_name,_toString(slice_number,"slice%d"),x_det.id());
    slice.setPlacement(slice_phv);
    // Increment x position for next slice.
    slice_pos_z += slice_thickness/2.;
    // Increment slice number.
    ++slice_number;
  }
  
  // chamber placements
  for(int stave_id = 1; stave_id <= 4; stave_id++){
    double angle = -pi/2.*(stave_id-1);
    RotationZYX lrot(angle,0,0);
    for (int layer_id = 1; layer_id <= number_of_chambers; layer_id++){
      double Zoff = -pDz + (layer_id-1)*layerThickness + Hcal_radiator_thickness + Hcal_chamber_thickness/2.;
      Position l_pos(0., 0., Zoff);
      Position l_new = lrot*l_pos;
      Transform3D ltran3D(lrot,l_new);
      PlacedVolume layer_phv = EndCapRingLogical.placeVolume(EndCapRingChamberLogical, ltran3D);
      layer_phv.addPhysVolID("layer",layer_id);
      layer_phv.addPhysVolID("stave",stave_id);
      
      //string l_name = _toString(layer_id,"layer%d");
      //string stave_name = _toString(stave_id,"stave%d");
      //DetElement layer(module_det, l_name+stave_name, det_id);
      layer.setPlacement(layer_phv);

      if(stave_id==1&&layer_id==1){
	cout << "Hcal_EndcapRing:  inner_thickness= " << inner_thickness << endl;
	cout << "Hcal_EndcapRing:  outer_thickness= " << thickness_sum << endl;
      }
      if(stave_id==1){//only one needed, according to wenxingfang's
	LayeredCalorimeterData::Layer caloLayer ;
	caloLayer.cellSize0 = cell_sizeX;
	caloLayer.cellSize1 = cell_sizeY;
	caloLayer.inner_nRadiationLengths   = nRadiationLengthsInside;
	caloLayer.inner_nInteractionLengths = nInteractionLengthsInside;
	caloLayer.inner_thickness           = inner_thickness;
	caloLayer.sensitive_thickness       = sensitive_thickness;
	caloLayer.outer_nRadiationLengths   = nRadiationLengths;
	caloLayer.outer_nInteractionLengths = nInteractionLengths;
	caloLayer.outer_thickness           = thickness_sum;
	
	caloLayer.distance = start_z + (layer_id-1)*layerThickness;
	caloLayer.absorberThickness = Hcal_radiator_thickness ;
      
	caloData->layers.push_back( caloLayer ) ;
      }
    }
  }
  
  // Placements
  double endcap_z_offset = start_z + pDz;
  for(int side = 0; side <= 1; side++){
    int module_id = (side==0) ? 6 : 0;
    double this_module_z_offset = (side==0) ? endcap_z_offset : -endcap_z_offset;
    // use reflect volume for z<0, therefore, same rotation
    // segmentation violation happen if EndCapRingLogical.reflect(), back to rotate Y
    // double this_module_rotY = (side==0) ? 0.0 : 0.0;
    double this_module_rotY = (side==0) ? 0.0 : pi;
    //double this_module_rotZ = (side==0) ? pi/8. : pi/8;
    RotationZYX rot(0,this_module_rotY,0);
    Transform3D tran3D(rot,Position(0,0,this_module_z_offset));
    
    PlacedVolume module_phv;
    //if(side==0) module_phv = envelope.placeVolume(EndCapRingLogical, tran3D);
    //else        module_phv = envelope.placeVolume(EndCapRingLogical.reflect(), tran3D);
    module_phv = envelope.placeVolume(EndCapRingLogical, tran3D);
    
    module_phv.addPhysVolID("module", module_id);
    //DetElement sd = (module_id==0) ? module_det : module_det.clone(_toString(side,"module%d"));
    module.setPlacement(module_phv);
  }

  det.addExtension<LayeredCalorimeterData>(caloData) ;

  return det;
}

DECLARE_DETELEMENT(SHcalRpc01_EndcapRing, create_detector)
