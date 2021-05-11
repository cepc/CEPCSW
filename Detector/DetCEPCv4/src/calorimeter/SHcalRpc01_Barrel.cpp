//====================================================================
//  SHcalRpc01 - Implementation from ILCSoft's Mokka version                              
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetType.h"
#include "DDSegmentation/TiledLayerGridXY.h"

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
using dd4hep::RotationZ;
using dd4hep::Transform3D;
using dd4hep::Box;
using dd4hep::Tube;
using dd4hep::PolyhedraRegular;
using dd4hep::SubtractionSolid;
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
  cout << "--------------------------" << endl;
  cout << "creating SHcalRpc01_Barrel" << endl;
  cout << "--------------------------" << endl;

  xml_det_t    x_det  = element;
  string       name   = x_det.nameStr();
  int          det_id = x_det.id();
  DetElement   det(name, det_id);

  Volume envelope = dd4hep::xml::createPlacedEnvelope(theDetector, element , det ) ;

  dd4hep::xml::setDetectorTypeFlag(element, det) ;

  if( theDetector.buildType() == BUILD_ENVELOPE ) return det ;

  xml_comp_t    x_staves          = x_det.staves();
  string   Hcal_radiator_material = x_staves.materialStr();
  Material      stavesMaterial    = theDetector.material(Hcal_radiator_material);
  Material      air               = theDetector.air();

  sens.setType("calorimeter");

  Readout readout = sens.readout();
  Segmentation seg = readout.segmentation();
  dd4hep::DDSegmentation::TiledLayerGridXY* tiledSeg = dynamic_cast<dd4hep::DDSegmentation::TiledLayerGridXY*> (seg.segmentation());
  assert(tiledSeg && "no TiledLayerGridXY found" );

  std::vector<double> cellSizeVector = seg.segmentation()->cellDimensions(0);
  double cell_sizeX      = cellSizeVector[0];
  double cell_sizeZ      = cellSizeVector[1];

  double Hcal_inner_radius            = theDetector.constant<double>("Hcal_inner_radius");
  double Hcal_outer_radius_set        = theDetector.constant<double>("Hcal_outer_radius");
  double Hcal_half_length             = theDetector.constant<double>("Hcal_half_length");
  int    Hcal_inner_symmetry          = theDetector.constant<int>("Hcal_inner_symmetry");
  int    Hcal_outer_symmetry          = 0;
  double Hcal_lateral_plate_thickness = theDetector.constant<double>("Hcal_lateral_structure_thickness");
  double Hcal_modules_gap             = theDetector.constant<double>("Hcal_modules_gap");
  double Ecal_outer_radius            = theDetector.constant<double>("Ecal_outer_radius");
  int    Hcal_barrel_number_modules   = theDetector.constant<int>("Hcal_barrel_number_modules");
  
  double hPrime   = Ecal_outer_radius + theDetector.constant<double>("Hcal_Ecal_gap");
  Hcal_inner_radius = hPrime / cos(pi/8.);
  
  double Hcal_normal_dim_z = (2*Hcal_half_length - (Hcal_barrel_number_modules-1)*Hcal_modules_gap)/Hcal_barrel_number_modules;

  xml_coll_t c(x_det,_U(layer));
  xml_comp_t x_layer = c;
  int         Hcal_nlayers = x_layer.repeat();

  double Hcal_radiator_thickness = 0;
  double layerThickness = 0.0;
  for(xml_coll_t k(x_layer,_U(slice)); k; ++k)  {
    xml_comp_t x_slice = k;
    layerThickness += x_slice.thickness();
    if(x_slice.materialStr()==Hcal_radiator_material) Hcal_radiator_thickness = x_slice.thickness();
  }
  cout << " cell size xy = " <<  cell_sizeX << " cell size z = " << cell_sizeZ << endl;
  cout << " layer_thickness (from slices) = " << layerThickness << " and radiator_thickness = " << Hcal_radiator_thickness << endl;
  double Hcal_chamber_thickness = layerThickness - Hcal_radiator_thickness; 

  int MinNumCellsInTransvPlane  = theDetector.constant<int>("Hcal_MinNumCellsInTransvPlane");
  double RPC_EdgeWidth          = theDetector.constant<double>("Hcal_gas_edge_width");
  double RPCGazInletInnerRadius = theDetector.constant<double>("Hcal_gasInlet_inner_radius");
  double RPCGazInletOuterRadius = theDetector.constant<double>("Hcal_gasInlet_outer_radius");
  double RPCGazInletLength      = theDetector.constant<double>("Hcal_gasInlet_length");
  double RPC_PadSeparation      = theDetector.constant<double>("Hcal_pad_separation");
  double Hcal_spacer_thickness  = theDetector.constant<double>("Hcal_spacer_thickness");
  double Hcal_spacer_separation = theDetector.constant<double>("Hcal_spacer_separation");

  //========== fill data for reconstruction ============================
  LayeredCalorimeterData* caloData = new LayeredCalorimeterData ;
  caloData->layoutType = LayeredCalorimeterData::BarrelLayout ;
  caloData->inner_symmetry = Hcal_inner_symmetry  ;
  caloData->outer_symmetry = Hcal_outer_symmetry  ;
  caloData->phi0 = 0 ; // fg: also hardcoded below

  // general calculated parameters
  double AngleRatio=0.76536686;//"k"
  double d_InnerOctoSize=AngleRatio*Hcal_inner_radius;//"d"
  double LMin = 2*RPC_EdgeWidth+cell_sizeX*MinNumCellsInTransvPlane+(MinNumCellsInTransvPlane+1)*RPC_PadSeparation;

  double Ynl = 0.5*d_InnerOctoSize - Hcal_nlayers*layerThickness;
  double Hcal_outer_radius = sqrt((LMin-Ynl)*(LMin-Ynl) + (hPrime + Hcal_nlayers*layerThickness)*(hPrime + Hcal_nlayers*layerThickness));
  if(Hcal_outer_radius!=Hcal_outer_radius_set){
    cout << "calculated Hcal_outer_radius != input, will impact HcalEndcap and HcalEndcapRing. Hcal_outer_radius = " << Hcal_outer_radius
	 << " but set as " << Hcal_outer_radius_set << " difference = " << Hcal_outer_radius-Hcal_outer_radius_set << endl;  
  }

  /// extent of the calorimeter in the r-z-plane [ rmin, rmax, zmin, zmax ] in cm.
  caloData->extent[0] = Hcal_inner_radius ;
  caloData->extent[1] = Hcal_outer_radius ;
  caloData->extent[2] = 0. ; // Barrel zmin is "0" by default.
  caloData->extent[3] = Hcal_half_length ;

  double Hcal_total_dim_y = Hcal_outer_radius - hPrime;
  
  // the  y_dim1_for_z kept as the original value in TDR
  double Hcal_regular_chamber_dim_z = Hcal_normal_dim_z - 2 *(Hcal_lateral_plate_thickness);
  //int N_cells_z =  static_cast <int> ( (Hcal_regular_chamber_dim_z - 2*RPC_EdgeWidth - RPC_PadSeparation) / (Hcal_cell_dim_x + RPC_PadSeparation) );
  //  Hcal_cell_dim_z=(Hcal_regular_chamber_dim_z-RPC_PadSeparation )/N_cells_z
  //                      - RPC_PadSeparation;
  Tube solidCaloTube(0, Hcal_outer_radius, Hcal_half_length);

  PolyhedraRegular solidOctogon(8, 0, hPrime, 4*Hcal_half_length);
  RotationZYX rotOctogon(dd4hep::twopi/16,0,0);
  SubtractionSolid solidCalo(solidCaloTube, solidOctogon, rotOctogon);
  Volume logicCalo(name+"_radiator", solidCalo, stavesMaterial);
  logicCalo.setAttributes(theDetector,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());
  PlacedVolume calo_pv = envelope.placeVolume(logicCalo, Position(0,0,0));
  DetElement calo(det, "envelope", det_id);
  calo.setPlacement(calo_pv);
  if(tiledSeg) tiledSeg->setOffsetY(-(Hcal_regular_chamber_dim_z/2.-RPC_EdgeWidth)+0.5*cell_sizeZ);
  for(int layer_id=1; layer_id<=Hcal_nlayers; layer_id++){
    double yn = sqrt(Hcal_outer_radius*Hcal_outer_radius - (hPrime + layer_id*layerThickness)*(hPrime + layer_id*layerThickness));
    double Yn = 0.5*d_InnerOctoSize - layer_id*layerThickness;

    double halfX = Hcal_chamber_thickness/2.;
    double halfY = (yn+Yn)/2.;
    
    LayeredCalorimeterData::Layer caloLayer ;
    caloLayer.cellSize0 = cell_sizeX;
    caloLayer.cellSize1 = cell_sizeZ;

    //double halfZ = Hcal_normal_dim_z / 2.;
    double halfZ = Hcal_regular_chamber_dim_z / 2.;
    
    double localXPos = hPrime + Hcal_radiator_thickness + Hcal_chamber_thickness/2. + (layer_id-1)*layerThickness;
    double localYPos = -Yn + 0.5*(Yn + yn);

    Box chamberSolid(halfY, halfZ, halfX);
    string chamberLogical_name      = name+_toString(layer_id,"_layer%d");
    Volume chamberLogical(chamberLogical_name, chamberSolid, air);
    chamberLogical.setAttributes(theDetector, x_layer.regionStr(), x_layer.limitsStr(), x_layer.visStr());

    if(tiledSeg) tiledSeg->setLayerOffsetX((-(halfY-RPC_EdgeWidth)+0.5*cell_sizeX)*2/cell_sizeX);

    string layer_name      = name+_toString(layer_id,"_layer%d");

    double nRadiationLengths=0.;
    double nInteractionLengths=0.;
    double thickness_sum=0;

    nRadiationLengths   = Hcal_radiator_thickness/(stavesMaterial.radLength());
    nInteractionLengths = Hcal_radiator_thickness/(stavesMaterial.intLength());

    double slice_pos_z = -halfX;
    int slice_number = 0;
    for(xml_coll_t k(x_layer,_U(slice)); k; ++k)  {
      xml_comp_t x_slice = k;
      if(x_slice.materialStr()==Hcal_radiator_material) continue;
      string   slice_name      = layer_name + _toString(slice_number,"_slice%d");
      double   slice_thickness = x_slice.thickness();
      Material slice_material  = theDetector.material(x_slice.materialStr());
      if(layer_id==1) cout<<"  Layer_slice:  "<<  slice_name<<" slice_thickness:  "<< slice_thickness<< endl;
      
      slice_pos_z += slice_thickness/2.;
      nRadiationLengths   += slice_thickness/(2.*slice_material.radLength());
      nInteractionLengths += slice_thickness/(2.*slice_material.intLength());
      thickness_sum       += slice_thickness/2;

      // Slice volume & box
      Box sliceSolid(halfY, halfZ, slice_thickness/2.);
      Volume sliceVol(slice_name, sliceSolid, slice_material);
      
      if ( x_slice.isSensitive() ) {
	sliceVol.setSensitiveDetector(sens);
	if(RPC_EdgeWidth>0){
	  double RPC_GazInlet_In_Z  = halfZ - RPC_EdgeWidth - RPCGazInletOuterRadius;
	  double RPC_GazInlet_In_Y  = halfY - RPC_EdgeWidth/2;
	  double RPC_GazInlet_Out_Z = -RPC_GazInlet_In_Z;
	  double RPC_GazInlet_Out_Y =  RPC_GazInlet_In_Y;

	  string mateialName = x_slice.attr<string>(_Unicode(edge_material));
	  Material edge_material = theDetector.material(mateialName);
	  Box solidRPCEdge1(halfY, halfZ, slice_thickness/2.);
	  Box solidRPCEdge2(halfY-RPC_EdgeWidth, halfZ-RPC_EdgeWidth, slice_thickness/2.);
	  SubtractionSolid solidRPCEdge(solidRPCEdge1, solidRPCEdge2, Position(0,0,0));
	  Volume logicRPCEdge(slice_name+"_edge", solidRPCEdge, edge_material);
	  logicRPCEdge.setAttributes(theDetector,x_slice.regionStr(),x_slice.limitsStr(),x_slice.visStr());
	  sliceVol.placeVolume(logicRPCEdge);

	  RotationZYX rotGaz(0, pi/2., 0);
	  Tube solidRPCGazInlet(RPCGazInletInnerRadius,RPCGazInletOuterRadius,RPC_EdgeWidth/*RPCGazInletLength*//2);
	  Volume logicRPCGazInlet(slice_name+"_GazInlet", solidRPCGazInlet, edge_material);
	  logicRPCGazInlet.setAttributes(theDetector,x_slice.regionStr(),x_slice.limitsStr(),x_slice.visStr());
	  logicRPCEdge.placeVolume(logicRPCGazInlet, Transform3D(rotGaz, Position(RPC_GazInlet_In_Y,RPC_GazInlet_In_Z, 0)));
	  logicRPCEdge.placeVolume(logicRPCGazInlet, Transform3D(rotGaz, Position(RPC_GazInlet_Out_Y,RPC_GazInlet_Out_Z, 0)));
	  
	  Tube solidRPCGazInsideInlet(0,RPCGazInletInnerRadius,RPC_EdgeWidth/*RPCGazInletLength*//2);
	  Volume logicRPCGazInsideInlet(slice_name+"_GazInsideInlet", solidRPCGazInsideInlet, slice_material);
	  logicRPCGazInsideInlet.setAttributes(theDetector,x_slice.regionStr(),x_slice.limitsStr(),"SeeThrough");
	  logicRPCEdge.placeVolume(logicRPCGazInsideInlet, Transform3D(rotGaz, Position(RPC_GazInlet_In_Y,RPC_GazInlet_In_Z, 0)));
	  logicRPCEdge.placeVolume(logicRPCGazInsideInlet, Transform3D(rotGaz,Position(RPC_GazInlet_Out_Y,RPC_GazInlet_Out_Z, 0)));
	}
	if(Hcal_spacer_thickness>0){
	  Tube solidRPCSpacer(0,Hcal_spacer_thickness/2,slice_thickness/2);
	  Material space_material = theDetector.material(x_slice.attr<string>(_Unicode(spacer_material)));
	  Volume logicRPCSpacer(slice_name+"_spacer", solidRPCSpacer, space_material);
	  logicRPCSpacer.setAttributes(theDetector,x_slice.regionStr(),x_slice.limitsStr(),x_slice.visStr());
	  RotationZYX rotSpacer(0, 0, 0);
	  
	  double gap_hZ = halfZ-RPC_EdgeWidth;
	  double gap_hY = halfY-RPC_EdgeWidth;
	  int y_number_of_separations = (int)(2*gap_hY/Hcal_spacer_separation);
	  int z_number_of_separations = (int)(2*gap_hZ/Hcal_spacer_separation);
	  double y_lateral_space = (2*gap_hY - y_number_of_separations*Hcal_spacer_separation)/2;
	  double z_lateral_space = (2*gap_hZ - z_number_of_separations*Hcal_spacer_separation)/2;
	  if(y_lateral_space < Hcal_spacer_thickness/2.){
	    y_number_of_separations = (int)((2*gap_hY-Hcal_spacer_thickness)/Hcal_spacer_separation);
	    y_lateral_space = (2*gap_hY - y_number_of_separations*Hcal_spacer_separation)/2;
	  }
	  if(z_lateral_space < Hcal_spacer_thickness/2.){
	    z_number_of_separations = (int)((2*gap_hZ-Hcal_spacer_thickness)/Hcal_spacer_separation);
	    z_lateral_space = (2*gap_hZ - z_number_of_separations*Hcal_spacer_separation)/2;
	  }
	  for(int y_counter = 0; y_counter <=y_number_of_separations; y_counter++){
	    double SpacerY = gap_hY - y_lateral_space - y_counter*Hcal_spacer_separation;
	    for(int z_counter = 0; z_counter <=z_number_of_separations; z_counter++){
	      double SpacerZ = gap_hZ - z_lateral_space - z_counter*Hcal_spacer_separation;
	      PlacedVolume space_pv = sliceVol.placeVolume(logicRPCSpacer, Transform3D(rotSpacer, Position(SpacerY,SpacerZ,0)));
	    }
	  }
	}

	caloLayer.inner_nRadiationLengths = nRadiationLengths;
	caloLayer.inner_nInteractionLengths = nInteractionLengths;
	caloLayer.inner_thickness = thickness_sum;
	if(layer_id==1) cout<<"Hcal_Barrel:  inner_thickness= "<<thickness_sum<<endl;
        //Store readout gasgap thickness
	caloLayer.sensitive_thickness = slice_thickness;
        //Reset counters to measure "outside" quantitites
	nRadiationLengths=0.;
	nInteractionLengths=0.;
	thickness_sum = 0.;
	
	sliceVol.setAttributes(theDetector,x_slice.regionStr(),x_slice.limitsStr(),"SeeThrough");
      }
      else{
	sliceVol.setAttributes(theDetector,x_slice.regionStr(),x_slice.limitsStr(),x_slice.visStr());
      }
      nRadiationLengths += slice_thickness/(2.*slice_material.radLength());
      nInteractionLengths += slice_thickness/(2.*slice_material.intLength());
      thickness_sum += slice_thickness/2;

      // slice PlacedVolume
      PlacedVolume slice_phv = chamberLogical.placeVolume(sliceVol,Position(0,0,slice_pos_z));
      if ( x_slice.isSensitive() ) {
	int slice_id  = (layer_id > Hcal_nlayers)? 1:-1;
	slice_phv.addPhysVolID("layer",layer_id).addPhysVolID("slice",slice_id);
      }
      DetElement sliceDetE(layer_name,_toString(slice_number,"slice%d"),x_det.id());
      sliceDetE.setPlacement(slice_phv);
      // Increment x position for next slice.
      slice_pos_z += slice_thickness/2.;
      // Increment slice number.
      ++slice_number;
    }
    caloLayer.outer_nRadiationLengths = nRadiationLengths;
    caloLayer.outer_nInteractionLengths = nInteractionLengths;
    caloLayer.outer_thickness = thickness_sum;
    if(layer_id==1) cout << "Hcal_Barrel:  outer_thickness= " << thickness_sum << endl;
    
    double chamber_y_offset = -(-Hcal_total_dim_y/2. + (layer_id-1)*layerThickness + layerThickness/2.);

    caloLayer.distance = Hcal_inner_radius + Hcal_total_dim_y/2.0 + chamber_y_offset ;
    caloLayer.absorberThickness = Hcal_radiator_thickness ;

    caloData->layers.push_back( caloLayer ) ;
    
    double stave_phi_offset, module_z_offset;
    
    stave_phi_offset = pi*0.5;
    for(int stave_id = 1; stave_id <= 8; stave_id++){
      double phirot = stave_phi_offset+(stave_id-1)*pi/4.;

      RotationZYX rot(pi/2, pi/2, 0); //phirot);
      RotationZ rotZ(phirot);
      RotationZYX rotAll = rotZ*rot;
      RotationZYX rotInverse(phirot, 0, 0);
      for(int module_id = 1; module_id <= Hcal_barrel_number_modules; module_id++){
        module_z_offset = - Hcal_half_length + Hcal_normal_dim_z/2. + (module_id-1)*(Hcal_normal_dim_z+Hcal_modules_gap);
	
        Position localPos(localXPos,localYPos,module_z_offset);
        Position newPos = rotInverse*localPos;

	Transform3D tran3D(rotAll, newPos);
	PlacedVolume pv = logicCalo.placeVolume(chamberLogical, tran3D);
	pv.addPhysVolID("stave",stave_id).addPhysVolID("module",module_id).addPhysVolID("layer",layer_id);
	DetElement layer(calo, name+_toString(stave_id,"_stave%d")+_toString(module_id,"_module%d")+_toString(layer_id,"_layer%d"), det_id);
	layer.setPlacement(pv);
      }
    }
  }

  det.addExtension< LayeredCalorimeterData >( caloData ) ;

  return det;
}

DECLARE_DETELEMENT(SHcalRpc01_Barrel, create_detector)
