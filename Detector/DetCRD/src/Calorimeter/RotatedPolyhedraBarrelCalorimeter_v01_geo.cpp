//==========================================================================
// RotatedPolyhedraBarrelCalorimeter_v01 implementation 
//--------------------------------------------------------------------------
// Author: FU Chengdong, IHEP
//==========================================================================

#include "DD4hep/DetFactoryHelper.h"
#include "XML/Utilities.h"
#include "XML/Layering.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens) {
  std::cout << "This is the RotatedPolyhedraBarrelCalorimeter_v01:"  << std::endl;

  xml_det_t x_det = e;
  Layering layering(x_det);
  xml_comp_t staves = x_det.staves();
  xml_dim_t dim = x_det.dimensions();
  string det_name = x_det.nameStr();
  Material air = description.air();
  double totalThickness = layering.totalThickness();
  double gap = dim.gap();//xml_dim_t(x_det).gap();
  int numSides = dim.numsides();
  double zhalf = dim.zhalf();
  double rmin = dim.rmin();
  double phi0 = dim.phi0();
  double zpos = dim.zpos();
  std::cout << "rmin    = " << rmin << std::endl
	    << "rmax    = " << rmin+totalThickness << std::endl
	    << "nstave  = " << numSides << std::endl
	    << "phi0    = " << phi0 << std::endl
	    << "zoffset = " << zpos << std::endl;

  DetElement cal(det_name, x_det.id());

  // --- create an envelope volume and position it into the world ---------------------
  Volume envelope = dd4hep::xml::createPlacedEnvelope( description, e, cal ) ;
  dd4hep::xml::setDetectorTypeFlag( e, cal ) ;
  if( description.buildType() == BUILD_ENVELOPE ) return cal ;
  envelope.setVisAttributes(description, "SeeThrough");

  //Volume motherVol = description.pickMotherVolume(sdet);
  //PolyhedraRegular polyhedra(numSides, rmin, rmin + totalThickness, zhalf*2);
  //Volume wholeVol(det_name+"_Whole", polyhedra, air);
  //wholeVol.setAttributes(description, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());

  // Add the subdetector envelope to the structure.
  DetElement stave0(cal, "stave1", x_det.id());
  Material matStave = air;
  std::string visStave = "SeeThrough";
  if ( staves )  {
    visStave = staves.visStr();
    matStave = description.material(staves.materialStr());
  }

  double innerAngle = 2*M_PI/numSides;
  double halfInnerAngle = innerAngle/2;
  double tan_inner = std::tan(halfInnerAngle)*2;
  double innerFaceLen = rmin*tan_inner + totalThickness/std::sin(innerAngle);
  double outerFaceLen = (rmin+totalThickness)*tan_inner - totalThickness/std::sin(innerAngle);
  double staveThickness = totalThickness;

  Trapezoid staveTrdOuter(innerFaceLen/2, outerFaceLen/2, zhalf, zhalf, staveThickness/2);
  Volume staveOuterVol("stave_outer", staveTrdOuter, air);
  staveOuterVol.setVisAttributes(description, "SeeThrough");

  Trapezoid staveTrdInner(innerFaceLen/2 - gap/2, outerFaceLen/2 - gap/2, zhalf, zhalf, staveThickness/2);
  Volume staveInnerVol("stave_inner", staveTrdInner, matStave);
  staveInnerVol.setVisAttributes(description, visStave);

  double layerOuterAngle = (M_PI - innerAngle)/2;
  double layerInnerAngle = (M_PI/2 - layerOuterAngle);
  double layerSideAngle  = M_PI - layerOuterAngle*2;
  double layer_pos_z = -staveThickness/2;
  double layer_dim_x = innerFaceLen/2 - gap/2;
  int layer_num = 1;

  //#### LayeringExtensionImpl* layeringExtension = new LayeringExtensionImpl();
  //#### Position layerNormal(0,0,1);

  for (xml_coll_t xc(x_det, _U(layer)); xc; ++xc) {
    xml_comp_t x_layer = xc;
    int repeat = x_layer.repeat();            // Get number of times to repeat this layer.
    const Layer* lay = layering.layer(layer_num - 1); // Get the layer from the layering engine.
    // Loop over repeats for this layer.
    for (int j = 0; j < repeat; j++) {
      string layer_name = _toString(layer_num, "layer%d");
      double layer_thickness = lay->thickness();
      DetElement layer(stave0, layer_name, layer_num);
      //### layeringExtension->setLayer(layer_num, layer, layerNormal);
      
      layer_dim_x -= layer_thickness/std::tan(layerSideAngle);
      // Layer position in Z within the stave.
      layer_pos_z += layer_thickness/2;
      // Layer box & volume
      Volume layer_vol(layer_name, Box(layer_dim_x, zhalf, layer_thickness/2), air);

      // Create the slices (sublayers) within the layer.
      double slice_pos_z = -layer_thickness/2;
      int slice_number = 1;
      for (xml_coll_t xk(x_layer, _U(slice)); xk; ++xk) {
        xml_comp_t x_slice = xk;
        string slice_name = _toString(slice_number, "slice%d");
        double slice_thickness = x_slice.thickness();
        Material slice_material = description.material(x_slice.materialStr());
        DetElement slice(layer, slice_name, slice_number);

        slice_pos_z += slice_thickness/2;
        // Slice volume & box
        Volume slice_vol(slice_name, Box(layer_dim_x, zhalf, slice_thickness/2), slice_material);

        if (x_slice.isSensitive()) {
          sens.setType("calorimeter");
          slice_vol.setSensitiveDetector(sens);
        }
        // Set region, limitset, and vis.
        slice_vol.setAttributes(description, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());
        // slice PlacedVolume
        PlacedVolume slice_phv = layer_vol.placeVolume(slice_vol, Position(0, 0, slice_pos_z));
        slice_phv.addPhysVolID("slice", slice_number);

        slice.setPlacement(slice_phv);
        // Increment Z position for next slice.
        slice_pos_z += slice_thickness / 2;
        // Increment slice number.
        ++slice_number;
      }
      // Set region, limitset, and vis.
      layer_vol.setAttributes(description, x_layer.regionStr(), x_layer.limitsStr(), x_layer.visStr());

      // Layer physical volume.
      PlacedVolume layer_phv = staveInnerVol.placeVolume(layer_vol, Position(0, 0, layer_pos_z));
      layer_phv.addPhysVolID("layer", layer_num);
      layer.setPlacement(layer_phv);

      // Increment the layer X dimension.
      //layer_dim_x += layer_thickness * std::tan(layerInnerAngle);    // * 2;
      // Increment the layer Z position.
      layer_pos_z += layer_thickness / 2;
      // Increment the layer number.
      ++layer_num;
    }
  }

  // Add stave inner physical volume to outer stave volume.
  PlacedVolume pv = staveOuterVol.placeVolume(staveInnerVol, dd4hep::Position(gap/2,0,0));
  if ( staves )  {
    stave0.setVisAttributes(description, staves.visStr(), staveInnerVol);
    stave0.setVisAttributes(description, staves.visStr(), staveOuterVol);
  }
  // Place the staves.
  double innerRotation = innerAngle;
  double offsetRotation = phi0 - M_PI;//-innerRotation/2;
  double sectCenterRadius = rmin + totalThickness/2;
  double offset = totalThickness/std::sin(innerAngle)/2;
  double rotX = M_PI / 2;
  //double rotY = -offsetRotation;
  //double posX = -sectCenterRadius * std::sin(rotY);
  //double posY = sectCenterRadius * std::cos(rotY);

  for (int istave = 1; istave <= numSides; istave++){
  //for (int istave = 1; istave <= 3; istave++) {
    DetElement stave = istave > 1 ? stave0.clone(_toString(istave,"stave%d")) : stave0;
    double rotY = offsetRotation - (istave-1)*innerRotation;
    double posX = sectCenterRadius*std::sin(rotY) + offset*std::cos(rotY);
    double posY = -sectCenterRadius*std::cos(rotY) + offset*std::sin(rotY);
    Transform3D trafo(RotationZYX(0, rotY, rotX), Translation3D(posX, posY, zpos));
    PlacedVolume pv = envelope.placeVolume(staveOuterVol, trafo);
    // Not a valid volID: pv.addPhysVolID("stave", 0);
    pv.addPhysVolID("stave", istave);
    stave.setPlacement(pv);
    //cal.add(stave);
  }
  
  //placeStaves(cal, stave, rmin, numSides, totalThickness, envelopeVol, innerAngle, staveOuterVol);
  // Set envelope volume attributes.
  //envelope.setAttributes(description, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());

  //double z_offset = dim.hasAttr(_U(z_offset)) ? dim.z_offset() : 0.0;
  //Transform3D transform(RotationZ(M_PI / numSides), Translation3D(0, 0, z_offset));
  //PlacedVolume pv_whole = envelope.placeVolume(wholeVol, transform);
  //pv_whole.addPhysVolID("system", cal.id());
  //pv_whole.addPhysVolID("barrel", 0);
  //cal.setPlacement(pv_whole);

  //#### cal.addExtension<SubdetectorExtension>(new SubdetectorExtensionImpl(cal));
  //#### cal.addExtension<LayeringExtension>(layeringExtension);
  return cal;
}

DECLARE_DETELEMENT(DD4hep_RotatedPolyhedraBarrelCalorimeter_v01, create_detector)
DECLARE_DEPRECATED_DETELEMENT(RotatedPolyhedraBarrelCalorimeter_v01, create_detector)
