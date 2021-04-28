//====================================================================
//  Detector description for CEPC Reference Detector ECal Barrel.
//--------------------------------------------------------------------
//
//  Author: Fangyi Guo    email: guofangyi@ihep.ac.cn
//  Ecal cosists of long crystal bar 1cm*1cm*~40cm 
//  8 parts cover 2pi phi range. Inner radius, height, z-length can change from xml.
//  In each part, crystal bar crosses in odd-even layer. 
//
//  v0r0: nothing but BGO crystal
//  v01: use double-layer as basic block, for cross-location in reconstruction. 
//		layout: 8 module 
//						->14 double-layer(dlayer) in R-direction 
//						->10 part along Z-direction
//						->4 block along phi-direction   (4*10*14 blocks in each module).
//						->2 sub-layer (bars in 0 sub-layer along phi direction, in 1 sub-layer along z direction).
//						->N bar. 
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"
#include "DDSegmentation/Segmentation.h"

#define MYDEBUG(x) std::cout << __FILE__ << ":" << __LINE__ << ": " << x << std::endl;
#define MYDEBUGVAL(x) std::cout << __FILE__ << ":" << __LINE__ << ": " << #x << ": " << x << std::endl;

using dd4hep::rec::LayeredCalorimeterData;
using namespace dd4hep;
using namespace std;
static dd4hep::Ref_t create_detector(dd4hep::Detector& theDetector,
                                     xml_h e,
                                     dd4hep::SensitiveDetector sens) {

  xml_det_t x_det = e;
  
  std::string det_name = x_det.nameStr();
  std::string det_type = x_det.typeStr();
  MYDEBUGVAL(det_name);
  MYDEBUGVAL(det_type);
  int detid = x_det.id();
  
  //Global geometry
  double R0 = theDetector.constant<double>("ecalbarrel_inner_radius");
  double h0 = theDetector.constant<double>("ecalbarrel_thickness");
  double Z0 = theDetector.constant<double>("ecalbarrel_zlength");
  double barx = theDetector.constant<double>("bar_x");   //Crystal size in R direction. 
  double bary = theDetector.constant<double>("bar_y");   //Crystal size in z/phi direction (z for odd layer, phi for even layer).
  
  double dim_x1 = R0*tan(22.5*degree) + sqrt(2)*h0/2.;
  double dim_x2 = dim_x1 - h0;
  double dim_y = Z0/2.;
  double dim_z = h0/2.;		                 
  double dx = dim_x1 - R0*tan(22.5*degree);  //transport distance in x-axis
  double r0 = R0+h0/2.;                      //rotation radius 
  
  //Crystal bar size
  int Nlayers = (int)h0/(2*barx);         //14 double-layers. 
  int Nblock_z   = 11;                    //block number in z direction
  int Nblock_phi = 4;                     //block number in phi direction
  double barz_s0;                         //Crystal bar lenghth in sub-layer 0(phi direction). Depends on layer number. 
  double barz_s1 = Z0/Nblock_z;           //Crystal bar lenghth in sub-layer 1(z direction, 46cm).
  int Nbar_phi;                           //Crystal bar number in each block, in phi direction.
  int Nbar_z = (int)barz_s1/bary;         //Crystal bar number in each block, in z direction.
  
  //Define detector and motherVolume(world)
  dd4hep::DetElement ECAL(det_name, detid);
  dd4hep::Volume motherVol = theDetector.pickMotherVolume(ECAL);
  
  // Create a Tube-like envelope representing the whole detector volume
  dd4hep::PolyhedraRegular envelope(8, 22.5*degree, R0, (R0+h0), Z0);
  dd4hep::Material	air(theDetector.material("Air"));
  dd4hep::Volume	   envelopeVol(det_name, envelope, air);
  dd4hep::PlacedVolume	envelopePlv = motherVol.placeVolume(envelopeVol, Position(0,0,0));
  envelopePlv.addPhysVolID("system",x_det.id());
  envelopeVol.setVisAttributes(theDetector, "InvisibleWithChildren" );
  ECAL.setPlacement(envelopePlv);

  //Define specific material and volumen for detElement
  dd4hep::Material mat_BGO(theDetector.material("G4_BGO"));
  dd4hep::Trapezoid trap(dim_x1, dim_x2, dim_y, dim_y, dim_z);
  dd4hep::Volume det_vol("trap_vol", trap, mat_BGO);
  det_vol.setVisAttributes(theDetector, "InvisibleWithChildren");

  dd4hep::Trapezoid subtrap(dim_x1, dim_x2, dim_y/Nblock_z, dim_y/Nblock_z, dim_z);
  dd4hep::Volume det_stave("stave_vol", subtrap, mat_BGO);
  det_stave.setVisAttributes(theDetector, "InvisibleWithChildren");  

  // Create extension objects for reconstruction
  LayeredCalorimeterData* caloData = new LayeredCalorimeterData ;
  for(int il=0;il<Nlayers; il++){
    //used for reconstruction, so write a 1*1*2 layer cell size. No absorber or dead-meaterial.
    dd4hep::rec::LayeredCalorimeterData::Layer _caloLayer;
    _caloLayer.distance                        = R0+il*2*barx;
    _caloLayer.phi0                               = 0;
    _caloLayer.absorberThickness            = 0;
    _caloLayer.inner_nRadiationLengths   = 0.01;
    _caloLayer.inner_nInteractionLengths = 0.01;
    _caloLayer.outer_nRadiationLengths   = 0.01;
    _caloLayer.outer_nInteractionLengths = 0.01;
    _caloLayer.inner_thickness              = barx;    //1cm
    _caloLayer.outer_thickness              = barx;    //1cm
    _caloLayer.sensitive_thickness       = 2*barx; //2cm
    _caloLayer.cellSize0                         = barx;    //1cm
    _caloLayer.cellSize1                         = barx;    //1cm
  caloData->layers.push_back(_caloLayer);
  }

  caloData->layoutType = LayeredCalorimeterData::BarrelLayout ;
  caloData->inner_symmetry = 8  ;
  caloData->outer_symmetry = 8  ;
  caloData->phi0 = 0 ; // hardcoded

  // extent of the calorimeter in the r-z-plane [ rmin, rmax, zmin, zmax ] in mm.
  caloData->extent[0] = R0 ;
  caloData->extent[1] = R0+h0;
  caloData->extent[2] = 0. ;
  caloData->extent[3] = Z0 ;


  //Loop to place crystalls in one part
  //Outer loop: layer (odd layer). 
  dd4hep::Volume bar_s1("bar_s1", dd4hep::Box(bary/2, barz_s1/2, barx/2), mat_BGO);
  bar_s1.setVisAttributes(theDetector, "VisibleRed");
  bar_s1.setSensitiveDetector(sens);

  dd4hep::DetElement stavedet(ECAL, "trap",detid);
  for(int ilayer=1; ilayer<=Nlayers; ilayer++){
    double lx = dim_x1 - ilayer*2*barx;
    barz_s0 = floor(2*lx/Nblock_phi);
    Nbar_phi = (int)barz_s0/bary;
    dd4hep::Volume bar_s0("bar_s0", dd4hep::Box(barz_s0/2, bary/2, barx/2), mat_BGO);
    bar_s0.setVisAttributes(theDetector, "VisibleGreen");
    bar_s0.setSensitiveDetector(sens);
    for(int iphi=1; iphi<=Nblock_phi; iphi++){
    	  dd4hep::Volume block("block", dd4hep::Box(barz_s0/2, barz_s1/2, barx), mat_BGO);
    	  block.setVisAttributes(theDetector, "VisibleGreen");
    	  //std::string blockname = "Block_"+std::to_string(ilayer)+"_"+std::to_string(iphi)+"_"+std::to_string(iz);
    	  std::string blockname = "Block_"+std::to_string(ilayer)+"_"+std::to_string(iphi);
    	  dd4hep::DetElement sd(stavedet, blockname, detid);
    
        //sub-layer 0: bars along phi. length=barz_s0. Bar num=Nbar_z
      for(int ibar0=1;ibar0<=Nbar_z;ibar0++){
    	    dd4hep::PlacedVolume plv_bar0 = block.placeVolume(bar_s0, Position(0,(2*ibar0-1)*bary/2-barz_s1/2, -barx/2));
    	    plv_bar0.addPhysVolID("slayer",0).addPhysVolID("bar",ibar0);
    	    std::string barname0 = "CrystalBar_s0_"+std::to_string(ibar0);	
    	    dd4hep::DetElement bardet0(sd, barname0, detid);
    	    bardet0.setPlacement(plv_bar0);
    	  }
    
    	  //sub-layer1 
    	  for(int ibar1=1;ibar1<=Nbar_phi;ibar1++){
    	    dd4hep::PlacedVolume plv_bar1 = block.placeVolume(bar_s1, Position((2*ibar1-1)*bary/2-barz_s0/2, 0, barx/2));
    	  	 plv_bar1.addPhysVolID("slayer",1).addPhysVolID("bar",ibar1);
    	  	 std::string barname1 = "CrystalBar_s1_"+std::to_string(ibar1);
    	  	 dd4hep::DetElement bardet1(sd, barname1, detid);
    	  	 bardet1.setPlacement(plv_bar1);
    	  }
    
    	  //dd4hep::PlacedVolume plv = det_vol.placeVolume(block, Position(lx-(2*iphi-1)*barz_s0/2, (2*iz-1)*barz_s1/2-dim_y, (2*ilayer-1)*bary-dim_z));
    	  dd4hep::PlacedVolume plv = det_stave.placeVolume(block, Position(lx-(2*iphi-1)*barz_s0/2, 0, (2*ilayer-1)*bary-dim_z));
    	  plv.addPhysVolID("dlayer", ilayer).addPhysVolID("part", iphi);
    	  sd.setPlacement(plv);
    
      }
    }

  for(int iz=1; iz<=Nblock_z; iz++){
    dd4hep::PlacedVolume plv = det_vol.placeVolume(det_stave, Position(0, (2*iz-1)*barz_s1/2-dim_y, 0) );
    plv.addPhysVolID("stave", iz);
    DetElement sd(stavedet, _toString(iz,"stave%3d"), detid);
    sd.setPlacement(plv);    
  }
  
  for(int i=0;i<8;i++){
    double rotAngle = 45*i*degree;
    double posx = -r0*sin(rotAngle) - dx*cos(rotAngle);
  	 double posy = r0*cos(rotAngle) - dx*sin(rotAngle);
  	 dd4hep::Transform3D transform(dd4hep::RotationZ(rotAngle)*dd4hep::RotationX(-90*degree),  dd4hep::Position(posx, posy, 0.));
  	 dd4hep::PlacedVolume plv = envelopeVol.placeVolume(det_vol, transform);
  	 plv.addPhysVolID("module", i);
  	 DetElement sd(ECAL, _toString(i,"trap%3d"), detid);
  	 sd.setPlacement(plv);
  }
  
  sens.setType("calorimeter");
  ECAL.addExtension< LayeredCalorimeterData >( caloData ) ; 

  MYDEBUG("create_detector DONE. ");
  return ECAL;
}

DECLARE_DETELEMENT(CRDEcalBarrel_v01, create_detector)
