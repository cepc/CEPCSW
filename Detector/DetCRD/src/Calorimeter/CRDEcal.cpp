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

	double dim_x1 = R0*tan(22.5*degree) + sqrt(2)*h0/2.;
	double dim_x2 = dim_x1 - h0;
	double dim_y = Z0/2.;
	double dim_z = h0/2.;		                 
	double dx = dim_x1 - R0*tan(22.5*degree);  //transport distance in x-axis
	double r0 = R0+h0/2.;							  //rotation radius 

	 //Crystal bar size
	double barx = 10*mm;             //Crystal size in R direction. 
	double bary = 10*mm;             //Crystal size in z/phi direction (z for odd layer, phi for even layer).
	int Nlayers = (int)h0/barx;
	int Nbarz_odd = (int)Z0/bary;
	int Nbarphi_odd; 	  				 //Depends on layer. Layer1~5: 5. Layer7~23: 4. Layer25~27: 3. 
	int Nbarz_even = 12; 
	int Nbarphi_even;                //Depends on layer
	double barz_even = Z0/Nbarz_even; //~38 cm
	double barz_odd; 					  //Depends on layer

	//Define detector and motherVolume(world)
	dd4hep::DetElement ECAL(det_name, detid);
	dd4hep::Volume motherVol = theDetector.pickMotherVolume(ECAL);

	// Create a Tube-like envelope representing the whole detector volume
	dd4hep::PolyhedraRegular envelope(8, 22.5*degree, R0, (R0+h0), Z0);
	dd4hep::Material	air(theDetector.material("Air"));
	dd4hep::Volume	envelopeVol(det_name, envelope, air);
	dd4hep::PlacedVolume	envelopePlv = motherVol.placeVolume(envelopeVol, Position(0,0,0));
	envelopeVol.setVisAttributes(theDetector, "InvisibleWithChildren" );
    envelopeVol.setRegion(theDetector,x_det.regionStr());
	ECAL.setPlacement(envelopePlv);

    //Define specific material and volumen for detElement
	dd4hep::Material mat_BGO(theDetector.material("G4_BGO"));
	dd4hep::Trapezoid trap(dim_x1, dim_x2, dim_y, dim_y, dim_z);
	dd4hep::Volume det_vol("trap_vol", trap, mat_BGO);
	det_vol.setVisAttributes(theDetector, "VisibleGreen");
	dd4hep::DetElement stavedet(ECAL, "trap",detid);

	//Loop to place crystalls in one part
	//Outer loop: layer (odd layer). 
	for(int ilayer=1; ilayer<=Nlayers; ilayer+=2){
		double lx = dim_x1 - ilayer*10*mm;
		//Loop in phi direction
		if(ilayer<=5) Nbarphi_odd=5; else if(ilayer<=23) Nbarphi_odd=4; else Nbarphi_odd=3;
		for(int iphi=1; iphi<=Nbarphi_odd; iphi++){
			barz_odd = 2*lx/Nbarphi_odd; 
			//Loop in Z direction
			dd4hep::Volume bar_odd("box_bar", dd4hep::Box(barz_odd/2, bary/2, barx/2), mat_BGO);
    		bar_odd.setSensitiveDetector(sens);
			for(int iz=1; iz<=Nbarz_odd;iz++){
				dd4hep::PlacedVolume plv = det_vol.placeVolume(bar_odd, Position(lx-(2*iphi-1)*barz_odd/2, (2*iz-1)*barx/2-dim_y, (2*ilayer-1)*bary/2-dim_z));
				plv.addPhysVolID("layer", ilayer).addPhysVolID("block", iphi).addPhysVolID("bar", iz);
				std::string barname = "CrystalBar_"+std::to_string(ilayer)+"_"+std::to_string(iphi)+"_"+std::to_string(iz);
				dd4hep::DetElement sd(stavedet, barname, detid);
				sd.setPlacement(plv);
			}
		}
	}

	//Loop in even layer
	dd4hep::Volume bar_even("box_bar", dd4hep::Box(bary/2, barz_even/2, barx/2), mat_BGO);
	bar_even.setSensitiveDetector(sens);
	for(int ilayer=2; ilayer<=Nlayers; ilayer+=2){
		double lx = dim_x1 - ilayer*10*mm;
		//Loop in phi direction
		Nbarphi_even = (int)floor(2*lx/bary);
		for(int iphi=1; iphi<=Nbarphi_even; iphi++){
			//Loop in Z direction
			for(int iz=1; iz<=Nbarz_even;iz++){
				dd4hep::PlacedVolume plv = det_vol.placeVolume(bar_even, Position(lx-(2*iphi-1)*bary/2, (2*iz-1)*barz_even/2-dim_y, (2*ilayer-1)*bary/2-dim_z));
				plv.addPhysVolID("layer", ilayer).addPhysVolID("block", iphi).addPhysVolID("bar", iz);
				std::string barname = "CrystalBar_"+std::to_string(ilayer)+"_"+std::to_string(iphi)+"_"+std::to_string(iz);
				dd4hep::DetElement sd(stavedet, barname, detid);
				sd.setPlacement(plv);
			}
		}
	}

	for(int i=0;i<8;i++){
		double rotAngle = 45*i*degree;
		double posx = -r0*sin(rotAngle) - dx*cos(rotAngle);
		double posy = r0*cos(rotAngle) - dx*sin(rotAngle);
		dd4hep::Transform3D transform(dd4hep::RotationZ(rotAngle)*dd4hep::RotationX(-90*degree),  dd4hep::Position(posx, posy, 0.));
		dd4hep::PlacedVolume plv = envelopeVol.placeVolume(det_vol, transform);
		plv.addPhysVolID("system", i);
		DetElement sd(ECAL, _toString(i,"trap%3d"), detid);
		sd.setPlacement(plv);
	}

	sens.setType("calorimeter");

	MYDEBUG("create_detector DONE. ");
	return ECAL;
}

DECLARE_DETELEMENT(CRDEcalBarrel, create_detector)
