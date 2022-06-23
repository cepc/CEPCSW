//==========================================================================
// RotatedCrystalCalorimeter implementation
//--------------------------------------------------------------------------
// rmin  : inner radius (minimum of 8 vertex of crystal)
// rmax  : outer radius (if gap=0, maximum of 8 vertex of crystal)
// alpha : rotate angle
// nphi  : number of crystals rotated in phi direction
// nz    : number of crystals in z direction
// gap   : space in tail of crystal for SiPM, electronics, etc.  
// Author: FU Chengdong, IHEP
//==========================================================================

#include "DD4hep/DetFactoryHelper.h"
#include "XML/Utilities.h"
#include "XML/Layering.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens) {
  std::cout << "This is the RotatedCrystalCalorimeter_v01:"  << std::endl;

  xml_det_t x_det = e;
  xml_dim_t dim = x_det.dimensions();
  string det_name = x_det.nameStr();
  Material air = description.air();
  double zhalf = dim.zhalf();
  double rmin  = dim.rmin();
  double rmax  = dim.rmax();
  double alpha = dim.alpha();
  int    nphi  = dim.nphi();
  int    nz    = dim.nz();
  double gap   = dim.gap();
  std::cout << "rmin    = " << rmin/mm << std::endl
	    << "rmax    = " << rmax/mm << std::endl
	    << "zhalf   = " << zhalf/mm << std::endl
	    << "alpha   = " << alpha/degree << std::endl
	    << "nphi    = " << nphi << std::endl
	    << "nz      = " << nz << std::endl;
  
  DetElement cal(det_name, x_det.id());
  if(alpha>90*dd4hep::degree){
    std::cout << "alpha>90 degree! return null cal!" << std::endl;
    return cal;
  }
  
  // --- create an envelope volume and position it into the world ---------------------
  Volume envelope = dd4hep::xml::createPlacedEnvelope( description, e, cal );
  dd4hep::xml::setDetectorTypeFlag( e, cal ) ;
  if( description.buildType() == BUILD_ENVELOPE ) return cal;
  envelope.setVisAttributes(description, "SeeThrough");

  DetElement module0(cal, "module1", x_det.id());
  
  double yhalf = zhalf/nz;
  Tube moduleTube(rmin, rmax, yhalf);
  Volume moduleVol("module", moduleTube, air);
  moduleVol.setVisAttributes(description, "SeeThrough");
  //
  //          |         __C
  //          |      __/  |
  //          |   __/     | 
  //          |__/        |
  //         B|___________|__E=C' of another crystal
  //         A|           D 
  //          |             OA=rmin,OC=OE=rmax
  //          |             /_BAD=alpha
  //          |O
  double dphi = 2*M_PI/nphi;
  double angleEAO = M_PI - alpha;
  double angleAEO = std::asin(sin(angleEAO)/rmax*rmin);
  double AE = rmax/sin(angleEAO)*sin(alpha-angleAEO);
  double angleOpen = dphi;
  double angleBottom = (M_PI-angleOpen)/2;
  double disNeighbor = rmin*sin(dphi/2)*2;
  //double x1 = disNeighbor/sin(angleBottom)*sin(angleEAO-angleBottom);
  double x1 = disNeighbor/sin(angleBottom)*sin(angleOpen+angleBottom-alpha);
  double edge = AE - disNeighbor/sin(angleBottom)*sin(alpha);
  double length = edge*cos(angleOpen/2) - gap;
  double x2 = x1 + length*tan(angleOpen/2)*2;
  double center2A = sqrt(x1*x1/4+length*length/4);
  double angleCenterAO = angleEAO+(M_PI-angleBottom-atan(length/x1));
  double center2O = sqrt(rmin*rmin+x1*x1/4+length*length/4-2*rmin*center2A*cos(angleCenterAO));
  double phi0Center = -asin(sin(angleCenterAO)/center2O*center2A);
  std::cout << "Single crystal's open angle = " << angleOpen/degree << std::endl;
  std::cout << "AE = " << AE/mm << ", /_AOE = " << (alpha-angleAEO)/degree << std::endl;
  std::cout << "OA = " << rmin/mm << ", /_AEO = " << angleAEO/degree << std::endl;
  std::cout << "OE = " << rmax/mm << ", /_EAO = " << angleEAO/degree << std::endl;
  std::cout << "x1 = " << x1/mm << " x2 = " << x2/mm << " y = " << yhalf*2/mm << " l = " << length/mm << std::endl; 

  xml_comp_t x_crystal = x_det.child( _Unicode(crystal) );
  Material crystalMat = description.material(x_crystal.materialStr());
  string crystalVis = x_crystal.visStr();

  Trapezoid crystalTrap(x1/2, x2/2, yhalf, yhalf, length/2);
  Volume crystalVol("crystal_whole", crystalTrap, crystalMat);
  crystalVol.setVisAttributes(description, crystalVis);

  sens.setType("calorimeter");
  // TODO: crystal pack
  crystalVol.setSensitiveDetector(sens);
  for (xml_coll_t xp(x_crystal, _U(slice)); xp; ++xp){
    //
    xml_comp_t x_slice = xp;
  }
  
  for(int crystal_id=1; crystal_id<=nphi; crystal_id++){
    double angleRot = -alpha + dphi/2 + (crystal_id-1)*dphi;
    double phiCenter = phi0Center + (crystal_id-1)*dphi;
    Transform3D trafo(RotationZYX(0, angleRot+M_PI/2, M_PI/2), Translation3D(center2O*cos(phiCenter), center2O*sin(phiCenter), 0));
    PlacedVolume pv = moduleVol.placeVolume(crystalVol, trafo);
    pv.addPhysVolID("crystal", crystal_id);
    DetElement crystal(module0, _toString(crystal_id,"crystal%d"), x_det.id());
    crystal.setPlacement(pv);
  }

  for (int module_id = 1; module_id <= nz; module_id++){
    DetElement module = module_id > 1 ? module0.clone(_toString(module_id,"module%d")) : module0;
    PlacedVolume pv = envelope.placeVolume(moduleVol, Position(0,0,-zhalf+(2*module_id-1)*yhalf));
    pv.addPhysVolID("module", module_id);
    module.setPlacement(pv);
  }

  return cal;
}

DECLARE_DETELEMENT(DD4hep_RotatedCrystalCalorimeter_v01, create_detector)
