#ifndef __ILDFTDDETECTOR__
#define __ILDFTDDETECTOR__

/** Petal based FTD to be used for ILD DBD studies 
 * WARNING: Still very experimental
 *
 * @author S.Aplin DESY, Robin Glattauer HEPHY
 */

#include "kaltest/TVKalDetector.h"

class TNode;
class TVector3;
class IGeomSvc;

namespace gear{
  class GearMgr ;
}


class ILDFTDKalDetector : public TVKalDetector {
public:
  
  /** Initialize the FTD from GEAR */
  ILDFTDKalDetector( const gear::GearMgr& gearMgr, IGeomSvc* geoSvc );
  
  
private:
  
  struct FTD_Petal {
    
    int    ipetal;
    double phi;
    double alpha;
    double rInner;
    double height;
    double innerBaseLength;
    double outerBaseLength;
    double senThickness;
    double supThickness;
    double senZPos;
    bool faces_ip;
    
  };
  
  
  struct FTD_Disk {
    int nPetals;
    double phi0;
    double dphi;
    
    double alpha;
    double rInner;
    double height;
    double innerBaseLength;
    double outerBaseLength;
    double senThickness;
    double supThickness;
    
    double stripAngle;
    
    double senZPos_even_front;
    double senZPos_odd_front;
    
    bool isDoubleSided;
    bool isStripReadout;
    
    int nSensors;
    
    
  };
  
 
  void build_staggered_design();
  
  //void create_petal(TVector3 measurement_plane_centre, FTD_Petal petal, int CellID);
  /**
   * @param zpos the z position of the front measurement surface (middle of front sensitive)
   */
  void create_segmented_disk_layers(int idisk, int nsegments, bool even_petals, double phi0, double zpos );
  
  
  void setupGearGeom( const gear::GearMgr& gearMgr ) ;
  
  int _nDisks ;
  double _bZ ;
  
   
  std::vector<FTD_Disk> _FTDgeo;
  
};

#endif
