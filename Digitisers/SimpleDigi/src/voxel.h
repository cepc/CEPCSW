#ifndef _voxel_included_
#define _voxel_included_

// A header file which defines a voxel class for the TPC
#include <vector>
//#include "ThreeVector.h"
#include <CLHEP/Vector/ThreeVector.h>

using namespace std;

class Voxel_tpc{

 public:
  Voxel_tpc();
  // the intialation in the constructor here would be preferable though I don't know how to intialise
  // the array xyz[3] here with pos[3], for the mean time the constructor will be put in the .cc file
  //  Voxel_tpc(int row, int phi, int z, double pos[3]) : row_index(row), phi_index(phi), z_index(z){}
  Voxel_tpc(int row, int phi, int z, double pos[3], double posRPhi[2], double edep, double rPhiRes, double zRes);
  Voxel_tpc(int row, int phi, int z, CLHEP::Hep3Vector coord, double edep, double rPhiRes, double zRes);
  ~Voxel_tpc();

  void setAdjacent(Voxel_tpc * p_voxel) { _adjacent_voxels.push_back(p_voxel);};
  void setIsClusterHit() { _isClusterHit = true;};
  void setIsMerged() { _isMerged = true;};
  bool IsClusterHit() { return _isClusterHit;};
  bool IsMerged() { return _isMerged;};
  int clusterFind(vector <Voxel_tpc*>* hitList);


  int getRowIndex() {return _row_index;};
  int getPhiIndex() {return _phi_index;};
  int getZIndex() {return _z_index;};
  Voxel_tpc * getFirstAdjacent() {return *(_adjacent_voxels.begin());};
  Voxel_tpc * getAdjacent(int i) {return _adjacent_voxels[i];};
  int getNumberOfAdjacent() {return _adjacent_voxels.size();};
  double getX() {return _coord.x();};
  double getY() {return _coord.y();};
  double getZ() {return _coord.z();};
  double getR() {return _coord.perp();};
  double getPhi() {return _coord.phi();};
  double getEDep() {return _edep;};
  double getRPhiRes() {return _rPhiRes;};
  double getZRes() {return _zRes;};
  const CLHEP::Hep3Vector getHep3Vector() {return _coord;};
  //  bool compare_phi( Voxel_tpc * & a, Voxel_tpc * & b);



 private:
  int _row_index;
  int _phi_index;
  int _z_index;
  vector <Voxel_tpc *> _adjacent_voxels;
  CLHEP::Hep3Vector _coord;
  double _edep;
  double _rPhiRes;
  double _zRes;
  bool _isMerged;
  bool _isClusterHit;
};
#endif
