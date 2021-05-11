#include "DataHelper/TrackerHitHelper.h"
#include "Identifier/CEPCConf.h"

#include "TMatrixF.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include <bitset>

std::array<float,6> CEPC::GetCovMatrix(edm4hep::TrackerHit& hit, bool useSpacePointBuilderMethod){
  if(hit.isAvailable()){
    int type = hit.getType();
    if(std::bitset<32>(type)[CEPCConf::TrkHitTypeBit::COMPOSITE_SPACEPOINT]){
      return hit.getCovMatrix();
    }
    else if(std::bitset<32>(type)[CEPCConf::TrkHitTypeBit::PLANAR]){
      float thetaU = hit.getCovMatrix(0);
      float phiU   = hit.getCovMatrix(1);
      float dU     = hit.getCovMatrix(2);
      float thetaV = hit.getCovMatrix(3);
      float phiV   = hit.getCovMatrix(4);
      float dV     = hit.getCovMatrix(5);
      return ConvertToCovXYZ(dU, thetaU, phiU, dV, thetaV, phiV, useSpacePointBuilderMethod);
    }
    else{
      std::cout << "Warning: not SpacePoint and Planar, return original cov matrix preliminaryly." << std::endl;
      return hit.getCovMatrix();
    }
  }
  std::array<float,6> cov{0.,0.,0.,0.,0.,0.};
  return cov;
}

float CEPC::GetResolutionRPhi(edm4hep::TrackerHit& hit){
  if(hit.isAvailable()){
    int type = hit.getType();
    if(std::bitset<32>(type)[CEPCConf::TrkHitTypeBit::COMPOSITE_SPACEPOINT]){
      return sqrt(hit.getCovMatrix(0)+hit.getCovMatrix(2));
    }
    else if(std::bitset<32>(type)[CEPCConf::TrkHitTypeBit::PLANAR]){
      return hit.getCovMatrix(2);
    }
    else{
      std::cout << "Warning: not SpacePoint and Planar, return value from original cov matrix preliminaryly." << std::endl;
      return sqrt(hit.getCovMatrix(0)+hit.getCovMatrix(2));
    }
  }
  return 0.;
}

float CEPC::GetResolutionZ(edm4hep::TrackerHit& hit){
  if(hit.isAvailable()){
    int type = hit.getType();
    if(std::bitset<32>(type)[CEPCConf::TrkHitTypeBit::COMPOSITE_SPACEPOINT]){
      return sqrt(hit.getCovMatrix(5));
    }
    else if(std::bitset<32>(type)[CEPCConf::TrkHitTypeBit::PLANAR]){
      return hit.getCovMatrix(5);
    }
    else{
      std::cout << "Warning: not SpacePoint and Planar, return value from original cov matrix preliminaryly." << std::endl;
      return sqrt(hit.getCovMatrix(5));
    }
  }
  return 0.;
}

std::array<float, 6> CEPC::ConvertToCovXYZ(float dU, float thetaU, float phiU, float dV, float thetaV, float phiV, bool useSpacePointBuilderMethod){
  std::array<float,6> cov{0.,0.,0.,0.,0.,0.};
  if(!useSpacePointBuilderMethod){
    TMatrixF diffs(2,3);
    TMatrixF diffsT(3,2);
    diffs(0,0) = sin(thetaU)*cos(phiU);
    diffs(0,1) = sin(thetaU)*sin(phiU);
    diffs(0,2) = cos(thetaU);
    diffs(1,0) = sin(thetaV)*cos(phiV);
    diffs(1,1) = sin(thetaV)*sin(phiV);
    diffs(1,2) = cos(thetaV);

    diffsT.Transpose(diffs);

    TMatrixF covMatrixUV(2,2);
    covMatrixUV(0,0) = dU*dU;
    covMatrixUV(0,1) = 0;
    covMatrixUV(1,0) = 0;
    covMatrixUV(1,1) = dV*dV;

    TMatrixF covMatrixXYZ(3,3);
    covMatrixXYZ = diffsT*covMatrixUV*diffs;
    cov[0] = covMatrixXYZ(0,0);
    cov[1] = covMatrixXYZ(1,0);
    cov[2] = covMatrixXYZ(1,1);
    cov[3] = covMatrixXYZ(2,0);
    cov[4] = covMatrixXYZ(2,1);
    cov[5] = covMatrixXYZ(2,2);
  }
  else{ // Method used in SpacePointBuilder, results are almost same with above
    CLHEP::Hep3Vector u_sensor(sin(thetaU)*cos(phiU), sin(thetaU)*sin(phiU), cos(thetaU));
    CLHEP::Hep3Vector v_sensor(sin(thetaV)*cos(phiV), sin(thetaV)*sin(phiV), cos(thetaV));
    CLHEP::Hep3Vector w_sensor = u_sensor.cross(v_sensor);
    CLHEP::HepRotation rot_sensor(u_sensor, v_sensor, w_sensor);
    CLHEP::HepMatrix rot_sensor_matrix;
    rot_sensor_matrix = rot_sensor;
    CLHEP::HepSymMatrix cov_plane(3,0);
    cov_plane(1,1) = dU*dU;
    cov_plane(2,2) = dV*dV;
    CLHEP::HepSymMatrix cov_xyz= cov_plane.similarity(rot_sensor_matrix);
    cov[0] = cov_xyz[0][0];
    cov[1] = cov_xyz[1][0];
    cov[2] = cov_xyz[1][1];
    cov[3] = cov_xyz[2][0];
    cov[4] = cov_xyz[2][1];
    cov[5] = cov_xyz[2][2];
  }
  return cov;
}
