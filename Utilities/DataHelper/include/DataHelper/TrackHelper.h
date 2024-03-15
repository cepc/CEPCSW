#ifndef TRACKHELPER_H
#define TRACKHELPER_H
#include "edm4hep/TrackState.h"
#include "TMatrixDSym.h"
#include "TVector3.h"

namespace CEPC{
    //get track position and momentum from TrackState
    void getPosMomFromTrackState(const edm4hep::TrackState& trackState,
            double Bz, TVector3& pos,TVector3& mom,double& charge,
            TMatrixDSym& covMatrix_6);

    //Set track state from position, momentum and charge
    void getTrackStateFromPosMom(edm4hep::TrackState& trackState,double Bz,
            TVector3 pos,TVector3 mom,double charge,TMatrixDSym covMatrix_6);

}

#endif
