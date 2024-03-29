/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert & Johannes Rauch

   This file is part of GENFIT.

   GENFIT is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   GENFIT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with GENFIT.  If not, see <http://www.gnu.org/licenses/>.
*/
/** @addtogroup genfit
 * @{
 */

#ifndef genfit_PlanarMeasurementSDT_h
#define genfit_PlanarMeasurementSDT_h

#include "AbsMeasurement.h"
#include "AbsHMatrix.h"
#include "MeasurementOnPlane.h"

namespace edm4hep{
    class TrackerHit;
    class SimTrackerHit;
}

namespace genfit {

class AbsTrackRep;

/** @brief Measurement class implementing a planar hit geometry (1 or 2D).
 *
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Johannes Rauch  (Technische Universit&auml;t M&uuml;nchen, original author)
 *
 * The main feature of this type of hit is, that the detector plane
 * is defined by the detector hardware. 
 */
class PlanarMeasurementSDT : public AbsMeasurement {

 public:
  PlanarMeasurementSDT(int nDim = 1);
  PlanarMeasurementSDT(const TVectorD& rawHitCoords, const TMatrixDSym& rawHitCov, int detId, int hitId, TrackPoint* trackPoint);

  virtual ~PlanarMeasurementSDT() {;}

  virtual AbsMeasurement* clone() const override {return new PlanarMeasurementSDT(*this);}

  int getPlaneId() const {return planeId_;}

  virtual SharedPlanePtr constructPlane(const StateOnPlane& state) const override;

  virtual std::vector<MeasurementOnPlane*> constructMeasurementsOnPlane(const StateOnPlane& state) const override;

  virtual const AbsHMatrix* constructHMatrix(const AbsTrackRep*) const override;

  virtual void setPlane(const SharedPlanePtr& physicalPlane, int planeId = -1) {physicalPlane_ = physicalPlane; planeId_ = planeId;}

  /** @brief Use if the coordinate for 1D hits measured in V direction.
   *
   * Per default for 1D planar hits, the coordinate is measured in U direction.
   * With this function you can set it to be measured in V direction.
   * This affects the outcoe of constructHMatrix().
   */
  void setStripV(bool v = true) {stripV_ = v;}

  void setTrackerHit(const edm4hep::TrackerHit* hit){trackerHit_=hit;}
  const edm4hep::TrackerHit* getTrackerHit(){return trackerHit_;}

 protected:
  SharedPlanePtr physicalPlane_;   //! This is persistent, but '!' makes ROOT shut up.
  int planeId_; // planeId id is -1 per default
  bool stripV_;
  const edm4hep::SimTrackerHit* simTrackerHit_;
  const edm4hep::TrackerHit* trackerHit_;

 public:

  ClassDefOverride(PlanarMeasurementSDT,1)
};

} /* End of namespace genfit */
/** @} */

#endif // genfit_PlanarMeasurementSDT_h
