#ifndef GEARTPCKALDETECTOR_H
#define GEARTPCKALDETECTOR_H

#include "kaltest/TVKalDetector.h"

#include "GearTPCMeasLayer.h"

#include <map>

namespace gear{
  class GearMgr ;
}

namespace kaldet{

  /**
   * The LCTPC implementation for a TPC which is completely instantiated from GEAR.
   * 
   */
class GearTPCKalDetector : public TVKalDetector {

public:
    /** 
     * The constructor. All information to initialise the TPC is taken from GEAR.
     *
     * As a pragmatic approach to avoid dealing with conditions data and material databases,
     * the information about the material budget and the resolution of the layers
     * is taken from the GEAR file as user parameters. If the parameters are not found in the
     * file the previously hard coded parameters are used as default, which ensures backward
     * compatibility.
     *
     * The gas properties for the matrial budget can be given as user parameters 
     * for the TPCParameters:
     * \param  TPCGas_A The mean atomic mass (default 36.2740552)
     * \param  TPCGas_Z The mean number of protons (default 16.4)
     * \param  TPCGas_density The density (default 0.749e-3 in which units?)
     * \param  TPCGas_radlen The radiation length (default 2.392e4 in which units?)
     *
     * The default gas parameters (are supposed to) correspond to Ar/CH4 90/10.
     * N.B.: KILLENB: I think there is a bug in the calculation, the mean A should be
     * 37.6 instead of 36.3 (see source code).
     * In addition the description as a single TMaterial is not good. 
     * Using TMixture would be better.
     *
     * The reslution is calculated as \f$\sigma_x = \sqrt{x_0^2 + x_1^2 \cdot z}\f$.
     * This requires z to be proportional to the drift distance, i.\ e. z=0 is at the readout.

     * The resolution of the layers can be given as user parameters in each TPCModule 
     * section of the GEAR xml file.
     * \param sigmax0 The constant part of the x resolution (default 38.3e-3 mm)
     * \param sigmax1 The drift distance dependent part of the x resolution 
     *                (default 6.74e-3 mm/sqrt(mm) )
     * \param sigmaz0 The constant part of the z resolution (default 0.5 mm)
     * \param sigmaz1 The drift distance dependent part the z resolution
     *                (default 10.2e-3 mm/sqrt(mm) )
     */
    GearTPCKalDetector(const gear::GearMgr& gearMgr);

    /// The destructor.
    virtual ~GearTPCKalDetector();

    /**
     * Get access to the measurement layers using moduleID and row.
     * Do not directly access the measurement layers using At() 
     * because the order depends on the order in the gear file.
     * Throws a gear::Exception if the row on the module is not defined.
     */
    virtual GearTPCMeasLayer const * GetMeasLayer(int moduleID, int row) const;

protected:
    /// Map which contains the information which measurement layer is stored
    /// at which position in the array.
    std::map< std::pair<int, int >, Int_t > moduleRowToMeasurementLayerMap;
};

}// namespace kaldet
#endif //GEARTPCKALDETECTOR_H
