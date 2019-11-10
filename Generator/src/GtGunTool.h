#ifndef GtGunTool_h
#define GtGunTool_h

/*
 * Description: 
 *   A particle gun to generate particles.
 *   User could specify following:
 *   * PDGID or Particle Name
 *   * Status: this is used for extension.
 *   * Momentum or TotalEnergy or KineticEnergy
 *   * Position and Time
 */

#include <GaudiKernel/AlgTool.h>
#include "IGenTool.h"

class GtGunTool: public extends<AlgTool, IGenTool> {
public:
    using extends::extends;

    // Overriding initialize and finalize
    StatusCode initialize() override;
    StatusCode finalize() override;

    // IGenTool
    bool mutate(MyHepMC::GenEvent& event) override;
    bool finish() override;
    bool configure_gentool() override;

};


#endif
