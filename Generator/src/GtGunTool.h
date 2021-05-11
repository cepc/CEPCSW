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
#include <Gaudi/Property.h>
#include "IGenTool.h"

#include <vector>

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

private:

    Gaudi::Property<std::vector<std::string>> m_particles{this, "Particles"};

    Gaudi::Property<std::vector<double>> m_positionXs{this, "PositionXs"};
    Gaudi::Property<std::vector<double>> m_positionYs{this, "PositionYs"};
    Gaudi::Property<std::vector<double>> m_positionZs{this, "PositionZs"};


    Gaudi::Property<std::vector<double>> m_energymins{this, "EnergyMins"};
    Gaudi::Property<std::vector<double>> m_energymaxs{this, "EnergyMaxs"};

    Gaudi::Property<std::vector<double>> m_thetamins{this, "ThetaMins"};
    Gaudi::Property<std::vector<double>> m_thetamaxs{this, "ThetaMaxs"};

    Gaudi::Property<std::vector<double>> m_phimins{this, "PhiMins"};
    Gaudi::Property<std::vector<double>> m_phimaxs{this, "PhiMaxs"};

};


#endif
