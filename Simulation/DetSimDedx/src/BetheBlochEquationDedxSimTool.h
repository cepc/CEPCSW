#ifndef BetheBlochEquationDedxSimTool_h
#define BetheBlochEquationDedxSimTool_h

#include "DetSimInterface/IDedxSimTool.h"
#include <GaudiKernel/AlgTool.h>
#include "edm4hep/MCParticle.h"

class BetheBlochEquationDedxSimTool: public extends<AlgTool, IDedxSimTool> {
    public:
        using extends::extends;

        StatusCode initialize() override;
        StatusCode finalize() override;
        double dedx(const G4Step* aStep) override;
        double dedx(const edm4hep::MCParticle& mc) override;
        double dndx(double betagamma) override;

    private:

        Gaudi::Property<float> m_material_Z{this, "material_Z", 2};//Default is Helium
        Gaudi::Property<float> m_material_A{this, "material_A", 4};
        Gaudi::Property<float> m_material_density{this, "material_density", 0.000178};//g/cm^3
        Gaudi::Property<float> m_scale{this, "scale", 1};
        Gaudi::Property<float> m_resolution{this, "resolution", 0};
        float m_me;// Here me is the electron rest mass
        float m_K; // K was set as a constant.
        float m_I; // Mean excitation energy
};

#endif
