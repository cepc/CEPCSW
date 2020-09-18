#ifndef BetheBlochEquationDedxSimTool_h
#define BetheBlochEquationDedxSimTool_h

#include "DetSimInterface/IDedxSimTool.h"
#include <GaudiKernel/AlgTool.h>
#include <random>

class BetheBlochEquationDedxSimTool: public extends<AlgTool, IDedxSimTool> {
    public:
        using extends::extends;

        StatusCode initialize() override;
        StatusCode finalize() override;
        double dedx(const G4Step* aStep) override;

    private:

        Gaudi::Property<float> m_material_Z{this, "material_Z", 2};//Default is Helium
        Gaudi::Property<float> m_material_A{this, "material_A", 4};
        Gaudi::Property<float> m_material_density{this, "material_density", 0.000178};//g/cm^3
        Gaudi::Property<float> m_scale{this, "scale", 1};
        Gaudi::Property<float> m_resolution{this, "resolution", 0.01};
        float m_me;// Here me is the electron rest mass
        float m_K; // K was set as a constant.
        float m_I; // Mean excitation energy

        std::default_random_engine m_generator;
        std::normal_distribution<double>* m_distribution;
};

#endif
