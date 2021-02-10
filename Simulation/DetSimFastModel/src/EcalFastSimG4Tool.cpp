#include "EcalFastSimG4Tool.h"

#include "G4Region.hh"
#include "G4RegionStore.hh"

#include "G4VFastSimulationModel.hh"
#include "EcalFastSimG4Model.h"

DECLARE_COMPONENT(EcalFastSimG4Tool);

StatusCode EcalFastSimG4Tool::initialize() {
    StatusCode sc;

    return sc;
}

StatusCode EcalFastSimG4Tool::finalize() {
    StatusCode sc;

    return sc;
}

bool EcalFastSimG4Tool::CreateFastSimulationModel() {
    // In this method:
    // * Retrieve the G4Region
    // * Create Model
    // * Associate model and region

    G4String model_name = "EcalFastSimG4Model";
    for (auto region_name: m_regions.value()) {
        G4Region* aEnvelope = G4RegionStore::GetInstance()->GetRegion(region_name);
        if (!aEnvelope) {
            error() << "Failed to find G4Region '" << region_name << "'" << endmsg;
            return false;
        }

        EcalFastSimG4Model* model = new EcalFastSimG4Model(model_name+region_name, aEnvelope);
        info() << "Create Model " << model_name << " for G4Region " << region_name << endmsg;
    }
    return true;
}
