#include "DummyFastSimG4Tool.h"

#include "G4Region.hh"
#include "G4RegionStore.hh"

#include "G4VFastSimulationModel.hh"
#include "DummyFastSimG4Model.h"

DECLARE_COMPONENT(DummyFastSimG4Tool)

StatusCode DummyFastSimG4Tool::initialize() {
    StatusCode sc;

    return sc;
}

StatusCode DummyFastSimG4Tool::finalize() {
    StatusCode sc;

    return sc;
}

bool DummyFastSimG4Tool::CreateFastSimulationModel() {
    // In this method:
    // * Retrieve the G4Region
    // * Create Model
    // * Associate model and region

    G4String model_name = "DummyFastSimG4Model";
    for (auto region_name: m_regions.value()) {
        G4Region* aEnvelope = G4RegionStore::GetInstance()->GetRegion(region_name);
        if (!aEnvelope) {
            error() << "Failed to find G4Region '" << region_name << "'" << endmsg;
            return false;
        }

        DummyFastSimG4Model* model = new DummyFastSimG4Model(model_name+region_name, aEnvelope);
        info() << "Create Model " << model_name << " for G4Region " << region_name << endmsg;
    }
    return true;
}
