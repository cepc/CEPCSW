#ifndef Edm4hepWriterAnaElemTool_h
#define Edm4hepWriterAnaElemTool_h

#include "GaudiKernel/AlgTool.h"
#include "FWCore/DataHandle.h"
#include "DetSimInterface/IAnaElemTool.h"

#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/CaloHitContributionCollection.h"

class Edm4hepWriterAnaElemTool: public extends<AlgTool, IAnaElemTool> {

public:

    using extends::extends;

    /// IAnaElemTool interface
    // Run
    virtual void BeginOfRunAction(const G4Run*) override;
    virtual void EndOfRunAction(const G4Run*) override;

    // Event
    virtual void BeginOfEventAction(const G4Event*) override;
    virtual void EndOfEventAction(const G4Event*) override;

    // Tracking
    virtual void PreUserTrackingAction(const G4Track*) override;
    virtual void PostUserTrackingAction(const G4Track*) override;

    // Stepping
    virtual void UserSteppingAction(const G4Step*) override;


    /// Overriding initialize and finalize
    StatusCode initialize() override;
    StatusCode finalize() override;

private:
    // In order to associate MCParticle with contribution, we need to access MC Particle.
    DataHandle<edm4hep::MCParticleCollection> m_mcParCol{"MCParticle", 
            Gaudi::DataHandle::Writer, this};

    // Generic collections for Tracker and Calorimeter
    DataHandle<edm4hep::SimTrackerHitCollection> m_trackerCol{"SimTrackerCol", 
            Gaudi::DataHandle::Writer, this};
    DataHandle<edm4hep::SimCalorimeterHitCollection> m_calorimeterCol{"SimCalorimeterCol", 
            Gaudi::DataHandle::Writer, this};
    DataHandle<edm4hep::CaloHitContributionCollection> m_caloContribCol{"SimCaloContributionCol", 
            Gaudi::DataHandle::Writer, this};

    // Dedicated collections for CEPC
    DataHandle<edm4hep::SimTrackerHitCollection> m_VXDCol{"VXDCollection", 
            Gaudi::DataHandle::Writer, this};
    DataHandle<edm4hep::SimTrackerHitCollection> m_FTDCol{"FTDCollection", 
            Gaudi::DataHandle::Writer, this};
    DataHandle<edm4hep::SimTrackerHitCollection> m_SITCol{"SITCollection", 
            Gaudi::DataHandle::Writer, this};
    DataHandle<edm4hep::SimTrackerHitCollection> m_TPCCol{"TPCCollection", 
            Gaudi::DataHandle::Writer, this};
    DataHandle<edm4hep::SimTrackerHitCollection> m_SETCol{"SETCollection", 
            Gaudi::DataHandle::Writer, this};

};

#endif
