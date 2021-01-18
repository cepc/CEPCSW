#ifndef Edm4hepWriterAnaElemTool_h
#define Edm4hepWriterAnaElemTool_h

#include <map>

#include "GaudiKernel/AlgTool.h"
#include "k4FWCore/DataHandle.h"
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

    // Ecal
    DataHandle<edm4hep::SimCalorimeterHitCollection> m_EcalBarrelCol{"EcalBarrelCollection", 
            Gaudi::DataHandle::Writer, this};
    DataHandle<edm4hep::CaloHitContributionCollection> m_EcalBarrelContributionCol{
            "EcalBarrelContributionCollection", 
            Gaudi::DataHandle::Writer, this};
    DataHandle<edm4hep::SimCalorimeterHitCollection> m_EcalEndcapsCol{"EcalEndcapsCollection", 
            Gaudi::DataHandle::Writer, this};
    DataHandle<edm4hep::CaloHitContributionCollection> m_EcalEndcapsContributionCol{
            "EcalEndcapsContributionCollection", 
            Gaudi::DataHandle::Writer, this};
    DataHandle<edm4hep::SimCalorimeterHitCollection> m_EcalEndcapRingCol{"EcalEndcapRingCollection", 
            Gaudi::DataHandle::Writer, this};
    DataHandle<edm4hep::CaloHitContributionCollection> m_EcalEndcapRingContributionCol{
            "EcalEndcapRingContributionCollection", 
            Gaudi::DataHandle::Writer, this};

    // Hcal
    DataHandle<edm4hep::SimCalorimeterHitCollection> m_HcalBarrelCol{"HcalBarrelCollection", 
            Gaudi::DataHandle::Writer, this};
    DataHandle<edm4hep::CaloHitContributionCollection> m_HcalBarrelContributionCol{
            "HcalBarrelContributionCollection", 
            Gaudi::DataHandle::Writer, this};
    DataHandle<edm4hep::SimCalorimeterHitCollection> m_HcalEndcapsCol{"HcalEndcapsCollection", 
            Gaudi::DataHandle::Writer, this};
    DataHandle<edm4hep::CaloHitContributionCollection> m_HcalEndcapsContributionCol{
            "HcalEndcapsContributionCollection", 
            Gaudi::DataHandle::Writer, this};
    DataHandle<edm4hep::SimCalorimeterHitCollection> m_HcalEndcapRingCol{"HcalEndcapRingCollection", 
            Gaudi::DataHandle::Writer, this};
    DataHandle<edm4hep::CaloHitContributionCollection> m_HcalEndcapRingContributionCol{
            "HcalEndcapRingContributionCollection", 
            Gaudi::DataHandle::Writer, this};

    // Coil
    DataHandle<edm4hep::SimTrackerHitCollection> m_COILCol{"COILCollection",
	Gaudi::DataHandle::Writer, this};

    // Muon
    DataHandle<edm4hep::SimCalorimeterHitCollection> m_MuonBarrelCol{"MuonBarrelCollection",
	Gaudi::DataHandle::Writer, this};
    DataHandle<edm4hep::CaloHitContributionCollection> m_MuonBarrelContributionCol{
      "MuonBarrelContributionCollection",
	Gaudi::DataHandle::Writer, this};
    DataHandle<edm4hep::SimCalorimeterHitCollection> m_MuonEndcapsCol{"MuonEndcapsCollection",
	Gaudi::DataHandle::Writer, this};
    DataHandle<edm4hep::CaloHitContributionCollection> m_MuonEndcapsContributionCol{
      "MuonEndcapsContributionCollection",
	Gaudi::DataHandle::Writer, this};

    // Drift Chamber
    // - DriftChamberHitsCollection
    DataHandle<edm4hep::SimTrackerHitCollection> m_DriftChamberHitsCol{
            "DriftChamberHitsCollection", 
            Gaudi::DataHandle::Writer, this};


private:
    // in order to associate the hit contribution with the primary track,
    // we have a bookkeeping of every track.
    // The primary track will assign the same key/value.

    // Following is an example:
    //    1 -> 1,
    //    2 -> 2,
    //    3 -> 1,
    // Now, if parent of trk #4 is trk #3, using the mapping {3->1} could 
    // locate the primary trk #1.

    std::map<int, int> m_track2primary;

};

#endif
