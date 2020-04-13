#include "ExampleAnaElemTool.h"

#include "G4Event.hh"
#include "G4THitsCollection.hh"

#include "DD4hep/Detector.h"
#include "DD4hep/Plugins.h"
#include "DDG4/Geant4Converter.h"
#include "DDG4/Geant4Mapping.h"
#include "DDG4/Geant4HitCollection.h"
#include "DDG4/Geant4Data.h"
#include "DDG4/Geant4Hits.h"

DECLARE_COMPONENT(ExampleAnaElemTool)

void
ExampleAnaElemTool::BeginOfRunAction(const G4Run*) {
    G4cout << "Begin Run of detector simultion..." << G4endl;
}

void
ExampleAnaElemTool::EndOfRunAction(const G4Run*) {
    G4cout << "End Run of detector simultion..." << G4endl;
}

void
ExampleAnaElemTool::BeginOfEventAction(const G4Event* anEvent) {
    msg() << "Event " << anEvent->GetEventID() << endmsg;
}

void
ExampleAnaElemTool::EndOfEventAction(const G4Event* anEvent) {
    auto mcCol = m_mcParCol.get();
    msg() << "mcCol size: " << mcCol->size() << endmsg;
    // save all data

    // create collections.
    auto trackercols = m_trackerCol.createAndPut();
    auto calorimetercols = m_calorimeterCol.createAndPut();
    auto calocontribcols = m_caloContribCol.createAndPut();

    auto vxdcols = m_VXDCol.createAndPut();
    auto ftdcols = m_FTDCol.createAndPut();
    auto sitcols = m_SITCol.createAndPut();
    auto tpccols = m_TPCCol.createAndPut();
    auto setcols = m_SETCol.createAndPut();

    // readout defined in DD4hep
    auto lcdd = &(dd4hep::Detector::getInstance());
    auto allReadouts = lcdd->readouts();

    for (auto& readout : allReadouts) {
        info() << "Readout " << readout.first << endmsg;
    }

    // retrieve the hit collections
    G4HCofThisEvent* collections = anEvent->GetHCofThisEvent();
    if (!collections) {
        warning() << "No collections found. " << endmsg;
        return;
    }
    int Ncol = collections->GetNumberOfCollections();
    for (int icol = 0; icol < Ncol; ++icol) {
        G4VHitsCollection* collect = collections->GetHC(icol);
        if (!collect) {
            warning() << "Collection iCol " << icol << " is missing" << endmsg;
            continue;
        }
        size_t nhits = collect->GetSize();
        info() << "Collection " << collect->GetName()
               << " #" << icol
               << " has " << nhits << " hits."
               << endmsg;
        if (nhits==0) {
            // just skip this collection.
            continue;
        }

        plcio::SimTrackerHitCollection* tracker_col_ptr = nullptr;
        plcio::SimCalorimeterHitCollection* calo_col_ptr = nullptr;

        // the mapping between hit collection and the data handler
        if (collect->GetName() == "VXDCollection") {
            tracker_col_ptr = vxdcols;
        } else if (collect->GetName() == "FTDCollection") {
            tracker_col_ptr = ftdcols;
        } else if (collect->GetName() == "SITCollection") {
            tracker_col_ptr = sitcols;
        } else if (collect->GetName() == "TPCCollection") {
            tracker_col_ptr = tpccols;
        } else if (collect->GetName() == "SETCollection") {
            tracker_col_ptr = setcols;
        } else if (collect->GetName() == "SETCollection") {
            tracker_col_ptr = setcols;
        } else if (collect->GetName() == "CaloHitsCollection") {
            calo_col_ptr = calorimetercols;
        } else {
            warning() << "Unknown collection name: " << collect->GetName()
                      << ". The SimTrackerCol will be used. " << endmsg;
            tracker_col_ptr = trackercols;
        }



        // There are different types (new and old)

        dd4hep::sim::Geant4HitCollection* coll = dynamic_cast<dd4hep::sim::Geant4HitCollection*>(collect);
        if (coll) {
            info() << " cast to dd4hep::sim::Geant4HitCollection. " << endmsg;
            for(size_t i=0; i<nhits; ++i)   {

                dd4hep::sim::Geant4HitData* h = coll->hit(i);

                dd4hep::sim::Geant4Tracker::Hit* trk_hit = dynamic_cast<dd4hep::sim::Geant4Tracker::Hit*>(h);
                if ( 0 != trk_hit )   {
                    dd4hep::sim::Geant4HitData::Contribution& t = trk_hit->truth;
                    int trackID = t.trackID;
                    // t.trackID = m_truth->particleID(trackID);
                }
                // Geant4Calorimeter::Hit* cal_hit = dynamic_cast<Geant4Calorimeter::Hit*>(h);
                // if ( 0 != cal_hit )   {
                //     Geant4HitData::Contributions& c = cal_hit->truth;
                //     for(Geant4HitData::Contributions::iterator j=c.begin(); j!=c.end(); ++j)  {
                //         Geant4HitData::Contribution& t = *j;
                //         int trackID = t.trackID;
                //         // t.trackID = m_truth->particleID(trackID);
                //     }
                // }
            }
            continue;
        }

        typedef G4THitsCollection<dd4hep::sim::Geant4Hit> HitCollection;
        HitCollection* coll2 = dynamic_cast<HitCollection*>(collect);

        if (coll2) {
            info() << " cast to G4THitsCollection<dd4hep::sim::Geant4Hit>. " << endmsg;

            int n_trk_hit = 0;
            int n_cal_hit = 0;

            for(size_t i=0; i<nhits; ++i)   {
                dd4hep::sim::Geant4Hit* h = dynamic_cast<dd4hep::sim::Geant4Hit*>(coll2->GetHit(i));
                if (!h) {
                    warning() << "Failed to cast to dd4hep::sim::Geant4Hit. " << endmsg;
                    continue;
                }

                dd4hep::sim::Geant4TrackerHit* trk_hit = dynamic_cast<dd4hep::sim::Geant4TrackerHit*>(h);
                if (trk_hit) {
                    ++n_trk_hit;
                    // auto edm_trk_hit = trackercols->create();
                    auto edm_trk_hit = tracker_col_ptr->create();

                    // Refer to: ./DDG4/lcio/LCIOConversions.cpp
                    edm_trk_hit->setCellID0((trk_hit->cellID >>    0         ) & 0xFFFFFFFF);
                    edm_trk_hit->setCellID1((trk_hit->cellID >> sizeof(int)*8) & 0xFFFFFFFF);
                    edm_trk_hit->setEDep(trk_hit->energyDeposit/CLHEP::GeV);
                    edm_trk_hit->setTime(trk_hit->truth.time/CLHEP::ns);
                    edm_trk_hit->setPathLength(trk_hit->length/CLHEP::mm);
                    // lc_hit->setMCParticle(lc_mcp);
                    double pos[3] = {trk_hit->position.x()/CLHEP::mm,
                                     trk_hit->position.y()/CLHEP::mm,
                                     trk_hit->position.z()/CLHEP::mm};
                    edm_trk_hit->setPosition(plcio::DoubleThree(pos));

                    float mom[3] = {trk_hit->momentum.x()/CLHEP::GeV,
                                    trk_hit->momentum.y()/CLHEP::GeV,
                                    trk_hit->momentum.z()/CLHEP::GeV};
                    edm_trk_hit->setMomentum(plcio::FloatThree(mom));
                }

                dd4hep::sim::Geant4CalorimeterHit* cal_hit = dynamic_cast<dd4hep::sim::Geant4CalorimeterHit*>(h);
                if (cal_hit) {
                    ++n_cal_hit;
                    auto edm_calo_hit = calo_col_ptr->create();
                    edm_calo_hit->setCellID0((cal_hit->cellID >> 0            ) & 0xFFFFFFFF);
                    edm_calo_hit->setCellID1((cal_hit->cellID >> sizeof(int)*8) & 0xFFFFFFFF);
                    edm_calo_hit->setEnergy(cal_hit->energyDeposit);
                    float pos[3] = {cal_hit->position.x()/CLHEP::mm,
                                    cal_hit->position.y()/CLHEP::mm,
                                    cal_hit->position.z()/CLHEP::mm};
                    edm_calo_hit->setPosition(plcio::FloatThree(pos));

                    // contribution
                    typedef dd4hep::sim::Geant4CalorimeterHit::Contributions Contributions;
                    typedef dd4hep::sim::Geant4CalorimeterHit::Contribution Contribution;
                    for (Contributions::const_iterator j = cal_hit->truth.begin();
                         j != cal_hit->truth.end(); ++j) {
                        const Contribution& c = *j;
                        // The legacy Hit object does not contains positions of contributions.
                        // float contrib_pos[] = {float(c.x/mm), float(c.y/mm), float(c.z/mm)};
                        auto edm_calo_contrib = calocontribcols->create();
                        edm_calo_contrib.setPDG(c.pdgID);
                        edm_calo_contrib.setEnergy(c.deposit/CLHEP::GeV);
                        edm_calo_contrib.setTime(c.time/CLHEP::ns);
                        edm_calo_contrib.setStepPosition(plcio::FloatThree(pos));
                        edm_calo_contrib.setParticle(mcCol->at(0));
                        edm_calo_hit->addContribution(edm_calo_contrib);
                    }
                }

            }

            info() << n_trk_hit << " hits cast to dd4hep::sim::Geant4TrackerHit. " << endmsg;
            info() << n_cal_hit << " hits cast to dd4hep::sim::Geant4CalorimeterHit. " << endmsg;


            continue;
        }

        warning() << "Failed to convert to collection "
                  << collect->GetName()
                  << endmsg;
        
    }
}

void
ExampleAnaElemTool::PreUserTrackingAction(const G4Track*) {

}

void
ExampleAnaElemTool::PostUserTrackingAction(const G4Track*) {

}

void
ExampleAnaElemTool::UserSteppingAction(const G4Step*) {

}

StatusCode
ExampleAnaElemTool::initialize() {
    StatusCode sc;

    return sc;
}

StatusCode
ExampleAnaElemTool::finalize() {
    StatusCode sc;

    return sc;
}


