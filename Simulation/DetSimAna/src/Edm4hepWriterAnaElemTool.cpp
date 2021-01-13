#include "Edm4hepWriterAnaElemTool.h"

#include "G4Event.hh"
#include "G4THitsCollection.hh"
#include "G4EventManager.hh"
#include "G4TrackingManager.hh"
#include "G4SteppingManager.hh"

#include "DD4hep/Detector.h"
#include "DD4hep/Plugins.h"
#include "DDG4/Geant4Converter.h"
#include "DDG4/Geant4Mapping.h"
#include "DDG4/Geant4HitCollection.h"
#include "DDG4/Geant4Data.h"
#include "DDG4/Geant4Hits.h"

DECLARE_COMPONENT(Edm4hepWriterAnaElemTool)

void
Edm4hepWriterAnaElemTool::BeginOfRunAction(const G4Run*) {
    G4cout << "Begin Run of detector simultion..." << G4endl;
}

void
Edm4hepWriterAnaElemTool::EndOfRunAction(const G4Run*) {
    G4cout << "End Run of detector simultion..." << G4endl;
}

void
Edm4hepWriterAnaElemTool::BeginOfEventAction(const G4Event* anEvent) {
    msg() << "Event " << anEvent->GetEventID() << endmsg;

    // reset
    m_track2primary.clear();

}

void
Edm4hepWriterAnaElemTool::EndOfEventAction(const G4Event* anEvent) {
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

    auto ecalbarrelcol            = m_EcalBarrelCol.createAndPut();
    auto ecalbarrelcontribcols    = m_EcalBarrelContributionCol.createAndPut();
    auto ecalendcapscol           = m_EcalEndcapsCol.createAndPut();
    auto ecalendcapscontribcols   = m_EcalEndcapsContributionCol.createAndPut();
    auto ecalendcapringcol        = m_EcalEndcapRingCol.createAndPut();
    auto ecalendcapringcontribcol = m_EcalEndcapRingContributionCol.createAndPut();

    auto hcalbarrelcol            = m_HcalBarrelCol.createAndPut();
    auto hcalbarrelcontribcols    = m_HcalBarrelContributionCol.createAndPut();
    auto hcalendcapscol           = m_HcalEndcapsCol.createAndPut();
    auto hcalendcapscontribcols   = m_HcalEndcapsContributionCol.createAndPut();
    auto hcalendcapringcol        = m_HcalEndcapRingCol.createAndPut();
    auto hcalendcapringcontribcol = m_HcalEndcapRingContributionCol.createAndPut();

    auto coilcols = m_COILCol.createAndPut();

    auto muonbarrelcol            = m_MuonBarrelCol.createAndPut();
    auto muonbarrelcontribcols    = m_MuonBarrelContributionCol.createAndPut();
    auto muonendcapscol           = m_MuonEndcapsCol.createAndPut();
    auto muonendcapscontribcols   = m_MuonEndcapsContributionCol.createAndPut();

    auto driftchamberhitscol = m_DriftChamberHitsCol.createAndPut();

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

        edm4hep::SimTrackerHitCollection* tracker_col_ptr = nullptr;
        edm4hep::SimCalorimeterHitCollection* calo_col_ptr = nullptr;
        edm4hep::CaloHitContributionCollection* calo_contrib_col_ptr = nullptr;

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
            calo_contrib_col_ptr = calocontribcols;
        } else if (collect->GetName() == "EcalBarrelCollection") {
            calo_col_ptr = ecalbarrelcol;
            calo_contrib_col_ptr = ecalbarrelcontribcols;
        } else if (collect->GetName() == "EcalEndcapsCollection") {
            calo_col_ptr = ecalendcapscol;
            calo_contrib_col_ptr = ecalendcapscontribcols;
        } else if (collect->GetName() == "EcalEndcapRingCollection") {
            calo_col_ptr = ecalendcapringcol;
            calo_contrib_col_ptr = ecalendcapringcontribcol;
        } else if (collect->GetName() == "HcalBarrelCollection") {
            calo_col_ptr = hcalbarrelcol;
            calo_contrib_col_ptr = hcalbarrelcontribcols;
        } else if (collect->GetName() == "HcalEndcapsCollection") {
            calo_col_ptr = hcalendcapscol;
            calo_contrib_col_ptr = hcalendcapscontribcols;
        } else if (collect->GetName() == "HcalEndcapRingCollection") {
            calo_col_ptr = hcalendcapringcol;
            calo_contrib_col_ptr = hcalendcapringcontribcol;
	} else if (collect->GetName() == "COILCollection") {
	    tracker_col_ptr = coilcols;
	} else if (collect->GetName() == "MuonBarrelCollection") {
	    calo_col_ptr = muonbarrelcol;
	    calo_contrib_col_ptr = muonbarrelcontribcols;
        } else if (collect->GetName() == "MuonEndcapsCollection") {
	    calo_col_ptr = muonendcapscol;
	    calo_contrib_col_ptr = muonendcapscontribcols;
        } else if (collect->GetName() == "DriftChamberHitsCollection") {
            tracker_col_ptr = driftchamberhitscol;
        } else {
            warning() << "Unknown collection name: " << collect->GetName()
                      << ". Please register in Edm4hepWriterAnaElemTool. " << endmsg;
            continue;
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
                if (trk_hit && tracker_col_ptr) {
                    ++n_trk_hit;
                    // auto edm_trk_hit = trackercols->create();
                    auto edm_trk_hit = tracker_col_ptr->create();

                    edm_trk_hit.setCellID(trk_hit->cellID);
                    edm_trk_hit.setEDep(trk_hit->energyDeposit/CLHEP::GeV);
                    edm_trk_hit.setTime(trk_hit->truth.time/CLHEP::ns);
                    edm_trk_hit.setPathLength(trk_hit->length/CLHEP::mm);
                    // lc_hit->setMCParticle(lc_mcp);
                    double pos[3] = {trk_hit->position.x()/CLHEP::mm,
                                     trk_hit->position.y()/CLHEP::mm,
                                     trk_hit->position.z()/CLHEP::mm};
                    edm_trk_hit.setPosition(edm4hep::Vector3d(pos));

                    float mom[3] = {trk_hit->momentum.x()/CLHEP::GeV,
                                    trk_hit->momentum.y()/CLHEP::GeV,
                                    trk_hit->momentum.z()/CLHEP::GeV};
                    edm_trk_hit.setMomentum(edm4hep::Vector3f(mom));

                    // get the truth or contribution
                    auto& truth = trk_hit->truth;
                    int trackID = truth.trackID;
                    
                    int pritrkid = m_track2primary[trackID];
                    if (pritrkid <= 0) {
                        error() << "Failed to find the primary track for trackID #" << trackID << endmsg;
                        pritrkid = 1;
                    }

                    edm_trk_hit.setMCParticle(mcCol->at(pritrkid-1));

                    if (pritrkid != trackID) {
                        // If the track is a secondary, then the primary track id and current track id is different
                        edm_trk_hit.setProducedBySecondary(true);
                    }
                }

                dd4hep::sim::Geant4CalorimeterHit* cal_hit = dynamic_cast<dd4hep::sim::Geant4CalorimeterHit*>(h);
                if (cal_hit && calo_col_ptr) {
                    ++n_cal_hit;
                    auto edm_calo_hit = calo_col_ptr->create();
                    edm_calo_hit.setCellID(cal_hit->cellID);
                    edm_calo_hit.setEnergy(cal_hit->energyDeposit/CLHEP::GeV);
                    float pos[3] = {cal_hit->position.x()/CLHEP::mm,
                                    cal_hit->position.y()/CLHEP::mm,
                                    cal_hit->position.z()/CLHEP::mm};
                    edm_calo_hit.setPosition(edm4hep::Vector3f(pos));

                    // contribution
                    typedef dd4hep::sim::Geant4CalorimeterHit::Contributions Contributions;
                    typedef dd4hep::sim::Geant4CalorimeterHit::Contribution Contribution;
                    for (Contributions::const_iterator j = cal_hit->truth.begin();
                         j != cal_hit->truth.end(); ++j) {
                        const Contribution& c = *j;
                        // The legacy Hit object does not contains positions of contributions.
                        // float contrib_pos[] = {float(c.x/mm), float(c.y/mm), float(c.z/mm)};
                        auto edm_calo_contrib = calo_contrib_col_ptr->create();
                        edm_calo_contrib.setPDG(c.pdgID);
                        edm_calo_contrib.setEnergy(c.deposit/CLHEP::GeV);
                        edm_calo_contrib.setTime(c.time/CLHEP::ns);
                        edm_calo_contrib.setStepPosition(edm4hep::Vector3f(pos));

                        // from the track id, get the primary track
                        int pritrkid = m_track2primary[c.trackID];
                        if (pritrkid<=0) {
                            error() << "Failed to find the primary track for trackID #" << c.trackID << endmsg;
                            pritrkid = 1;
                        }

                        edm_calo_contrib.setParticle(mcCol->at(pritrkid-1)); // todo
                        edm_calo_hit.addToContributions(edm_calo_contrib);
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
Edm4hepWriterAnaElemTool::PreUserTrackingAction(const G4Track* track) {
    int curtrkid = track->GetTrackID();
    int curparid = track->GetParentID();
    int pritrkid = curparid;

    // try to find the primary track id from the parent track id.
    if (curparid) {
        auto it = m_track2primary.find(curparid);
        if (it == m_track2primary.end()) {
            error() << "Failed to find primary track for track id " << curparid << endmsg;
        } else {
            pritrkid = it->second;
        }
    } else {
        // curparid is 0, it is primary
        pritrkid = curtrkid;
    }


    m_track2primary[curtrkid] = pritrkid;
}

void
Edm4hepWriterAnaElemTool::PostUserTrackingAction(const G4Track* track) {
    int curtrkid = track->GetTrackID(); // starts from 1
    int curparid = track->GetParentID();

    if (curparid == 0) {
        // select the primary tracks (parentID == 0)
        auto mcCol = m_mcParCol.get();

        if (curtrkid-1>=mcCol->size()) {
            error() << "out of range: curtrkid is " << curtrkid
                    << " while the MCParticle size is " << mcCol->size() << endmsg;
            return;
        }
        auto primary_particle = mcCol->at(curtrkid-1);

        const G4ThreeVector& stop_pos  = track->GetPosition();
        edm4hep::Vector3d endpoint(stop_pos.x()/CLHEP::mm,
                                   stop_pos.y()/CLHEP::mm,
                                   stop_pos.z()/CLHEP::mm);
        primary_particle.setEndpoint(endpoint);

        const G4ThreeVector& stop_mom = track->GetMomentum();

        edm4hep::Vector3f mom_endpoint(stop_mom.x()/CLHEP::GeV,
                                       stop_mom.y()/CLHEP::GeV,
                                       stop_mom.z()/CLHEP::GeV);
        primary_particle.setMomentumAtEndpoint(mom_endpoint);

        // ===================================================================
        // Update the flags of the primary particles
        // ===================================================================

        // processes
        static G4VProcess* proc_decay = nullptr;
        if (!proc_decay) {
            G4ProcessManager* pm 
                = track->GetDefinition()->GetProcessManager();
            G4int nprocesses = pm->GetProcessListLength();
            G4ProcessVector* pv = pm->GetProcessList();
            
            for(G4int i=0; i<nprocesses; ++i){
                if((*pv)[i]->GetProcessName()=="Decay"){
                    proc_decay = (*pv)[i];
                }
            }

        }

        // flags
        bool is_decay = false;

        G4TrackingManager* tm = G4EventManager::GetEventManager() 
            -> GetTrackingManager();
        G4TrackVector* secondaries = tm->GimmeSecondaries();
        if(secondaries) {


            size_t nSeco = secondaries->size();
            for (size_t i = 0; i < nSeco; ++i) {
                G4Track* sectrk = (*secondaries)[i];
                G4ParticleDefinition* secparticle = sectrk->GetDefinition();
                const G4VProcess* creatorProcess = sectrk->GetCreatorProcess();

                // select the necessary processes
                if (creatorProcess==proc_decay) {
                    info() << "Creator Process is Decay for secondary particle: "
                           << " idx: " << i
                           << " trkid: " << sectrk->GetTrackID() // not valid until track
                           << " particle: " << secparticle->GetParticleName()
                           << " pdg: " << secparticle->GetPDGEncoding()
                           << " at position: " << sectrk->GetPosition() //
                           << " time: " << sectrk->GetGlobalTime()
                           << " momentum: " << sectrk->GetMomentum() // 
                           << endmsg;
                    is_decay = true;

                    // create secondaries in MC particles
                    // todo: convert the const collection to non-const
                    auto mcCol = const_cast<edm4hep::MCParticleCollection*>(m_mcParCol.get());
                    edm4hep::MCParticle mcp = mcCol->create();
                    mcp.setPDG(secparticle->GetPDGEncoding());
                    mcp.setGeneratorStatus(0); // not created by Generator
                    mcp.setCreatedInSimulation(1);
                    mcp.setCharge(secparticle->GetPDGCharge());
                    mcp.setTime(sectrk->GetGlobalTime()/CLHEP::ns); // todo
                    mcp.setMass(secparticle->GetPDGMass());

                    const G4ThreeVector& sec_init_pos = sectrk->GetPosition();
                    double x=sec_init_pos.x()/CLHEP::mm;
                    double y=sec_init_pos.y()/CLHEP::mm;
                    double z=sec_init_pos.z()/CLHEP::mm;

                    const G4ThreeVector& sec_init_mom = sectrk->GetMomentum();
                    double px=sec_init_mom.x()/CLHEP::GeV;
                    double py=sec_init_mom.y()/CLHEP::GeV;
                    double pz=sec_init_mom.z()/CLHEP::GeV;
                    mcp.setVertex(edm4hep::Vector3d(x,y,z)); // todo
                    mcp.setEndpoint(edm4hep::Vector3d(x,y,z)); // todo
                    mcp.setMomentum(edm4hep::Vector3f(px,py,pz)); // todo
                    mcp.setMomentumAtEndpoint(edm4hep::Vector3f(px,py,pz)); //todo

                    mcp.addToParents(primary_particle);
                    primary_particle.addToDaughters(mcp);
                }
            }
        }

        // now update
        if (is_decay) {
            primary_particle.setDecayedInTracker(is_decay);
            primary_particle.setDecayedInCalorimeter(is_decay);
        }

    } else {
        // TODO: select other interested tracks
    }

}

void
Edm4hepWriterAnaElemTool::UserSteppingAction(const G4Step*) {

}

StatusCode
Edm4hepWriterAnaElemTool::initialize() {
    StatusCode sc;

    return sc;
}

StatusCode
Edm4hepWriterAnaElemTool::finalize() {
    StatusCode sc;

    return sc;
}


