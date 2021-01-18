/***********************************************************************************\
* (c) Copyright 1998-2019 CERN for the benefit of the LHCb and ATLAS collaborations *
*                                                                                   *
* This software is distributed under the terms of the Apache version 2 licence,     *
* copied verbatim in the file "LICENSE".                                            *
*                                                                                   *
* In applying this licence, CERN does not waive the privileges and immunities       *
* granted to it by virtue of its status as an Intergovernmental Organization        *
* or submit itself to any jurisdiction.                                             *
\***********************************************************************************/
/* The algorithm is created by Mingrui ZHAO.
 * It is transplanted from Marlin to Gaudi by Shunan ZHANG and Mingrui ZHAO.
 * It is now maintained by Mingrui ZHAO (mingrui.zhao@mail.labz0.org)
 *  please contact if you have any question. 
 * 
 * ------------------------------------
 * The algorithm inspects the MCParticles and the correponding tracks.
 * The output of the algorithm is a tuple.
 * In each element of the tuple, it contains:
 *    The information of the MCParticle,
 *    The information of the corresponding tracks.
 *        If the number of reconstructed track candidates is larger than 1, 
 *        each reconstructed track will occupy an element of the tuple, and they
 *        will be differred through the "nCandidates",
 *
 * Establishment of the correspondance:
 * MCParticle -(1)-> SimTrackerHit -(2)-> TrackerHits -(3)-> Track
 * 
 *
 *
 */
// Include files
#include "TrackInspectAlg.h"

#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "DataHelper/HelixClass.h"

#include "CLHEP/Units/SystemOfUnits.h"
#include <math.h>
#include <TTree.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <vector>
#include <iostream>

DECLARE_COMPONENT( TrackInspectAlg )

//------------------------------------------------------------------------------
TrackInspectAlg::TrackInspectAlg( const std::string& name, ISvcLocator* pSvcLocator )
    : Algorithm( name, pSvcLocator ) {
        declareProperty("TrackCollection", _inTrackColHdl, "Handle of the Input Track collection");
        declareProperty("MCParticleCollection", _inMCParticleColHdl, "Handle of the Input MC particle collection");

        declareProperty("TPCTrackerHitRelations", _TPCRelColHdl, "Handle of the TPC Tracker Hit relations");
        declareProperty("VXDTrackerHitRelations", _VXDRelColHdl, "Handle of the TPC Tracker Hit relations");
        declareProperty("SITTrackerHitRelations", _SITRelColHdl, "Handle of the TPC Tracker Hit relations");
        declareProperty("SETTrackerHitRelations", _SETRelColHdl, "Handle of the TPC Tracker Hit relations");
        declareProperty("FTDTrackerHitRelations", _FTDRelColHdl, "Handle of the TPC Tracker Hit relations");

        declareProperty("useTPC", _useTPC, "flag whether to use TPC hits");
        declareProperty("useVXD", _useVXD, "flag whether to use VXD hits");
        declareProperty("useSIT", _useSIT, "flag whether to use SIT hits");
        declareProperty("useSET", _useSET, "flag whether to use SET hits");
        declareProperty("useFTD", _useFTD, "flag whether to use FTD hits");
       
        m_thisName = name;
    }

//------------------------------------------------------------------------------
StatusCode TrackInspectAlg::initialize(){
    info() << "Booking Ntuple" << endmsg;
    
    NTuplePtr nt1(ntupleSvc(), "MyTuples/Track"+m_thisName);
    if ( !nt1 ) {
        m_tuple = ntupleSvc()->book("MyTuples/Track"+m_thisName,CLID_ColumnWiseTuple,"Tracking result");
        if ( 0 != m_tuple ) {
            m_tuple->addItem        ("nmc",                 m_nParticles, 0, 1000 ).ignore();
            m_tuple->addIndexedItem ("vx",                  m_nParticles, vx                 ).ignore();
            m_tuple->addIndexedItem ("vy",                  m_nParticles, vy                 ).ignore();
            m_tuple->addIndexedItem ("vz",                  m_nParticles, vz                 ).ignore();
            m_tuple->addIndexedItem ("ex",                  m_nParticles, ex                 ).ignore();
            m_tuple->addIndexedItem ("ey",                  m_nParticles, ey                 ).ignore();
            m_tuple->addIndexedItem ("ez",                  m_nParticles, ez                 ).ignore();
            m_tuple->addIndexedItem ("Omega",               m_nParticles, Omega              ).ignore();
            m_tuple->addIndexedItem ("D0",                  m_nParticles, D0                 ).ignore();
            m_tuple->addIndexedItem ("Z0",                  m_nParticles, Z0                 ).ignore();
            m_tuple->addIndexedItem ("Phi",                 m_nParticles, Phi                ).ignore();
            m_tuple->addIndexedItem ("TanLambda",           m_nParticles, TanLambda          ).ignore();
            m_tuple->addIndexedItem ("TRUEPX",              m_nParticles, TRUEPX             ).ignore();
            m_tuple->addIndexedItem ("TRUEPY",              m_nParticles, TRUEPY             ).ignore();
            m_tuple->addIndexedItem ("TRUEPZ",              m_nParticles, TRUEPZ             ).ignore();
            m_tuple->addIndexedItem ("TRUEPE",              m_nParticles, TRUEPE             ).ignore();
            m_tuple->addIndexedItem ("TRUEPT",              m_nParticles, TRUEPT             ).ignore();
            m_tuple->addIndexedItem ("TRUEP",               m_nParticles, TRUEP              ).ignore();
            m_tuple->addIndexedItem ("TRUEETA",             m_nParticles, TRUEETA            ).ignore();
            m_tuple->addIndexedItem ("TRUEY",               m_nParticles, TRUEY              ).ignore();
            m_tuple->addIndexedItem ("TRUETHETA",           m_nParticles, TRUETHETA          ).ignore();
            m_tuple->addIndexedItem ("eventNumber",         m_nParticles, eventNumber        ).ignore();
            m_tuple->addIndexedItem ("particleNumber",      m_nParticles, particleNumber     ).ignore();
            m_tuple->addIndexedItem ("totalCandidates",     m_nParticles, totalCandidates    ).ignore();
            m_tuple->addIndexedItem ("nCandidate",          m_nParticles, nCandidate         ).ignore();
            m_tuple->addIndexedItem ("nHits",               m_nParticles, nHits              ).ignore();
            m_tuple->addIndexedItem ("pid",                 m_nParticles, pid                ).ignore();
            m_tuple->addIndexedItem ("WMSelectionVariable", m_nParticles, WMSelectionVariable).ignore();
        }
        else { // did not manage to book the N tuple....
            fatal() << "Cannot book MyTuples/Track"+m_thisName <<endmsg; 
            return StatusCode::FAILURE;
        }
    }
    else{
        m_tuple = nt1;
    }

    _nEvt = 0;
    return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------
StatusCode TrackInspectAlg::execute(){
    debug() << "TrackInspectAlg::execute() ------ start ------"  << endmsg;
    hitmap.clear();
    mcpHitMap.clear();
    matchvec.clear();

    std::vector<const edm4hep::MCRecoTrackerAssociationCollection*> relCols;


    for (auto relCol: relCols) {
    	if (relCol){
	    for (auto rel: *relCol){
		    std::pair<edm4hep::ConstTrackerHit, edm4hep::ConstMCParticle> p = std::make_pair(rel.getRec(), rel.getSim().getMCParticle());
		    if (hitmap.find(p) == hitmap.end()) hitmap[p] = 0.;
		    hitmap[p] += rel.getWeight();
	    }
	}
    }


    // Establish the relation of MCParticle --> Track
    // Put the relation of MCParticle --> Track to matchvec
    const edm4hep::TrackCollection* trkCol = nullptr;
    try {
        trkCol = _inTrackColHdl.get();
    }
    catch ( GaudiException &e ) {
        debug() << "Collection " << _inTrackColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
    }
    const edm4hep::MCParticleCollection* mcpCol = nullptr;
    try {
        mcpCol = _inMCParticleColHdl.get();
    }
    catch ( GaudiException &e ) {
        debug() << "Collection " << _inMCParticleColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
    }
    if (trkCol && mcpCol){
        for (auto track: *trkCol){
            for (auto particle: *mcpCol){
                double match_weight = match(particle, track);
                if (match_weight > 0.2){
                    std::tuple<edm4hep::ConstMCParticle, edm4hep::ConstTrack, double> tuple = std::make_tuple(particle, track, match_weight);
                    matchvec.push_back(tuple);
                }
            }
        }
    }

    if (mcpCol){
        // MCParticleHitAssociator(mcpCol);
        m_nParticles = 0;
        for (auto particle: *mcpCol) {
            std::vector<edm4hep::ConstTrack> theTracks = MCParticleTrackAssociator(particle);

            if (theTracks.size() == 0) {
                totalCandidates[m_nParticles] = 0;
                nCandidate[m_nParticles] = -1;
                Fill(particle, edm4hep::Track(nullptr));
                m_nParticles++;
            }
            else {
                for (unsigned j = 0; j < theTracks.size(); j++) {
                    totalCandidates[m_nParticles] = theTracks.size();
                    nCandidate[m_nParticles] = j;
                    Fill(particle, theTracks[j]);
                    m_nParticles++;
                }   
            }   
        }
        debug() << "MCParticle: " << m_nParticles << endmsg;
    }
    m_tuple->write();
    _nEvt++;
    debug() << "TrackInspectAlg::execute() ------ end ------"  << endmsg;
    return StatusCode::SUCCESS;
}

double TrackInspectAlg::match(edm4hep::ConstMCParticle particle, edm4hep::ConstTrack track){

    int NHits = track.trackerHits_size();

    double matchedHits = 0;
    double usedHits = 0;
    for (int i = 0; i < NHits; i++) {
        edm4hep::ConstTrackerHit hit = track.getTrackerHits(i);
        usedHits++;
        std::pair<edm4hep::ConstTrackerHit, edm4hep::ConstMCParticle> ele = std::make_pair(hit, particle);
        //std::cout << "lookup --> " << ele.first << std::endl;
        //if (hitmap.find(ele) != hitmap.end() ) {
        //std::cout << "find --> " << hitmap[ele] << std::endl;
        //} 
        if (hitmap.find(ele) != hitmap.end() && hitmap[ele] > 0.2) {
            matchedHits++;
        }   
    }   

    // UTIL::BitField64* encoder = new UTIL::BitField64(lcio::ILDCellID0::encoder_string);
    //     encoder->setValue(hit->getCellID0());
    //     int detID = (*encoder)[lcio::ILDCellID0::subdet];
    //     if (detID < 0 || !usedDetectorsArray[detID]) continue;
    // delete encoder;

    return matchedHits / usedHits;
}

void TrackInspectAlg::Fill(edm4hep::ConstMCParticle particle, edm4hep::ConstTrack theTrack) {
    pid[m_nParticles] = particle.getPDG();

    vx[m_nParticles] = particle.getVertex().x;
    vy[m_nParticles] = particle.getVertex().y;
    vz[m_nParticles] = particle.getVertex().z;

    ex[m_nParticles] = particle.getEndpoint().x;
    ey[m_nParticles] = particle.getEndpoint().y;
    ez[m_nParticles] = particle.getEndpoint().z;

    TLorentzVector v2;
    v2.SetXYZM(
        particle.getMomentum().x, 
        particle.getMomentum().y,
        particle.getMomentum().z,
        particle.getMass() );
    TRUEPX[m_nParticles] = v2.X();
    TRUEPY[m_nParticles] = v2.Y();
    TRUEPZ[m_nParticles] = v2.Z();
    TRUEPE[m_nParticles] = v2.E();
    TRUEPT[m_nParticles] = v2.Pt();
    TRUEP[m_nParticles]  = v2.P();
    TRUEETA[m_nParticles] = v2.PseudoRapidity();
    TRUEY[m_nParticles] = v2.Rapidity();
    TRUETHETA[m_nParticles] = v2.Theta();

    // TVector3 theMomentum = TVector3(
    //         particle->getMomentum().x,
    //         particle->getMomentum().y,
    //         particle->getMomentum().z );

    // WMSelectionVariable = (
    //         theParticle->isDecayedInTracker()==0 &&
    //         theMomentum.Perp()>1.0 &&
    //         theParticle->isDecayedInCalorimeter() &&
    //         theParticle->getGeneratorStatus() == 1 &&
    //         sin(theMomentum.Theta()) > 0.18 );

    if (theTrack.isAvailable()) {
        for (std::vector<edm4hep::TrackState>::const_iterator it = theTrack.trackStates_end() - 1; it != theTrack.trackStates_begin() - 1; it--){
            edm4hep::TrackState trackState = *it;
            Omega[m_nParticles] = trackState.omega;
            TanLambda[m_nParticles] = trackState.tanLambda;
            Phi[m_nParticles] = trackState.phi;
            D0[m_nParticles] = trackState.D0;
            Z0[m_nParticles] = trackState.Z0;
            nHits[m_nParticles] = theTrack.trackerHits_size();
        }
    }   
    else {
        Omega[m_nParticles] = -1;
        TanLambda[m_nParticles] = -1;
        Phi[m_nParticles] = -1;
        D0[m_nParticles] = -1;
        Z0[m_nParticles] = -1; 
        nHits[m_nParticles] = 0;
    }   
}

std::vector<edm4hep::ConstTrack> TrackInspectAlg::MCParticleTrackAssociator(edm4hep::ConstMCParticle theParticle) {
    std::vector<edm4hep::ConstTrack> theTracks;
    // std::cout << "The particle: " << theParticle.getPDG() << " " << theParticle << std::endl;
    for (auto matchtuple: matchvec){
        if (std::get<0>(matchtuple) == theParticle){
            if (std::get<2>(matchtuple) > _weight){
                theTracks.push_back(std::get<1>(matchtuple));
            }
        }
    }   
    return theTracks;
}

void TrackInspectAlg::initializeRelationCollections(std::vector<const edm4hep::MCRecoTrackerAssociationCollection*> &relCols) {

    // Use TPC 
    if (_useTPC) {
	    const edm4hep::MCRecoTrackerAssociationCollection* relCol = nullptr;
	    // Establish the relation of MCParticle --> TrackerHit
	    try {
		relCol = _TPCRelColHdl.get();
	    }
	    catch ( GaudiException &e ) {
		debug() << "Collection " << _TPCRelColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
	    }
            relCols.push_back(relCol);
    }

    // Use VXD
    if (_useVXD) {
	    const edm4hep::MCRecoTrackerAssociationCollection* relCol = nullptr;
	    // Establish the relation of MCParticle --> TrackerHit
	    try {
		relCol = _VXDRelColHdl.get();
	    }
	    catch ( GaudiException &e ) {
		debug() << "Collection " << _VXDRelColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
	    }
            relCols.push_back(relCol);
    }

    // Use SIT
    if (_useSIT) {
	    const edm4hep::MCRecoTrackerAssociationCollection* relCol = nullptr;
	    // Establish the relation of MCParticle --> TrackerHit
	    try {
		relCol = _SITRelColHdl.get();
	    }
	    catch ( GaudiException &e ) {
		debug() << "Collection " << _SITRelColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
	    }
            relCols.push_back(relCol);
    }

    // Use SET
    if (_useSET) {
	    const edm4hep::MCRecoTrackerAssociationCollection* relCol = nullptr;
	    // Establish the relation of MCParticle --> TrackerHit
	    try {
		relCol = _SETRelColHdl.get();
	    }
	    catch ( GaudiException &e ) {
		debug() << "Collection " << _SETRelColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
	    }
            relCols.push_back(relCol);
    }

    // Use FTD
    if (_useFTD) {
	    const edm4hep::MCRecoTrackerAssociationCollection* relCol = nullptr;
	    // Establish the relation of MCParticle --> TrackerHit
	    try {
		relCol = _FTDRelColHdl.get();
	    }
	    catch ( GaudiException &e ) {
		debug() << "Collection " << _FTDRelColHdl.fullKey() << " is unavailable in event " << _nEvt << endmsg;
	    }
            relCols.push_back(relCol);
    }
}

//------------------------------------------------------------------------------
StatusCode TrackInspectAlg::finalize(){
    debug() << "Finalizing..." << endmsg;

    return StatusCode::SUCCESS;
}
