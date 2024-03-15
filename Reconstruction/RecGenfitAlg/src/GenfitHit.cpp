#include "GenfitHit.h"
#include "DataHelper/TrackerHitHelper.h"
#include "DetSegmentation/GridDriftChamber.h"
#include "GenfitUnit.h"

#include "DD4hep/Detector.h"
#include "DD4hep/DetElement.h"
#include "DD4hep/Segmentations.h"
#include "DD4hep/DD4hepUnits.h"
#include "edm4hep/TrackerHit.h"
#include "edm4hep/SimTrackerHit.h"
#include "TRandom.h"

#include <iostream>

GenfitHit::GenfitHit(const edm4hep::TrackerHit* trackerHit,
        const edm4hep::SimTrackerHit* simTrackerHit,
        const dd4hep::DDSegmentation::BitFieldCoder* decoder,
        const dd4hep::DDSegmentation::GridDriftChamber* gridDriftChamber,
        double driftVelocity,double driftDistanceErr){
    m_trackerHit=trackerHit;
    m_simTrackerHit=simTrackerHit;
    m_decoder=decoder;
    m_gridDriftChamber=gridDriftChamber;
    m_driftVelocity=driftVelocity;
    m_driftDistanceErr=driftDistanceErr*GenfitUnit::cm;
    //driftVelocity um/ns
    m_driftDistanceTruth=m_trackerHit->getTime()*driftVelocity*GenfitUnit::um;
    m_driftDistance=m_driftDistanceTruth*GenfitUnit::cm;
//    if(driftDistanceErr>0) m_driftDistance+=gRandom->Gaus(0,fabs(driftDistanceErr*GenfitUnit::cm));//FIXME
}

unsigned long long GenfitHit::getCellID() const {
    return m_trackerHit->getCellID();
}

int GenfitHit::getLayer() const {
    return m_decoder->get(getCellID(),"layer");
}

int GenfitHit::getCell() const {
    return m_decoder->get(getCellID(),"cellID");
}

int GenfitHit::getLeftRightAmbig() const {
    TVector3 momTruth=getTruthMom();
    TVector3 pocaOnTrack=getTruthPos();//FIXME, not poca on track
    TVector3 trackDir=momTruth.Unit();
    TVector3 wireDir=(getEnd1()-getEnd0()).Unit();
    TVector3 pocaOnWire=
        m_gridDriftChamber->wirePos_vs_z(getCellID(),pocaOnTrack.Z());
    TVector3 pocaDir=(pocaOnWire-pocaOnTrack).Unit();
    //TVector3 a=pocaDir.Cross(trackDir);
    int lrAmbig=(pocaDir.Cross(trackDir))*wireDir;
    return fabs(lrAmbig)/lrAmbig;
}

TVector3 GenfitHit::getEnd0() const {
    TVector3 end0;
    TVector3 end1;
    m_gridDriftChamber->cellposition(m_trackerHit->getCellID(),end0,end1);//dd4hep unit
    end0*=(1./dd4hep::cm)*GenfitUnit::cm;
    return end0;
}

TVector3 GenfitHit::getEnd1() const {
    TVector3 end0;
    TVector3 end1;
    m_gridDriftChamber->cellposition(m_trackerHit->getCellID(),end0,end1);//dd4hep unit
    end1*=(1./dd4hep::cm)*GenfitUnit::cm;
    return end1;
}

TVector3 GenfitHit::getTruthPos()const{
    edm4hep::Vector3d pos=m_simTrackerHit->getPosition();
    return TVector3(pos.x/dd4hep::cm*GenfitUnit::cm,
            pos.y/dd4hep::cm*GenfitUnit::cm,
            pos.z/dd4hep::cm*GenfitUnit::cm);
}

TVector3 GenfitHit::getTruthMom()const{
    edm4hep::Vector3f mom=m_simTrackerHit->getMomentum();
    return TVector3(mom.x/dd4hep::GeV*GenfitUnit::GeV,
            mom.y/dd4hep::GeV*GenfitUnit::GeV,
            mom.z/dd4hep::GeV*GenfitUnit::GeV);
}

void GenfitHit::print()const{
    //TODO
    std::cout<<"driftDistanceTruth cm "<<m_driftDistanceTruth;
    std::cout<<"driftDistance cm "<<m_driftDistance;
    std::cout<<"driftDistanceErr cm "<<m_driftDistanceErr<<std::endl;
}

