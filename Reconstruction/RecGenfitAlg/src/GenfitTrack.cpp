#include "GenfitTrack.h"
#include "GenfitMsg.h"
#include "GenfitField.h"

//CEPCSW
#include "DataHelper/HelixClass.h"

//Externals
#include "DD4hep/DD4hepUnits.h"
#include "edm4hep/MCParticle.h"
#include "edm4hep/Track.h"
#include "edm4hep/ReconstructedParticle.h"
#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/Vector3d.h"
#include "DetSegmentation/GridDriftChamber.h"

//genfit
#include "Track.h"
#include "MeasuredStateOnPlane.h"
#include "RKTrackRep.h"
#include "TrackPoint.h"
#include "StateOnPlane.h"
#include "KalmanFitterInfo.h"
#include "KalmanFittedStateOnPlane.h"
#include "AbsTrackRep.h"
#include "FitStatus.h"
#include "SpacepointMeasurement.h"
#include "WireMeasurementNew.h"

//ROOT
#include "TRandom.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMatrixDSym.h"

const int GenfitTrack::s_PDG[2][5]
={{-11,-13,211,321,2212},{11,13,-211,-321,-2212}};

    bool
sortDCHit(edm4hep::ConstTrackerHit hit1,edm4hep::ConstTrackerHit hit2)
{
    bool isEarly=hit1.getTime()<hit2.getTime();
    return isEarly;
}

    GenfitTrack::GenfitTrack(const GenfitField* genfitField, const dd4hep::DDSegmentation::GridDriftChamber* seg, const char* name)
:m_genfitField(genfitField),m_gridDriftChamber(seg),m_name(name), m_track(nullptr) ,m_reps() ,m_debug(0)
{

}

GenfitTrack::~GenfitTrack()
{
    ///Note: track reps and points will be deleted in the destructor of track
    ///implemented in genfit::Track::Clear()
    delete m_track;
}

/// create a Genfit track from track state, without trackRep
/// Initialize track with seed states
/// NO unit conversion here
bool GenfitTrack::createGenfitTrack(int pdgType,int charge,
        TLorentzVector posInit, TVector3 momInit, TMatrixDSym covMInit)
{
    TVectorD seedState(6);
    TMatrixDSym seedCov(6);

    //for(int i = 0; i < 6; ++i) {
    //  for(int j = 0; j < 6; ++j) {
    //    seedCov(i,j)=covMInit(i,j);
    //  }
    //}
    //yzhang FIXME
    //seed position
    for(int i = 0; i < 3; ++i) {
        seedState(i)=posInit[i];
        //yzhang TODO from covMInit to seedCov
        double resolution = 0.1;//*dd4hep::mm/dd4hep::cm;
        seedCov(i,i)=resolution*resolution;
        if(i==2) seedCov(i,i)=0.5*0.5;
    }
    //seed momentum
    for(int i = 3; i < 6; ++i){
        //seedState(i)=momInit[i-3]*(dd4hep::GeV);
        seedState(i)=momInit[i-3];
        //yzhang TODO from covMInit to seedCov
        seedCov(i,i)=0.01;//pow(resolution / sqrt(3),2);
    }

    if(nullptr==m_track) m_track=new genfit::Track();
    m_track->setStateSeed(seedState);
    m_track->setCovSeed(seedCov);

    /// new a track representation and add to the track
    int chargeId=0;
    charge<0 ? chargeId=0 : chargeId=1;

    GenfitMsg::get()<<MSG::DEBUG<<m_name<<" CreateGenfitTrack seed pos("
        <<seedState[0]<<" " <<seedState[1]<<" " <<seedState[2]<<")cm ("
        <<seedState[3]<<" " <<seedState[4]<<" " <<seedState[5]<<")GeV pdg "
        <<s_PDG[chargeId][pdgType]<<endmsg;
    GenfitMsg::get()<<MSG::DEBUG<<"seedCov "<<endmsg;

    addTrackRep(s_PDG[chargeId][pdgType]);

    if(m_debug>0) seedCov.Print();
    return true;
}

///Create a Genfit track with MCParticle, unit conversion here
bool GenfitTrack::createGenfitTrackFromMCParticle(int pidType,
        const edm4hep::MCParticle& mcParticle, double eventStartTime)
{
    ///get track parameters from McParticle
    edm4hep::Vector3d mcPocaPos = mcParticle.getVertex();//mm
    edm4hep::Vector3f mcPocaMom = mcParticle.getMomentum();//GeV
    GenfitMsg::get()<<MSG::DEBUG<<"seedPos poca "<< mcPocaPos.x
        <<" "<<mcPocaPos.y<<" "<<mcPocaPos.z<<" mm "<<endmsg;
    GenfitMsg::get()<<MSG::DEBUG<<"seedMom poca "<< mcPocaMom.x
        <<" "<<mcPocaMom.y<<" "<<mcPocaMom.z<<" GeV "<<endmsg;

    ///Pivot to first layer to avoid correction of beam pipe
    edm4hep::Vector3d firstLayerPos(1e9,1e9,1e9);
    edm4hep::Vector3f firstLayerMom(1e9,1e9,1e9);
    pivotToFirstLayer(mcPocaPos,mcPocaMom,firstLayerPos,firstLayerMom);

    //TODO convert unit
    ///Get seed position and momentum
    TLorentzVector seedPos(firstLayerPos.x,firstLayerPos.y,firstLayerPos.z,
            eventStartTime);
    TVector3 seedMom(firstLayerMom.x,firstLayerMom.y,firstLayerMom.z);
    GenfitMsg::get()<<MSG::DEBUG<<"seedPos "<< firstLayerPos.x
        <<" "<<firstLayerPos.y<<" "<<firstLayerPos.z<<endmsg;
    GenfitMsg::get()<<MSG::DEBUG<<"seedMom "<< firstLayerMom.x
        <<" "<<firstLayerMom.y<<" "<<firstLayerMom.z<<endmsg;

    ///Get error matrix of seed track
    TMatrixDSym covM(5);//FIXME, TODO

    ///Create a genfit track with seed
    GenfitMsg::get()<<MSG::DEBUG<<"createGenfitTrack " ;
    if(!GenfitTrack::createGenfitTrack(pidType,mcParticle.getCharge(),
                seedPos,seedMom,covM)){
        GenfitMsg::get()<<MSG::ERROR<<"GenfitTrack"
            <<" Error in createGenfitTrack" <<endmsg;
        return false;
    }else{
        GenfitMsg::get()<<MSG::DEBUG<<"GenfitTrack "
            <<"createGenfitTrackFromMCParticle track created" <<endmsg;
    }
    return true;
}//end of createGenfitTrackFromMCParticle

///Create a Genfit track with MCParticle, unit conversion here
bool GenfitTrack::createGenfitTrackFromEDM4HepTrack(int pidType,
        const edm4hep::Track& track, double eventStartTime)
{
    //std::cout<<__FILE__<<"   "<<__LINE__<<" bz kilogauss "<<m_genfitField->getBz({0.,0.,0.})/dd4hep::kilogauss<<std::endl;
    //std::cout<<__FILE__<<"   "<<__LINE__<<" bz tesla "<<m_genfitField->getBz({0.,0.,0.})/dd4hep::tesla<<std::endl;
    //std::cout<<__FILE__<<"   "<<__LINE__<<" bz "<<m_genfitField->getBz({0.,0.,0.})<< dd4hep::kilogauss <<" "<<dd4hep::tesla<<std::endl;
    //TODO
    //pivotToFirstLayer(mcPocaPos,mcPocaMom,firstLayerPos,firstLayerMom);
    //Get track parameters
    edm4hep::TrackState trackStat=track.getTrackStates(0);//FIXME?
    HelixClass helixClass;
    helixClass.Initialize_Canonical(trackStat.phi,trackStat.D0,
            trackStat.Z0,trackStat.omega,trackStat.tanLambda,
            m_genfitField->getBz({0.,0.,0.})*dd4hep::kilogauss/dd4hep::tesla);
    TLorentzVector posInit(helixClass.getReferencePoint()[0],
            helixClass.getReferencePoint()[1],
            helixClass.getReferencePoint()[2],eventStartTime);
    posInit.SetX(posInit.X()*dd4hep::mm);
    posInit.SetY(posInit.Y()*dd4hep::mm);
    posInit.SetZ(posInit.Z()*dd4hep::mm);
    TVector3 momInit(helixClass.getMomentum()[0],
            helixClass.getMomentum()[1],helixClass.getMomentum()[2]);
    momInit.SetX(momInit.x()*dd4hep::GeV);
    momInit.SetY(momInit.y()*dd4hep::GeV);
    momInit.SetZ(momInit.z()*dd4hep::GeV);
    TMatrixDSym covMInit;
    if(!createGenfitTrack(pidType,
                int(trackStat.omega/fabs(trackStat.omega)),//charge,FIXME?
                posInit,momInit,covMInit)){
        return false;
    }
    return true;
}

/// Add a 3d SpacepointMeasurement with MC truth position smeared by sigma
bool GenfitTrack::addSpacePointMeasurement(const TVectorD& pos,
        double sigma, int detID, int hitID)
{

    /// Convert from CEPCSW unit to genfit unit, cm
    TVectorD pos_t(3);
    pos_t(0) = pos(0)*dd4hep::mm;
    pos_t(1) = pos(1)*dd4hep::mm;
    pos_t(2) = pos(2)*dd4hep::mm;
    /// smear hit position with same weight
    TVectorD pos_smeared(pos_t);
    for (int ii=0;ii<3;ii++){
        pos_smeared[ii] += gRandom->Gaus(0, sigma/TMath::Sqrt(3.));
    }

    /// New a SpacepointMeasurement
    TMatrixDSym hitCov(3);
    hitCov(0,0) = sigma*sigma;
    hitCov(1,1) = sigma*sigma;
    hitCov(2,2) = sigma*sigma;

    GenfitMsg::get()<< MSG::DEBUG<<m_name<<" addSpacePointMeasurement"
        <<hitID<<" " <<pos_t[0]<<" "<<pos_t[1]<<" "<<pos_t[2]<<"cm smeared "
        <<pos_smeared[0]<<" "<<pos_smeared[1]<<" "<<pos_smeared[2]
        <<" sigma "<<sigma<<" cm"<<endmsg;

    genfit::SpacepointMeasurement* sMeas =
        new genfit::SpacepointMeasurement(pos_smeared,hitCov,detID,hitID,nullptr);
    genfit::TrackPoint* trackPoint = new genfit::TrackPoint(sMeas,m_track);
    m_track->insertPoint(trackPoint);

    return true;
}


/// Add a WireMeasurement, no Unit conversion here
void GenfitTrack::addWireMeasurement(double driftDistance,
        double driftDistanceError, const TVector3& endPoint1,
        const TVector3& endPoint2, int lrAmbig, int detID, int hitID)
{
    try{
        /// New a WireMeasurement
        genfit::WireMeasurementNew* wireMeas = new genfit::WireMeasurementNew(
                driftDistance, driftDistanceError, endPoint1, endPoint2, detID,
                hitID, nullptr);
        wireMeas->setMaxDistance(0.5);//0.5 cm FIXME
        wireMeas->setLeftRightResolution(lrAmbig);

        GenfitMsg::get()<<MSG::DEBUG<<m_name<<" Add wire measurement(cm) "<<hitID
            <<" ep1("<<endPoint1[0]<<" "<<endPoint1[1]<<" "<<endPoint1[2]
            <<") ep2("<<endPoint2[0]<<" "<<endPoint2[1]<<" "<<endPoint2[2]
            <<") drift "<<driftDistance<<" driftErr "<<driftDistanceError
            <<" lr "<<lrAmbig<<" detId "<<detID << " hitId "<< hitID
            <<endmsg;

        ///New a TrackPoint,create connection between meas. and trackPoint
        genfit::TrackPoint* trackPoint=new genfit::TrackPoint(wireMeas,m_track);
        wireMeas->setTrackPoint(trackPoint);

        m_track->insertPoint(trackPoint);

    }catch(genfit::Exception e){
        GenfitMsg::get()<<MSG::DEBUG<<m_name
            <<"Add wire measurement exception"<<endmsg;
        e.what();
    }
}//end of addWireMeasurementOnTrack

//Add wire measurement on wire, unit conversion here
bool GenfitTrack::addWireMeasurementOnTrack(edm4hep::Track& track,double sigma)
{
    for(int iHit=0;iHit<track.trackerHits_size();iHit++){
        edm4hep::ConstTrackerHit hit=track.getTrackerHits(iHit);

        double driftVelocity=40;//FIXME, TODO, um/ns
        double driftDistance=hit.getTime()*driftVelocity*dd4hep::um;
        TVector3 endPointStart(0,0,0);
        TVector3 endPointEnd(0,0,0);
        m_gridDriftChamber->cellposition(hit.getCellID(),endPointStart,
                endPointEnd);
        int lrAmbig=0;
        GenfitMsg::get()<<MSG::DEBUG<<m_name<<" time  "<<hit.getTime()<<" driftVelocity "
            <<driftVelocity<<std::endl;
        GenfitMsg::get()<<MSG::DEBUG<<m_name<<" wire pos " <<endPointStart.X()<<" "
            <<endPointStart.Y()<<" " <<endPointStart.Z()<<" "
            <<endPointEnd.X()<<" " <<endPointEnd.Y()<<" "
            <<endPointEnd.Z()<<" " <<std::endl;
        endPointStart.SetX(endPointStart.x()*dd4hep::cm);
        endPointStart.SetY(endPointStart.y()*dd4hep::cm);
        endPointStart.SetZ(endPointStart.z()*dd4hep::cm);
        endPointEnd.SetX(endPointEnd.x()*dd4hep::cm);
        endPointEnd.SetY(endPointEnd.y()*dd4hep::cm);
        endPointEnd.SetZ(endPointEnd.z()*dd4hep::cm);
        addWireMeasurement(driftDistance,sigma*dd4hep::cm,endPointStart,
                endPointEnd,lrAmbig,hit.getCellID(),iHit);
    }
    return true;
}//end of addWireMeasurementOnTrack of Track

/// Get MOP
bool GenfitTrack::getMOP(int hitID,
        genfit::MeasuredStateOnPlane& mop, genfit::AbsTrackRep* trackRep) const
{
    if(nullptr == trackRep) trackRep = getRep();
    try{
        mop = m_track->getFittedState(hitID,trackRep);
    }catch(genfit::Exception e){
        e.what();
        return false;
    }
    return true;
}

/// New and add a track representation to track
genfit::RKTrackRep* GenfitTrack::addTrackRep(int pdg)
{
    /// create a new track representation
    genfit::RKTrackRep* rep = new genfit::RKTrackRep(pdg);
    m_reps.push_back(rep);
    m_track->addTrackRep(rep);
    //try{
    //  genfit::MeasuredStateOnPlane stateInit(rep);
    //  rep->setPosMomCov(stateInit, pos, mom, covM);
    //}catch(genfit::Exception e){
    //  GenfitMsg::get() << MSG::DEBUG<<m_name<<" Exception in set track status"
    //  <<endmsg     ;
    //  std::cout<<e.what()<<std::endl;
    //  return false;
    //}
    return rep;
}

/// Get the position from genfit::Track::getStateSeed
const TLorentzVector GenfitTrack::getSeedStatePos()const
{
    TVectorD seedStat(6);
    seedStat = m_track->getStateSeed();
    TVector3 p(seedStat[0],seedStat[1],seedStat[2]);
    p = p*dd4hep::cm;
    TLorentzVector pos(p[0],p[1],p[2],9999);//FIXME
    return pos;
}

/// Get the momentum from genfit::Track::getStateSeed
const TVector3 GenfitTrack::getSeedStateMom() const
{
    TVectorD seedStat(6); seedStat = m_track->getStateSeed();
    TVector3 mom(seedStat[3],seedStat[4],seedStat[5]);
    return mom*dd4hep::GeV;
}

/// Get the seed states of momentum and position
void GenfitTrack::getSeedStateMom(TLorentzVector& pos, TVector3& mom) const
{
    TVectorD seedStat(6); seedStat = m_track->getStateSeed();
    mom = TVector3(seedStat[3],seedStat[4],seedStat[5])*dd4hep::GeV;
    seedStat = m_track->getStateSeed();
    TVector3 p = TVector3(seedStat[0],seedStat[1],seedStat[2])*dd4hep::cm;
    pos.SetXYZT(p[0],p[1],p[2],9999);//FIXME time
}

unsigned int GenfitTrack::getNumPoints() const
{
    return m_track->getNumPoints();
}

/// Test the fit result FIXME
bool GenfitTrack::fitSuccess(int repID) const
{

    genfit::FitStatus* fitStatus = m_track->getFitStatus(getRep(repID));

    /// Test fitting converged
    if (!fitStatus->isFitted()||!fitStatus->isFitConverged()
            ||fitStatus->isFitConvergedFully()) {
        GenfitMsg::get() << MSG::DEBUG<<m_name<< "Fitting is failed... isFitted"
            <<fitStatus->isFitted()<<" , isFitConverged "
            <<fitStatus->isFitConverged()<<", isFitConvergedFully "
            <<fitStatus->isFitConvergedFully()<<endmsg;
        return false;
    }

    double chi2 = fitStatus->getChi2();
    double ndf  = fitStatus->getNdf();
    GenfitMsg::get() << MSG::INFO << "Fit result: chi2 "<<chi2 <<" ndf "<<ndf
        << " chi2/ndf = " << chi2/ndf<<endmsg;

    /// Test fitting chi2
    if (chi2<= 0) {
        GenfitMsg::get() << MSG::DEBUG<<m_name<< "Fit chi2<0 (chi2,ndf) = (" <<
            chi2 << "," << ndf  << ")"<<endmsg;
        return false;
    }
    return true;
}

void GenfitTrack::setDebug(int debug)
{
    m_debug = debug;
    for(unsigned int i=0;i<m_reps.size();i++){ m_reps[i]->setDebugLvl(debug); }
}

void GenfitTrack::printSeed() const
{
    TLorentzVector pos = getSeedStatePos();
    TVector3 mom = getSeedStateMom();
    print(pos,mom);
    GenfitMsg::get()<<MSG::DEBUG<<m_name<<" NumPoints "<<getNumPoints()<<endmsg;
}

void GenfitTrack::printFitted(int repID) const
{
    TLorentzVector fittedPos;
    TVector3 fittedMom;
    TMatrixDSym cov;

    GenfitMsg::get()<<MSG::DEBUG<<m_name<< "printFitted nHit="
        <<m_track->getNumPoints()<<endmsg;
    for(unsigned int iHit=0; iHit<m_track->getNumPoints(); iHit++){
        if (getPosMomCovMOP((int) iHit, fittedPos, fittedMom, cov, repID)){
            //print(fittedPos,fittedMom,to_string(iHit).c_str());//TODO
        }else{
            GenfitMsg::get()<<MSG::DEBUG<<m_name<<"Hit "<<iHit
                <<" have no fitter info"<<endmsg;
        }
    }
}

/// Print track information
void GenfitTrack::print( TLorentzVector pos, TVector3 mom,
        const char* comment) const
{
    TVector3 pglo = pos.Vect();
    TVector3 mglo = mom;

    // TODO
    GenfitMsg::get()<<MSG::DEBUG<<m_name<<" "<<comment<<endmsg;

    if(m_debug>1){
        for(unsigned int i=0;i<m_reps.size();i++){ m_reps[i]->Print(); }
    }
    //for(unsigned int i=0; i<m_track->getNumPoints(); i++){
    //  m_track->getPoint(i)->print();
    //}
}

/// Get position, momentum, cov on plane of hitID-th hit
bool GenfitTrack::getPosMomCovMOP(int hitID, TLorentzVector& pos,
        TVector3& mom, TMatrixDSym& cov, int repID) const
{
    TVector3 p;
    genfit::MeasuredStateOnPlane mop;
    if(!getMOP(hitID,mop,getRep(repID))) return false;
    mop.getPosMomCov(p,mom,cov);
    pos.SetVect(p*dd4hep::cm);
    pos.SetT(9999);//FIXME
    mom = mom*(dd4hep::GeV);
    //FIXME
    //TrackingUtils::CovConvertUnit(cov, dd4hep::cm, dd4hep::GeV);
    return true;
}

int GenfitTrack::getNumPointsWithFittedInfo(int repID) const
{
    int nHitWithFittedInfo = 0;
    int nHit = m_track->getNumPointsWithMeasurement();
    for(int i=0; i<nHit; i++){
        if(nullptr != m_track->getPointWithFitterInfo(i,getRep(repID))){
            nHitWithFittedInfo++;
        }
    }
    return nHitWithFittedInfo;
}

int GenfitTrack::getFittedState(TLorentzVector& pos, TVector3& mom,
        TMatrixDSym& cov, int repID, bool biased) const
{
    //check number of hit with fitter info
    if(getNumPointsWithFittedInfo(repID)<=2) return 1;

    //get track rep
    genfit::AbsTrackRep* rep = getRep(repID);
    if(nullptr == rep) return 2;

    //get first or last measured state on plane
    genfit::MeasuredStateOnPlane mop;
    try{
        mop = m_track->getFittedState(biased);
    }catch(genfit::Exception e){
        std::cout<<" getNumPointsWithFittedInfo "
            <<getNumPointsWithFittedInfo(repID)
            <<" no TrackPoint with fitted info "<<std::endl;
        GenfitMsg::get()<<MSG::DEBUG<<m_name
            <<"Exception in getFittedState"<<endmsg;
        std::cout<<e.what()<<std::endl;
        return 3;
    }

    //get state
    TVector3 p;
    mop.getPosMomCov(p,mom,cov);
    pos.SetVect(p*dd4hep::cm);
    pos.SetT(9999);//FIXME
    mom = mom*(dd4hep::GeV);

    return 0;//success
}

// Get point with fitter info
int GenfitTrack::getDetIDWithFitterInfo(int hitID, int idRaw) const
{
    return m_track->getPointWithFitterInfo(hitID)->
        getRawMeasurement(idRaw)->getDetId();
}

int GenfitTrack::getPDG(int id) const
{
    return m_reps[id]->getPDG();
}

int GenfitTrack::getPDGCharge(int id) const
{
    return m_reps[id]->getPDGCharge();
}

const genfit::FitStatus*
GenfitTrack::getFitStatus(int repID) const
{
    return m_track->getFitStatus(getRep(repID));
}

/// Extrapolate track to the cloest point of approach(POCA) to the wire of hit,
/// return StateOnPlane of this POCA
/// inputs
///  pos,mom ... position & momentum at starting point, unit= [mm]&[GeV/c]
///              (converted to [cm]&[GeV/c] inside this function)
///  hit ... destination
/// outputs poca [mm] (converted from [cm] in this function) ,global
double GenfitTrack::extrapolateToHit( TVector3& poca, TVector3& pocaDir,
        TVector3& pocaOnWire, double& doca, TVector3 pos, TVector3 mom,
        TVector3 end0,//one end of the hit wire
        TVector3 end1,//the orhter end of the hit wire
        int debug,
        int repID,
        bool stopAtBoundary,
        bool calcJacobianNoise)const
{

    //genfit::MeasuredStateOnPlane state = getMOP(iPoint); // better?
    //genfit::MeasuredStateOnPlane state = getMOP(0);      // better?
    //To do the extrapolation without IHitSelection,above 2 lines are not used.
    pos = pos*dd4hep::cm;//FIXME
    mom = mom*dd4hep::GeV;

    //std::cout<<__LINE__<<" extrapolate to Hit pos and mom"<<std::endl;
    //pos.Print();
    //mom.Print();

    genfit::AbsTrackRep* rep = new genfit::RKTrackRep(getRep(repID)->getPDG());
    rep->setDebugLvl(debug);
    genfit::MeasuredStateOnPlane state(rep);
    rep->setPosMom(state, pos, mom);


    //genfit::MeasuredStateOnPlane state(m_track->getRep(repID));
    //m_track->getRep(repID)->setPosMom(state, pos, mom);

    //m_track->PrintSeed();
    double extrapoLen(0);
    //std::cout<<" wire1 "<<std::endl;
    //end0.Print();
    //std::cout<<" wire0 "<<std::endl;
    //end1.Print();
    //std::cout<<" state "<<std::endl;
    //state.Print();


    // forth
    try {
        //genfit::RKTrackRep* rkRep =
        //  dynamic_cast<genfit::RKTrackRep*> (m_track->getRep(repID));
        //extrapoLen = rkRep->extrapolateToLine(state, end0, end1, poca,
        //    pocaDir, pocaOnWire, stopAtBoundary, calcJacobianNoise);
        extrapoLen = rep->extrapolateToLine(state, end0, end1, poca,
                pocaDir, pocaOnWire, stopAtBoundary, calcJacobianNoise);
    } catch (genfit::Exception& e) {
        GenfitMsg::get() << MSG::ERROR
            <<"Exception in GenfigFitter::ExtrapolateToHit"<<e.what()<<endmsg;
        return extrapoLen;
    }

    //poca = poca*(dd4hep::cm);
    //pocaOnWire = pocaOnWire*(dd4hep::cm);
    pocaOnWire = pocaOnWire;
    doca = (pocaOnWire-poca).Mag();
    //TODO: debug pocaOnWire
    GenfitMsg::get() << MSG::DEBUG << " poca "<<poca.x()<<","<<poca.y()
        <<" "<<poca.z()<<" doca "<<doca<<endmsg;
    GenfitMsg::get() << MSG::DEBUG << " pocaOnWire "<<pocaOnWire.x()
        <<","<<pocaOnWire.y()<<" "<<pocaOnWire.z()<<" doca "<<doca<<endmsg;

    return extrapoLen*(dd4hep::cm);
}//end of ExtrapolateToHit


///Add space point measurement from edm4hep::Track to genfit track
int GenfitTrack::addSpacePointMeasurementOnTrack(const edm4hep::Track& track,
        double sigma) {

    //Sort sim DC hit by time
    std::vector<edm4hep::ConstTrackerHit> sortedDCTrackHitCol;
    for(int i=0;i<track.trackerHits_size();i++){
        sortedDCTrackHitCol.push_back(track.getTrackerHits(i));
    }
    std::sort(sortedDCTrackHitCol.begin(),sortedDCTrackHitCol.end(),sortDCHit);

    int hitID=0;
    for(auto dCTrackerHit: sortedDCTrackHitCol){
        int rand=gRandom->Integer(sortedDCTrackHitCol.size());
        if(rand>200){ continue; }
        edm4hep::Vector3d pos=dCTrackerHit.getPosition();
        TVectorD p(3);
        p[0]=pos.x*dd4hep::mm;
        p[1]=pos.y*dd4hep::mm;
        p[2]=pos.z*dd4hep::mm;
        unsigned long long detID = dCTrackerHit.getCellID();
        if(addSpacePointMeasurement(p,sigma,detID,hitID)){
            hitID++;
        }else{
            GenfitMsg::get()<<MSG::ERROR<<"addSpacePointMeasurementOnTrack "
                <<detID<<" faieled" <<endmsg;
        }
    }
    GenfitMsg::get()<<MSG::DEBUG<<"GenfitTrack nHitAdded "<<hitID<<endmsg;
    return hitID;
}

bool GenfitTrack::storeTrack(edm4hep::ReconstructedParticle& recParticle,
        int pidType, int ndfCut, double chi2Cut)
{

    /// Get fit status
    const genfit::FitStatus* fitState = m_track->getFitStatus();
    double chi2 = fitState->getChi2();
    double ndf = fitState->getNdf();
    double charge = fitState->getCharge();
    int isFitted = fitState->isFitted();
    int isConverged = fitState->isFitConverged();
    int isConvergedFully = fitState->isFitConvergedFully();
    TMatrixDSym fittedCov;
    TLorentzVector fittedPos;
    TVector3 fittedMom;
    int fittedState=getFittedState(fittedPos,fittedMom,fittedCov);
    GenfitMsg::get()<<MSG::DEBUG<<"fit result: get status OK? pidType "
        <<pidType<<" fittedState==0 " <<(0==fittedState)<<" isFitted "<<isFitted
        <<" isConverged "<<isConverged<<" ndf "<<ndf<<endmsg;
    if((0!=fittedState) || (!isFitted) || !isConvergedFully || (ndf<ndfCut)){
        GenfitMsg::get()<<MSG::DEBUG<<" fitting failed"<<endmsg;
    }else{
        GenfitMsg::get()<<MSG::DEBUG<<"fit result: Pos("<<
            fittedPos.X()<<" "<<
            fittedPos.Y()<<" "<<
            fittedPos.Z()<<") mom("<<
            fittedMom.X()<<" "<<
            fittedMom.Y()<<" "<<
            fittedMom.Z()<<") p_tot "<<
            fittedMom.Mag()<<" pt "<<
            fittedMom.Perp()<<
            " chi2 "<<chi2<<
            " ndf "<<ndf
            <<endmsg;
    }
    float pos[3]={fittedPos.X(),fittedPos.Y(),fittedPos.Z()};
    float mom[3]={fittedMom.X(),fittedMom.Y(),fittedMom.Z()};
    HelixClass helix;
    helix.Initialize_VP(pos,mom,charge,m_genfitField->getBz(fittedPos.Vect()));


    /////track status at POCA to origin
    //TVector3 origin(pivot);
    //TVector3 pocaToOrigin_pos(1e9*dd4hep::cm,1e9*dd4hep::cm,1e9*dd4hep::cm);
    //TVector3 pocaToOrigin_mom(1e9*dd4hep::GeV,1e9*dd4hep::GeV,1e9*dd4hep::GeV);

    //if(ExtrapolateToPoint(pocaToOrigin_pos,pocaToOrigin_mom,
    //            m_track,origin) > 1e6*dd4hep::cm){
    //    log<<"extrapolate to origin failed";
    //    return false;
    //}


    //new TrackState
    edm4hep::TrackState* trackState = new edm4hep::TrackState();
    trackState->location=0;//edm4hep::AtIP;
    trackState->D0=helix.getD0();
    trackState->phi=helix.getPhi0();
    trackState->omega=helix.getOmega();
    trackState->Z0=helix.getZ0();
    trackState->tanLambda=helix.getTanLambda();
    trackState->referencePoint=helix.getReferencePoint();

    //    std::array<float,15> covMatrix;
    //    int k=0;
    //    for(int i=0;i<5;i++){
    //        for(int j=0;j<5;j++){
    //            if(i<=j) covMatrix[k]=;//FIXME
    //        }
    //    }
    //    trackState.covMatrix=

    //new Track
    edm4hep::Track* track = new edm4hep::Track();
    //track->setType();
    track->setChi2(fitState->getChi2());
    track->setNdf(fitState->getNdf());
    //track->setDEdx();
    //track->setRadiusOfInnermostHit();//FIXME
    //track->addToTrackerHits();

    //new ReconstructedParticle
    //recParticle->setType();
    //dcRecParticle->setEnergy();

    edm4hep::Vector3f momVec3(helix.getMomentum()[0],
            helix.getMomentum()[1],helix.getMomentum()[2]);
    recParticle.setMomentum(momVec3);
    //recParticle->setReferencePoint(referencePoint);
    recParticle.setCharge(helix.getCharge());
    //    recParticle->setMass();
    //    recParticle->setCovMatrix();
    //    rcRecParticle->setStartVertex();
    //recParticle->addToTracks(track);

    return true;
}

void GenfitTrack::pivotToFirstLayer(edm4hep::Vector3d& pos,
        edm4hep::Vector3f& mom, edm4hep::Vector3d& firstPos,
        edm4hep::Vector3f& firstMom)
{
    //FIXME, TODO
    firstPos=pos;
    firstMom=mom;
}

