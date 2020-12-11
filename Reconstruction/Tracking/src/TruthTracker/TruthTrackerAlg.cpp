#include "TruthTrackerAlg.h"
#include "DataHelper/HelixClass.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/MCRecoTrackerAssociationCollection.h"
#include "edm4hep/MCRecoParticleAssociationCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "DDSegmentation/Segmentation.h"
#include "DetSegmentation/GridDriftChamber.h"
#include "DetInterface/IGeomSvc.h"
//#include "edm4hep/SimCalorimeterHitCollection.h"
//#include "edm4hep/CaloHitContributionCollection.h"

#include "DD4hep/Detector.h"
#include "DD4hep/IDDescriptor.h"
#include "DD4hep/Plugins.h"
#include "DD4hep/DD4hepUnits.h"
#include "CLHEP/Random/RandGauss.h"

DECLARE_COMPONENT(TruthTrackerAlg)

TruthTrackerAlg::TruthTrackerAlg(const std::string& name, ISvcLocator* svcLoc)
: GaudiAlgorithm(name, svcLoc), m_dd4hep(nullptr), m_gridDriftChamber(nullptr),m_decoder(nullptr)
{
    declareProperty("MCParticle", m_mcParticleCol,
            "Handle of the input MCParticle collection");
    declareProperty("DigiDCHitCollection", m_digiDCHitsCol,
            "Handle of digi DCHit collection");
    declareProperty("DCHitAssociationCollection", m_dcHitAssociationCol,
        "Handle of association collection");
    declareProperty("DCTrackCollection", m_dcTrackCol,
        "Handle of association collection");
    declareProperty("DCRecParticleCollection", m_dcRecParticleCol,
            "Handle of drift chamber reconstructed particle collection");
    declareProperty("DCRecParticleAssociationCollection", m_dcRecParticleAssociationCol,
            "Handle of drift chamber reconstructed particle collection");
    declareProperty("debug", m_debug=false);
}

StatusCode TruthTrackerAlg::initialize()
{
    m_me = 0.511*pow(10,6);//0.511 MeV to eV
    m_K = 0.307075;//const
    ///Get geometry
    m_geomSvc=service<IGeomSvc>("GeomSvc");
    if (!m_geomSvc) {
        error() << "Failed to get GeomSvc." << endmsg;
        return StatusCode::FAILURE;
    }
    ///Get Detector
    m_dd4hep=m_geomSvc->lcdd();
    if (nullptr==m_dd4hep) {
        error() << "Failed to get dd4hep::Detector." << endmsg;
        return StatusCode::FAILURE;
    }
    //Get Field
    m_dd4hepField=m_geomSvc->lcdd()->field();
    ///Get Readout
    dd4hep::Readout readout=m_dd4hep->readout(m_readout_name);
    ///Get Segmentation
    m_gridDriftChamber=dynamic_cast<dd4hep::DDSegmentation::GridDriftChamber*>
        (readout.segmentation().segmentation());
    if(nullptr==m_gridDriftChamber){
        error() << "Failed to get the GridDriftChamber" << endmsg;
        return StatusCode::FAILURE;
    }
    ///Get Decoder
    m_decoder = m_geomSvc->getDecoder(m_readout_name);
    if (nullptr==m_decoder) {
        error() << "Failed to get the decoder" << endmsg;
        return StatusCode::FAILURE;
    }
    if(m_WriteAna){

        NTuplePtr nt( ntupleSvc(), "MyTuples/TruthTrack_evt" );
        if ( nt ) m_tuple = nt;
        else {
            m_tuple = ntupleSvc()->book( "MyTuples/TruthTrack_evt", CLID_ColumnWiseTuple, "TruthTrack" );
            if ( m_tuple ) {
              m_tuple->addItem( "N_hit", m_hit, 0, 100000 ).ignore();
              m_tuple->addItem( "hit_x" , m_hit, m_hit_x  ).ignore();
              m_tuple->addItem( "hit_y" , m_hit, m_hit_y  ).ignore();
              m_tuple->addItem( "hit_z" , m_hit, m_hit_z  ).ignore();
              m_tuple->addItem( "hit_dedx" , m_hit, m_hit_dedx ).ignore();
              m_tuple->addItem( "track_dedx" , m_track_dedx  ).ignore();
              m_tuple->addItem( "track_dedx_BB" , m_track_dedx_BB  ).ignore();
              m_tuple->addItem( "track_px" , m_track_px  ).ignore();
              m_tuple->addItem( "track_py" , m_track_py  ).ignore();
              m_tuple->addItem( "track_pz" , m_track_pz  ).ignore();
              m_tuple->addItem( "track_mass" , m_track_mass  ).ignore();
              m_tuple->addItem( "track_pid" , m_track_pid  ).ignore();
              m_tuple->addItem( "mc_px" , m_mc_px  ).ignore();
              m_tuple->addItem( "mc_py" , m_mc_py  ).ignore();
              m_tuple->addItem( "mc_pz" , m_mc_pz  ).ignore();
            } 
            else { // did not manage to book the N tuple....
              error() << "    Cannot book N-tuple:" << long( m_tuple ) << endmsg;
              return StatusCode::FAILURE;
            }
        }
    }

    return GaudiAlgorithm::initialize();
}

StatusCode TruthTrackerAlg::execute()
{
    ///Retrieve MC particle(s)
    const edm4hep::MCParticleCollection* mcParticleCol=nullptr;
    mcParticleCol=m_mcParticleCol.get();
    if(nullptr==mcParticleCol){
        debug()<<"MCParticleCollection not found"<<endmsg;
        return StatusCode::SUCCESS;
    }
    ///Retrieve DC digi
    const edm4hep::TrackerHitCollection* digiDCHitsCol=nullptr;
    digiDCHitsCol=m_digiDCHitsCol.get();//FIXME DEBUG
    if(nullptr==digiDCHitsCol){
        debug()<<"TrackerHitCollection not found"<<endmsg;
        //return StatusCode::SUCCESS;
    }

    ///Output Track collection
    edm4hep::TrackCollection* dcTrackCol=m_dcTrackCol.createAndPut();
    ////TODO
    /////Output MCRecoTrackerAssociationCollection collection
    //const edm4hep::MCRecoTrackerAssociationCollection*
    //    mcRecoTrackerAssociationCol=nullptr;
    //if(nullptr==mcRecoTrackerAssociationCol){
    //    log<<MSG::DEBUG<<"MCRecoTrackerAssociationCollection not found"
    //        <<endmsg;
    //    return StatusCode::SUCCESS;
    //}
    //mcRecoTrackerAssociationCol=m_mcRecoParticleAssociation.get();

    ///Convert MCParticle to Track and ReconstructedParticle
    if(m_debug){
        debug()<<"MCParticleCol size="<<mcParticleCol->size()<<endmsg;
    }
    if(mcParticleCol->size() != 1) throw ("current only support one track each event!!!");
    m_hit = 0;
    for(auto mcParticle : *mcParticleCol){
        //if(fabs(mcParticle.getCharge()<1e-6) continue;//Skip neutral particles
        edm4hep::Vector3d posV= mcParticle.getVertex();
        edm4hep::Vector3f momV= mcParticle.getMomentum();//GeV
        float pos[3]={posV.x, posV.y, posV.z};
        float mom_truth[3]={momV.x, momV.y, momV.z};
        float mom[3]      ={momV.x, momV.y, momV.z};
        float mom_smear  = CLHEP::RandGauss::shoot(1, m_mom_resolution);
        mom[0]=mom[0]*mom_smear;
        mom[1]=mom[1]*mom_smear;
        mom[2]=mom[2]*mom_smear;
        //FIXME
        //pivotToFirstLayer(mcPocaPos,mcPocaMom,firstLayerPos,firstLayerMom);
        double B[3]={1e9,1e9,1e9};
        m_dd4hepField.magneticField({0.,0.,0.},B);
        if(m_debug)std::cout<<" PDG "<<mcParticle.getPDG()<<" charge "
            << mcParticle.getCharge()
            <<" pos "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]
            <<" mom "<<mom[0]<<" "<<mom[1]<<" "<<mom[2]
            <<" Bxyz "<<B[0]/dd4hep::tesla<<" "<<B[1]/dd4hep::tesla
            <<" "<<B[2]/dd4hep::tesla<<" tesla"<<std::endl;

        HelixClass helix;
        helix.Initialize_VP(pos,mom,(int) mcParticle.getCharge(),B[2]/dd4hep::tesla);
        edm4hep::TrackState trackState;
        trackState.D0=helix.getD0();
        trackState.phi=helix.getPhi0();
        trackState.omega=helix.getOmega();
        trackState.Z0=helix.getZ0();
        trackState.tanLambda=helix.getTanLambda();
        trackState.referencePoint=helix.getReferencePoint();
        std::array<float,15> covMatrix;
        for(int i=0;i<15;i++){covMatrix[i]=999.;}//FIXME
        trackState.covMatrix=covMatrix;

        if(m_debug){
            std::cout<<" helix  d0 "<<trackState.D0<<
                " phi "<<trackState.phi<<
                " omega "<<trackState.omega<<
                " z0 "<<trackState.Z0<<
                " tanLambda "<<trackState.tanLambda<<
                " referencePoint "<<trackState.referencePoint<< std::endl;
        }
        //new Track
        edm4hep::Track track = dcTrackCol->create();
        track.addToTrackStates(trackState);
        std::vector<double> tmp_dedx_vec;
        if(nullptr!=digiDCHitsCol) {
            //track.setType();//TODO
            //track.setChi2(trackState->getChi2());
            track.setNdf(digiDCHitsCol->size()-5);
            //track.setDEdx();//TODO
            //track.setRadiusOfInnermostHit();//TODO
            for(int i=0; i<digiDCHitsCol->size(); i++ ){
                edm4hep::TrackerHit digiDC=digiDCHitsCol->at(i);
                //if(Sim->MCParti!=current) continue;//TODO
                track.addToTrackerHits(digiDC);
                tmp_dedx_vec.push_back(digiDC.getEdx());
                if(m_WriteAna){
                    m_hit_x[m_hit] = digiDC.getPosition()[0]; 
                    m_hit_y[m_hit] = digiDC.getPosition()[1]; 
                    m_hit_z[m_hit] = digiDC.getPosition()[2]; 
                    m_hit_dedx[m_hit] = digiDC.getEdx(); 
                    m_hit ++ ;
                }
            }
        }
        int usedhit=-1;
        double track_dedx = cal_dedx_bitrunc(m_truncate, tmp_dedx_vec, usedhit);
        double track_dedx_BB = BetheBlochEquationDedx(mcParticle) ;
        float dedx_smear  = CLHEP::RandGauss::shoot(1, m_track_dedx_resolution);
        track.setDEdx(track_dedx*dedx_smear);
        if(m_WriteAna){
            m_track_dedx = track_dedx*dedx_smear ;
            m_track_dedx_BB = track_dedx_BB ;
            m_track_px    = mom[0] ;
            m_track_py    = mom[1] ;
            m_track_pz    = mom[2] ;
            //std::cout<<"m_track_mass="<<mcParticle.getMass()<<",m_track_pid="<<mcParticle.getPDG()<<std::endl;
            m_track_mass  = mcParticle.getMass() ;
            m_track_pid   = mcParticle.getPDG() ;
            m_mc_px    = mom_truth[0] ;
            m_mc_py    = mom_truth[1] ;
            m_mc_pz    = mom_truth[2] ;
        }
        //TODO
        //new ReconstructedParticle
        //recParticle->setType();
        //dcRecParticle->setEnergy();
        //edm4hep::Vector3f momVec3(helix.getMomentum()[0],
        //        helix.getMomentum()[1],helix.getMomentum()[2]);
        //recParticle->setMomentum(momVec3);
    }
    if(m_WriteAna){
        StatusCode status = m_tuple->write();
        if ( status.isFailure() ) {
          error() << "    Cannot fill N-tuple:" << long( m_tuple ) << endmsg;
          return StatusCode::FAILURE;
        }
    }

    if(m_debug) std::cout<<"Output DCTrack size="<<dcTrackCol->size()<<std::endl;
    return StatusCode::SUCCESS;
}

StatusCode TruthTrackerAlg::finalize()
{
    return GaudiAlgorithm::finalize();
}

double TruthTrackerAlg::cal_dedx_bitrunc(float truncate, std::vector<double> phlist, int & usedhit )
{   
    sort(phlist.begin(),phlist.end());
    int nsampl = (int)( phlist.size()*truncate );
    int smpl = (int)(phlist.size()*(truncate+0.05));
    int min_cut = (int)( phlist.size()*0.05 + 0.5 );
    double qSum = 0;
    unsigned i = 0;
    for(std::vector<double>::iterator ql= phlist.begin();ql!=phlist.end();ql++)
    {   
        i++;
        if(i<= smpl && i>=min_cut ) qSum += (*ql);
    }
    double trunc=-99;
    usedhit = smpl-min_cut+1;
    if(usedhit>0)  trunc=qSum/usedhit;
    
    return trunc;
}

double TruthTrackerAlg::BetheBlochEquationDedx(const edm4hep::MCParticle& mcp)
{
    int z = mcp.getCharge();
    if (z == 0) return 0;
    float m_I = m_material_Z*10;  // Approximate

    double M = mcp.getMass();
    double gammabeta=sqrt(mcp.getMomentum()[0]*mcp.getMomentum()[0]+mcp.getMomentum()[1]*mcp.getMomentum()[1]+mcp.getMomentum()[2]*mcp.getMomentum()[2])/M;
    if(gammabeta<0.01)return 0;//too low momentum
    float beta = gammabeta/sqrt(1.0+pow(gammabeta,2));
    float gamma = gammabeta/beta;
    M = M*pow(10,9);//to eV
    float Tmax = 2*m_me*pow(gammabeta,2)/(1+(2*gamma*m_me/M)+pow(m_me/M,2));
    float dedx = m_K*pow(z,2)*m_material_Z*(0.5*log(2*m_me*pow(gammabeta,2)*Tmax/pow(m_I,2))-pow(beta,2))/(m_material_A*pow(beta,2));    
    dedx = dedx*CLHEP::RandGauss::shoot(m_dedx_scale, m_dedx_resolution);
    return dedx; // (CLHEP::MeV/CLHEP::cm) / (CLHEP::g/CLHEP::cm3)
}
