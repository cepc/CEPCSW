#include "RecDCHDedxAlg.h"
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

DECLARE_COMPONENT(RecDCHDedxAlg)

RecDCHDedxAlg::RecDCHDedxAlg(const std::string& name, ISvcLocator* svcLoc)
: GaudiAlgorithm(name, svcLoc)
{
    declareProperty("DCHitAssociationCollection", m_dcHitAssociationCol, "Handle of association collection");
    declareProperty("DCTrackCollection", m_dcTrackCol,"Handle of input Track collection");
}

StatusCode RecDCHDedxAlg::initialize()
{
    if(m_WriteAna){

        NTuplePtr nt( ntupleSvc(), "MyTuples/Dedx_evt" );
        if ( nt ) m_tuple = nt;
        else {
            m_tuple = ntupleSvc()->book( "MyTuples/Dedx_evt", CLID_ColumnWiseTuple, "Dedx" );
            if ( m_tuple ) {
              m_tuple->addItem( "N_track"       , m_n_track, 0, 1000 ).ignore();
              m_tuple->addItem( "track_dedx"    , m_n_track, m_track_dedx     ).ignore();
              m_tuple->addItem( "track_dedx_BB" , m_n_track, m_track_dedx_BB  ).ignore();
              m_tuple->addItem( "track_px"      , m_n_track, m_track_px       ).ignore();
              m_tuple->addItem( "track_py"      , m_n_track, m_track_py       ).ignore();
              m_tuple->addItem( "track_pz"      , m_n_track, m_track_pz       ).ignore();
              m_tuple->addItem( "track_mass"    , m_n_track, m_track_mass     ).ignore();
              m_tuple->addItem( "track_pid"     , m_n_track, m_track_pid      ).ignore();
            }
            else { // did not manage to book the N tuple....
              error() << "    Cannot book N-tuple:" << long( m_tuple ) << endmsg;
              return StatusCode::FAILURE;
            }
        }
    }

    return GaudiAlgorithm::initialize();
}


StatusCode RecDCHDedxAlg::execute()
{

    if(m_WriteAna){
        m_n_track = 0;
    }
    if(m_method==1){// get the primary mc particle of the reconstructed track and do dedx sampling using beta-gamma curve  

        m_dedx_simtool = ToolHandle<IDedxSimTool>(m_dedx_sim_option.value());
        if (!m_dedx_simtool) {
            error() << "Failed to find dedx sampling tool :"<<m_dedx_sim_option.value()<< endmsg;
            return StatusCode::FAILURE;
        }
        
        const edm4hep::TrackCollection* trkCol=nullptr;
        trkCol=m_dcTrackCol.get();
        if(nullptr==trkCol){
            throw ("Error: TrackCollection not found.");
            return StatusCode::SUCCESS;
        }
        edm4hep::TrackCollection* ptrkCol = const_cast<edm4hep::TrackCollection*>(trkCol);

        const edm4hep::MCRecoTrackerAssociationCollection* assoCol=nullptr;
        assoCol = m_dcHitAssociationCol.get();
        if(nullptr==assoCol){
            throw ("Error: MCRecoTrackerAssociationCollection not found.");
            return StatusCode::SUCCESS;
        }
        for(unsigned j=0; j<ptrkCol->size(); j++){
            std::map< edm4hep::ConstMCParticle, int > map_mc_count;
            auto tmp_track = ptrkCol->at(j);
            for(unsigned k=0; k< tmp_track.trackerHits_size(); k++){
                for(unsigned z=0; z< assoCol->size(); z++){
                    if(assoCol->at(z).getRec().id() != tmp_track.getTrackerHits(k).id()) continue;
                    
                    auto tmp_mc = assoCol->at(z).getSim().getMCParticle();
                    if(tmp_mc.parents_size()!=0) continue;
                    if(map_mc_count.find(tmp_mc) != map_mc_count.end()) map_mc_count[tmp_mc] += 1;
                    else                                                map_mc_count[tmp_mc] = 1;
                    
                }
            }
            edm4hep::MCParticle tmp_mc;
            int max_cout = 0;
            for(auto iter = map_mc_count.begin(); iter != map_mc_count.end(); iter++ ){
                if(iter->second > max_cout){
                    max_cout = iter->second;
                    //tmp_mc = iter->first;
                    tmp_mc.setMass(iter->first.getMass());
                    tmp_mc.setCharge(iter->first.getCharge());
                    tmp_mc.setMomentum(iter->first.getMomentum());
                    if(m_debug){
                        std::cout<<"mass="<<iter->first.getMass()<<",charge="<<iter->first.getCharge()<<",mom_x"<<iter->first.getMomentum()[0]<<",mom_y"<<iter->first.getMomentum()[1]<<",mom_z"<<iter->first.getMomentum()[2]<<std::endl;
                    }
                }
            }
            double dedx_ori = tmp_track.getDEdx();
            double dedx = 0.0;
            dedx = m_dedx_simtool->dedx(tmp_mc);
            float dedx_smear  = CLHEP::RandGauss::shoot(m_scale, m_resolution);
            if(m_debug) std::cout<<"ori dedx="<<dedx_ori<<", dedx from sampling ="<<dedx<<",smear="<<dedx_smear<<",final="<<dedx*dedx_smear<<",scale="<<m_scale<<",res="<<m_resolution<<std::endl;
            ptrkCol->at(j).setDEdx(dedx*dedx_smear);
            if(m_WriteAna){
                m_track_dedx[m_n_track]= dedx*dedx_smear;
                m_track_dedx_BB[m_n_track]= dedx;
                m_track_px     [m_n_track]= tmp_mc.getMomentum()[0];
                m_track_py     [m_n_track]= tmp_mc.getMomentum()[1];
                m_track_pz     [m_n_track]= tmp_mc.getMomentum()[2];
                m_track_mass   [m_n_track]= tmp_mc.getMass();
                m_track_pid    [m_n_track]= tmp_mc.getPDG();
                m_n_track ++;
            }
        }
   }//method 1
   if(m_method==2){// get the primary mc particle of the reconstructed track and do dndx sampling according the mc beta-gamma 

       m_dedx_simtool = ToolHandle<IDedxSimTool>(m_dedx_sim_option.value());
       if (!m_dedx_simtool) {
           error() << "Failed to find dedx sampling tool :"<<m_dedx_sim_option.value()<< endmsg;
           return StatusCode::FAILURE;
       }
       
       const edm4hep::TrackCollection* trkCol=nullptr;
       trkCol=m_dcTrackCol.get();
       if(nullptr==trkCol){
           throw ("Error: TrackCollection not found.");
           return StatusCode::SUCCESS;
       }
       edm4hep::TrackCollection* ptrkCol = const_cast<edm4hep::TrackCollection*>(trkCol);

       const edm4hep::MCRecoTrackerAssociationCollection* assoCol=nullptr;
       assoCol = m_dcHitAssociationCol.get();
       if(nullptr==assoCol){
           throw ("Error: MCRecoTrackerAssociationCollection not found.");
           return StatusCode::SUCCESS;
       }
       for(unsigned j=0; j<ptrkCol->size(); j++){
           std::map< edm4hep::ConstMCParticle, int > map_mc_count;
           auto tmp_track = ptrkCol->at(j);
           for(unsigned k=0; k< tmp_track.trackerHits_size(); k++){
               for(unsigned z=0; z< assoCol->size(); z++){
                   if(assoCol->at(z).getRec().id() != tmp_track.getTrackerHits(k).id()) continue;
                   
                   auto tmp_mc = assoCol->at(z).getSim().getMCParticle();
                   if(tmp_mc.parents_size()!=0) continue;
                   if(map_mc_count.find(tmp_mc) != map_mc_count.end()) map_mc_count[tmp_mc] += 1;
                   else                                                map_mc_count[tmp_mc] = 1;
                   
               }
           }
           edm4hep::MCParticle tmp_mc;
           int max_cout = 0;
           for(auto iter = map_mc_count.begin(); iter != map_mc_count.end(); iter++ ){
               if(iter->second > max_cout){
                   max_cout = iter->second;
                   //tmp_mc = iter->first;
                   tmp_mc.setMass(iter->first.getMass());
                   tmp_mc.setCharge(iter->first.getCharge());
                   tmp_mc.setMomentum(iter->first.getMomentum());
                   if(m_debug){
                       std::cout<<"mass="<<iter->first.getMass()<<",charge="<<iter->first.getCharge()<<",mom_x"<<iter->first.getMomentum()[0]<<",mom_y"<<iter->first.getMomentum()[1]<<",mom_z"<<iter->first.getMomentum()[2]<<std::endl;
                   }
               }
           }
           double dndx = 0.0;
           double betagamma = tmp_mc.getMass() !=0 ? sqrt(tmp_mc.getMomentum()[0]*tmp_mc.getMomentum()[0] + tmp_mc.getMomentum()[1]*tmp_mc.getMomentum()[1] + tmp_mc.getMomentum()[2]*tmp_mc.getMomentum()[2] )/tmp_mc.getMass() : -1 ;
           dndx = m_dedx_simtool->dndx(betagamma);
           float dndx_smear  = CLHEP::RandGauss::shoot(m_scale, m_resolution);
           if(m_debug) std::cout<<"dndx from sampling ="<<dndx<<",smear="<<dndx_smear<<",final="<<dndx*dndx_smear<<",scale="<<m_scale<<",res="<<m_resolution<<std::endl;
           ptrkCol->at(j).setDEdx(dndx*dndx_smear);//update the dndx
           if(m_WriteAna){
               m_track_dedx[m_n_track]= dndx*dndx_smear;
               m_track_dedx_BB[m_n_track]= dndx;
               m_track_px     [m_n_track]= tmp_mc.getMomentum()[0];
               m_track_py     [m_n_track]= tmp_mc.getMomentum()[1];
               m_track_pz     [m_n_track]= tmp_mc.getMomentum()[2];
               m_track_mass   [m_n_track]= tmp_mc.getMass();
               m_track_pid    [m_n_track]= tmp_mc.getPDG();
               m_n_track ++;
           }
       }
   }//method 2
   if(m_WriteAna){
       StatusCode status = m_tuple->write();
       if ( status.isFailure() ) {
         error() << "    Cannot fill N-tuple:" << long( m_tuple ) << endmsg;
         return StatusCode::FAILURE;
       }
   }

   return StatusCode::SUCCESS;
}

StatusCode RecDCHDedxAlg::finalize()
{
    return GaudiAlgorithm::finalize();
}

double RecDCHDedxAlg::cal_dedx_bitrunc(float truncate, std::vector<double> phlist, int & usedhit )
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
/*
double RecDCHDedxAlg::BetheBlochEquationDedx(const edm4hep::MCParticle& mcp)
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
    dedx = dedx*CLHEP::RandGauss::shoot(m_scale, m_resolution);
    return dedx; // (CLHEP::MeV/CLHEP::cm) / (CLHEP::g/CLHEP::cm3)
}
*/
