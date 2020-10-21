#include "DumpIDAlg.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/CaloHitContributionCollection.h"

#include "DD4hep/Detector.h"
#include "DD4hep/IDDescriptor.h"
#include "DD4hep/Plugins.h"

DECLARE_COMPONENT(DumpIDAlg)

DumpIDAlg::DumpIDAlg(const std::string& name, ISvcLocator* svcLoc)
: GaudiAlgorithm(name, svcLoc), m_dd4hep_geo(nullptr), m_decoder(nullptr)
{

}

StatusCode DumpIDAlg::initialize()
{
    m_geosvc = service<IGeomSvc>("GeomSvc");
    if (!m_geosvc) {
        error() << "Failed to find GeomSvc." << endmsg;
        return StatusCode::FAILURE;
    }
    m_dd4hep_geo = m_geosvc->lcdd();
    if (!m_dd4hep_geo) {
        error() << "failed to retrieve dd4hep_geo: " << m_dd4hep_geo << endmsg;
        return StatusCode::FAILURE;
    }

    // get the DD4hep readout
    const std::string name_readout = "EcalBarrelCollection";
    m_decoder = m_geosvc->getDecoder(name_readout);
    if (!m_decoder) {
        error() << "Failed to get the decoder. " << endmsg;
        return StatusCode::FAILURE;
    }

    //  Book N-tuple 1
    NTuplePtr nt1( ntupleSvc(), "MyTuples/1" );
    if ( nt1 ) {
        m_tuple_id = nt1;
    } else {
        m_tuple_id = ntupleSvc()->book( "MyTuples/1", CLID_RowWiseTuple, "Row-wise N-Tuple example" );
        if ( m_tuple_id ) {
            m_tuple_id->addItem( "system", m_id_system ).ignore();
            m_tuple_id->addItem( "module", m_id_module ).ignore();
            m_tuple_id->addItem( "stave", m_id_stave ).ignore();
            m_tuple_id->addItem( "tower", m_id_tower ).ignore();
            m_tuple_id->addItem( "layer", m_id_layer ).ignore();
            m_tuple_id->addItem( "wafer", m_id_wafer ).ignore();
            m_tuple_id->addItem( "cellX", m_id_cellX ).ignore();
            m_tuple_id->addItem( "cellY", m_id_cellY ).ignore();

        } else { // did not manage to book the N tuple....
            error() << "    Cannot book N-tuple:" << long( m_tuple_id ) << endmsg;
            return StatusCode::FAILURE;
        }
    }


    return GaudiAlgorithm::initialize();
}

StatusCode DumpIDAlg::execute()
{


    auto ecalBarrelCol = m_EcalBarrelCol.get();
    for (auto calohit: *ecalBarrelCol) {
        auto cellid = calohit.getCellID();

        m_id_system = m_decoder->get(cellid, "system");
        m_id_module = m_decoder->get(cellid, "module");
        m_id_stave = m_decoder->get(cellid, "stave");
        m_id_tower = m_decoder->get(cellid, "tower");
        m_id_layer = m_decoder->get(cellid, "layer");
        m_id_wafer = m_decoder->get(cellid, "wafer");
        m_id_cellX = m_decoder->get(cellid, "cellX");
        m_id_cellY = m_decoder->get(cellid, "cellY");

        info() << "Calo hit cellid: " << cellid
               << " system: " << m_id_system
               << " module: " << m_id_module
               << " stave: " << m_id_stave
               << " tower: " << m_id_tower
               << " layer: " << m_id_layer
               << " wafer: " << m_id_wafer
               << " cellX: " << m_id_cellX
               << " cellY: " << m_id_cellY
               << endmsg;

        // calculate I/J/K 

        m_tuple_id->write();

    }
    return StatusCode::SUCCESS;
}

StatusCode DumpIDAlg::finalize()
{

    return GaudiAlgorithm::finalize();
}




