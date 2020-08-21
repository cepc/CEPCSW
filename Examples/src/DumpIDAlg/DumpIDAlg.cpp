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
    m_geosvc = service<IGeoSvc>("GeoSvc");
    if (!m_geosvc) {
        error() << "Failed to find GeoSvc." << endmsg;
        return StatusCode::FAILURE;
    }
    m_dd4hep_geo = m_geosvc->lcdd();
    if (!m_dd4hep_geo) {
        error() << "failed to retrieve dd4hep_geo: " << m_dd4hep_geo << endmsg;
        return StatusCode::FAILURE;
    }

    // get the DD4hep readout
    auto readouts = m_dd4hep_geo->readouts();
    const std::string name_readout = "EcalBarrelCollection";
    if (readouts.find(name_readout) != readouts.end()) {
        dd4hep::Readout readout = m_dd4hep_geo->readout(name_readout);

        auto m_idspec = readout.idSpec();

        info() << "The idspec is " << m_idspec.fieldDescription() << " for " << name_readout << endmsg;

        dd4hep::DDSegmentation::BitFieldCoder* decoder = m_idspec.decoder();
        
        m_decoder = decoder;
    }

    if (!m_decoder) {
        error() << "Failed to get the decoder. " << endmsg;
        return StatusCode::FAILURE;
    }

    return GaudiAlgorithm::initialize();
}

StatusCode DumpIDAlg::execute()
{


    auto ecalBarrelCol = m_EcalBarrelCol.get();
    for (auto calohit: *ecalBarrelCol) {
        auto cellid = calohit.getCellID();

        int id_system = m_decoder->get(cellid, "system");
        int id_module = m_decoder->get(cellid, "module");
        int id_stave = m_decoder->get(cellid, "stave");
        int id_tower = m_decoder->get(cellid, "tower");
        int id_layer = m_decoder->get(cellid, "layer");
        int id_wafer = m_decoder->get(cellid, "wafer");
        int id_cellX = m_decoder->get(cellid, "cellX");
        int id_cellY = m_decoder->get(cellid, "cellY");

        info() << "Calo hit cellid: " << cellid
               << " system: " << id_system
               << " module: " << id_module
               << " stave: " << id_stave
               << " tower: " << id_tower
               << " layer: " << id_layer
               << " wafer: " << id_wafer
               << " cellX: " << id_cellX
               << " cellY: " << id_cellY
               << endmsg;
    }
    return StatusCode::SUCCESS;
}

StatusCode DumpIDAlg::finalize()
{

    return GaudiAlgorithm::finalize();
}




