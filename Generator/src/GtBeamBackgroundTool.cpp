#include "GtBeamBackgroundTool.h"
#include "IBeamBackgroundFileParser.h"

#include "BeamBackgroundFileParserV0.h"
#include "GuineaPigPairsFileParser.h"

#include "TVector3.h" // for rotation
DECLARE_COMPONENT(GtBeamBackgroundTool)

StatusCode GtBeamBackgroundTool::initialize() {
    // check the consistency of the properties

    // create the instances of the background parsers

    for (auto& [label, inputfn]: m_inputmaps) {
        std::string format = "BeamBackgroundFileParserV0";

        auto itFormat = m_formatmaps.find(label);
        if (itFormat != m_formatmaps.end()) {
            format = itFormat->second;
        }

        if (format == "BeamBackgroundFileParserV0") {
            init_BeamBackgroundFileParserV0(label, inputfn);
        } else if (format == "GuineaPigPairsFileParser") {
            init_GuineaPigPairsFileParser(label, inputfn);
        } else {
            init_BeamBackgroundFileParserV0(label, inputfn);
        }
        
    }

    // check the size
    if (m_beaminputs.empty()) {
        error() << "Empty Beam Background File Parser. " << endmsg;
        return StatusCode::FAILURE;
    }

    return StatusCode::SUCCESS;
}

StatusCode GtBeamBackgroundTool::finalize() {
    return StatusCode::SUCCESS;
}


bool GtBeamBackgroundTool::mutate(MyHepMC::GenEvent& event) {
    if (m_beaminputs.empty()) {
        return false;
    }
    // TODO: should sample according to the rates
    // dummy: get one and stop to generate
    for (auto& [label, parser]: m_beaminputs) {
        IBeamBackgroundFileParser::BeamBackgroundData beamdata;
        auto isok = parser->load(beamdata);
        if (not isok) {
            error() << "Failed to load beam background data from the parser " << label << endmsg;
            return false;
        }
        // fill the value
        float charge = beamdata.pdgid == 11 ? -1: 1;

        TVector3 pos(beamdata.x,beamdata.y,beamdata.z);
        TVector3 mom(beamdata.px,beamdata.py,beamdata.pz);

        auto itrot = m_rotYmaps.find(label);
        if (itrot  != m_rotYmaps.end() ) {
            info() << "Apply rotation along Y " << itrot->second << endmsg;

            pos.RotateY(itrot->second);
            mom.RotateY(itrot->second);
        }

        // create the MC particle
        auto mcp = event.m_mc_vec.create();
        mcp.setPDG(beamdata.pdgid);
        mcp.setGeneratorStatus(1);
        mcp.setSimulatorStatus(1);
        mcp.setCharge(static_cast<float>(charge));
        mcp.setTime(beamdata.t);
        mcp.setMass(beamdata.mass);
        mcp.setVertex(edm4hep::Vector3d(pos.X(), pos.Y(), pos.Z())); 
        mcp.setMomentum(edm4hep::Vector3f(mom.X(), mom.Y(), mom.Z()));

    }
    
    return true;
}

bool GtBeamBackgroundTool::finish() {
    return true;
}

bool GtBeamBackgroundTool::configure_gentool() {

    return true;
}

bool GtBeamBackgroundTool::init_BeamBackgroundFileParserV0(const std::string& label,
                                                           const std::string& inputfn) {
    double beamE = 120.;
    auto itBeamE = m_Ebeammaps.find(label);
    if (itBeamE != m_Ebeammaps.end()) {
        beamE = itBeamE->second;
    }
    info() << "Initializing beam background ... "
           << label << " "
           << beamE << " "
           << inputfn
           << endmsg;
    m_beaminputs[label] = std::make_shared<BeamBackgroundFileParserV0>(inputfn, 11, beamE);

    return true;
}

bool GtBeamBackgroundTool::init_GuineaPigPairsFileParser(const std::string& label,
                                                         const std::string& inputfn) {

    m_beaminputs[label] = std::make_shared<GuineaPigPairsFileParser>(inputfn);

    return true;
}
