#ifndef GtBeamBackgroundTool_h
#define GtBeamBackgroundTool_h

/*
 * Description:
 *   This tool is used to simulation the non-collision beam backgrounds.
 *
 *   The properties:
 *     - InputFileMap
 *         this is a map to store the label and the input filename
 *     - InputFormatMap
 *         this is a map to store the label and the input format
 *     - InputRateMap
 *         this is a map to store the label and the rate
 * 
 *     Note: the label (key) should be consistent
 *
 * About the design:
 *   IBeamBackgroundFileParser is the interface to load the next event.
 *   Different file formats should be implemented in the corresponding parsers. 
 *   The format will be used to create the corresponding instance.
 *
 * Author:
 *   Tao Lin <lintao AT ihep.ac.cn>
 */

#include <GaudiKernel/AlgTool.h>
#include <Gaudi/Property.h>
#include "IGenTool.h"
#include "IBeamBackgroundFileParser.h"

#include <vector>
#include <map>


class GtBeamBackgroundTool: public extends<AlgTool, IGenTool> {
public:
    using extends::extends;

    // Overriding initialize and finalize
    StatusCode initialize() override;
    StatusCode finalize() override;

    // IGenTool
    bool mutate(MyHepMC::GenEvent& event) override;
    bool finish() override;
    bool configure_gentool() override;


private:
    bool init_BeamBackgroundFileParserV0(const std::string& label, const std::string& inputfn);
    bool init_GuineaPigPairsFileParser(const std::string& label, const std::string& inputfn);

private:
    Gaudi::Property<std::map<std::string, std::string>> m_inputmaps{this, "InputFileMap"};
    Gaudi::Property<std::map<std::string, std::string>> m_formatmaps{this, "InputFormatMap"};
    Gaudi::Property<std::map<std::string, double>>      m_ratemaps {this, "InputRateMap"};

    // unit of beam energy: GeV
    Gaudi::Property<std::map<std::string, double>>      m_Ebeammaps{this, "InputBeamEnergyMap"};

    // unit of the rotation along Y: rad
    Gaudi::Property<std::map<std::string, double>>      m_rotYmaps {this, "RotationAlongYMap"};

private:
    std::map<std::string, std::shared_ptr<IBeamBackgroundFileParser>> m_beaminputs;

};

#endif
