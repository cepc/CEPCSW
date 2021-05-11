#ifndef ISensDetTool_h
#define ISensDetTool_h
/*
 * ISensDetTool is a tool to configure and create a Geant4's sensitive detector.
 * After create the SD, the Geant4 will take ownership of the SD.
 *
 * This tool is used to replace the DDG4's Geant4SensitiveDetector.
 * It will be invoked in ConstructSDandField().
 *
 * -- 12 June 2020, Tao Lin <lintao@ihep.ac.cn>
 */

#include "GaudiKernel/IAlgTool.h"

class G4VSensitiveDetector;

class ISensDetTool: virtual public IAlgTool {
public:

    DeclareInterfaceID(ISensDetTool, 0, 1);

    virtual ~ISensDetTool() {};

    virtual G4VSensitiveDetector* createSD(const std::string& name) = 0;

};

#endif
