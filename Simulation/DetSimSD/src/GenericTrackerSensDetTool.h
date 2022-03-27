#ifndef GenericTrackerSensDetTool_h
#define GenericTrackerSensDetTool_h

/*
 * GenericTrackerSensDetTool is used to create Time Projection Chamber SD.
 */

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolHandle.h"
#include "DetSimInterface/ISensDetTool.h"
#include "DetInterface/IGeomSvc.h"

#include "DD4hep/DD4hepUnits.h"

class GenericTrackerSensDetTool: public extends<AlgTool, ISensDetTool> {
  
 public:

  using extends::extends;

  /// Overriding initialize and finalize
  StatusCode initialize() override;
  StatusCode finalize() override;

  /// Override ISensDetTool
  virtual G4VSensitiveDetector* createSD(const std::string& name) override;

private:

  // in order to initialize SD, we need to get the lcdd()
  SmartIF<IGeomSvc> m_geosvc;

};

#endif
