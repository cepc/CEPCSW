#ifndef TimeProjectionChamberSensDetTool_h
#define TimeProjectionChamberSensDetTool_h

/*
 * TimeProjectionChamberSensDetTool is used to create Time Projection Chamber SD.
 */

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolHandle.h"
#include "DetSimInterface/ISensDetTool.h"
#include "DetInterface/IGeomSvc.h"

#include "DD4hep/DD4hepUnits.h"

class TimeProjectionChamberSensDetTool: public extends<AlgTool, ISensDetTool> {
  
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
  Gaudi::Property<int>    m_sdTypeOption{this, "TypeOption", 0};
  Gaudi::Property<double> m_threshold{this, "ThresholdEnergyDeposit", 32*dd4hep::eV};
  Gaudi::Property<double> m_lowPtCut{this, "LowPtCut", 10*dd4hep::MeV};
  Gaudi::Property<double> m_lowPtMaxHitSeparation{this, "LowPtMaxHitSeparation", 5*dd4hep::mm};
  Gaudi::Property<bool>   m_sameStepLimit{this, "SameStepLimit", true};
  Gaudi::Property<bool>   m_writeMCTruthForLowPtHits{this, "WriteMCTruthForLowPtHits", false};
};

#endif
