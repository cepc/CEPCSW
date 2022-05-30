#ifndef CylinderDigiAlg_H
#define CylinderDigiAlg_H

#include "k4FWCore/DataHandle.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/MCRecoTrackerAssociationCollection.h"

#include <DDRec/DetectorData.h>
#include "DetInterface/IGeomSvc.h"

class CylinderDigiAlg : public GaudiAlgorithm{
 public:
  
  CylinderDigiAlg(const std::string& name, ISvcLocator* svcLoc);
 
  virtual StatusCode initialize() ;
  virtual StatusCode execute() ; 
  virtual StatusCode finalize() ;
 
protected:
  
  SmartIF<IGeomSvc> m_geosvc;
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;

  // Input collections
  DataHandle<edm4hep::SimTrackerHitCollection>            m_inputColHdls{"DriftChamberHitsCollection", Gaudi::DataHandle::Reader, this};
  // Output collections
  DataHandle<edm4hep::TrackerHitCollection>               m_outputColHdls{"DCTrackerHits", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::MCRecoTrackerAssociationCollection> m_assColHdls{"DCTrackerHitAssociationCollection", Gaudi::DataHandle::Writer, this};

  Gaudi::Property<float> m_resRPhi{this, "ResRPhi", 0.1};
  Gaudi::Property<float> m_resZ   {this, "ResZ", 2.828};

  int m_nEvt=0;
};
#endif
