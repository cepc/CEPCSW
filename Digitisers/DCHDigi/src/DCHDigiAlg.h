#ifndef DCH_DIGI_ALG_H
#define DCH_DIGI_ALG_H

#include "FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/MCRecoTrackerAssociationCollection.h"

#include <DDRec/DetectorData.h>
#include <DDRec/CellIDPositionConverter.h>
#include "DetInterface/IGeoSvc.h"




class DCHDigiAlg : public GaudiAlgorithm
{
 
public:
 
  DCHDigiAlg(const std::string& name, ISvcLocator* svcLoc);
 
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual StatusCode initialize() ;
 
  /** Called for every event - the working horse.
   */
  virtual StatusCode execute() ; 
 
  /** Called after data processing for clean up.
   */
  virtual StatusCode finalize() ;
 
protected:

  SmartIF<IGeoSvc> m_geosvc;
  typedef std::vector<float> FloatVec;
  int _nEvt ;

  //float m_length;
  dd4hep::rec::CellIDPositionConverter* m_cellIDConverter;
 
  //Gaudi::Property<float> m_scale     { this, "Scale", 1 };
  //Gaudi::Property<float> m_resolution{ this, "Res", 0.01 };

  // Input collections
  DataHandle<edm4hep::SimTrackerHitCollection> r_SimDCHCol{"DriftChamberHitsCollection", Gaudi::DataHandle::Reader, this};
  // Output collections
  DataHandle<edm4hep::TrackerHitCollection>    w_DigiDCHCol{"DigiDCHitsCollection", Gaudi::DataHandle::Writer, this};
  DataHandle<edm4hep::MCRecoTrackerAssociationCollection>    w_AssociationCol{"DCHitAssociationCollection", Gaudi::DataHandle::Writer, this};
};

#endif
