#ifndef Calo_DIGI_ALG_H
#define Calo_DIGI_ALG_H

#include "FWCore/DataHandle.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "edm4hep/SimCalorimeterHitConst.h"
#include "edm4hep/SimCalorimeterHit.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"

#include <DDRec/DetectorData.h>
#include <DDRec/CellIDPositionConverter.h>
#include "DetInterface/IGeoSvc.h"

/** ======= PlanarDigiProcessor / CaloDigiAlg ========== <br>
 * Creates TrackerHits from SimTrackerHits, smearing them according to the input parameters. 
 * The SimTrackerHits should come from a planar detector like VXD, SIT, SET or FTD.
 * 
 * WARNING: this processor depends on correctly set CellID0s and is NOT backwards compatible to 
 * SimTrackerHit output with wrong CellID0s!!!
 * 
 * The positions of "digitized" TrackerHits are obtained by gaussian smearing positions
 * of SimTrackerHits in u and v direction. 
 * <h4>Input collections and prerequisites</h4> 
 * Processor requires a collection of SimTrackerHits <br>
 * <h4>Output</h4>
 * Processor produces collection of smeared TrackerHits<br>
 * @param SimTrackHitCollectionName The name of input collection of SimTrackerHits <br>
 * (default name VXDCollection) <br>
 * @param TrackerHitCollectionName The name of output collection of smeared TrackerHits <br>
 * (default name VTXTrackerHits) <br>
 * @param SimTrkHitRelCollection The name of the TrackerHit SimTrackerHit relation collection <br>
 * (default name VTXTrackerHitRelations) <br>
 * @param ResolutionU resolution in direction of u (in mm) <br>
 * (default value 0.004) <br>
 * @param ResolutionV Resolution in direction of v (in mm) <br>
 * (default value 0.004) <br>
 * @param IsStrip whether the hits are 1 dimensional strip measurements <br>
 * (default value false)<br>
 * <br>
 * 
 */



class CaloDigiAlg : public GaudiAlgorithm
{
 
public:
 
  CaloDigiAlg(const std::string& name, ISvcLocator* svcLoc);
 
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
 
  std::string GetCode(const edm4hep::SimCalorimeterHit *const pHit) const;
 
protected:

  const std::vector<dd4hep::rec::LayeredCalorimeterData::Layer>* m_layerLayout;
 
  SmartIF<IGeoSvc> m_geosvc;
  typedef std::vector<float> FloatVec;

  int _nEvt ;
  float m_length;
  dd4hep::rec::CellIDPositionConverter* m_cellIDConverter;
 
  Gaudi::Property<float> m_scale{ this, "Scale", 1 };

  // Input collections
  DataHandle<edm4hep::SimCalorimeterHitCollection> r_SimCaloCol{"SimCaloCol", Gaudi::DataHandle::Reader, this};
  // Output collections
  DataHandle<edm4hep::CalorimeterHitCollection>    w_DigiCaloCol{"DigiCaloCol", Gaudi::DataHandle::Writer, this};
};

#endif
