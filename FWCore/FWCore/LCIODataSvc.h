#ifndef FWCORE_LCIODATASVC_H
#define FWCORE_LCIODATASVC_H

#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/IConversionSvc.h"
// LCIO
#include "podio/CollectionBase.h"
#include "podio/CollectionIDTable.h"
#include "podio/EventStore.h"
#include "podio/ROOTReader.h"

#include "IO/LCReader.h"
#include "EVENT/LCCollection.h"
#include "EVENT/MCParticle.h"
#include "plcio/MCParticleCollection.h"
#include "plcio/MCParticle.h"
#include "plcio/EventHeaderCollection.h"

#include "src/components/LCIO2Plcio.h"
#include <utility>
// Forward declarations

/** @class LCIOEvtSvc EvtDataSvc.h
 *
 *   An EvtDataSvc for LCIO classes
 *
 *  @author B. Hegner
 */
class LCIODataSvc : public DataSvc {
public:

  friend class SvcFactory<LCIODataSvc>;

  typedef std::vector<std::pair<std::string, podio::CollectionBase*>> CollRegistry;

  virtual StatusCode initialize();
  virtual StatusCode reinitialize();
  virtual StatusCode finalize();
  virtual StatusCode clearStore();

  /// Standard Constructor
  LCIODataSvc(const std::string& name, ISvcLocator* svc);

  /// Standard Destructor
  virtual ~LCIODataSvc();

  // Use DataSvc functionality except where we override
  using DataSvc::registerObject;
  /// Overriding standard behaviour of evt service
  /// Register object with the data store.
  virtual StatusCode registerObject(const std::string& fullPath, DataObject* pObject) final;

  StatusCode readCollection(const std::string& collectionName, int collectionID);

  virtual const CollRegistry& getCollections() const { return m_collections; }
  virtual const CollRegistry& getReadCollections() const { return m_readCollections; }
  virtual podio::CollectionIDTable* getCollectionIDs() { return m_collectionIDs; }

  /// Set the collection IDs (if reading a file)
  void setCollectionIDs(podio::CollectionIDTable* collectionIds);
  /// Resets caches of reader and event store, increases event counter
  void endOfRead();


  TTree* eventDataTree() {return m_eventDataTree;}

private:

  EVENT::LCEvent* evt = nullptr;
  // eventDataTree
  TTree* m_eventDataTree;
  /// LCIO reader for ROOT files
  IO::LCReader* m_reader;
  /// LCIO reader for ROOT files
  plcio::EventHeaderCollection* pl_evtcol;
  /// podio::ROOTReader m_reader;
  /// LCIO EventStore, used to initialise collections
  /// podio::EventStore m_provider;
  /// Counter of the event number
  int m_eventNum{0};
  /// Number of events in the file / to process
  int m_eventMax{-1};
  /// the current file index in the m_filenames vector
  int m_fileIndex{0};


  SmartIF<IConversionSvc> m_cnvSvc;

  // special members for podio handling
  std::vector<std::pair<std::string, podio::CollectionBase*>> m_collections;
  std::vector<std::pair<std::string, podio::CollectionBase*>> m_readCollections;
  podio::CollectionIDTable* m_collectionIDs;

protected:
//  bool exist_MCP = false;
  LCIO2Plcio cvtor;
//  EVENT::LCCollection* mcpcol_lc;
//  plcio::MCParticleCollection* mcpcol_pl;

  /// ROOT file name the input is read from. Set by option filename
  std::vector<std::string> m_filenames;
  std::string m_filename;
};
#endif  // CORE_LCIODATASVC_H
