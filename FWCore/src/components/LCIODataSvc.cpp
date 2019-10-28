#include "FWCore/LCIODataSvc.h"
#include "LCIO2Plcio.h"
#include "GaudiKernel/IConversionSvc.h"
#include "GaudiKernel/IEventProcessor.h"
#include "GaudiKernel/ISvcLocator.h"

#include "IOIMPL/LCFactory.h"
#include "FWCore/DataWrapper.h"
#include "lcio.h"
#include "plcio/MCParticleCollection.h"
#include "plcio/MCParticle.h"
#include "plcio/EventHeaderCollection.h"
#include "EVENT/MCParticle.h"

#include "TTree.h"

typedef std::vector<lcio::MCParticle*> MCParticleVec ;  

DECLARE_SERVICE_FACTORY(LCIODataSvc)
/// Service initialisation
StatusCode LCIODataSvc::initialize() {
  // Nothing to do: just call base class initialisation
  StatusCode status = DataSvc::initialize();
  ISvcLocator* svc_loc = serviceLocator();

  // Attach data loader facility
  m_cnvSvc = svc_loc->service("EventPersistencySvc");
  status = setDataLoader(m_cnvSvc);

  if ( name() != "EventDataSvc" ) {
    service("EventDataSvc", m_pIDP, true);
    if ( m_pIDP == nullptr ) {
      error() << "Could not get the EventDataSvc instance" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  m_reader = IOIMPL::LCFactory::getInstance()->createLCReader();

  if (m_filename != "") {
    m_filenames.push_back(m_filename);
  }

  if (m_filenames.size() > 0) {
    m_reader->open(m_filenames[0]);
    m_eventMax = m_reader->getNumberOfEvents();
  }

  return status;
}
/// Service reinitialisation
StatusCode LCIODataSvc::reinitialize() {
  // Do nothing for this service
  return StatusCode::SUCCESS;
}
/// Service finalization
StatusCode LCIODataSvc::finalize() {
  m_reader->close();
  delete m_reader;
  m_reader = nullptr;
  m_cnvSvc = 0;  // release
  DataSvc::finalize().ignore();
  return StatusCode::SUCCESS;
}

StatusCode LCIODataSvc::clearStore() {
  for (auto& collNamePair : m_collections) {
    if (collNamePair.second != nullptr) {
      collNamePair.second->clear();
    }
  }
  for (auto& collNamePair : m_readCollections) {
    if (collNamePair.second != nullptr) {
      collNamePair.second->clear();
    }
  }
  DataSvc::clearStore().ignore();
  m_collections.clear();
  m_readCollections.clear();
  return StatusCode::SUCCESS;
}

void LCIODataSvc::endOfRead() {
  if (m_eventMax != -1) {
    // m_provider.clearCaches();
    // m_reader.endOfEvent();
    if ( ++m_eventNum >= m_eventMax ) {
      if ( ++m_fileIndex < m_filenames.size() ) {  // move to next file
        m_reader->close();
        m_reader->open( m_filenames[m_fileIndex] );
        m_eventMax += m_reader->getNumberOfEvents();
      }
      else {  // reach to the end of the file list
        info() << "Reached end of file with event " << m_eventMax << endmsg;
        IEventProcessor* eventProcessor;
        service("ApplicationMgr", eventProcessor);
        eventProcessor->stopRun();
      }
    }
  }
  evt = nullptr;
}

void LCIODataSvc::setCollectionIDs(podio::CollectionIDTable* collectionIds) {
  if (m_collectionIDs != nullptr) {
    delete m_collectionIDs;
  }
  m_collectionIDs = collectionIds;
}

/// Standard Constructor
LCIODataSvc::LCIODataSvc(const std::string& name, ISvcLocator* svc)
    : DataSvc(name, svc), m_collectionIDs(new podio::CollectionIDTable()) {

  m_eventDataTree = new TTree("events", "Events tree");
  declareProperty("inputs", m_filenames = {}, "Names of the files to read");
  declareProperty("input", m_filename = "", "Name of the file to read");

    }

/// Standard Destructor
LCIODataSvc::~LCIODataSvc() {}


StatusCode LCIODataSvc::readCollection(const std::string& collName, int collectionID) {

  StatusCode stat = StatusCode::SUCCESS;
  podio::CollectionBase* collection(nullptr);

  if( evt == nullptr ){
    evt = m_reader->readNextEvent();
    cvtor.clear();

    // basicly set EventHeader;
    pl_evtcol = new plcio::EventHeaderCollection();
    plcio::EventHeader evt_header;
    evt_header = (plcio::EventHeader) pl_evtcol->create();
    evt_header.setEventNumber( evt->getEventNumber() );
    evt_header.setRunNumber( evt->getRunNumber() );
    evt_header.setTimeStamp( evt->getTimeStamp() );
    evt_header.setDetectorName( evt->getDetectorName() );

    // wrap event header collection into Data service;
    auto wrapper = new DataWrapper<podio::CollectionBase>;
    int id = m_collectionIDs->add("EventHeader");
    pl_evtcol->setID(id);
    wrapper->setData(pl_evtcol);

    if ( m_pIDP ) {
      m_pIDP->registerObject("EventHeader", wrapper);
    }
    else {
      m_readCollections.emplace_back(std::make_pair("EventHeader", pl_evtcol));
      DataSvc::registerObject("EventHeader", wrapper);
    }
  }

  debug() << "reading collection name: " << collName  << "." << endmsg;
  EVENT::LCCollection* lc_col;
  std::vector<std::string> vec_colns = *evt->getCollectionNames();
  std::vector<std::string>::iterator it = find(vec_colns.begin(), vec_colns.end(), collName);
  if( it != vec_colns.end() ){
    lc_col = evt->getCollection(collName); 
  }
  else
    return stat;
//  debug() << "Got collection: " << collName  << "." << endmsg;

  std::string TypeName = lc_col->getTypeName();
//  if( !exist_MCP && (TypeName == "MCParticle") ) exist_MCP = true;
//  if( TypeName == "MCParticle" ){
//    if( !exist_MCP ) exist_MCP = true;
//    mcpcol_pl = LCIO2Plcio::Core_MCParticle(lc_col);
//    LCIO2Plcio::setLCIOMCParticleCollection(mcpcol_lc);
//    LCIO2Plcio::setPlcioMCParticleCollection(mcpcol_pl);
//  }
  cvtor.setCollName(collName);
  collection = cvtor.Convertor_getPlcio( lc_col );
  pl_evtcol->at(0)->addCollectionName(collName);
  pl_evtcol->at(0)->addCollectionType(TypeName);

  auto wrapper = new DataWrapper<podio::CollectionBase>;
  int id = m_collectionIDs->add(collName);
  collection->setID(id);
  wrapper->setData(collection);

//  info() << "readCollection completed." << endmsg;

  if ( m_pIDP ) {
    stat = m_pIDP->registerObject(collName, wrapper);
  }
  else {
    m_readCollections.emplace_back(std::make_pair(collName, collection));
    stat = DataSvc::registerObject(collName, wrapper);
  }
  return stat;
}

StatusCode LCIODataSvc::registerObject(const std::string& fullPath, DataObject* pObject) {
  DataWrapperBase* wrapper = dynamic_cast<DataWrapperBase*>(pObject);
  if (wrapper != nullptr) {
    podio::CollectionBase* coll = wrapper->collectionBase();
    if (coll != nullptr) {
      size_t pos = fullPath.find_last_of("/");
      std::string shortPath(fullPath.substr(pos + 1, fullPath.length()));
      int id = m_collectionIDs->add(shortPath);
      coll->setID(id);
      m_collections.emplace_back(std::make_pair(shortPath, coll));
    }
  }
  return DataSvc::registerObject(fullPath, pObject);
}
