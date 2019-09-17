#include "LCIOInput.h"

#include "TFile.h"
#include "TROOT.h"

#include "FWCore/DataWrapper.h"
#include "FWCore/LCIODataSvc.h"

DECLARE_COMPONENT(LCIOInput)

LCIOInput::LCIOInput(const std::string& name, ISvcLocator* svcLoc) : GaudiAlgorithm(name, svcLoc) {}

StatusCode LCIOInput::initialize() {
  if (GaudiAlgorithm::initialize().isFailure()) return StatusCode::FAILURE;

  // check whether we have the LCIOEvtSvc active
  m_LCIODataSvc = dynamic_cast<LCIODataSvc*>(evtSvc().get());
  if (nullptr == m_LCIODataSvc) return StatusCode::FAILURE;

  auto idTable = m_LCIODataSvc->getCollectionIDs();
  for (auto& name : m_collectionNames) {
    debug() << "Finding collection " << name << " in collection registry." << endmsg;
/*
    if (!idTable->present(name)) {
      error() << "Requested product " << name << " not found." << endmsg;
      return StatusCode::FAILURE;
    }
*/
    m_collectionIDs.push_back(idTable->add(name));
  }
  return StatusCode::SUCCESS;
}

StatusCode LCIOInput::execute() {
  size_t cntr = 0;
  // Re-create the collections from ROOT file
  for (auto& id : m_collectionIDs) {
    const std::string& collName = m_collectionNames.value().at(cntr++);
    debug() << "Registering collection to read " << collName << " with id " << id << endmsg;
    if (m_LCIODataSvc->readCollection(collName, id).isFailure()) {
      return StatusCode::FAILURE;
    }
  }
  // Tell data service that we are done with requested collections
  m_LCIODataSvc->endOfRead();
  return StatusCode::SUCCESS;
}

StatusCode LCIOInput::finalize() {
  if (GaudiAlgorithm::finalize().isFailure()) return StatusCode::FAILURE;
  return StatusCode::SUCCESS;
}
