#ifndef FWCORE_CONVERTOR_H
#define FWCORE_CONVERTOR_H

// LCIO

#include <map>
#include <iostream>
#include <string>
#include "lcio.h"
// #include "IO/LCReader.h"
// #include "EVENT/LCCollection.h"
#include "EVENT/MCParticle.h"
#include "EVENT/SimTrackerHit.h"
#include "EVENT/SimCalorimeterHit.h"
#include "plcio/SimTrackerHit.h"
#include "plcio/SimTrackerHitCollection.h"
#include "plcio/SimCalorimeterHit.h"
#include "plcio/SimCalorimeterHitCollection.h"
#include "plcio/MCParticle.h"
#include "plcio/MCParticleCollection.h"

#include <utility>
// Forward declarations

/** @class LCIO2Plcio LCIO2Plcio.h
 *
 *   An LCIO2Plcio for Converting from LCCollection to plcio collection;
 *
 *  @author jhZou, gjCao 
 */

// typedef plcio::MCParticleCollection* (*fptr) (EVENT::LCCollection*);
typedef podio::CollectionBase* (*fptr) (EVENT::LCCollection*);
typedef std::vector<std::pair<EVENT::LCCollection*, podio::CollectionBase*>> CollectionsVec;
typedef std::map<std::string, CollectionsVec> CollectionsMap;

class LCIO2Plcio{
public:

  /// Standard Constructor
  LCIO2Plcio();
  LCIO2Plcio(EVENT::LCCollection*);

  /// Standard Destructor
  virtual ~LCIO2Plcio(){}

  void test(){ printf("MYTESTFUC\n"); }
  void clear(){ map_cols.clear(); };

//  plcio::MCParticleCollection* Convertor_getPlcio(EVENT::LCCollection*);
  podio::CollectionBase* Convertor_getPlcio(EVENT::LCCollection*);
  static podio::CollectionBase* Convertor_MCParticle(EVENT::LCCollection*);
  static podio::CollectionBase* Convertor_LCRunHeader(EVENT::LCCollection*);
  static podio::CollectionBase* Convertor_SimTrackerHit(EVENT::LCCollection*);
  static podio::CollectionBase* Convertor_SimCalorimeterHit(EVENT::LCCollection*);

  static void void_Core_MCParticle(EVENT::MCParticle*, plcio::MCParticle&);
  static plcio::MCParticleCollection* Core_MCParticle(EVENT::LCCollection*);

  bool isReady(const std::string&);
  void setPlcioMCParticleCollection(plcio::MCParticleCollection*);
  void setLCIOMCParticleCollection(EVENT::LCCollection*);

private:
  std::string TypeName;
  // maintain a log vec about data read;
  std::vector<std::string> vec_Types;

  // maintain a map from keyword to function pointer.
  std::map<std::string, fptr> map_cvt;
  static CollectionsMap map_cols;

//  plcio::MCParticleCollection* hitcol_pl;
//  EVENT::LCCollection* hitcol_lc;
};
#endif  // CORE_CONVERTOR_H
