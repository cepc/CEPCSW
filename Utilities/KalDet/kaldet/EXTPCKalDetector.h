#ifndef LCTPC_EXTPCKALDETECTOR_H
#define LCTPC_EXTPCKALDETECTOR_H

#include "GearTPCKalDetector.h"

/**
 *  A backward compatibility class for GearTPCKalDetector.
 *  It basically provides a static instance of the detector which can be
 *  accessed via the GetInstance() method.
 *  In addition it provides the static GetVDrift() and GetBField(), which are used
 *  in the old code. The use of this class is highly depreciated.
 *
 *  \deprecated EXTPCKalDetector
 */

class EXTPCKalDetector: public kaldet::GearTPCKalDetector
{
  private:
  /// As this this a singleton class the constructor is private.
  EXTPCKalDetector(const gear::GearMgr& gearMgr);

public:
  /// The destructor.
  virtual ~EXTPCKalDetector();
  
  /// Static access function to the singleton instance.
  static EXTPCKalDetector * GetInstance();

  /// Returns the hard coded drift velocity of 76.e-3 mm/ns. 
  static Double_t GetVdrift() { return fgVdrift; }

  /// Static function to access the magnetic field.
  static Double_t GetBfield();

 private:
  static Double_t           fgVdrift;   //< The drift velocity.
  static EXTPCKalDetector * fgInstance; //< The singleton pointer.

  Double_t fBField; //< The magnetic field

  //  ClassDef(EXTPCKalDetector, 1)  // User defined detector class

};

#endif // LCTPC_EXTPCKALDETECTOR_H
