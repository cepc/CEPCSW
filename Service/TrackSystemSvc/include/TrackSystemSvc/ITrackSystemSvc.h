#ifndef ITrackSystemSvc_h
#define ITrackSystemSvc_h

#include "GaudiKernel/IService.h"
#include "IMarlinTrkSystem.h"

// ITrackSystemSvc is the interface between Gaudi and Track.

class ITrackSystemSvc: virtual public IService {
 public:
  DeclareInterfaceID(ITrackSystemSvc, 0, 1); // major/minor version
  
  virtual ~ITrackSystemSvc() = default;
  
  //Get the track manager
  virtual MarlinTrk::IMarlinTrkSystem* getTrackSystem(void* address=0) = 0;
  
  virtual void removeTrackSystem(void* address=0) = 0;
};

#endif
